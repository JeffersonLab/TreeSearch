//*-- Author :    Ole Hansen, Jefferson Lab   31-Jan-2013

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Tracker                                                       //
//                                                                           //
// Abstract base class for set of tracker planes analyzed using TreeSearch.  //
// TreeSearch is the recusrsive template matching algorithm described in     //
// DellOrso et al., NIM A 287, 436 (1990).                                   //
//                                                                           //
// Examples of actual implementations are TreeSearch::MWDC for sets of       //
// horizontal drift chambers and TreeSearch::GEMTracker for sets of GEM      //
// tracker planes.                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Tracker.h"
#include "Plane.h"
#include "Projection.h"
#include "Hit.h"
#include "Road.h"
#include "Helper.h"

#include "THaDetMap.h"
#include "THaEvData.h"
#include "THaTrack.h"

#include "TString.h"
#include "TMath.h"
#include "THashTable.h"
#include "TVector2.h"
#include "TDecompChol.h"
#include "TSystem.h"
#include "TThread.h"
#include "TCondition.h"
#include "TMutex.h"
#include "TBits.h"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <string>
#include <stdexcept>

#ifdef TESTCODE
#include "TStopwatch.h"
#endif

using namespace std;
typedef string::size_type ssiz_t;

namespace {

// Helper classes for describing a DAQ hardware module and using it with a
// THashTable. CSpair can be used to find elements in the THashTable w/o
// creating large dummy objects.
//TODO: replace in favor or retrieving the information from db_cratemap
class CSpair : public TObject {
public:
  CSpair( UShort_t crate, UShort_t slot ) : fCrate(crate), fSlot(slot) {}
  virtual ULong_t Hash() const {
    UInt_t cs = (static_cast<UInt_t>(fCrate)<<16) + static_cast<UInt_t>(fSlot);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
    return TString::Hash( &cs, sizeof(cs) );
#else
    return TMath::Hash( &cs, sizeof(cs) );
#endif
  }
  virtual Bool_t IsEqual( const TObject* obj ) const {
    const CSpair* m;
    if( !obj or !(m = dynamic_cast<const CSpair*>(obj)) ) return kFALSE;
    return ( fCrate == m->fCrate and fSlot == m->fSlot );
  }
  UShort_t  fCrate;
  UShort_t  fSlot;
};
class DAQmodule : public CSpair {
public:
  DAQmodule( UShort_t crate, UShort_t slot, UInt_t model, UInt_t nchan )
    : CSpair(crate, slot), fModel(model), fNchan(nchan) {}
  virtual ~DAQmodule() {}
  virtual void Copy( TObject& obj ) const {
    TObject::Copy(obj);
    DAQmodule* m = dynamic_cast<DAQmodule*>(&obj);
    if( !m ) return;
    m->fCrate = fCrate; m->fSlot = fSlot; m->fModel = fModel;
    m->fNchan = fNchan;
  }
  virtual void Print( Option_t* ) const {
    cout << "DAQmodule: "
	 << " crate = " << fCrate
	 << " slot = "  << fSlot
	 << " model = " << fModel
	 << " nchan = "   << fNchan
	 << endl;
  }
  UInt_t    fModel;
  UInt_t    fNchan;
};

///////////////////////////////////////////////////////////////////////////////

} // end namespace

namespace TreeSearch {

typedef vector<Plane*>::size_type vrsiz_t;
typedef vector<Plane*>::iterator  vriter_t;
typedef vector<Projection*>::size_type vpsiz_t;
typedef vector<Projection*>::iterator  vpiter_t;
typedef vector<vector<Int_t> >::iterator vviter_t;

// Global constant indicating invalid/uninitialized data
const Double_t kBig = 1e38;

#define ALL(c) (c).begin(), (c).end()

static void DoTrack( void* ptr );

//_____________________________________________________________________________
// Support classes for data-processing threads

struct TrackThread {
  Projection*  proj;    // Projection to be processed
  Int_t*       status;  // Return status (shared)
  UInt_t*      running; // Bitfield indicating threads to wait for (shared)
  TMutex*      start_m; // Mutex for start condition
  TMutex*      done_m;  // Mutex for done condition
  TCondition*  start;   // Start condition of tracking threads
  TCondition*  done;    // Condition indicating all tracking threads done
  TThread*     thread;  // The actual thread running with these arguments
  TrackThread()
    : proj(0), status(0), running(0), start(0), done(0), thread(0) {}
  ~TrackThread() { 
    if( thread ) {
      TThread::Delete(thread);
      delete thread;
    }
  }
};
  
//_____________________________________________________________________________
class ThreadCtrl {
  //TODO: better encapsulation
  friend THaAnalysisObject::EStatus Tracker::Init(const TDatime&);
  friend Int_t Tracker::CoarseTrack(TClonesArray&);
public:
  ThreadCtrl() : fTrackStatus(0), fTrackToDo(0),
		 fTrackStartM(new TMutex), fTrackDoneM(new TMutex),
		 fTrackStart(new TCondition(fTrackStartM)),
		 fTrackDone(new TCondition(fTrackDoneM)) {}
  ~ThreadCtrl() {
    // Terminate tracking threads
    if( !fTrack.empty() ) {
      fTrackStartM->Lock();
      fTrackToDo = 0;
      for( vector<TrackThread>::iterator it = fTrack.begin(); it !=
	     fTrack.end(); ++it ) {
	TThread* th = (*it).thread;
	if( th and th->GetState() == TThread::kRunningState )
	  SETBIT(fTrackToDo, (*it).proj->GetType());
      }
      // negative sign is termination flag
      SETBIT(fTrackToDo, 31);
      fTrackStart->Broadcast();
      fTrackDoneM->Lock();
      fTrackStartM->UnLock();
      while( true ) {
	Int_t ret = fTrackDone->Wait();
	if( ret == 0 and fTrackToDo == BIT(31) )
	  break;
      }
      fTrackDoneM->UnLock();
    }
    // At this point, all threads should have released their condition
    // waits and mutexes, so clean up is safe
    delete fTrackDone;
    delete fTrackStart;
    delete fTrackDoneM;
    delete fTrackStartM;
  }
  TThread* AddTrackThread( Projection* proj ) {
    assert( proj );
    fTrack.push_back( TrackThread() );
    TrackThread* t = &fTrack.back();
    t->proj    = proj;
    t->status  = &fTrackStatus;
    t->running = &fTrackToDo;
    t->start_m = fTrackStartM;
    t->done_m  = fTrackDoneM;
    t->start   = fTrackStart;
    t->done    = fTrackDone;
    string tn;
    tn = "trk_"; tn.append( proj->GetName() );
    t->thread  = new TThread( tn.c_str(), DoTrack, (void*)t );
    return t->thread;
  }
private:
  vector<TrackThread>  fTrack;        // Tracking thread descriptors
  Int_t                fTrackStatus;  // Common status variable
  UInt_t               fTrackToDo;    // Bitfield of projections to wait for
  TMutex*              fTrackStartM;  // Mutex for start condition
  TMutex*              fTrackDoneM;   // Mutex for done condition
  TCondition*          fTrackStart;   // Start condition
  TCondition*          fTrackDone;    // Finish condition
};

//_____________________________________________________________________________
Tracker::Tracker( const char* name, const char* desc, THaApparatus* app )
  : THaTrackingDetector(name,desc,app), fCrateMap(0), 
    fMaxThreads(1), fThreads(0),
    f3dMatchvalScalefact(1), f3dMatchCut(0),
    fMaxCorrMismatches(0), fMaxCorrNsigma(1.0),
    fMinNdof(1), fFailNhits(0), fFailNpat(0),
    fNcombos(0), fN3dFits(0), fEvNum(0),
    t_track(0), t_3dmatch(0), t_3dfit(0), t_coarse(0)
{ 
  // Constructor
}

//_____________________________________________________________________________
Tracker::~Tracker()
{
  // Destructor. Delete objects & subdetectors and unregister variables
  if (fIsSetup)
    RemoveVariables();
  
  delete fThreads;
  if( fMaxThreads > 1 )
    gSystem->Unload("libThread");

  DeleteContainer( fPlanes );
  DeleteContainer( fProj );
}

//_____________________________________________________________________________
Int_t Tracker::Begin( THaRunBase* run )
{
#ifdef TESTCODE
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->Begin(run);
#endif
  return 0;
}

//_____________________________________________________________________________
void Tracker::Clear( Option_t* opt )
{
  // Clear event-by-event data, including those of the planes and projections
  THaTrackingDetector::Clear(opt);
  
  // Clear the planes and projections, but only if we're not called from Init()
  if( !opt or *opt != 'I' ) {
    for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
      fPlanes[iplane]->Clear(opt);

    for( vpsiz_t k = 0; k < fProj.size(); ++k )
      fProj[k]->Clear(opt);
  }
  fFailNpat = fFailNhits = 0;
#ifdef TESTCODE
  size_t nbytes = (char*)&t_coarse - (char*)&fNcombos + sizeof(t_coarse);
  memset( &fNcombos, 0, nbytes );
#endif
}

//_____________________________________________________________________________
Int_t Tracker::Decode( const THaEvData& evdata )
{
  // Decode all planes and fill hitpatterns per projection
  
  //  static const char* const here = "Decode";

#ifdef TESTCODE
  // Save current event number for diagnostics
  fEvNum = evdata.GetEvNum();
#endif

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << "============== Event ";
#ifdef TESTCODE
    cout << "# " << fEvNum << " ";
#endif
    cout << "==============" << endl;
  }
#endif

  // Decode the planes, then fill the hitpatterns in the projections
  //TODO: multithread?
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    Projection* theProj = *it;
    Int_t nhits = theProj->Decode( evdata );
    // Sanity cut on overfull planes. nhits < 0 indicates overflow
    if( nhits < 0 ) {
      fFailNhits = 1;
      continue;
    }
    // No need to fill the hitpattern if no tracking requested
    if( TestBit(kDoCoarse) )
      theProj->FillHitpattern();
  }

  return 0;
}

//_____________________________________________________________________________
Int_t Tracker::End( THaRunBase* run )
{
#ifdef TESTCODE
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->End(run);
#endif
  return 0;
}

//_____________________________________________________________________________
class AddFitCoord
{
public:
  AddFitCoord() {}
  void operator() ( Road* rd, Road::Point* p, const vector<Double_t>& coef )
  {
    const Projection* proj = rd->GetProjection();
    Double_t slope2d  = rd->GetSlope();
    Double_t trkpos2d = rd->GetPos() + p->z * slope2d;
    Double_t cosa     = proj->GetCosAngle();
    Double_t sina     = proj->GetSinAngle();
    Double_t slope    = coef[1]*cosa + coef[3]*sina;
    Double_t pos      = coef[0]*cosa + coef[2]*sina + slope * p->z;
    p->hit->GetPlane()->AddFitCoord
      ( FitCoord( p->hit, rd, p->x, trkpos2d, slope2d, pos, slope ));
  }
};

//_____________________________________________________________________________
#ifdef VERBOSE
class PrintFitPoint
{
public:
  PrintFitPoint() {}
  void operator() ( Road*, Road::Point* p, const vector<Double_t>& )
  {
    cout << p->hit->GetPlane()->GetName() << " " 
	 << "z = " << p->z << " x = " << p->x << endl;
  }
};
#endif

//_____________________________________________________________________________
class FillFitMatrix
{
public:
  FillFitMatrix( TMatrixDSym& AtA, TVectorD& Aty ) 
    : fAtA(AtA), fAty(Aty), fNpoints(0) 
  {
    assert( fAtA.GetNrows() == 4 && fAty.GetNrows() == 4 );
  }

  void operator() ( Road* rd, Road::Point* p, const vector<Double_t>& )
  {
    const Projection* proj = rd->GetProjection();
    Double_t cosa = proj->GetCosAngle();
    Double_t sina = proj->GetSinAngle();
    Double_t Ai[4] = { cosa, cosa * p->z, sina, sina * p->z };
    Double_t s2 = 1.0/(p->res()*p->res());
    for( int j = 0; j<4; ++j ) {
      for( int k = j; k<4; ++k ) {
	fAtA(j,k) += Ai[j] * s2 * Ai[k];
      }
      fAty(j) += Ai[j] * s2 * p->x;
    }
    ++fNpoints;
  }
  Int_t GetNpoints() const { return fNpoints; }
private:
  TMatrixDSym& fAtA;
  TVectorD&    fAty;
  Int_t        fNpoints;
};

//_____________________________________________________________________________
class CalcChisquare
{
public:
  CalcChisquare() : fChi2(0) {}
  void operator() ( Road* rd, Road::Point* p, const vector<Double_t>& coef )
  {
    const Projection* proj = rd->GetProjection();
    Double_t cosa  = proj->GetCosAngle();
    Double_t sina  = proj->GetSinAngle();
    Double_t slope = coef[1]*cosa + coef[3]*sina;
    Double_t pos   = coef[0]*cosa + coef[2]*sina + slope*p->z;
    Double_t diff  = (pos - p->x) / p->res();
    fChi2 += diff*diff;
  }
  void     Clear() { fChi2 = 0; }
  Double_t GetResult() const { return fChi2; }
private:
  Double_t fChi2;
};

//_____________________________________________________________________________
template< typename Action >
Action Tracker::ForAllTrackPoints( const Rvec_t& roads, 
				   const vector<Double_t>& coef, Action action )
{
  // Apply action to all points from roads and the track given by coef

  for( Rvec_t::const_iterator it = roads.begin(); it != roads.end(); ++it ) {
    Road* rd = *it;
    const Road::Pvec_t& points = rd->GetPoints();
    for( Road::Pvec_t::const_iterator it2 = points.begin();
	 it2 != points.end(); ++it2 ) {
      Road::Point* p = *it2;
      action( rd, p, coef );
    }
  }
  return action;
}

//_____________________________________________________________________________
void Tracker::FitErrPrint( Int_t err ) const
{
  static const char* const here = "FitTrack";

#ifdef TESTCODE
  Error( Here(here), "Event %d: Failure fitting 3D track, err = %d. "
	 "Should never happen. Call expert.", fEvNum, err );
#else
  Error( Here(here), "Failure fitting 3D track, err = %d. "
	 "Should never happen. Call expert.", err );
#endif
}

//_____________________________________________________________________________
Int_t Tracker::FitTrack( const Rvec_t& roads, vector<Double_t>& coef,
			 Double_t& chi2, TMatrixDSym* coef_covar ) const
{
  // Linear minimization routine to fit a physical straight line track through
  // the hits registered in the different projection planes.
  //
  // This is a much streamlined version of ROOT's TLinearFitter that solves
  // the normal equations with weights, (At W A) a = (At W) y, using Cholesky
  // decomposition, TDecompChol. The model used is
  //   y_i = P_i * T_i
  //       = ( x + z_i * mx, y + z_i * my) * ( cos(a_i), sin(a_i) )
  // where
  //   y_i is the measured coordinate in the i-th plane at z_i
  //   P_i is the physical track intersection point with the z_i plane
  //   T_i is the axis unit vector of the i-th plane
  //   a_i is the angle of the coordinate axis of the i-th plane w.r.t. x
  //   x,mx,y,my are the track parameters to be fitted, origin x,y and
  //       slopes mx,my.
  //   
  // "roads" contains a set of Roads that successfully combine in 3-d, one
  // Road* per projection. Each road, in turn, contains a set of (y_i,z_i)
  // coordinates, at most one per plane of the projection type, that
  // give the best 2D track fit within the road.
  //
  // "coef" are the fitted track parameters, x, x'(=mx), y, y'(=my).
  // The reference system is the z=0 plane (usually the first chamber).
  //
  // "chi2" is the chi2 of the fit (not normalized).
  //
  // "coef_covar" is the covariance matrix of the fit results. Since it is
  // expensive to calculate, it is only filled if the argument is supplied.
  //
  // The return value is the number of degrees of freedom of the fit, i.e.
  // npoints-4 > 0, or negative if too few points or matrix inversion error

  // Fill the (At W A) matrix and (At W y) vector with the measured points
  TMatrixDSym AtA(4);
  TVectorD Aty(4);
  FillFitMatrix f = ForAllTrackPoints(roads, coef, FillFitMatrix(AtA, Aty));

  Int_t npoints = f.GetNpoints();
  assert( npoints > 4 );
  if( npoints <=4 ) return -1; // Meaningful fit not possible

  // FillFitMatrix only fills the upper triangle for efficiency
  //TODO: make a function of FillFitMatrix?
  for( int j = 0; j<4; ++j )
    for( int k = j+1; k<4; ++k )
      AtA(k,j) = AtA(j,k);

  // Invert the characteristic matrix and solve the normal equations.
  // As in ROOT's TLinearFitter, we use a Cholesky decomposition
  // to do this (since AtA is symmetric and positive-definite).
  // For more speed but less accuracy, one could use TMatrixDSymCramerInv.
  TDecompChol chol(AtA);
  Bool_t ok = chol.Solve(Aty);
  assert(ok);
  if( !ok ) return -2; //Urgh, decomposition failed. Should never happen

  // Copy results to output vector in order x, x', y, y'
  assert( Aty.GetNrows() == 4 );
  coef.assign( Aty.GetMatrixArray(), Aty.GetMatrixArray()+Aty.GetNrows() );

#ifdef VERBOSE
  if( fDebug > 2 ) {
    cout << "Points in 3D fit:" << endl;
    ForAllTrackPoints( roads, coef, PrintFitPoint() );
  }
#endif

  // Calculate chi2
  chi2 = ForAllTrackPoints(roads, coef, CalcChisquare()).GetResult();

  // Calculate covariance matrix of the parameters, if requested
  if( coef_covar ) {
    if( coef_covar->GetNrows() != 4 )
      coef_covar->ResizeTo(4,4);
    Bool_t ok = chol.Invert(*coef_covar);
    assert(ok);
    if( !ok ) return -3; // Urgh, inversion failed. Should never happen
  }

#ifdef VERBOSE
  if( fDebug > 1 ) 
    cout << "3D fit:  x/y = " << coef[0] << "/" << coef[2] << " "
	 << "mx/my = " << coef[1] << "/" << coef[3] << " "
	 << "ndof = " << npoints-4 << " rchi2 = " << chi2/(double)(npoints-4)
	 << endl;
#endif
  
  return npoints-4;
}

//_____________________________________________________________________________
void Tracker::FindNearestHits( Plane* rp, const THaTrack* track, 
			       const Rvec_t& roads ) const
{
  // For the given plane, find the hit nearest to the given track
  // and register it in the plane's fit coordinates array.
  // The given roads are the ones generating the given track.
  // This routine is used for efficiency and alignment studies and testing.

  assert( !roads.empty() );

  //TODO: evaluate how much is common between MWDC/GEM

  // // Search for the hit with position closest to the track crossing
  // // position in this plane. The hits are sorted by position, so
  // // the search can be made fast.
  // Double_t z     = rp->GetZ();
  // Double_t cosa  = rp->GetProjection()->GetCosAngle();
  // Double_t sina  = rp->GetProjection()->GetSinAngle();
  // Double_t slope = track->GetDTheta()*cosa + track->GetDPhi()*sina;
  // Double_t x     = track->GetDX()*cosa + track->GetDY()*sina + slope*z;
  // Double_t pmin  = kBig;
  // Hit* hmin = 0;
  // // Binary search the hits array for the track crossing position x, similar
  // // to std::lower_bound(). This finds the first hit with strip position >= x.
  // const TSeqCollection* hits = rp->GetHits();
  // Int_t first = 0, last = hits->GetSize();
  // // GetSize() is incorrect for TObjArrays and TClonesArrays
  // const TObjArray* arr = dynamic_cast<const TObjArray*>(hits);
  // if( arr )
  //   last = arr->GetLast()+1;
  // Int_t len = last - first;
  // Int_t half, middle;
  // while( len > 0 ) {
  //   half = len >> 1;
  //   middle = first + half;
  //   if( static_cast<Hit*>(hits->At(middle))->GetPos() < x ) {
  //     first = middle + 1;
  //     len -= half + 1;
  //   } else
  //     len = half;
  // }
  // // Decide whether the hit >= x or the first one < x are closest.
  // if( last > 0 ) {
  //   assert( first <= last );
  //   Hit *hnext = 0, *hprev = 0;
  //   if( first < last ) {
  //     hnext = static_cast<Hit*>(hits->At(first));
  //     assert( hnext->GetPos() >= x );
  //   }
  //   if( first > 0 ) {
  //     hprev = static_cast<Hit*>(hits->At(first-1));
  //     assert( hprev->GetPos() < x );
  //     if( hnext ) {
  // 	assert( hprev->GetPos() < hnext->GetPos() );
  // 	if( x - hprev->GetPos() < hnext->GetPos() - x )
  // 	  hnext = 0;
  // 	else
  // 	  hprev = 0;
  //     }
  //   }
  //   assert( (hprev != 0) xor (hnext != 0) );
  //   if( hnext )
  //     hmin = hnext;
  //   else
  //     hmin = hprev;
  //   pmin = hmin->GetPos();
  // }
  // // The road vector does not necessarily contain all projections, so
  // // search for the road of the type of this readout plane, taking advantage
  // // of the fact that the road vector is sorted by type
  // Road* rd = 0;
  // Int_t k = min( roads.size(), (Rvec_t::size_type)rp->GetType() );
  // do {
  //   if( roads[k]->GetProjection()->GetType() > rp->GetType() )
  //     --k;
  //   else {
  //     if( roads[k]->GetProjection() == rp->GetProjection() )
  // 	rd = roads[k];
  //     break;
  //   }
  // } while( k>=0 );
  // Double_t slope2d = rd ? rd->GetSlope() : kBig;
  // Double_t pos2d   = rd ? rd->GetPos() + z * slope2d  : kBig;

  // // Finally, record the hit info in the readout plane
  // rp->AddFitCoord( FitCoord(hmin, rd, pmin, pos2d, slope2d, x, slope) );
}

//_____________________________________________________________________________
THaTrack* Tracker::NewTrack( TClonesArray& tracks, const FitRes_t& fit_par )
{
  // Make new track with given parameters. Used by CoarseTrack to generate
  // all GEM tracks

  Double_t d_x  = fit_par.coef[0];
  Double_t d_xp = fit_par.coef[1];
  Double_t d_y  = fit_par.coef[2];
  Double_t d_yp = fit_par.coef[3];

  // Correct the track coordinates for the position offset of the detector
  // system. The track's coordinates are given in reference to the following
  // coordinate systems:
  //
  // "fp" coordinates ("tr.x" etc.): Tracker system
  // "detector" coordinates ("tr.d_x" etc.): Plane system
  //
  // NB: The track origin is always given at z=0 of the Tracker system.
  //
  // Currently assumes that the Plane and Tracker systems have parallel
  // axes, i.e. there is no rotation.
  Double_t xp = d_xp;
  Double_t yp = d_yp;
  Double_t x  = d_x + fOrigin.X() - d_xp*fOrigin.Z();
  Double_t y  = d_y + fOrigin.Y() - d_yp*fOrigin.Z();

  THaTrack* newTrack = AddTrack( tracks, x, y, xp, yp );
  assert( newTrack );
  newTrack->SetD( d_x, d_y, d_xp, d_yp );
  newTrack->SetChi2( fit_par.chi2, fit_par.ndof );
  //TODO: make a TrackID?

  ForAllTrackPoints( *fit_par.roads, fit_par.coef, AddFitCoord() );
  for( Rpvec_t::const_iterator it = fCalibPlanes.begin(); it !=
	 fCalibPlanes.end(); ++it ) {
    FindNearestHits( (*it), newTrack, *fit_par.roads );
  }

  return newTrack;
}

//_____________________________________________________________________________
class CheckTypes : public unary_function<Road*,void>
{
  // Function object for testing the projection occupancy of roads.
  // Applied to a road, it saves the road's plane type (u,v,x..).
  // Testing the object returns true if all types specified in the 'req'
  // argument to the constructor have been seen in all the roads tested.
public:
  CheckTypes( UInt_t req ) : fReq(req), fActive(0) {}
  void operator() ( const Road* rd )
  { fActive |= 1U << rd->GetProjection()->GetType(); }
  operator bool()       const { return (fActive & fReq) == fReq; }
  bool     operator!()  const { return (fActive & fReq) != fReq; }
  UInt_t   GetTypes()   const { return fActive; }
private:
  UInt_t fReq;
  UInt_t fActive;
};

//_____________________________________________________________________________
class AnySharedHits : public unary_function<Road*,bool>
{
public:
  AnySharedHits() {}
  bool operator() ( const Road* rd )
  {
    const Hset_t& test_hits = rd->GetHits();
//     cout << "--- testing: " << endl;
//     PrintHits(test_hits);
    for( Hset_t::const_iterator it = test_hits.begin(); it != test_hits.end();
	 ++it ) {
      if( fHits.find(*it) != fHits.end() )
	return true;
    }
    return false;
  }
  const AnySharedHits& use( const Rset_t& tuple )
  {
    fHits.clear();
    for( Rset_t::const_iterator it = tuple.begin(); it != tuple.end(); ++it ) {
      const Road* rd = *it;
      fHits.insert( ALL(rd->GetHits()) );
    }
//     PrintHits(fHits);
    return *this;
  }
private:
  Hset_t fHits;
};

//_____________________________________________________________________________
class Tracker::TrackFitWeight
{
  // Object for sorting fits. The smaller a track's TrackFitWeight, the
  // "better" a track is considered to be.
public:
  TrackFitWeight( const Tracker::FitRes_t& fit_par ) :
    fNdof(fit_par.ndof), fChi2(fit_par.chi2) { assert(fNdof); }

  bool operator<( const TrackFitWeight& rhs ) const
  {
    // The "best" tracks have the largest number of hits and the smallest chi2
    //TODO: devalue tracks with very large chi2?
    if( fNdof > rhs.fNdof ) return true;
    if( fNdof < rhs.fNdof ) return false;
    return ( fChi2 < rhs.fChi2 );
  }
  operator double() const { return fChi2/(double)fNdof; }
private:
  UInt_t   fNdof;
  Double_t fChi2;
};

//_____________________________________________________________________________
template< typename Container, typename Weight, typename TestF,
	  typename QuitF > Double_t
OptimalN( const Container& choices, const multimap<Weight,Container>& weights,
	  vector<Container>& picks, TestF testf, QuitF quitf )
{
  // This is the second-level de-ghosting algorithm, operating on fitted
  // 3D tracks.  It selects the best set of tracks if multiple roads
  // combinations are present in one or more projection.
  //
  // Requirements:
  //  "Container": STL container (set or sorted vector) supporting begin(),
  //               insert(), end(), swap(), and containing "Element"s sorted
  //               by Element's operator<.  Here, Element = Road*
  //  "Weight":    sorting functor, supporting operator< and operator double()
  //  "TestF":     Functor for eliminating additional leftover elements,
  //               must support:
  //                  operator() (const Element&): test this element,
  //                                         eliminate if true
  //                  use(const Container&): test against elements in this 
  //                                         container
  //               (use() is in lieu of a constructor because the object
  //                is constructed by the caller (possibly with parameters)).
  //  "QuitF":     Functor for terminating search early. Applied to all 
  //               leftover elements. Must support
  //                  operator() (const Element&)
  //                  operator bool(): End search if false after running on
  //                                   all elements

  //TODO: how efficient? O(N^2)? (N iterations on include)

  picks.clear();
  Container choices_left(choices);

  //TODO: try minimizing chi2 sum
  double wsum = 0;
  typename multimap<Weight,Container>::const_iterator it = weights.begin();
  for( ; it != weights.end(); ++it ) {
    const Container& tuple = (*it).second;
    if( includes(ALL(choices_left), ALL(tuple)) ) {
      // Pick next set of still-available elements in order of ascending weight
      picks.push_back( tuple );
      // Sum of weights of chosen sets, converted to double
      wsum += (*it).first;
      // Remove chosen elements from the remaining choices
      Container choices_less_tuple;
      set_difference( ALL(choices_left), ALL(tuple),
		      inserter(choices_less_tuple, choices_less_tuple.end()) );
      // Test each remaining element against the chosen set of elements
      // and remove it if it compares true (=remove elems with common traits)
      Container new_choices_left;
      remove_copy_if( ALL(choices_less_tuple),
		      inserter(new_choices_left, new_choices_left.end()),
		      testf.use(tuple) );
      // Quit if the quit function, applied to each element still remaining
      // now, is no longer true
      if( !for_each(ALL(new_choices_left), quitf) )
	break;
      // Update the set of remaining elements
      choices_left.swap( new_choices_left );
    }
  }

  return wsum;
}

//_____________________________________________________________________________
inline
void Tracker::Add3dMatch( const Rvec_t& selected, Double_t matchval,
			  list< pair<Double_t,Rvec_t> >& combos_found,
			  Rset_t& unique_found ) const
{
  // Save road combination with good matchvalue

  combos_found.push_back( make_pair(matchval,selected) );

  // Since not all of the roads in 'roads' may make a match,
  // keep track of unique roads found for each projection type

  unique_found.insert( ALL(selected) );

#ifdef VERBOSE
  if( fDebug > 3 )
    cout << "ACCEPTED" << endl;
#endif
}

//_____________________________________________________________________________
UInt_t Tracker::MatchRoads( vector<Rvec_t>& roads,
			    list< pair<Double_t,Rvec_t> >& combos_found,
			    Rset_t& unique_found )
{
  // Match roads from different projections. This works differently
  // depending on the number of projections and the geometry:
  //
  // 2 projections:
  //  Match amplitudes of hits in shared readout planes (if available)
  //
  // 3 or more projections:
  //  Match via geometric closeness of intersection points in both the front
  //  and back planes of the detector

  // ===> TODO: extract common code, split off amplitude matching into GEM

  // The number of projections that we work with (must be >= 2)
//   vector<Rvec_t>::size_type nproj = roads.size();
//   assert( nproj >= 2 );

//   combos_found.clear();
//   unique_found.clear();

//   // Number of all possible combinations of the input roads
//   // May overflow for extremely busy events
//   UInt_t ncombos;
//   bool inrange = true;
//   try {
//     ncombos = accumulate( ALL(roads), (UInt_t)1, SizeMul<Rvec_t>() );
//   }
//   catch( overflow_error ) {
//     ncombos = 0;
//     inrange = false;
//   }
//   bool fast_3d = TestBit(k3dFastMatch);
//   bool correlate_amplitudes = TestBit(k3dCorrAmpl);

//   // Limitation: amplitude correlation algo only implemented for 2 projections
//   assert( nproj == 2 or not correlate_amplitudes );

//   // Fast 3D matching and amplitude correlation are mutually exclusive
//   assert( not (fast_3d and correlate_amplitudes) );

// #ifdef VERBOSE
//   if( fDebug > 0 ) {
//     if( inrange )
//       cout << "Matching ";
//     else
//       cout << "Too many combinations trying to match ";
//     for( vector<Rvec_t>::size_type i = 0; i < nproj; ++i ) {
//       cout << roads[i].size();
//       if( i+1 < nproj ) cout << "x";
//     }
//     cout << " track projections in 3D";
//     if( inrange ) {
//       cout << " (";
//       if( fast_3d )
// 	cout << "fast";
//       else if( correlate_amplitudes )
// 	cout << "amplitude correlation";
//       else
// 	cout << "generic";
//       cout << " algo, " << ncombos << " combination";
//       if( ncombos != 1 ) cout << "s";
//       cout << "):";
//     } else {
//       cout << ". Giving up.";
//     }
//     cout << endl;
//   }
// #endif
// #ifdef TESTCODE
//   fNcombos = ncombos;
// #endif

//   if( ncombos == 0 or (nproj == 2 and not correlate_amplitudes) )
//     return 0;

//   // Vector holding a combination of roads to test. One road from each 
//   // projection 
//   Rvec_t selected;
//   selected.reserve(nproj);

  UInt_t nfound = 0;
//   Double_t zback = fPlanes.back()->GetZ();

//   if( correlate_amplitudes ) {
//     assert( nproj == 2 );
//     assert( not (roads[0].empty() or roads[1].empty()) );

//     const Rpvec_t& xplanes =
//       roads[0].front()->GetProjection()->GetListOfPlanes();
//     const Rpvec_t& yplanes =
//       roads[1].front()->GetProjection()->GetListOfPlanes();

//     assert( not xplanes.empty() and (xplanes.size() == yplanes.size()) );
//     assert( xplanes.front()->GetType() != yplanes.front()->GetType() );

//     selected.assign( nproj, 0 );

//     UInt_t nplanes = xplanes.size();
//     TBits xybits(nplanes), ybits(nplanes); 
//     // Look at all possible combinations of x-roads and y-roads.
//     // Try to match them via the ADC amplitudes of the hits in the shared
//     // readout planes. Amplitudes should correlate well in the absence
//     // of pileup. 
//     for( Rvec_t::iterator itx = roads[0].begin(); itx != roads[0].end();
// 	 ++itx ) {
//       Road* xroad = *itx;
//       const Road::Pvec_t& xpoints = xroad->GetPoints();
//       UInt_t xpat = xroad->GetPlanePattern();
//       for( Rvec_t::iterator ity = roads[1].begin(); ity != roads[1].end();
// 	   ++ity ) {
// 	Road* yroad = *ity;

// 	// xpat and ypat are bitpatterns of the plane numbers that have hits.
// 	// The AND of these patters is the pattern of planes where both read-
// 	// out directions have hits (matching or not). Check here if there
// 	// are enough such common active planes, else all following work can
// 	// be skipped.
// 	UInt_t ypat = yroad->GetPlanePattern();
// 	xybits.Set( nplanes, &xpat );
// 	ybits.Set( nplanes, &ypat );
// 	xybits &= ybits;
// 	assert( xybits.CountBits() <= nplanes );
// 	if( xybits.CountBits() + fMaxCorrMismatches < nplanes )
// 	  continue;

// 	// For all points (=hits that yield the best fit) of this xroad, 
// 	// get the yroad point in the same plane (if available), then check
// 	// if their ADC amplitudes approximately match (within a hard cut)
// 	const Road::Pvec_t& ypoints = yroad->GetPoints();
// 	Road::Pvec_t::const_iterator ityp = ypoints.begin();
// 	UInt_t nmatches = 0;
// 	Double_t matchval = 0.0;
// 	for( Road::Pvec_t::const_iterator itxp = xpoints.begin();
// 	     itxp != xpoints.end(); ++itxp ) {
// 	  const Road::Point* xp = *itxp;
// 	  assert( xp and xp->hit );
// 	  const Plane* xplane = xp->hit->GetPlane();
// 	  assert( xplane );
// 	  UInt_t xnum = xplane->GetPlaneNum();
// 	  assert( xnum < xplanes.size() );
// 	  // No hit in the other readout direction of this plane?
// 	  if( not ybits.TestBitNumber(xnum) )
// 	    continue;
// 	  // Move y-iterator forward until it gets to the same plane
// 	  while( ityp != ypoints.end() and
// 		 (*ityp)->hit->GetPlaneNum() != xnum ) {
// 	    ++ityp;
// 	  }
// 	  // Must find a corresponding hit here, else bug in ybits test
// 	  assert( ityp != ypoints.end() );
// 	  if( ityp == ypoints.end() ) break;
// 	  const Road::Point* yp = *ityp;
// 	  assert( yp and yp->hit );
// 	  const Plane* yplane = yp->hit->GetPlane();
// 	  assert( yplane );
// 	  assert( yplane->GetPartner() == xplane );

// 	  // Get ratio of hit amplitudes
// 	  Double_t xampl = xp->hit->GetADCsum();
// 	  Double_t yampl = yp->hit->GetADCsum();
// 	  assert( xampl > 0.0 and yampl > 0.0 ); // ensured in Decoder
// 	  if( xampl < 1.0 or yampl < 1.0 )
// 	    continue;
// 	  Double_t ratio = xampl/yampl;
// 	  Double_t asym = (xampl - yampl)/(xampl + yampl);
// 	  // Compute the cutoff. In general, the sigma of the distribution can
// 	  // be amplitude-dependent, so we get it via a calibration function
// 	  // that is a property of each readout plane
// // 	  Double_t xsigma = xplane->GetAmplSigma( xampl );
// // 	  Double_t ysigma = yplane->GetAmplSigma( yampl );
// 	  // Apply overall scale factor ("number of sigmas") from the database
// // 	  Double_t cutoff = fMaxCorrNsigma / yampl *
// // 	    TMath::Sqrt( xsigma*xsigma + ysigma*ysigma*ratio*ratio );
// 	  Double_t cutoff = fMaxCorrNsigma;
// #ifdef VERBOSE
// 	  if( fDebug > 3 ) {
// 	    cout << xplane->GetName() << yplane->GetName()
// 		 << " ampl = (" << xampl << ", " << yampl << ")"
// 		 << ",\tratio = " << ratio
// 		 << ", asym = " << asym << endl;
// 	  }
// #endif
// 	  // Count readout planes whose x/y-hit amplitudes match
// 	  if( TMath::Abs(asym) < cutoff ) {
// 	    matchval += TMath::Abs(asym);  // not really used later (yet)
// 	    ++nmatches;
// 	  }
// 	} //xpoints
	  
// #ifdef VERBOSE
// 	if( fDebug > 3 ) {
// 	  cout <<   "nmatches = " << nmatches 
// 	       << ", matchval = " << matchval << "   ";
// 	}
// #endif
// 	// If enough planes have correlated x/y hit amplitudes, save this
// 	// x/y road pair as a possible 3D track candidate. The 3D track fit
// 	// will separate the wheat from the chaff later
// 	if( nmatches + fMaxCorrMismatches >= nplanes ) {
// 	  selected[0] = xroad;
// 	  selected[1] = yroad;
// 	  Add3dMatch( selected, matchval, combos_found, unique_found );
// 	  ++nfound;
//         }
// #ifdef VERBOSE
// 	else if( fDebug > 3 ) { cout << endl; }
// #endif
//       } //yroads
//     }   //xroads

//   } else if( fast_3d ) {
//     // special case n==3 and symmetric angles of planes 0 and 1:
//     //  - intersect first two proj in front and back
//     //  - calc perp distances to 3rd, add in quadrature -> matchval
//     assert( nproj == 3 && fProj.size() == 3 );

//     // Fetch coefficients for coordinate transformations
//     vpiter_t ip = fProj.begin();
//     assert( (*ip)->GetType() == kUPlane );
//     Double_t su = (*ip)->GetSinAngle();
//     Double_t cu = (*ip)->GetCosAngle();
//     ++ip;
//     assert( (*ip)->GetType() == kVPlane );
//     Double_t sv = (*ip)->GetSinAngle();
//     Double_t cv = (*ip)->GetCosAngle();
//     Double_t inv_denom = 1.0/(sv*cu-su*cv);
//     // Only the scaled coefficients are needed (cf. Road::Intersect)
//     su *= inv_denom; cu *= inv_denom; sv *= inv_denom; cv *= inv_denom;
//     // Components of the 3rd projection's axis
//     ++ip;
//     assert( (*ip)->GetType() == kXPlane or (*ip)->GetType() == kYPlane );
//     Double_t xax_x = ((*ip)->GetAxis()).X();
//     Double_t xax_y = ((*ip)->GetAxis()).Y();

//     // The selected roads from each of the three projections
//     Road* tuple[3];
//     // Number of roads in u/v projections
//     UInt_t nrd[2] = { roads[0].size(), roads[1].size() };
//     // Indices of currently selected u/v pair
//     UInt_t ird[2];
//     Plane *front_plane = fPlanes.front(), *back_plane = fPlanes.back();

//     // For fast access to the relevant position range, sort the 3rd projection
//     // by ascending front position
//     sort( ALL(roads[2]), Road::PosIsLess() );
//     Road::PosIsNear pos_near( TMath::Sqrt(f3dMatchCut) );

//     Double_t matchval = 0.0;
//     // Time-critical loop, may run O(1e5) times per event with noisy input
//     ird[0] = 0;
//     while( ird[0] != nrd[0] ) {
//       tuple[0] = roads[0][ird[0]];
//       Double_t uf = tuple[0]->GetPos();
//       Double_t ub = tuple[0]->GetPos(zback);
//       Double_t usf = uf * sv;
//       Double_t ucf = uf * cv;
//       Double_t usb = ub * sv;
//       Double_t ucb = ub * cv;
//       ird[1] = 0;
//       while( ird[1] != nrd[1] ) {
// 	tuple[1] = roads[1][ird[1]];
// 	Double_t v = tuple[1]->GetPos();
// 	Double_t xf = usf - v * su;
// 	Double_t yf = v * cu - ucf;
// 	if( front_plane->Contains(xf,yf) ) {
// 	  v = tuple[1]->GetPos(zback);
// 	  Double_t xb = usb - v * su;
// 	  Double_t yb = v * cu - ucb;
// 	  if( back_plane->Contains(xb,yb) ) {
// 	    Double_t pf = xf*xax_x + yf*xax_y; // front x from u/v
// 	    Double_t pb = xb*xax_x + yb*xax_y; // back x from u/v
// 	    // Find range of roads in 3rd projection near this front x
// 	    pair<Rvec_t::iterator,Rvec_t::iterator> range = 
// 	      equal_range( ALL(roads[2]), pf, pos_near );
// 	    // Test the candidate x-roads for complete matches
// 	    for( Rvec_t::iterator it=range.first; it != range.second; ++it ) {
// 	      tuple[2] = *it;
// 	      Double_t d1 = tuple[2]->GetPos()      - pf;
// 	      Double_t d2 = tuple[2]->GetPos(zback) - pb;
// 	      //TODO; weigh with uncertainties?
// 	      matchval = d1*d1 + d2*d2;
// #ifdef VERBOSE
// 	      if( fDebug > 3 ) {
// 		if( matchval < f3dMatchCut || fDebug > 4 ) {
// 		  cout << tuple[0]->GetProjection()->GetName()
// 		       << tuple[1]->GetProjection()->GetName()
// 		       << " front = " << "(" << xf << "," << yf << ")" << endl;
// 		  cout << "front x  = (" << tuple[2]->GetPos() * xax_x << ","
// 		       << tuple[2]->GetPos() * xax_y << ")" << endl;
// 		  cout << "front dist = " << d1 << endl;
// 		  cout << tuple[0]->GetProjection()->GetName()
// 		       << tuple[1]->GetProjection()->GetName()
// 		       << " back = " << "(" << xb << "," << yb << ")" << endl;
// 		  cout << "back x =  (" 
// 		       << tuple[2]->GetPos(zback) * xax_x << ","
// 		       << tuple[2]->GetPos(zback) * xax_y << ")" << endl;
// 		  cout << "back dist = " << d2 << endl;
// 		  cout << "matchval = " << matchval*f3dMatchvalScalefact
// 		       << endl;
// 		}
// 	      }
// #endif
// 	      // Check if match, if so then add
// 	      if( matchval < f3dMatchCut ) {
// 		++nfound;
// 		selected.assign( tuple, tuple+3 );
// 		Add3dMatch( selected, matchval, combos_found, unique_found );
// 	      }
// 	    }
// 	  }
// 	}
// 	++ird[1];
//       }
//       ++ird[0];
//     }

//   } else {
//     // general algorithm:
//     //  - find all front and back intersections [ n(n-1)/2 each ] 
//     //  - compute weighted center of gravity of intersection points
//     //  - sum dist^2 of points to center of gravity -> matchval

//     assert( nproj >= 3 );

//     vector<TVector2> fxpts, bxpts;
//     fxpts.reserve( nproj*(nproj-1)/2 );
//     bxpts.reserve( nproj*(nproj-1)/2 );

//     for( UInt_t i = 0; i < ncombos; ++i ) {
//       Double_t matchval = 0.0;

//       NthCombination( i, roads, selected );
//       assert( selected.size() == nproj );

//       fxpts.clear();
//       bxpts.clear();
//       TVector2 fctr, bctr;
//       for( Rvec_t::iterator it1 = selected.begin(); it1 != selected.end();
// 	   ++it1 ) {
// 	for( Rvec_t::iterator it2 = it1+1; it2 != selected.end(); ++it2 ) {
// 	  //TODO: weigh with uncertainties of coordinates?
// 	  fxpts.push_back( (*it1)->Intersect(*it2, 0.0) );
// 	  bxpts.push_back( (*it1)->Intersect(*it2, zback) );
// 	  fctr += fxpts.back();
// 	  bctr += bxpts.back();
// #ifdef VERBOSE
// 	  if( fDebug > 3 ) {
// 	    cout << (*it1)->GetProjection()->GetName()
// 		 << (*it2)->GetProjection()->GetName()
// 		 << " front(" << fxpts.size() << ") = ";
// 	    fxpts.back().Print();
// 	    cout << (*it1)->GetProjection()->GetName()
// 		 << (*it2)->GetProjection()->GetName()
// 		 << " back (" << bxpts.size() << ") = ";
// 	    bxpts.back().Print();
// 	  }
// #endif
// 	}
//       }
//       assert( fxpts.size() <= nproj*(nproj-1)/2 );
//       assert( bxpts.size() == fxpts.size() );
//       fctr /= static_cast<Double_t>( fxpts.size() );
//       if( !fPlanes.front()->Contains(fctr) )
// 	continue;
//       bctr /= static_cast<Double_t>( fxpts.size() );
//       if( !fPlanes.back()->Contains(bctr) )
// 	continue;
//       for( vector<TVector2>::size_type k = 0; k < fxpts.size(); ++k ) {
// 	matchval += (fxpts[k]-fctr).Mod2() + (bxpts[k]-bctr).Mod2();
//       }
// #ifdef VERBOSE
//       if( fDebug > 3 ) {
// 	cout << "fctr = "; fctr.Print();
// 	cout << "bctr = "; bctr.Print();
// 	cout << "matchval = " << matchval << endl;
//       }
// #endif
//       // We could just connect fctr and bctr here to get an approximate
//       // 3D track. But the linear minimization below is the right way
//       // to do this.

//       if( matchval < f3dMatchCut ) {
// 	++nfound;
// 	Add3dMatch( selected, matchval, combos_found, unique_found );
//       }
//     } //for(ncombos)

//   } //matching methods

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << nfound << " match";
    if( nfound != 1 )
      cout << "es";
    cout << " found" << endl;
  }
#endif
  return nfound;
}

//_____________________________________________________________________________
inline
Bool_t Tracker::PassTrackCuts( const FitRes_t& fit_par ) const
{
  // Test results of 3D track fit (in fit_par) against hard cuts from
  // database

  if( fit_par.ndof < fMinNdof )
    return false;

  pdbl_t chi2_interval = GetChisqLimits(fit_par.ndof);
  if( fit_par.chi2 < chi2_interval.first or
      fit_par.chi2 > chi2_interval.second )
    return false;

  return true;
}

//_____________________________________________________________________________
static void DoTrack( void* ptr )
{
  TrackThread* arg = reinterpret_cast<TrackThread*>(ptr);

  //  TThread::SetCancelOn();
  bool terminate = false;
  arg->start_m->Lock();
  while( !terminate ) {

    // Wait for start condition
    while( true ) {
      Int_t ret = arg->start->Wait();
      if( ret == 0 and TESTBIT(*arg->running, arg->proj->GetType()) )
	break;
    }
    terminate = TESTBIT(*arg->running, 31);
    arg->start_m->UnLock();

    Int_t nrd = 0;
    if( !terminate )
      // Process this event
      nrd = arg->proj->Track();

    // Set error flag, if necessary, and decrement run counter
    arg->done_m->Lock();
    if( nrd < 0 ) {
      *arg->status = 1;
    }
    // Clear the bit for this projection in the status bitfield
    assert( TESTBIT(*arg->running, arg->proj->GetType()) );
    CLRBIT( *arg->running, arg->proj->GetType() );
    // If all bits are zero, all threads have finished processing
    if( (*arg->running & ~BIT(31)) == 0 )
      arg->done->Signal();

    if( !terminate )
      // Ensure that we enter Wait() before the main thread can send the
      // next Broadcast()
      arg->start_m->Lock();  

    arg->done_m->UnLock();
  }
}

//_____________________________________________________________________________
Int_t Tracker::CoarseTrack( TClonesArray& tracks )
{
  // Find tracks from the hitpatterns, using the coarse hit drift times
  // uncorrected for track slope, timing offset, fringe field effects etc.

  //  static const char* const here = "CoarseTrack";

  if( fFailNhits )
    return -1;

  if( !TestBit(kDoCoarse) )
    return 0;

#ifdef TESTCODE
  TStopwatch timer, timer_tot;
#endif

  Int_t err = 0;

  bool do_thread = (fMaxThreads > 1);
  if( do_thread ) {
    // Wake the tracking threads
    fThreads->fTrackDoneM->Lock();
    fThreads->fTrackStatus = 0;
    EProjType type_beg = fProj.front()->GetType();
    EProjType type_end = fProj.back()->GetType();
    UInt_t all_todo = (BIT(type_end)-1) & ~(BIT(type_beg)-1);
    UInt_t max_todo = BIT(fMaxThreads)-1;
    Int_t start = type_beg;
    while( start < type_end ) {
      UInt_t now_todo = (max_todo << start) & all_todo;
      fThreads->fTrackStartM->Lock();
      fThreads->fTrackToDo = now_todo;
      fThreads->fTrackStart->Broadcast();
      fThreads->fTrackStartM->UnLock();
      // Wait for end of processing
      while( true ) {
	Int_t ret = fThreads->fTrackDone->Wait();
	if( ret == 0 and fThreads->fTrackToDo == 0 )
	  break;
      }
      start += fMaxThreads;
    }
    // Retrieve status
    err = fThreads->fTrackStatus;
    fThreads->fTrackDoneM->UnLock();
  } else {
    // Single-threaded execution: Track() each projection in turn
    for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
      Int_t ret = (*it)->Track();
      if( ret < 0 )
	err = 1;
    }
  }

  // Abort on error (e.g. too many patterns)
  if( err != 0 ) {
    fFailNpat = 1;
    return -1;
  }
  // Copy pointers to roads from each projection into local 2D vector
  // (projections, roads).
  // NB: subsequent code assumes that the first index runs in order
  // of ascending projection type, so this can't easily be done in the 
  // tracking threads.
  vector<Rvec_t>::size_type nproj = 0;
  vector<Rvec_t> roads;
  roads.reserve( fProj.size() );
  UInt_t found_types = 0;
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    Projection* proj = *it;
    EProjType type = proj->GetType();

    Int_t nrd = proj->GetNgoodRoads();
    if( nrd > 0 ) {
      // Count number of projections with at least one road
      ++nproj;
      SETBIT( found_types, type );
      // Add pointers to the good roads to the appropriate vector
      roads.push_back( Rvec_t() );
      roads.back().reserve(nrd);
      for( UInt_t i = 0; i < proj->GetNroads(); ++i ) {
	Road* rd = proj->GetRoad(i);
	assert(rd);
	if( rd->IsGood() )
	  roads.back().push_back(rd);
      }
    }
  }
  assert( roads.size() == nproj );
#ifdef TESTCODE
  t_track = 1e6*timer.RealTime();
  timer.Start();
#endif

  // Combine track projections to 3D tracks
  //===> TODO: parameterize minimum number of projections required
  if( nproj >= 2 ) {
    // Vector holding the results (vectors of roads with good matchval)
    list< pair<Double_t,Rvec_t> > road_combos;
    // Set of the unique roads occurring in the road_combos elements
    Rset_t unique_found;

    // Find matching combinations of roads
    UInt_t nfits = MatchRoads( roads, road_combos, unique_found );

#ifdef TESTCODE
    t_3dmatch = 1e6*timer.RealTime();
    timer.Start();
#endif
    // Fit each set of matched roads using linear least squares, yielding 
    // the 3D track parameters, x, x'(=mx), y, y'(=my)
    FitRes_t fit_par;
    fit_par.coef.reserve(4);
    if( nfits == 1 ) {
      // If there is only one combo (typical case), life is simple:
      fit_par.matchval    = road_combos.front().first;
      Rvec_t& these_roads = road_combos.front().second;
      fit_par.roads       = &these_roads;
      fit_par.ndof = FitTrack( these_roads, fit_par.coef, fit_par.chi2 );
      if( fit_par.ndof > 0 ) {
	if( PassTrackCuts(fit_par) )
	  NewTrack( tracks, fit_par );
      }
      else
	FitErrPrint( fit_par.ndof );
    }
    else if( nfits > 0 ) {
      // For multiple road combinations, find the set of tracks with the
      // lowest chi2s that uses each road at most once
      roads.clear();
      typedef map<Rset_t, FitRes_t> FitResMap_t;
      FitResMap_t fit_results;
      multimap< TrackFitWeight, Rset_t > fit_chi2;
      // Fit all combinations and sort the results by ascending chi2
      for( list< pair<Double_t,Rvec_t> >::iterator it = road_combos.begin();
	   it != road_combos.end(); ++it ) {
	fit_par.matchval    = (*it).first;
	Rvec_t& these_roads = (*it).second;
	fit_par.roads       = &these_roads;
	fit_par.ndof = FitTrack( these_roads, fit_par.coef, fit_par.chi2 );
	if( fit_par.ndof > 0 ) {
	  if( PassTrackCuts(fit_par) ) {
	    Rset_t road_tuple( ALL(these_roads) );
	    pair<FitResMap_t::iterator,bool> ins =
	      fit_results.insert( make_pair(road_tuple, fit_par) );
	    assert( ins.second );
	    fit_chi2.insert( make_pair( TrackFitWeight(fit_par), road_tuple) );
	  }
	} else
	  FitErrPrint( fit_par.ndof );
      }

#ifdef VERBOSE
      if( fDebug > 2 ) {
	cout << "Track candidates:" << endl;
	for( multimap<TrackFitWeight,Rset_t>::iterator it = fit_chi2.begin();
	     it != fit_chi2.end(); ++it ) {
	  FitResMap_t::iterator found = fit_results.find((*it).second);
	  FitRes_t& r = (*found).second;
	  cout 	 << "ndof = " << r.ndof
		 << " rchi2 = " << r.chi2/(double)r.ndof
		 << " x/y = " << r.coef[0] << "/" << r.coef[2]
		 << " mx/my = " << r.coef[1] << "/" << r.coef[3]
		 << endl;
	}
      }
#endif
      // Select "optimal" set of roads, minimizing sum of chi2s
      vector<Rset_t> best_roads;
      OptimalN( unique_found, fit_chi2, best_roads, 
		AnySharedHits(), CheckTypes(found_types) );

      // Now each selected road tuple corresponds to a new track
      for( vector<Rset_t>::iterator it = best_roads.begin(); it !=
	     best_roads.end(); ++it ) {
	// Retrieve the fit results for this tuple
	FitResMap_t::iterator found = fit_results.find(*it);
	assert( found != fit_results.end() );
	NewTrack( tracks, (*found).second );
      }
    } //if(nfits)
#ifdef TESTCODE
    fN3dFits = nfits;
    t_3dfit = 1e6*timer.RealTime();
#endif
#ifdef VERBOSE
    if( fDebug > 0 ) {
      Int_t ntr = tracks.GetLast()+1;
      cout << ntr << " track";
      if( ntr != 1 ) cout << "s";
      if( nfits > 1 ) cout << " after chi2 optimization";
      if( fDebug > 1 and ntr > 0 ) {
	cout << ":" << endl;
	for( Int_t i = 0; i < ntr; ++i ) {
	  THaTrack* tr = (THaTrack*)tracks.UncheckedAt(i);
	  cout << "3D track:  x/y = " << tr->GetX() << "/" << tr->GetY() 
	       << " mx/my = " << tr->GetTheta() << "/" << tr->GetPhi()
	       << " ndof = "  << tr->GetNDoF()
	       << " rchi2 = " << tr->GetChi2()/(double)tr->GetNDoF()
	       << endl;
	}
      } else
	cout << endl;
    }
  }
  else if( fDebug > 0 ) {
    cout << "No tracks, number of projections with roads = " << nproj;
    if( nproj > 0 ) {
      // Write out the names of the active projections
      cout << " (";
      UInt_t ndone = 0;
      for( EProjType k = kTypeBegin; k < kTypeEnd; ++k ) {
	if( TESTBIT(found_types,k) ) {
	  cout << kProjParam[k].name;
	  if( ++ndone != nproj ) cout << ",";
	}
      }
      assert( ndone == nproj );
      cout << ")";
    }
    cout << ", >=2 required" << endl;
#endif
  } //if(nproj>=2)

#ifdef TESTCODE
    t_coarse = 1e6*timer_tot.RealTime();
#endif
  // Quit here to let detectors CoarseProcess() the approximate tracks,
  // so that they can determine the corrections that we need when we
  // continue in FineTrack
  return 0;
}

//_____________________________________________________________________________
Int_t Tracker::FineTrack( TClonesArray& /* tracks */ )
{
  // Second-level tracking, applying fine corrections.
  //
  // - Correct the hit drift distances using timing, track slope,
  //   track momentum, magnet fringe field, etc., as applicable.
  // - Re-collect hits of roads and re-fit roads
  // - Re-combine track projections in 3D and re-fit 3D track

  if( !TestBit(kDoFine) )
    return 0;

  // TODO

  return 0;;
}

//_____________________________________________________________________________
Int_t Tracker::DefineVariables( EMode mode )
{
  // Initialize global variables

  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list
  
  RVarDef vars[] = {
    { "fail.nhits", "Too many hits in readout plane(s)",  "fFailNhits" },
    { 0 },
  };
  DefineVarsFromList( vars, mode );

  // Define tracking-related variables only if doing tracking
  if( TestBit(kDoCoarse) ) {
    RVarDef vars_tracking[] = {
      { "fail.npat",  "Too many treesearch patterns",       "fFailNpat" },
#ifdef TESTCODE
      { "ncombos",    "Number of road combinations",        "fNcombos" },
      { "nfits",      "Number of 3D track fits done",       "fN3dFits" },
      { "t_track",    "Time in 1st stage tracking (us)",    "t_track" },
      { "t_3dmatch",  "Time in MatchRoads (us)",            "t_3dmatch" },
      { "t_3dfit",    "Time fitting/selecting 3D tracks (us)", "t_3dfit" },
      { "t_coarse",   "Total time in CoarseTrack (us)",     "t_coarse" },
#endif
      { 0 }
    };
    DefineVarsFromList( vars_tracking, mode );
  }
  return 0;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Tracker::Init( const TDatime& date )
{
  // Initialize GEM. Calls standard Init(), then initializes subdetectors.

  static const char* const here = "Init";

  // The "cratemap" is only needed during Init of GEM and Plane
  assert( fCrateMap == 0 );
  fCrateMap = new THashTable(100);
  fCrateMap->SetOwner();

  // Initialize ourselves. This calls our ReadDatabase() and DefineVariables()
  EStatus status = THaTrackingDetector::Init(date);

  // Initialize the planes
  if( !status ) {
    for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      status = fPlanes[iplane]->Init(date);
      if( status )
	break;
    }
  }
  delete fCrateMap; fCrateMap = 0;
  if( status )
    return fStatus = status;

  // Sort planes by increasing z-position
  sort( ALL(fPlanes), Plane::ZIsLess() );

  // Associate planes and partners
  bool all_partnered = true;
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    Plane* thePlane = fPlanes[iplane];
    TString other( thePlane->GetPartnerName() );
    if( other.IsNull() ) {
      thePlane->SetPartner( 0 );
      all_partnered = false;
      continue;
    }
    vriter_t it = find_if( ALL(fPlanes), Plane::NameEquals(other) );
    if( it != fPlanes.end() ) {
      Plane* partner = *it;
      // Sanity checks
      if( thePlane == partner ) {
	Error( Here(here), "Plane %s: cannot partner a plane with itself. "
	       "Fix database.", other.Data() );
	return fStatus = kInitError;
      }
      // Partner planes must be of different types (multi-dimensional readout)
      if( thePlane->GetType() == partner->GetType() ) {
	Error( Here(here), "Partner planes %s and %s have the same type!"
	       " Fix database.", thePlane->GetName(), partner->GetName() );
	return fStatus = kInitError;
      }
      // 2D readouts must have essentially the same z-position
      if( TMath::Abs( thePlane->GetZ() - partner->GetZ() ) > 1e-3 ) {
	Error( Here(here), "Partner planes %s and %s must have the same "
	       "z-position within 1 mm. Fix database.", 
	       thePlane->GetName(), partner->GetName() );
	return fStatus = kInitError;
      }
      // Check for consistency
      TString ppname( partner->GetPartnerName() );
      if( ppname.IsNull() or ppname != thePlane->GetName() ) {
	Error( Here(here), "Inconsistent plane partnering. Partner(%s) "
	       "= %s, but partner(%s) = %s. Fix database.",
	       thePlane->GetName(), other.Data(), 
	       partner->GetName(), ppname.IsNull() ? "(none)":ppname.Data() );
	return fStatus = kInitError;
      }
      // Nothing to do if partner already set (prevents duplicate printouts)
      if( thePlane->GetPartner() == partner ) {
	assert( partner->GetPartner() == thePlane );
	continue;
      }
      // Mutually associate thePlane and partner
      if( fDebug > 0 )
	Info( Here(here), "Partnering plane %s with %s",
	      thePlane->GetName(), partner->GetName() );

      partner->SetPartner( thePlane );

    } else {
      Error( Here(here), "Partner plane %s of %s is not defined!"
	     " Fix database.", other.Data(), thePlane->GetName() );
      return fStatus = kInitError;
    }
  }

  // Set up the projections based on which plane types are defined
  assert( fProj.empty() );
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    Plane* thePlane = fPlanes[iplane];
    assert( !thePlane->GetProjection() );
    EProjType type = thePlane->GetType();
    vpiter_t it = find_if( ALL(fProj), Projection::TypeEquals(type) );
    Projection* proj;
    if( it != fProj.end() ) {
      proj = *it;
    } else {
      try {
	proj = new Projection( type,
			       kProjParam[type].name, 
			       kProjParam[type].angle*TMath::DegToRad(),
			       this 
			       );
      }
      catch( bad_alloc ) { proj = 0; }
      if( !proj or proj->IsZombie() ) {
	// Urgh. Something very bad is going on
	Error( Here(here), "Error creating projection %s. Call expert.", 
	       kProjParam[type].name );
	return fStatus = kInitError;
      }
      fProj.push_back( proj );
    }
    proj->AddPlane( thePlane );
  }

  // Sort projections by ascending EProjType
  sort( ALL(fProj), Projection::ByType() );

  // If exactly two projections are defined, switch on the 3D amplitude
  // correlation matching mode
  if( fProj.size() == 2 ) {
    SetBit( k3dCorrAmpl );
    if( fDebug > 0 )
      Info( Here(here), "2 projections defined: Enabling amplitude "
	    "correlation matching.");

    if( not all_partnered ) {
      Error( Here(here), "Amplitude correlation mode enabled, but one or "
	     "more readout planes are 1D (=no partner). Fix database." );
      return fStatus = kInitError;
    }
  }

  // Initialize the projections. This will read the database and set
  // the projections' angle and maxslope, which we need in the following
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    status = (*it)->Init(date);
    if( status )
      return fStatus = status;
  }

  // Sanity checks of U and V angles which the projections just read via Init
  vpiter_t iu = find_if( ALL(fProj), Projection::TypeEquals(kUPlane) );
  vpiter_t iv = find_if( ALL(fProj), Projection::TypeEquals(kVPlane) );
  if( iu != fProj.end() and iv != fProj.end() ) {
    // Both u and v planes are defined
    //TODO: isn't this an unnecessary limitation?
    //TODO: instead, check if any projections have (nearly) identical angles
    Double_t u_angle = (*iu)->GetAngle()*TMath::RadToDeg();
    Double_t v_angle = (*iv)->GetAngle()*TMath::RadToDeg();
    Int_t qu = TMath::FloorNint( u_angle/90.0 );
    Int_t qv = TMath::FloorNint( v_angle/90.0 );
    if( (qu&1) == (qv&1) ) {
      Error( Here(here), "Plane misconfiguration: uangle (%6.2lf) and vangle "
	     "(%6.2lf) are in equivalent quadrants. Fix database.",
	     u_angle, v_angle );
      return fStatus = kInitError;
    }
    //TODO: put this check into Projection::Init?
    //FIXME: this does not require both u and v
    Double_t du = u_angle - 90.0*qu;
    Double_t dv = v_angle - 90.0*qv;
    if( TMath::Abs(TMath::Abs(du)-45.0) > 44.0 or
	TMath::Abs(TMath::Abs(dv)-45.0) > 44.0 ) {
      Error( Here(here), "uangle (%6.2lf) or vangle (%6.2lf) too close "
	     "to 0 or 90 degrees. Fix database.", u_angle, v_angle );
      return fStatus = kInitError;
    }

    // Check if we can use the simplified 3D matching algorithm
    //FIXME: fast_3d should work for x,y,u(45) and similar, too
    //TODO: generalize to 4 projections with symmetry
    if( fProj.size() == 3 ) {
      assert( !TestBit(k3dCorrAmpl) );
      // This assumes that u and v are the first two defined projections
      assert( kUPlane < kVPlane && kVPlane < 2 );
      // The abs(angle) of the two rotated planes must be (nearly) the same
      Double_t uang = TMath::Abs( TVector2::Phi_mpi_pi( (*iu)->GetAngle() ));
      if( (TMath::Abs( TVector2::Phi_mpi_pi( (*iv)->GetAngle() ))-uang )
	  <  0.5*TMath::DegToRad() ) {
	SetBit(k3dFastMatch);
	Double_t tan = TMath::Tan( 0.5*TMath::Pi()-uang );
	// The scale factor converts the fast_3d matchvalue to the one computed
	// by the general algorithm
	f3dMatchvalScalefact = 2.0 * (1.0/3.0 + tan*tan );
	// Avoid scaling for every event
	f3dMatchCut /= f3dMatchvalScalefact;
      }
    }
  }

  // Determine width and maxslope of the projections
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    Projection* theProj = *it;
    EProjType type = theProj->GetType();
    Double_t width = 0.0;
    //TODO: can this be part of Projection::Init?
    for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      Plane* thePlane = fPlanes[iplane];
      //FIXME: loop over the projection's planes only -> do in Projection::Init
      if( thePlane->GetType() == type ) {
	// Determine the "width" of this projection plane (=width along the
	// projection coordinate).
	// The idea is that all possible hit positions in all planes of a 
	// given projection must fall within the range [-W/2,W/2]. It is normal
	// that some planes cover less than this width (e.g. because they are 
	// smaller or offset) - it is meant to be the enclosing range for all
	// planes. In this way, the tree search can use one fixed bin width.
	// The total width found here divided by the number of bins used in
	// tree search is the resolution of the roads.
	// Of course, this will become inefficient for grotesque geometries.
	Double_t s = thePlane->GetStripStart();
	Double_t d = thePlane->GetStripPitch();
	Double_t n = static_cast<Double_t>( thePlane->GetNelem() );
	Double_t lo = s - 0.5*d;
	Double_t hi = s + (n-0.5)*d;
	Double_t w = max( TMath::Abs(hi), TMath::Abs(lo) );
	if( w > width )
	  width = w;
      }
    }

    //TODO: make this work if projection does not have enough planes (>=2)
    // Set width of this projection
    width *= 2.0;
    if( width > 0.01 )
      theProj->SetWidth( width );
    else {
      Error( Here(here), "Error calculating width of projection \"%s\". "
	     "Strip pitch too small? Fix database.", theProj->GetName() );
      return fStatus = kInitError;
    }
    // maxslope is the maximum expected track slope in the projection.
    // width/depth is the maximum geometrically possible slope. It may be
    // further limited by the trigger acceptance, optics, etc.
    Double_t dz = TMath::Abs(theProj->GetZsize());
    if( dz > 0.01 ) {
      Double_t maxslope = width/dz;
      if( theProj->GetMaxSlope() < 0.01 ) {  // Consider unset
	theProj->SetMaxSlope( maxslope );
      } else if( theProj->GetMaxSlope() > maxslope ) {
	if( fDebug > 0 ) {
	  Warning( Here(here), "For projection \"%s\", maxslope from "
		   "database = %lf exceeds geometric maximum = %lf. "
		   "Using smaller value.",
		   theProj->GetName(), theProj->GetMaxSlope(), maxslope );
	}
	theProj->SetMaxSlope( maxslope );
      }
    } else {
      Error( Here(here), "Error calculating geometric maxslope for plane "
	     "type \"%s\". z-range of planes too small. Fix database.",
	     theProj->GetName() );
	return fStatus = kInitError;
    }

    // Now that the projection's list of planes, width, and maxslope is known,
    // do the level-2 initialization of the projections - load the pattern
    // database and initialize the hitpattern
    //FIXME: no longer need 2-part Init
    status = theProj->InitLevel2(date);
    if( status )
      return fStatus = status;

  }

  // If threading requested, load thread library and start up threads
  if( fMaxThreads > 1 ) {
    gSystem->Load("libThread");
    //FIXME: check for successful load, warn and degrade if not
    delete fThreads;
    fThreads = new ThreadCtrl;
    fThreads->fTrack.reserve( fProj.size() );
    fThreads->fTrackStartM->Lock();
    for( vpsiz_t k = 0; k < fProj.size(); ++k ) {
      TThread* t = fThreads->AddTrackThread( fProj[k] );
      t->Run();
    }
    fThreads->fTrackStartM->UnLock();
  }

  return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t Tracker::ReadDatabase( const TDatime& date )
{
  // Read GEM database

  static const char* const here = "ReadDatabase";

  fIsInit = kFALSE;
  // Delete existing configuration (in case we are re-initializing)
  DeleteContainer( fProj );
  DeleteContainer( fPlanes );
  fCalibPlanes.clear();

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Read fOrigin (detector position) and fSize. fOrigin is the position of
  // the GEM detector relative to some superior coordinate system
  // (typically the spectrometer detector stack reference frame). 
  // fOrigin will be added to all tracks generated; if fOrigin.Z() is not
  // zero, tracks will be projected into the z=0 plane.
  fOrigin.SetXYZ(0,0,0);
  Int_t err = ReadGeometry( file, date );
  if( err )
    return err;

  // Putting this container on the stack may cause strange stack overflows!
  vector<vector<Int_t> > *cmap = new vector<vector<Int_t> >;

  string planeconfig, calibconfig;
  Int_t do_adc_cut = 1;
  f3dMatchCut = 1e-4;
  Int_t event_display = 0, disable_tracking = 0, disable_finetrack = 0;
  Int_t maxmiss = -1, maxthreads = -1, ampcorr_maxmiss = -1;
  fMaxCorrNsigma = 1.0;
  Double_t conf_level = 1e-9;
  ResetBit( k3dCorrAmpl );  // Set in Init()
  ResetBit( k3dFastMatch );
  DBRequest request[] = {
    { "planeconfig",       &planeconfig,       kString },
    { "cratemap",          cmap,               kIntM,   5 },
    { "do_adc_cut",        &do_adc_cut,        kInt,    0, 1 },
    { "3d_matchcut",       &f3dMatchCut,       kDouble, 0, 1 },
    { "event_display",     &event_display,     kInt,    0, 1 },
    { "disable_tracking",  &disable_tracking,  kInt,    0, 1 },
    { "disable_finetrack", &disable_finetrack, kInt,    0, 1 },
    { "calibrate",         &calibconfig,       kString, 0, 1 },
    { "3d_maxmiss",        &maxmiss,           kInt,    0, 1 },
    { "3d_chi2_conflevel", &conf_level,        kDouble, 0, 1 },
    { "3d_ampcorr_maxmiss",&ampcorr_maxmiss,   kInt,    0, 1 },
    { "3d_ampcorr_nsigma", &fMaxCorrNsigma,    kDouble, 0, 1 },
    { "maxthreads",        &maxthreads,        kInt,    0, 1 },
    { 0 }
  };

  Int_t status = kInitError;
  err = LoadDB( file, date, request, fPrefix );
  fclose(file);

  if( !err ) {
    if( cmap->empty() ) {
      Error(Here(here), "No cratemap defined. Set \"cratemap\" in database.");
    } else {
      // Build the list of crate map elements
      for( vviter_t it = cmap->begin(); it != cmap->end(); ++it ) {
	vector<int>& row = *it;
	for( Int_t slot = row[1]; slot <= row[2]; ++slot ) {
	  DAQmodule* m =
	    new DAQmodule( row[0], slot, row[3], row[4] );
	  DAQmodule* found = static_cast<DAQmodule*>(fCrateMap->FindObject(m));
	  if( found ) { 
	    m->Copy(*found);  // Later entries override earlier ones
	    delete m;
	  }
	  else
	    fCrateMap->Add(m);
	}
      }
      status = kOK;
    }
  }
  delete cmap; cmap = 0;
  if( status != kOK )
    return status;

  // Set analysis control flags
  SetBit( kDoADCCut,      do_adc_cut );
  SetBit( kEventDisplay,  event_display );
  SetBit( kDoCoarse,      !disable_tracking );
  SetBit( kDoFine,        !(disable_tracking or disable_finetrack) );
  bool doing_tracking = TestBit(kDoCoarse);

  vector<string> planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Set \"planeconfig\" in database." );
    return kInitError;
  }
  vector<string> calibplanes = vsplit(calibconfig);

  // Set up the readout planes
  for( ssiz_t i=0; i<planes.size(); ++i ) {
    assert( !planes[i].empty() );
    const char* name = planes[i].c_str();
    vriter_t it = find_if( ALL(fPlanes), Plane::NameEquals( name ) );
    if( it != fPlanes.end() ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", name );
      return kInitError;
    }
    Plane* newplane = 0;
    try { newplane = new Plane( name, name, this ); }
    catch( bad_alloc ) { newplane = 0; }
    if( !newplane or newplane->IsZombie() ) {
      // Urgh. Something is very bad
      Error( Here(here), "Error creating readout plane %s. Call expert.",
	     name );
      return kInitError;
    }
    fPlanes.push_back( newplane );
    newplane->SetDebug( fDebug );
    newplane->SetDefinedNum( i );

    // Set calibration mode if requested for any planes
    if( !calibplanes.empty() ) {
      vector<string>::iterator its =
	find_if( ALL(calibplanes), bind2nd(equal_to<string>(), planes[i]) );
      if( its != calibplanes.end() ) {
	Info( Here(here), "Plane %s in calibration mode", name );
	newplane->EnableCalibration();
	fCalibPlanes.push_back( newplane );
	calibplanes.erase( its );
      }
    }
  }

  // Warn if any requested calibration plane(s) do not exist
  if( !calibplanes.empty() ) {
    string s("plane");
    if( calibplanes.size() > 1 )
      s.append("s");
    for( vector<string>::size_type i = 0; i < calibplanes.size(); ++i ) {
      s.append(" ");
      s.append(calibplanes[i]);
    }
    Warning( Here(here), "Requested calibration for undefined plane %s. "
	     "Error in database?", s.c_str() );
  }

  UInt_t nplanes = fPlanes.size();
  if( nplanes < 5 ) {
    Error( Here(here), "Insufficient number of planes = %u. Need at least 5. "
	   "Fix database.", nplanes );
    return kInitError;
  }
  if( fDebug > 0 )
    Info( Here(here), "Loaded %u planes", nplanes );

  if( doing_tracking ) {
    // Convert maximum number of missing hits to ndof of fits
    fMinNdof = ( maxmiss >= 0 ) ? nplanes - 4 - maxmiss : 1;

    // Set parameter for maximum number of amplitude correlation mismatches
    fMaxCorrMismatches = ( ampcorr_maxmiss >= 0 ) ? ampcorr_maxmiss : 0;

    // Determine Chi2 confidence interval limits for the selected CL and the
    // possible degrees of freedom of the 3D track fit
    if( conf_level < 0.0 or conf_level > 1.0 ) {
      Error( Here(here), "Illegal fit confidence level = %lf. "
	     "Must be 0-1. Fix database.", conf_level );
      return kInitError;
    }
    fChisqLimits.clear();
    fChisqLimits.resize( nplanes-3, make_pair<Double_t,Double_t>(0,0) );
    for( vec_pdbl_t::size_type dof = fMinNdof; dof < fChisqLimits.size();
	 ++dof ) {
      fChisqLimits[dof].first  = TMath::ChisquareQuantile( conf_level, dof );
      fChisqLimits[dof].second = 
	TMath::ChisquareQuantile( 1.0-conf_level, dof );
    }
  }
  // If maxthreads set, use it
  if( maxthreads > 0 ) {
    fMaxThreads = max(maxthreads,1);
    if( fMaxThreads > 20 ) { // Sanity limit
      fMaxThreads = 20;
      Warning( Here(here), "Excessive value of maxthreads = %d, "
	       "limited to %u", maxthreads, fMaxThreads );
    }
  } else {
    // If maxthreads not given, automatically use the number of cpus reported
    SysInfo_t sysifo;
    gSystem->GetSysInfo( &sysifo );
    if( sysifo.fCpus > 0 )
      fMaxThreads = sysifo.fCpus;
    else
      fMaxThreads = 1;
  }
  
  fIsInit = kTRUE;
  return kOK;
}

//_____________________________________________________________________________
void Tracker::Print(const Option_t* opt) const
{
  THaTrackingDetector::Print(opt);

  Int_t verbose = 0;
  if( opt ) {
    TString opt_s(opt);
    opt_s.ToLower();
    verbose = opt_s.CountChar('v');
  }

  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    fPlanes[iplane]->Print(opt);
    if( iplane+1 == fPlanes.size() )
      cout << endl;
  }

  for( vpsiz_t k = 0; k < fProj.size(); ++k ) {
    fProj[k]->Print(opt);
    if( verbose > 0 and k+1 != fProj.size() )
      cout << endl;
  }

}


//_____________________________________________________________________________
void Tracker::SetDebug( Int_t level )
{
  // Set debug level of this detector, including all readout planes
  // (subdetectors) and projections

  THaTrackingDetector::SetDebug( level );

  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->SetDebug( level );

  for( vpsiz_t k = 0; k < fProj.size(); ++k )
    fProj[k]->SetDebug( level );
}

//_____________________________________________________________________________
void Tracker::EnableEventDisplay( Bool_t b )
{
  // Enable event display support. Can only be called before initialization.

  if( !fIsInit ) {
    Error( Here("EnableEventDisplay"), "Cannot enable/disable event display "
	   "support after initialization." );
    return;
  }

  SetBit( kEventDisplay, b );
}

//_____________________________________________________________________________
inline
static DAQmodule* FindDAQmodule( UShort_t crate, UShort_t slot, 
				 const THashTable* table )
{
  if( !table ) return 0;
  CSpair m( crate, slot );
  return static_cast<DAQmodule*>( table->FindObject(&m) );
}

//_____________________________________________________________________________
UInt_t Tracker::LoadDAQmodel( THaDetMap::Module* mod ) const
{
  // Update detector map module 'mod' with the model number from the cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  UInt_t num = found ? found->fModel : 0;
  mod->SetModel( num );
  return num;
}

//_____________________________________________________________________________
UInt_t Tracker::GetDAQnchan( THaDetMap::Module* mod ) const
{
  // Return number of channels for detector map module 'mod' from cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  return found ? found->fNchan : 0;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

ClassImp(TreeSearch::Tracker)

