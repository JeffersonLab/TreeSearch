//*-- Author :    Ole Hansen, Jefferson Lab   06-Jun-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
// Reconstruction class for horizontal drift chambers used in BigBite.       //
// Tracking is done using a tree search algorithm (DellOrso et al.,          //
// NIM A 287, 436 (1990)), in essence a recursive template matching method.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "MWDC.h"
#include "WirePlane.h"
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

#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <string>

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
    UInt_t cs = static_cast<UInt_t>(fCrate)<<16 + static_cast<UInt_t>(fSlot);
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
  DAQmodule( UShort_t crate, UShort_t slot, UInt_t model, Double_t res,
	     UInt_t nchan )
    : CSpair(crate, slot), fModel(model), fResolution(res), fNchan(nchan) {}
  virtual ~DAQmodule() {}
  virtual void Copy( TObject& obj ) const {
    TObject::Copy(obj);
    DAQmodule* m = dynamic_cast<DAQmodule*>(&obj);
    if( !m ) return;
    m->fCrate = fCrate; m->fSlot = fSlot; m->fModel = fModel;
    m->fResolution = fResolution, m->fNchan = fNchan;
  }
  virtual void Print( Option_t* ) const {
    cout << "DAQmodule: "
	 << " crate = " << fCrate
	 << " slot = "  << fSlot
	 << " model = " << fModel
	 << " res = "   << fResolution
	 << " nchan = "   << fNchan
	 << endl;
  }
  UInt_t    fModel;
  Double_t  fResolution;
  UInt_t    fNchan;
};

///////////////////////////////////////////////////////////////////////////////

} // end namespace

namespace TreeSearch {

typedef vector<WirePlane*>::size_type vwsiz_t;
typedef vector<WirePlane*>::iterator  vwiter_t;
typedef vector<Projection*>::size_type vpsiz_t;
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
  Int_t*       running; // Bitfield indicating threads to wait for (shared)
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
  friend THaAnalysisObject::EStatus MWDC::Init(const TDatime&);
  friend Int_t MWDC::CoarseTrack(TClonesArray&);
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
  Int_t                fTrackToDo;    // Bitfield of projections to wait for
  TMutex*              fTrackStartM;  // Mutex for start condition
  TMutex*              fTrackDoneM;   // Mutex for done condition
  TCondition*          fTrackStart;   // Start condition
  TCondition*          fTrackDone;    // Finish condition
};

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* desc, THaApparatus* app )
  : THaTrackingDetector(name,desc,app), fCrateMap(0),
    fMaxThreads(1), fThreads(0),
    f3dMatchvalScalefact(1), f3dMatchCut(0), fMinNdof(1),
    fFailNhits(0), fFailNpat(0)
#ifdef TESTCODE
  , fNcombos(0), fN3dFits(0), fEvNum(0)
#endif
{ 
  // Constructor

  fRefMap = new THaDetMap;

  // Set up the projection objects (describing wire planes of same type)
  fProj.resize( kTypeEnd, 0 );
  // Default plane angles and names. Angles can be overridden via the database
  Double_t p_angle[]   = { -60.0, 60.0, 0.0, 90.0 };
  //FIXME: the names should be part of the enum
  const char* p_name[] = { "u", "v", "x", "y" };
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Projection* proj = new Projection( type, 
				       p_name[type], 
				       p_angle[type]*TMath::DegToRad(),
				       this 
				       );
    if( !proj or proj->IsZombie() ) {
      // Urgh. Something very bad is going on
      Error( Here("MWDC"), "Error creating projection %s. Call expert.", 
	     p_name[type] );
      MakeZombie();
      return;
    }
    fProj[type] = proj;
  }
}

//_____________________________________________________________________________
MWDC::~MWDC()
{
  // Destructor. Delete objects & subdetectors and unregister variables
  if (fIsSetup)
    RemoveVariables();
  
  delete fThreads;
  if( fMaxThreads > 1 )
    gSystem->Unload("libThread");

  DeleteContainer( fPlanes );
  DeleteContainer( fProj );
  delete fRefMap;
}

//_____________________________________________________________________________
void MWDC::Clear( Option_t* opt )
{
  // Clear event-by-event data, including those of wire planes and projections
  THaTrackingDetector::Clear(opt);
  
  // Clear the planes and projections, but only if we're not called from Init()
  if( !opt or *opt != 'I' ) {
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
      fPlanes[iplane]->Clear(opt);

    for( EProjType type = kTypeBegin; type < kTypeEnd; ++type )
      fProj[type]->Clear(opt);

    if( fRefMap->GetSize() > 0 )
      fRefTime.assign( fRefMap->GetSize(), kBig );  // not strictly necessary
  }
  fFailNpat = fFailNhits = 0;
#ifdef TESTCODE
  size_t nbytes = (char*)&t_coarse - (char*)&fNcombos + sizeof(t_coarse);
  memset( &fNcombos, 0, nbytes );
#endif
}

//_____________________________________________________________________________
Int_t MWDC::Decode( const THaEvData& evdata )
{
  // Decode all planes and fill hitpatterns per projection
  
  static const char* const here = "Decode";

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

  // Decode reference channels of the VME readout (if any)
  for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fRefMap->GetModule(imod);
    // By construction, this map has one channel per module
    Int_t chan = d->lo;
    Int_t nhits = evdata.GetNumHits( d->crate, d->slot, chan );
    if( nhits > 0 ) {
      Int_t data = evdata.GetData( d->crate, d->slot, chan, nhits-1 );
      if( nhits > 1 ) {
	Warning( Here(here), "%d hits on reference channel %d module %d", 
		 nhits, chan, imod );
      }
      fRefTime[imod] = d->resolution * data;
    } else {
      // TODO: At this point, one could look for backup reference channels
      // to recover the data. Left for later, if needed.
      Warning( Here(here), "No hits on reference channel %d module %d.",
	       chan, imod );
      fRefTime[imod] = kBig;
    }
  }   // modules

  // Decode the planes, then fill the hitpatterns in the projections
  //TODO: multithread
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Int_t nhits = fProj[type]->Decode( evdata );
    // Sanity cut on overfull wire planes. nhits < 0 indicates overflow
    if( nhits < 0 ) {
      fFailNhits = 1;
      continue;
    }
    // No need to fill the hitpattern if no tracking requested
    if( TestBit(kDoCoarse) )
      fProj[type]->FillHitpattern();
  }

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
    p->hit->GetWirePlane()->AddFitCoord
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
    cout << p->hit->GetWirePlane()->GetName() << " " 
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
Action MWDC::ForAllTrackPoints( const Rvec_t& roads, 
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
void MWDC::FitErrPrint( Int_t err ) const
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
Int_t MWDC::FitTrack( const Rvec_t& roads, vector<Double_t>& coef,
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
  //   y_i is the measured coordinate in the i-th wire plane at z_i
  //   P_i is the physical track intersection point with the z_i plane
  //   T_i is the axis unit vector of the i-th plane
  //   a_i is the angle of the coordinate axis of the i-th plane w.r.t. x
  //   x,mx,y,my are the track parameters to be fitted, origin x,y and
  //       slopes mx,my.
  //   
  // "roads" contains a set of Roads that successfully combine in 3-d, one
  // Road* per projection. Each road, in turn, contains a set of (y_i,z_i)
  // coordinates, at most one per wire plane of the projection type, that
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
void MWDC::FindNearestHits( WirePlane* wp, const THaTrack* track, 
			    const Rvec_t& roads ) const
{
  // For the given wire plane, find the hit nearest to the given track
  // and register it in the plane's fit coordinates array.
  // The given roads are the ones generating the given track.
  // This routine is used for efficiency and alignment studies and testing.

  assert( !roads.empty() );

  // Search for the hit with wire position closest to the track crossing
  // position in this plane. The hits are sorted by wire position, so
  // the search can be made fast.
  Double_t z     = wp->GetZ();
  Double_t cosa  = wp->GetProjection()->GetCosAngle();
  Double_t sina  = wp->GetProjection()->GetSinAngle();
  Double_t slope = track->GetTheta()*cosa + track->GetPhi()*sina;
  Double_t x     = track->GetX()*cosa + track->GetY()*sina + slope*z;
  Double_t dmin  = kBig;
  Double_t pmin  = kBig;
  Hit* hmin = 0;
  // Binary search the hits array for the track crossing position x, similar
  // to std::lower_bound(). This finds the first hit with wire position >= x.
  const TSeqCollection* hits = wp->GetHits();
  Int_t first = 0, last = hits->GetSize();
  // GetSize() is incorrect for TObjArrays and TClonesArrays
  const TObjArray* arr = dynamic_cast<const TObjArray*>(hits);
  if( arr )
    last = arr->GetLast()+1;
  Int_t len = last - first;
  Int_t half, middle;
  while( len > 0 ) {
    half = len >> 1;
    middle = first + half;
    if( static_cast<Hit*>(hits->At(middle))->GetWirePos() < x ) {
      first = middle + 1;
      len -= half + 1;
    } else
      len = half;
  }
  // Decide whether the wire >= x or the first one < x are closest.
  // If the track crosses between two adjacent wires, keep both.
  if( last > 0 ) {
    assert( first <= last );
    Hit *hnext = 0, *hprev = 0;
    if( first < last ) {
      hnext = static_cast<Hit*>(hits->At(first));
      assert( hnext->GetWirePos() >= x );
    }
    if( first > 0 ) {
      hprev = static_cast<Hit*>(hits->At(first-1));
      assert( hprev->GetWirePos() < x );
      if( hnext ) {
	assert( hprev->GetWireNum() < hnext->GetWireNum() );
	if( hprev->GetWireNum() + 1 < hnext->GetWireNum() ) {
	  if( x - hprev->GetWirePos() < hnext->GetWirePos() - x )
	    hnext = 0;
	  else
	    hprev = 0;
	}
      }
    }
    // Of the closest wire(s) found, find the closest drift distance.
    // If there are multiple hits one a wire, test all hits - without
    // making assumptions about the order of drift distances
    if( hnext ) {
      hmin = hnext;
      pmin = hnext->GetPosL();
      dmin = TMath::Abs(pmin-x);
      Int_t i = first;
      Hit* h;
      while( ++i < last and (h = static_cast<Hit*>(hits->At(i)))->GetWireNum()
	     == hnext->GetWireNum() ) {
	Double_t d = TMath::Abs(h->GetPosL()-x);
	if( d < dmin ) {
	  dmin = d;
	  hmin = h;
	  pmin = h->GetPosL();
	}
      }
    }
    if( hprev ) {
      Double_t d = TMath::Abs(hprev->GetPosR()-x);
      if( !hmin or d < dmin ) {
	dmin = d;
	hmin = hprev;
	pmin = hprev->GetPosR();
      }
      Int_t i = first-1;
      Hit* h;
      while( --i >= 0 and (h = static_cast<Hit*>(hits->At(i)))->GetWireNum()
	     == hprev->GetWireNum() ) {
	d = TMath::Abs(h->GetPosR()-x);
	if( d < dmin ) {
	  dmin = d;
	  hmin = h;
	  pmin = h->GetPosR();
	}
      }
    }
  }
  // The road vector does not necessarily contain all projections, so
  // search for the road of the type of this wire plane, taking advantage
  // of the fact that the road vector is sorted by type
  Road* rd = 0;
  Int_t k = min( roads.size(), (Rvec_t::size_type)wp->GetType() );
  do {
    if( roads[k]->GetProjection()->GetType() > wp->GetType() )
      --k;
    else {
      if( roads[k]->GetProjection() == wp->GetProjection() )
	rd = roads[k];
      break;
    }
  } while( k>=0 );
  Double_t slope2d = rd ? rd->GetSlope() : kBig;
  Double_t pos2d   = rd ? rd->GetPos() + z * slope2d  : kBig;

  // Finally, record the hit info in the wire plane
  wp->AddFitCoord( FitCoord(hmin, rd, pmin, pos2d, slope2d, x, slope) );
}

//_____________________________________________________________________________
THaTrack* MWDC::NewTrack( TClonesArray& tracks, const FitRes_t& fit_par )
{
  // Make new track with given parameters. Used by CoarseTrack to generate
  // all MWDC tracks

  Double_t x  = fit_par.coef[0];
  Double_t xp = fit_par.coef[1];
  Double_t y  = fit_par.coef[2];
  Double_t yp = fit_par.coef[3];

  THaTrack* newTrack = AddTrack( tracks, x, y, xp, yp );
  //TODO: make a TrackID?
  assert( newTrack );
  // The "detector" and "fp" TRANSPORT systems are the same for BigBite
  newTrack->SetD( x, y, xp, yp );
  newTrack->SetChi2( fit_par.chi2, fit_par.ndof );

  ForAllTrackPoints( *fit_par.roads, fit_par.coef, AddFitCoord() );
  for( Wpvec_t::const_iterator it = fCalibPlanes.begin(); it !=
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
class MWDC::TrackFitWeight
{
  // Object for sorting fits. The smaller a track's TrackFitWeight, the
  // "better" a track is considered to be.
public:
  TrackFitWeight( const MWDC::FitRes_t& fit_par ) :
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
void MWDC::Add3dMatch( const Rvec_t& selected, Double_t matchval,
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
UInt_t MWDC::MatchRoads( vector<Rvec_t>& roads,
			 list< pair<Double_t,Rvec_t> >& combos_found,
			 Rset_t& unique_found )
{
  // Match roads from different projections that intersect in the front
  // and back plane of the detector

  // The number of projections that we work with (must be >= 3)
  vector<Rvec_t>::size_type nproj = roads.size();
  assert( nproj >= 2 );

  combos_found.clear();
  unique_found.clear();

  // Number of all possible combinations of the input roads
  //TODO: protect against overflow?
  UInt_t ncombos = accumulate( ALL(roads), (UInt_t)1, SizeMul<Rvec_t>() );
  UInt_t nfound  = 0;
  bool   fast_3d = TestBit(k3dFastMatch);

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << "Matching ";
    for( vector<Rvec_t>::size_type i = 0; i < nproj; ++i ) {
      cout << roads[i].size();
      if( i+1 < nproj ) cout << "x";
    }
    cout << " track projections in 3D (";
    if( fast_3d ) cout << "fast";
    else          cout << "generic";
    cout << " algo, " << ncombos << " combinations):" << endl;
  }
#endif
#ifdef TESTCODE
  fNcombos = ncombos;
#endif

  if( ncombos == 0 )
    return 0;

  // TODO: separate code sections
  // Vector holding a combination of roads to test. One road from each 
  // projection 
  Rvec_t selected;
  selected.reserve(nproj);

  Double_t zback = fPlanes.back()->GetZ();
  if( fast_3d ) {
    // special case n==3 and symmetric angles of planes 0 and 1:
    //  - intersect first two proj in front and back
    //  - calc perp distances to 3rd, add in quadrature -> matchval
    assert( nproj == 3 && fProj.size() == 3 );

    // Fetch coefficients for coordinate transformations
    Prvec_t::iterator ip = fProj.begin();
    Double_t su = (*ip)->GetSinAngle();
    Double_t cu = (*ip)->GetCosAngle();
    ++ip;
    Double_t sv = (*ip)->GetSinAngle();
    Double_t cv = (*ip)->GetCosAngle();
    Double_t inv_denom = 1.0/(sv*cu-su*cv);
    // Only the scaled coefficients are needed (cf. Road::Intersect)
    su *= inv_denom; cu *= inv_denom; sv *= inv_denom; cv *= inv_denom;
    // Components of the 3rd projection's axis
    ++ip;
    Double_t xax_x = ((*ip)->GetAxis()).X();
    Double_t xax_y = ((*ip)->GetAxis()).Y();

    // The selected roads from each of the three projections
    Road* tuple[3];
    // Number of roads in u/v projections
    UInt_t nrd[2] = { roads[0].size(), roads[1].size() };
    // Indices of currently selected u/v pair
    UInt_t ird[2];
    WirePlane *front_plane = fPlanes.front(), *back_plane = fPlanes.back();

    // For fast access to the relevant position range, sort the 3rd projection
    // by ascending front position
    sort( ALL(roads[2]), Road::PosIsLess() );
    Road::PosIsNear pos_near( TMath::Sqrt(f3dMatchCut) );

    Double_t matchval = 0.0;
    // Time-critical loop, may run O(1e5) times per event with noisy input
    ird[0] = 0;
    while( ird[0] != nrd[0] ) {
      tuple[0] = roads[0][ird[0]];
      Double_t uf = tuple[0]->GetPos();
      Double_t ub = tuple[0]->GetPos(zback);
      Double_t usf = uf * sv;
      Double_t ucf = uf * cv;
      Double_t usb = ub * sv;
      Double_t ucb = ub * cv;
      ird[1] = 0;
      while( ird[1] != nrd[1] ) {
	tuple[1] = roads[1][ird[1]];
	Double_t v = tuple[1]->GetPos();
	Double_t xf = usf - v * su;
	Double_t yf = v * cu - ucf;
	if( front_plane->Contains(xf,yf) ) {
	  v = tuple[1]->GetPos(zback);
	  Double_t xb = usb - v * su;
	  Double_t yb = v * cu - ucb;
	  if( back_plane->Contains(xb,yb) ) {
	    Double_t pf = xf*xax_x + yf*xax_y; // front x from u/v
	    Double_t pb = xb*xax_x + yb*xax_y; // back x from u/v
	    // Find range of roads in 3rd projection near this front x
	    pair<Rvec_t::iterator,Rvec_t::iterator> range = 
	      equal_range( ALL(roads[2]), pf, pos_near );
	    // Test the candidate x-roads for complete matches
	    for( Rvec_t::iterator it=range.first; it != range.second; ++it ) {
	      tuple[2] = *it;
	      Double_t d1 = tuple[2]->GetPos()      - pf;
	      Double_t d2 = tuple[2]->GetPos(zback) - pb;
	      //TODO; weigh with uncertainties?
	      matchval = d1*d1 + d2*d2;
#ifdef VERBOSE
	      if( fDebug > 3 ) {
		if( matchval < f3dMatchCut || fDebug > 4 ) {
		  cout << tuple[0]->GetProjection()->GetName()
		       << tuple[1]->GetProjection()->GetName()
		       << " front = " << "(" << xf << "," << yf << ")" << endl;
		  cout << "front x  = (" << tuple[2]->GetPos() * xax_x << ","
		       << tuple[2]->GetPos() * xax_y << ")" << endl;
		  cout << "front dist = " << d1 << endl;
		  cout << tuple[0]->GetProjection()->GetName()
		       << tuple[1]->GetProjection()->GetName()
		       << " back = " << "(" << xb << "," << yb << ")" << endl;
		  cout << "back x =  (" 
		       << tuple[2]->GetPos(zback) * xax_x << ","
		       << tuple[2]->GetPos(zback) * xax_y << ")" << endl;
		  cout << "back dist = " << d2 << endl;
		  cout << "matchval = " << matchval*f3dMatchvalScalefact
		       << endl;
		}
	      }
#endif
	      // Check if match, if so then add
	      if( matchval < f3dMatchCut ) {
		++nfound;
		selected.assign( tuple, tuple+3 );
		Add3dMatch( selected, matchval, combos_found, unique_found );
	      }
	    }
	  }
	}
	++ird[1];
      }
      ++ird[0];
    }

  } else {
    // general algorithm:
    //  - find all front and back intersections [ n(n-1)/2 each ] 
    //  - compute weighted center of gravity of intersection points
    //  - sum dist^2 of points to center of gravity -> matchval

    vector<TVector2> fxpts, bxpts;
    fxpts.reserve( nproj*(nproj-1)/2 );
    bxpts.reserve( nproj*(nproj-1)/2 );

    for( UInt_t i = 0; i < ncombos; ++i ) {
      Double_t matchval = 0.0;

      NthCombination( i, roads, selected );
      assert( selected.size() == nproj );

      fxpts.clear();
      bxpts.clear();
      TVector2 fctr, bctr;
      for( Rvec_t::iterator it1 = selected.begin();
	   it1 != selected.end(); ++it1 ) {
	for( Rvec_t::iterator it2 = it1+1; it2 != selected.end(); ++it2 ) {
	  //TODO: weigh with uncertainties of coordinates?
	  fxpts.push_back( (*it1)->Intersect(*it2, 0.0) );
	  bxpts.push_back( (*it1)->Intersect(*it2, zback) );
	  fctr += fxpts.back();
	  bctr += bxpts.back();
#ifdef VERBOSE
	  if( fDebug > 3 ) {
	    cout << (*it1)->GetProjection()->GetName()
		 << (*it2)->GetProjection()->GetName()
		 << " front(" << fxpts.size() << ") = ";
	    fxpts.back().Print();
	    cout << (*it1)->GetProjection()->GetName()
		 << (*it2)->GetProjection()->GetName()
		 << " back (" << bxpts.size() << ") = ";
	    bxpts.back().Print();
	  }
#endif
	}
      }
      assert( fxpts.size() <= nproj*(nproj-1)/2 );
      assert( bxpts.size() == fxpts.size() );
      fctr /= static_cast<Double_t>( fxpts.size() );
      if( !fPlanes.front()->Contains(fctr) )
	continue;
      bctr /= static_cast<Double_t>( fxpts.size() );
      if( !fPlanes.back()->Contains(bctr) )
	continue;
      for( vector<TVector2>::size_type k = 0; k < fxpts.size(); ++k ) {
	matchval += (fxpts[k]-fctr).Mod2() + (bxpts[k]-bctr).Mod2();
      }
#ifdef VERBOSE
      if( fDebug > 3 ) {
	cout << "fctr = "; fctr.Print();
	cout << "bctr = "; bctr.Print();
	cout << "matchval = " << matchval << endl;
      }
#endif
      // We could just connect fctr and bctr here to get an approximate
      // 3D track. But the linear minimization below is the right way
      // to do this.

      if( matchval < f3dMatchCut ) {
	++nfound;
	Add3dMatch( selected, matchval, combos_found, unique_found );
      }
    } //for(ncombos)

  } //if(fast_3d) else

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
Bool_t MWDC::PassTrackCuts( const FitRes_t& fit_par ) const
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
Int_t MWDC::CoarseTrack( TClonesArray& tracks )
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
    Int_t all_todo = BIT(kTypeEnd)-1 & ~(BIT(kTypeBegin)-1);
    Int_t max_todo = BIT(fMaxThreads)-1;
    Int_t start = kTypeBegin;
    while( start < kTypeEnd ) {
      Int_t now_todo = (max_todo << start) & all_todo;
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
    for( Prvec_t::iterator it = fProj.begin(); it != fProj.end(); ++it ) {
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
  roads.reserve( kTypeEnd-kTypeBegin );
  UInt_t found_types = 0;
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Projection* proj = fProj[type];

    Int_t nrd = proj->GetNgoodRoads();
    if( nrd > 0 ) {
      // Count number of projections with at least one road
      ++nproj;
      found_types |= 1U << type;
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
  // The standard method requires at least 3 projections. This helps reject
  // noise and, if there is more than one u or v road, correlates u and v
  if( nproj >= 3 ) {
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
#endif
  }
  else if( nproj == 2 ) {
    // TODO: optional recovery code for missing projection

  } //if(nproj>=3)

#ifdef TESTCODE
    t_coarse = 1e6*timer_tot.RealTime();
#endif
  // Quit here to let detectors CoarseProcess() the approximate tracks,
  // so that they can determine the corrections that we need when we
  // continue in FineTrack
  return 0;
}

//_____________________________________________________________________________
Int_t MWDC::FineTrack( TClonesArray& tracks )
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
Int_t MWDC::DefineVariables( EMode mode )
{
  // Initialize global variables

  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list
  
  RVarDef vars[] = {
    { "fail.nhits", "Too many hits in wire plane(s)",  "fFailNhits" },
    { "fail.npat",  "Too many treesearch patterns",    "fFailNpat" },
#ifdef TESTCODE
    { "ncombos",    "Number of road combinations",     "fNcombos" },
    { "nfits",      "Number of 3D track fits done",    "fN3dFits" },
    { "t_track",    "Time in 1st stage tracking (us)", "t_track" },
    { "t_3dmatch",  "Time in MatchRoads (us)",         "t_3dmatch" },
    { "t_3dfit",    "Time fitting/selecting 3D tracks (us)", "t_3dfit" },
    { "t_coarse",   "Total time in CoarseTrack (us)",  "t_coarse" },
#endif
    { 0 }
  };
  DefineVarsFromList( vars, mode );
  return 0;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus MWDC::Init( const TDatime& date )
{
  // Initialize MWDC. Calls standard Init(), then initializes subdetectors.

  static const char* const here = "Init";

  fRefMap->Reset();
  fRefTime.clear();

  // The "cratemap" is only needed during Init of MWDC and WirePlane
  fCrateMap = new THashTable(100);
  fCrateMap->SetOwner();

  // Initialize ourselves. This calls our ReadDatabase() and DefineVariables()
  EStatus status = THaTrackingDetector::Init(date);

  // Initialize the wire planes
  if( !status ) {
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      status = fPlanes[iplane]->Init(date);
      if( status )
	break;
    }
  }
  delete fCrateMap; fCrateMap = 0;
  if( status )
    return fStatus = status;

  // Automatically create a separate detector map for the reference channels
  // used by the data channels. By construction, each "module" in this map
  // has exactly one channel.
  // This map is decoded in our Decode() for efficiency.
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    THaDetMap* planeMap = fPlanes[iplane]->GetDetMap();
    for( Int_t imod = 0; imod < planeMap->GetSize(); ++imod ) {
      THaDetMap::Module* d = planeMap->GetModule(imod);
      if( d->refchan < 0 ) continue;
      THaDetMap::Module* refmod = fRefMap->Find( d->crate, d->slot,
						 d->refchan );
      // To allow the data channels to look up their reference channel data
      // quickly, we save the index into the refmap in d->refindex.
      if( refmod ) {
	d->refindex = refmod->first;
      }	else {
	d->refindex = fRefMap->GetSize();
	fRefMap->AddModule( d->crate, d->slot, d->refchan, 
			    d->refchan, d->refindex, d->model );
	refmod = fRefMap->GetModule(d->refindex);
	if( !refmod ) {
	  Error( Here(here), "Error when adding reference channel cr/sl/ch="
		 "%u/%u/%u. Call expert.", d->crate, d->slot, d->refchan );
	  return fStatus = kInitError;
	}
	// The reference channel module's model and resolution have already
	// been determined when the plane's detector map was filled. Also,
	// refindex was checked to be a vaid channel number for this model.
	refmod->SetResolution( d->resolution );
	refmod->MakeTDC();
      }
    }
  }
  // Check if any reference channels are also defined as data channels
  for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
    THaDetMap::Module* r = fRefMap->GetModule(imod);
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      THaDetMap* planeMap = fPlanes[iplane]->GetDetMap();
      if( planeMap->Find( r->crate, r->slot, r->lo ) ) {
	Error( Here(here), "Reference channel cr/sl/ch=%u/%u/%u also defined "
	       "as data channel in plane \"%s\". Fix database.", 
	       r->crate, r->slot, r->lo, fPlanes[iplane]->GetName() );
	return fStatus = kInitError;
      }
    }
  }
  // Reference map set up correctly and without conflicts
  fRefTime.assign( fRefMap->GetSize(), kBig );

  // Sort planes by increasing z-position
  sort( ALL(fPlanes), WirePlane::ZIsLess() );

  // Associate planes and partners (names identical except trailing "p")
  if( !TestBit(kNoPartner) ) {
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      WirePlane* thePlane = fPlanes[iplane];
      TString name( thePlane->GetName() );
      if( name.EndsWith("p") ) {
	TString other = name.Chop();
	if( other.IsNull() )
	  continue;
	vwiter_t it = find_if( ALL(fPlanes),WirePlane::NameEquals( other ) );
	if( it != fPlanes.end() ) {
	  WirePlane* partner = *it;
	  // Partner planes must be of the same type!
	  if( thePlane->GetType() != partner->GetType() ) {
	    Error( Here(here), "Partner planes %s and %s have different types!"
		   " Fix database.", thePlane->GetName(), partner->GetName() );
	    return fStatus = kInitError;
	  }
	  if( fDebug > 0 )
	    Info( Here(here), "Partnering plane %s with %s",
		  thePlane->GetName(), partner->GetName() );
	  partner->SetPartner( thePlane );
	}
      }
    }
  }

  // Initialize the projections. This will read the database and set
  // the projections' angle and maxslope, which we need in the following
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    status = fProj[type]->Init(date);
    if( status )
      return fStatus = status;
  }

  // Sanity checks of U and V angles which the projections just read via Init
  Double_t u_angle = fProj[kUPlane]->GetAngle()*TMath::RadToDeg();
  Double_t v_angle = fProj[kVPlane]->GetAngle()*TMath::RadToDeg();
  Int_t qu = TMath::FloorNint( u_angle/90.0 );
  Int_t qv = TMath::FloorNint( v_angle/90.0 );
  if( qu&1 == qv&1 ) {
    Error( Here(here), "Plane misconfiguration: uangle (%6.2lf) and vangle "
	   "(%6.2lf) are in equivalent quadrants. Fix database.",
	   u_angle, v_angle );
    return fStatus = kInitError;
  }
  Double_t du = u_angle - 90.0*qu;
  Double_t dv = v_angle - 90.0*qv;
  if( TMath::Abs(TMath::Abs(du)-45.0) > 44.0 or
      TMath::Abs(TMath::Abs(dv)-45.0) > 44.0 ) {
    Error( Here(here), "uangle (%6.2lf) or vangle (%6.2lf) too close "
	   "to 0 or 90 degrees. Fix database.", u_angle, v_angle );
    return fStatus = kInitError;
  }

  // Check if we can use the simplified 3D matching algorithm
  ResetBit(k3dFastMatch);
  if( fProj.size() == 3 ) {
    // This assumes that the first two planes are the rotated ones
    assert( kUPlane < 2 && kVPlane < 2 );
    // The abs(angle) of the two rotated planes should be (nearly) the same
    Double_t uang =
      TMath::Abs( TVector2::Phi_mpi_pi(fProj[kUPlane]->GetAngle()) );
    if( (TMath::Abs( TVector2::Phi_mpi_pi(fProj[kVPlane]->GetAngle()) )-uang )
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

  // Determine the width of and add the wire planes to the projections
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Projection* theProj = fProj[type];
    Double_t width = 0.0;
    // Associate planes with plane types
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      WirePlane* thePlane = fPlanes[iplane];
      if( thePlane->GetType() == type ) {
	// Add planes not yet used to the corresponding projection
	if( !thePlane->GetProjection() ) {
	  WirePlane* partner = thePlane->GetPartner();
	  assert( !partner or partner->GetProjection() == 0 );
	  // AddPlane() sets the plane number, fProjection and fCoordOffset
	  // in the WirePlane objects
	  theProj->AddPlane( thePlane, partner );
	}
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
	Double_t s = thePlane->GetWireStart();
	Double_t d = thePlane->GetWireSpacing();
	Double_t n = static_cast<Double_t>( thePlane->GetNelem() );
	Double_t lo = s - 0.5*d;
	Double_t hi = s + (n-0.5)*d;
	Double_t w = max( TMath::Abs(hi), TMath::Abs(lo) );
	if( w > width )
	  width = w;
      }
    }
    UInt_t n = theProj->GetNplanes();
    // Require at least 3 planes per projection
    if( n < 3 ) {
      Error( Here(here), "Too few planes of type \"%s\" defined. "
	     "Need >= 3, have %u. Fix database.", theProj->GetName(), n );
      return fStatus = kInitError;
    }
    // Set width of this projection
    width *= 2.0;
    if( width > 0.01 )
      theProj->SetWidth( width );
    else {
      Error( Here(here), "Error calculating width of projection \"%s\". "
	     "Wire spacing too small? Fix database.", theProj->GetName() );
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

    // Now that the projection's list of planes, width, and maxslope is know,
    // do the level-2 initialization of the projections - load the pattern
    // database and initialize the hitpattern
    status = fProj[type]->InitLevel2(date);
    if( status )
      return fStatus = status;

  }

  // If threading requested, load thread library and start up threads
  if( fMaxThreads > 1 ) {
    gSystem->Load("libThread");
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
Int_t MWDC::ReadDatabase( const TDatime& date )
{
  // Read MWDC database

  static const char* const here = "ReadDatabase";

  fIsInit = kFALSE;
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    delete fPlanes[iplane];
  }
  // Delete existing configuration (in case we are re-initializing)
  fPlanes.clear();
  fCalibPlanes.clear();

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Read fOrigin (detector position) and fSize. fOrigin is the position of
  // the MWDC detector relative to some superior coordinate system
  // (typically the spectrometer detector stack reference frame). 
  // fOrigin will be added to all tracks generated; if fOrigin.Z() is not
  // zero, tracks will be projected into the z=0 plane.
  Int_t err = ReadGeometry( file, date );
  if( err )
    return err;

  // Putting this container on the stack may cause strange stack overflows!
  vector<vector<Int_t> > *cmap = new vector<vector<Int_t> >;

  string planeconfig, calibconfig;
  Int_t time_cut = 1, pairs_only = 0, mc_data = 0, nopartner = 0;
  f3dMatchCut = 1e-4;
  Int_t event_display = 0, disable_tracking = 0, disable_finetrack = 0;
  Int_t maxmiss = -1, maxthreads = -1;
  Double_t conf_level = 1e-9;
  DBRequest request[] = {
    { "planeconfig",       &planeconfig,       kString },
    { "cratemap",          cmap,               kIntM,   6 },
    { "timecut",           &time_cut,          kInt,    0, 1 },
    { "pairsonly",         &pairs_only,        kInt,    0, 1 },
    { "MCdata",            &mc_data,           kInt,    0, 1 },
    { "nopartner",         &nopartner,         kInt,    0, 1 },
    { "3d_matchcut",       &f3dMatchCut,       kDouble, 0, 1 },
    { "event_display",     &event_display,     kInt,    0, 1 },
    { "disable_tracking",  &disable_tracking,  kInt,    0, 1 },
    { "disable_finetrack", &disable_finetrack, kInt,    0, 1 },
    { "calibrate",         &calibconfig,       kString, 0, 1 },
    { "3d_maxmiss",        &maxmiss,           kInt,    0, 1 },
    { "3d_chi2_conflevel", &conf_level,        kDouble, 0, 1 },
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
	    new DAQmodule( row[0], slot, row[3], row[4]*1e-12, row[5] );
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
  SetBit( kDoTimeCut,     time_cut );
  SetBit( kPairsOnly,     pairs_only );
  SetBit( kMCdata,        mc_data );
  SetBit( kNoPartner,     nopartner );
  SetBit( kEventDisplay,  event_display );
  SetBit( kDoCoarse,      !disable_tracking );
  SetBit( kDoFine,        !(disable_tracking or disable_finetrack) );

  vector<string> planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Set \"planeconfig\" in database." );
    return kInitError;
  }
  vector<string> calibplanes = vsplit(calibconfig);

  // Set up the wire planes
  for( ssiz_t i=0; i<planes.size(); ++i ) {
    assert( !planes[i].empty() );
    const char* name = planes[i].c_str();
    vwiter_t it = find_if( ALL(fPlanes), WirePlane::NameEquals( name ) );
    if( it != fPlanes.end() ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", name );
      return kInitError;
    }
    WirePlane* newplane = new WirePlane( name, name, this );
    if( !newplane or newplane->IsZombie() ) {
      // Urgh. Something is very bad
      Error( Here(here), "Error creating wire plane %s. Call expert.", name );
      return kInitError;
    }
    fPlanes.push_back( newplane );
    newplane->SetDebug( fDebug );

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
    Warning( Here(here), "Requested calibration for undefined %s. "
	     "Typo in database?", s.c_str() );
  }

  UInt_t nplanes = fPlanes.size();
  if( fDebug > 0 )
    Info( Here(here), "Loaded %u planes", nplanes );

  // Convert maximum number of missing hits to ndof of fits
  if( maxmiss >= 0 )
    fMinNdof = nplanes - 4 - maxmiss;
  else
    fMinNdof = 1;

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
    fChisqLimits[dof].second = TMath::ChisquareQuantile( 1.0-conf_level, dof );
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
void MWDC::Print(const Option_t* opt) const
{
  THaTrackingDetector::Print(opt);

  Int_t verbose = 0;
  if( opt ) {
    TString opt_s(opt);
    opt_s.ToLower();
    verbose = opt_s.CountChar('v');
  }

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->Print(opt);

  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    if( type == kTypeBegin or verbose > 0 )
      cout << endl;
    fProj[type]->Print(opt);
  }

}


//_____________________________________________________________________________
void MWDC::SetDebug( Int_t level )
{
  // Set debug level of this detector, including all wire planes (subdetectors)
  // and projections

  THaTrackingDetector::SetDebug( level );

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->SetDebug( level );

  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type )
    fProj[type]->SetDebug( level );
}

//_____________________________________________________________________________
void MWDC::EnableEventDisplay( Bool_t b )
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
EProjType MWDC::NameToType( const char* name )
{
  // Return the index corresponding to the given plane name.
  // The comparison is not case-sensitive.

  if( name and *name ) {
    TString s(name);
    s.ToLower();
    for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
      if( !fProj[type] )
	continue;
      TString ps( fProj[type]->GetName() );
      ps.ToLower();
      if( s == ps )
	return type;
    }
  }
  return kUndefinedType;
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
UInt_t MWDC::LoadDAQmodel( THaDetMap::Module* mod ) const
{
  // Update detector map module 'mod' with the model number from the cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  UInt_t num = found ? found->fModel : 0;
  mod->SetModel( num );
  return num;
}

//_____________________________________________________________________________
Double_t MWDC::LoadDAQresolution( THaDetMap::Module* mod ) const
{
  // Update detector map module 'mod' with the resolution from the cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  Double_t res = found ? found->fResolution : 0.0;
  mod->SetResolution( res );
  return res;
}

//_____________________________________________________________________________
UInt_t MWDC::GetDAQnchan( THaDetMap::Module* mod ) const
{
  // Return number of channels for detector map module 'mod' from cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  return found ? found->fNchan : 0;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

ClassImp(TreeSearch::MWDC)

