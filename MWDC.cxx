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

#include <iostream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <map>

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
    if( !obj || !(m = dynamic_cast<const CSpair*>(obj)) ) return kFALSE;
    return ( fCrate == m->fCrate && fSlot == m->fSlot );
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
typedef vector<Projection*>::size_type vpsiz_t;
typedef vector<WirePlane*>::iterator vwiter_t;
typedef vector<vector<Int_t> >::iterator vviter_t;
typedef set<Road*> Rset_t;

struct MWDC::FitRes_t {
  vector<Double_t> coef;
  Double_t matchval;
  Double_t chi2;
  Int_t ndof;
};

// Global constant indicating invalid/uninitialized data
const Double_t kBig = 1e38;

#define ALL(c) (c).begin(), (c).end()

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* desc, THaApparatus* app )
  : THaTrackingDetector(name,desc,app), fCrateMap(0),
    f3dMatchvalScalefact(1), f3dMatchCut(0), fFailNhits(0), fFailNpat(0),
    fNcombos(0), fN3dFits(0)
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
    if( !proj || proj->IsZombie() ) {
      // Urgh. Something very bad is going on
      Error( Here("MWDC"), "Error creating projection %s. Call expert.", 
	     p_name[type] );
      MakeZombie();
      return;
    }
    fProj[type] = proj;
  }

  //  fBench = new THaBenchmark;
}

//_____________________________________________________________________________
MWDC::~MWDC()
{
  // Destructor. Delete objects & subdetectors and unregister variables
  if (fIsSetup)
    RemoveVariables();
  
  //  delete fBench;

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    delete fPlanes[iplane];
  for( vpsiz_t type = 0; type < fProj.size(); ++type )
    delete fProj[type];

  delete fRefMap;
}

//_____________________________________________________________________________
void MWDC::Clear( Option_t* opt )
{
  // Clear event-by-event data, including those of wire planes and projections
  THaTrackingDetector::Clear(opt);
  
  // Clear the planes and projections, but only if we're not called from Init()
  if( !opt || *opt != 'I' ) {
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
      fPlanes[iplane]->Clear(opt);

    for( EProjType type = kTypeBegin; type < kTypeEnd; ++type )
      fProj[type]->Clear(opt);

    if( fRefMap->GetSize() > 0 )
      fRefTime.assign( fRefMap->GetSize(), kBig );  // not strictly necessary
  }
  fN3dFits = fNcombos = fFailNpat = fFailNhits = 0;
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

    fProj[type]->FillHitpattern();
  }

  return 0;
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
  Int_t npoints = 0;
  for( Rvec_t::const_iterator it = roads.begin(); it != roads.end(); ++it ) {
    const Projection* proj = (*it)->GetProjection();
    Double_t cosa = proj->GetCosAngle();
    Double_t sina = proj->GetSinAngle();

    // Fit the 3D track using the points from the best fit in each Road
    const Road::Pvec_t& points = (*it)->GetPoints();
    for( Road::Pvec_t::const_iterator it2 = points.begin();
	 it2 != points.end(); ++it2 ) {
      const Road::Point* p = (*it2);
      ++npoints;
      Double_t Ai[4] = { cosa, cosa * p->z, sina, sina * p->z };
      Double_t s2 = 1.0/(p->res()*p->res());
      for( int j = 0; j<4; ++j ) {
	for( int k = j; k<4; ++k ) {
	  AtA(j,k) += Ai[j] * s2 * Ai[k];
	}
	Aty(j) += Ai[j] * s2 * p->x;
      }
    }
  }
  for( int j = 0; j<4; ++j )
    for( int k = j+1; k<4; ++k )
      AtA(k,j) = AtA(j,k);

  assert( npoints > 4 );
  if( npoints <=4 ) return -1; // Meaningful fit not possible

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

  // Calculate chi2 and update FitCoord data in the WirePlanes
#ifdef VERBOSE
  if( fDebug > 2 ) cout << "Points in 3D fit:" << endl;
#endif
  chi2 = 0;
  for( Rvec_t::const_iterator it = roads.begin(); it != roads.end(); ++it ) {
    Road* rd = *it;
    const Projection* proj = rd->GetProjection();
    Double_t cosa = proj->GetCosAngle();
    Double_t sina = proj->GetSinAngle();
    
    const Road::Pvec_t& points = rd->GetPoints();
    for( Road::Pvec_t::const_iterator it2 = points.begin();
	 it2 != points.end(); ++it2 ) {
      const Road::Point* p = *it2;
      Double_t slope = coef[1]*cosa + coef[3]*sina;
      Double_t x = coef[0]*cosa + coef[2]*sina + slope*p->z;
      Double_t diff = (x - p->x) / p->res();
      chi2 += diff*diff;
#ifdef VERBOSE
      if( fDebug > 2 ) cout << p->hit->GetWirePlane()->GetName() << " " 
			    << "z = " << p->z << " x = " << p->x << endl;
#endif
#ifdef TESTCODE
      // TESTCODE adds 3D results to existing rank 0 hit coords in WirePlanes
      assert( !p->coord.empty() && p->coord.front()->GetRank() == 0 );
      p->coord.front()->Set3DTrkInfo( x, slope );
#else
      // Production code creates the hit coordinates in the WirePlanes
      // using only the points from the 3D track fit(s)
      Double_t slope2d  = rd->GetSlope();
      Double_t trkpos2d = rd->GetPos() + p->z * slope2d;
      p->hit->GetWirePlane()->AddFitCoord
	( FitCoord( p->hit, rd, p->x, trkpos2d, slope2d, x, slope ));
#endif
    }
  }

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
	 << "mx/my = " << coef[1] << " " << coef[3] << " "
	 << "ndof = " << npoints-4 << " chi2 = " << chi2
	 << endl;
#endif
  
  return npoints-4;
}

//_____________________________________________________________________________
THaTrack* MWDC::NewTrack( TClonesArray& tracks, const FitRes_t& fit_par )
{
  // Make new track with given parameters. Used by CoarseTrack.

  THaTrack* newTrack = AddTrack( tracks,
				 fit_par.coef[0], fit_par.coef[2],
				 fit_par.coef[1], fit_par.coef[3] );
  //TODO: make a TrackID?
  assert( newTrack );
  // The "detector" and "fp" TRANSPORT systems are the same for BigBite
  newTrack->SetD( newTrack->GetX(), newTrack->GetY(),
		  newTrack->GetTheta(), newTrack->GetPhi() );
  newTrack->SetChi2( fit_par.chi2, fit_par.ndof );
  return newTrack;
}

//_____________________________________________________________________________
class CheckTypes : public unary_function<Road*,void>
{
public:
  CheckTypes( UInt_t req ) : fReq(req), fActive(0) {}
  void operator() ( const Road* rd )
  { fActive |= 1U << rd->GetProjection()->GetType(); }
  operator bool()       const { return (fActive & fReq) == fReq; }
  bool     operator!()  const { return !((bool)*this); }
  UInt_t   GetTypes()   const { return fActive; }
private:
  UInt_t fReq;
  UInt_t fActive;
};

//_____________________________________________________________________________
template< typename T, typename TestF > void
OptimalN( const set<T>& choices, const multimap< double, set<T> >& weights,
	  vector< set<T> >& picks, TestF testf )
{
  // This is the second-level de-ghosting algorithm, operating on fitted
  // 3D tracks.  It selects the best set of tracks if multiple roads
  // combinations are present in one or more projection.
  // The key value of 'weights' is the chi2/dof of the fit corresponding to
  // the tuple.

  picks.clear();
  set<T> choices_left(choices);

  //TODO: try minimizing chi2 sum
  //  double wsum = 0;
  typename multimap<double,set<T> >::const_iterator it = weights.begin();
  for( ; it != weights.end(); ++it ) {
    const set<T>& tuple = (*it).second;
    if( includes(ALL(choices_left), ALL(tuple)) ) {
      picks.push_back( tuple );
      //      wsum += (*it).first;
      set<T> new_choices_left;
      set_difference( ALL(choices_left), ALL(tuple),
		      inserter(new_choices_left, new_choices_left.end()) );
      if( !for_each(ALL(new_choices_left), testf) )
	break;
      choices_left.swap( new_choices_left );
    }
  }

  return;
}

//_____________________________________________________________________________
Int_t MWDC::CoarseTrack( TClonesArray& tracks )
{
  // Find tracks from the hitpatterns, using the coarse hit drift times
  // uncorrected for track slope, timing offset, fringe field effects etc.

  //  static const char* const here = "CoarseTrack";

  if( fFailNhits )
    return -1;

  vector<Rvec_t>::size_type nproj = 0;
  vector<Rvec_t> roads;
  roads.reserve( kTypeEnd-kTypeBegin );
  bool err = false;
  UInt_t found_types = 0;
  //TODO: multithread this loop
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Projection* proj = fProj[type];

    // For each projection separately, do TreeSearch, build roads, and
    // fit tracks to uncorrected hit positions
    Int_t nrd = proj->Track();

    if( nrd < 0 )
      err = true;

    if( nrd > 0 ) {
      // Count number of projections with at least one road
      ++nproj;
      found_types |= 1U << type;
      // Add pointers to the good roads to the appropriate vector
      roads.push_back( Rvec_t() );
      roads.back().reserve( proj->GetNgoodRoads() );
      for( UInt_t i = 0; i < proj->GetNroads(); ++i ) {
	Road* rd = proj->GetRoad(i);
	assert(rd);
	if( rd->IsGood() )
	  roads.back().push_back(rd);
      }
    }
  }
  // Abort on error (e.g. too many patterns)
  if( err ) {
    fFailNpat = 1;
    return -1;
  }
  assert( roads.size() == nproj );

  // Combine track projections to 3D tracks
  // The standard method requires at least 3 projections. This helps reject
  // noise and, if there is more than one u or v road, correlates u and v
  if( nproj >= 3 ) {
    Rvec_t selected;
    selected.reserve(nproj);
    //TODO: protect against overflow?
    fNcombos = accumulate( ALL(roads), (UInt_t)1, SizeMul<Rvec_t>() );
    //TODO: cut on excessive number of fits
    // Vector holding the results (vectors of roads with good matchval)
    vector< pair<Double_t,Rvec_t> > road_combos;
    road_combos.reserve(fNcombos);
    // Set of the unique roads occurring in the road_combos elements
    Rset_t unique_found;
#ifdef VERBOSE
    if( fDebug > 2 ) 
      cout << "Matching " << nproj << " track projections in 3D (trying "
	   << fNcombos << " combinations):" << endl;
#endif
    Double_t zback = fPlanes.back()->GetZ();
    vector<TVector2> fxpts, bxpts;
    bool fast_3d = TestBit(k3dFastMatch);
    if( fast_3d ) {
      assert( nproj == 3 );
    } else {
      fxpts.reserve( nproj*(nproj-1)/2 );
      bxpts.reserve( nproj*(nproj-1)/2 );
    }
    for( UInt_t i = 0; i < fNcombos; ++i ) {
      NthCombination( i, roads, selected );
      assert( selected.size() == nproj );

      Double_t matchval = 0.0;
      if( fast_3d ) {
	// special case n==3 and symmetric angles of planes 0 and 1:
	//  - intersect first two proj in front and back
	//  - calc perp distances to 3rd, add in quadrature -> matchval
	TVector2 front = selected[0]->Intersect( selected[1], 0.0 );
	if( !fPlanes.front()->Contains(front) )
	  continue;
	TVector2 back = selected[0]->Intersect( selected[1], zback );
	if( !fPlanes.back()->Contains(back) )
	  continue;
	const TVector2& xax = selected[2]->GetProjection()->GetAxis();
	TVector2 xf = selected[2]->GetPos(0.0)   * xax;
	TVector2 xb = selected[2]->GetPos(zback) * xax;
 	TVector2 d1 = front.Proj(xax) - xf;
 	TVector2 d2 = back.Proj(xax)  - xb;
	// The scale factor converts this matchvalue to the one computed by
	// the general algorithm below
	//TODO; weigh with uncertainties?
 	matchval = ( d1.Mod2() + d2.Mod2() ) * f3dMatchvalScalefact;
#ifdef VERBOSE
	if( fDebug > 2 ) {
	  cout << selected[0]->GetProjection()->GetName()
	       << selected[1]->GetProjection()->GetName()
	       << " front = "; front.Print();
	  cout << "front x  = "; xf.Print();
	  cout << "front dist = " << d1.Mod() << endl;
	  cout << selected[0]->GetProjection()->GetName()
	       << selected[1]->GetProjection()->GetName()
	       << " back = "; back.Print();
	  cout << "back x =  "; xb.Print();
	  cout << "back dist = " << d2.Mod() << endl;
	}
#endif
      } else {
	// general algorithm:
	//  - find all front and back intersections [ n(n-1)/2 each ] 
	//  - compute weighted center of gravity of intersection points
	//  - sum dist^2 of points to center of gravity -> matchval
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
	    if( fDebug > 2 ) {
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
	bctr /= static_cast<Double_t>( fxpts.size() );
	if( !fPlanes.front()->Contains(fctr) or
	    !fPlanes.back()->Contains(bctr) )
	  continue;
	for( vector<TVector2>::size_type k = 0; k < fxpts.size(); ++k ) {
	  matchval += (fxpts[k]-fctr).Mod2() + (bxpts[k]-bctr).Mod2();
	}
#ifdef VERBOSE
	if( fDebug > 2 ) {
	  cout << "fctr = "; fctr.Print();
	  cout << "bctr = "; bctr.Print();
	  cout << "matchval = " << matchval;
	}
#endif
	// We could just connect fctr and bctr here to get an approximate
	// 3D track. But the linear minimization below is the right way
	// to do this.

      } //if k3dFastMatch else

      if( matchval < f3dMatchCut ) {
	// Save road combination with good matchvalue
	road_combos.push_back( make_pair(matchval,selected) );
	// Since not all of the roads in 'roads' may make a match,
	// keep track of unique roads found for each projection type
	unique_found.insert( ALL(selected) );
#ifdef VERBOSE
	if( fDebug > 2 ) cout << "ACCEPTED";
#endif
      }
#ifdef VERBOSE
      if( fDebug > 2 ) cout << endl;;
#endif
    } //for fNcombos

    fN3dFits = road_combos.size();

    // Fit each set of matched roads using linear least squares, yielding 
    // the 3D track parameters, x, x'(=mx), y, y'(=my)
    FitRes_t fit_par;
    fit_par.coef.reserve(4);
    if( fN3dFits == 1 ) {
      // If there is only one combo (typical case), life is simple:
      fit_par.matchval    = road_combos.front().first;
      Rvec_t& these_roads = road_combos.front().second;
      fit_par.ndof = FitTrack( these_roads, fit_par.coef, fit_par.chi2 );
      if( fit_par.ndof > 0 )
	NewTrack( tracks, fit_par );
      else
	FitErrPrint( fit_par.ndof );
    }
    else if( fN3dFits > 0 ) {
      // For multiple road combinations, find the set of tracks with the
      // lowest chi2s that uses each road at most once
      roads.clear();
      typedef map<Rset_t, FitRes_t> FitResMap_t;
      FitResMap_t fit_results;
      multimap< double, Rset_t > fit_chi2;
      // Fit all combinations and sort the results by ascending chi2
      for( vector< pair<Double_t,Rvec_t> >::iterator it = road_combos.begin();
	   it != road_combos.end(); ++it ) {
	fit_par.matchval    = (*it).first;
	Rvec_t& these_roads = (*it).second;
	fit_par.ndof = FitTrack( these_roads, fit_par.coef, fit_par.chi2 );
	if( fit_par.ndof > 0 ) {
	  Rset_t road_tuple( ALL(these_roads) );
	  pair<FitResMap_t::iterator,bool> ins =
	    fit_results.insert( make_pair(road_tuple, fit_par) );
	  assert( ins.second );
	  double rchi2 = fit_par.chi2/(double)fit_par.ndof;
	  fit_chi2.insert( make_pair(rchi2, road_tuple) );
	} else
	  FitErrPrint( fit_par.ndof );
      }

      // Select "optimal" set of roads, minimizing sum of chi2s
      vector<Rset_t> best_roads;
      OptimalN( unique_found, fit_chi2, best_roads, CheckTypes(found_types) );

      // Now each selected road tuple corresponds to a new track
      for( vector<Rset_t>::iterator it = best_roads.begin(); it !=
	     best_roads.end(); ++it ) {
	// Retrieve the fit results for this tuple
	FitResMap_t::iterator found = fit_results.find(*it);
	assert( found != fit_results.end() );
	NewTrack( tracks, (*found).second );
      }
    }//if(fN3dFits)
  }// if(nproj>=3)

  // TODO: optional recovery code for missing projection

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

  //TODO

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << "=========== end of event ";
#ifdef TESTCODE
    cout << "# " << fEvNum << " ";
#endif
    cout << "==============" << endl;
  }
#endif
  return 0;;
}

//_____________________________________________________________________________
Int_t MWDC::DefineVariables( EMode mode )
{
  // Initialize global variables and lookup table for decoder


  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list
  
  RVarDef vars[] = {
    { "fail.nhits",   "Too many hits in wire plane(s)",  "fFailNhits" },
    { "fail.npat",    "Too many treesearch patterns",    "fFailNpat" },
    { "ncombos",      "Number of road combinations",     "fNcombos" },
    { "nfits",        "Number of 3D track fits done",    "fN3dFits" },
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
  if( TMath::Abs(TMath::Abs(du)-45.0) > 44.0 ||
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
      f3dMatchvalScalefact = 2.0*(1.0/3.0 + tan*tan );
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
	  assert( !partner || partner->GetProjection() == 0 );
	  // AddPlane() sets the plane number in the WirePlane objects
	  theProj->AddPlane( thePlane, partner );
	  // Save pointer to the projection object with each plane and partner
	  thePlane->SetProjection(theProj);
	  if( partner ) 
	    partner->SetProjection(theProj);
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
	Double_t sina = theProj->GetSinAngle(), cosa = theProj->GetCosAngle();
	TVector3 wp_off = thePlane->GetOrigin() - fOrigin;
	Double_t off = wp_off.X()*sina - wp_off.Y()*cosa;
	Double_t s = thePlane->GetWireStart();
	Double_t d = thePlane->GetWireSpacing();
	Double_t n = static_cast<Double_t>( thePlane->GetNelem() );
	Double_t lo = s - 0.5*d + off;
	Double_t hi = s + (n-0.5)*d + off;
	Double_t w = TMath::Max( TMath::Abs(hi), TMath::Abs(lo) );
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
    // Set width of this projection plane
    width *= 2.0;
    if( width > 0.01 )
      theProj->SetWidth( width );
    else {
      Error( Here(here), "Error calculating width of projection plane \"%s\". "
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
// 	Warning( Here(here), "For plane type \"%s\", maxslope from database = "
// 		 "%lf exceeds geometric maximum = %lf. Using smaller value.",
// 		 theProj->GetName(), theProj->GetMaxSlope(), maxslope );
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

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Putting this container on the stack may cause a stack overflow!
  vector<vector<Int_t> > *cmap = new vector<vector<Int_t> >;

  string planeconfig;
  Int_t time_cut = 1, pairs_only = 0, mc_data = 0, nopartner = 0;
  f3dMatchCut = 1e-4;
  Int_t event_display = false;
  DBRequest request[] = {
    { "planeconfig",  &planeconfig, kString },
    { "cratemap",     cmap,         kIntM,    6 },
    { "timecut",      &time_cut,    kInt,     0, 1 },
    { "pairsonly",    &pairs_only,  kInt,     0, 1 },
    { "MCdata",       &mc_data,     kInt,     0, 1 },
    { "nopartner",    &nopartner,   kInt,     0, 1 },
    { "3d_matchcut",  &f3dMatchCut, kDouble,  0, 1 },
    { "event_display",&event_display, kInt,   0, 1 },
    { 0 }
  };

  Int_t status = kInitError;
  Int_t err = LoadDB( file, date, request, fPrefix );
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

  vector<string> planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Set \"planeconfig\" in database." );
    return kInitError;
  }

  // Set analysis control flags
  SetBit( kDoTimeCut, time_cut );
  SetBit( kPairsOnly, pairs_only );
  SetBit( kMCdata,    mc_data );
  SetBit( kNoPartner, nopartner );

  // Set up the wire planes
  for( ssiz_t i=0; i<planes.size(); ++i ) {
    TString name(planes[i].c_str());
    if( name.IsNull() )
      continue;
    vwiter_t it = find_if( ALL(fPlanes), WirePlane::NameEquals( name ) );
    if( it != fPlanes.end() ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", 
	     name.Data() );
      return kInitError;
    }
    WirePlane* newplane = new WirePlane( name, name, this );
    if( !newplane || newplane->IsZombie() ) {
      // Urgh. Something is very bad
      Error( Here(here), "Error creating wire plane %s. Call expert.", 
	     name.Data() );
      return kInitError;
    }
    fPlanes.push_back( newplane );
  }

  if( fDebug > 0 )
    Info( Here(here), "Loaded %u planes", fPlanes.size() );

  SetBit( kEventDisplay, (event_display != 0) );

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
    if( type == kTypeBegin || verbose > 0 )
      cout << endl;
    fProj[type]->Print(opt);
  }

//   if( verbose > 1 ) {
//     fBench->Print();
//   }
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
// void MWDC::EnableBenchmarks( Bool_t b )
// {
//   fDoBench = b;
// }

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

  if( name && *name ) {
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

//_____________________________________________________________________________
// Int_t MWDC::End(THaRunBase* run)
// {
//     //    fBench->Print();
//   return THaTrackingDetector::End(run);

// }



//_____________________________________________________________________________
#if 0
// TEST
  {
    // Print all defined detector map items
    string prefix = "B.mwdc.";
    THaDetMap* map = new THaDetMap;
    vector<Int_t>* vi = new vector<Int_t>;
    UInt_t flags = THaDetMap::kDoNotClear|THaDetMap::kFillRefChan;
    for( ssiz_t i=0; i<planes.size(); ++i ) {
      string key = prefix + planes[i];
      key += ".detmap";
      if( LoadDBarray(file,date,key.c_str(),*vi) == 0 ) {
	map->Fill(*vi,flags);
      }
      vi->clear();
    }
    map->Sort();
    map->Print();
    delete map;
    delete vi;
  }
#endif	  


} // end namespace TreeSearch

ClassImp(TreeSearch::MWDC)

///////////////////////////////////////////////////////////////////////////////
