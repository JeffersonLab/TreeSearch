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


#define ALL(c) (c).begin(), (c).end()

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* desc, THaApparatus* app )
  : Tracker(name,desc,app)
    
{ 
  // Constructor

  fRefMap = new THaDetMap;
}

//_____________________________________________________________________________
MWDC::~MWDC()
{
  // Destructor

  delete fRefMap;
}

//_____________________________________________________________________________
Int_t MWDC::Decode( const THaEvData& evdata )
{
  // Decode all planes and fill hitpatterns per projection
  
  static const char* const here = "Decode";

  // Decode reference channels of the VME readout (if any)
  for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fRefMap->GetModule(imod);
    // By construction, this map has one channel per module
    Int_t chan = d->lo;
    UInt_t nhits = evdata.GetNumHits( d->crate, d->slot, chan );
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

  // Standard decoding: decode all projectionsa and fill hitpatterns
  Tracker::Decode( evdata );

  return 0;
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
  if( (qu&1) == (qv&1) ) {
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
TClass* MWDC::GetPlaneClass() const
{
  // Return ROOT class of the Planes used by this implementation

  return WirePlane::Class();
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

ClassImp(TreeSearch::MWDC)

