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

#include "TString.h"
#include "TMath.h"
#include "THaDetMap.h"
#include "THaEvData.h"
#include "THashTable.h"
#include "TVector2.h"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <utility>

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
typedef vector<Road*> Rvec_t;

#define ALL(c) (c).begin(), (c).end()

// Global constant indicating invalid/uninitialized data
const Double_t kBig = 1e38;

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* desc, THaApparatus* app )
  : THaTrackingDetector(name,desc,app), fCrateMap(0),
    f3dSimpleMatch(false), f3dMatchvalScalefact(1), f3dMatchCut(0)
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
    fProj[type]->Decode( evdata );

  //TODO: check for excessive plane occupancy here and abort if too full
  // (via "max occupancy" parameter or similar)

    fProj[type]->FillHitpattern();

  }

  return 0;
}

//_____________________________________________________________________________
Int_t MWDC::CoarseTrack( TClonesArray& tracks )
{
  // Find tracks from the hitpatterns, using the coarse hit drift times
  // uncorrected for track slope, timing offset, fringe field effects etc.

  vector<Rvec_t>::size_type nproj = 0;
  vector<Rvec_t> roads;
  roads.reserve( kTypeEnd-kTypeBegin );
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Projection* proj = fProj[type];

    // For each projection separately, do TreeSearch, build roads, and
    // fit tracks to uncorrected hit positions
    proj->Track();

    if( proj->GetNgoodRoads() > 0 ) {
      // Count number of projections with at least one road
      ++nproj;
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
  assert( roads.size() == nproj );

  // Combine track projections to 3D tracks
  // The standard method requires at least 3 projections. This helps reject
  // noise and, if there is more than one u or v road, correlates u and v
  if( nproj >= 3 ) {
    multimap<Double_t,Rvec_t> proj_combos;
    TVector2 xax( fProj[kTypeBegin+2]->GetCosAngle(),
		  fProj[kTypeBegin+2]->GetSinAngle() ); // for n==3 algo
    Rvec_t selected;
    selected.reserve(nproj);
    UInt_t n_combos = accumulate( ALL(roads), (UInt_t)1, SizeMul<Rvec_t>() );
    for( UInt_t i = 0; i < n_combos; ++i ) {
      NthCombination( i, roads, selected );
      assert( selected.size() == nproj );

      Double_t matchval = 0.0;
      Double_t zback = fPlanes.back()->GetZ();
      if( f3dSimpleMatch ) {
	assert( nproj == 3 );
	// special case n==3 and symmetric angles of planes 0 and 1:
	//  - intersect first two proj in front and back
	//  - calc perp distances to 3rd, add in quadrature -> matchval

	TVector2 front = selected[0]->Intersect( selected[1], 0.0 );
	TVector2 back  = selected[0]->Intersect( selected[1], zback );
	TVector2 xf    = selected[2]->GetPos(0.0)   * xax;
	TVector2 xb    = selected[2]->GetPos(zback) * xax;
 	TVector2 d1    = front.Proj(xax) - xf;
 	TVector2 d2    = back.Proj(xax)  - xb;
 	matchval = ( d1.Mod2() + d2.Mod2() ) * f3dMatchvalScalefact;
      } else {
	// general algorithm:
	//  - find all front and back intersections [ n(n-1)/2 each ] 
	//  - compute center of gravity of intersection points
	//  - sum dist^2 of points to center of gravity -> matchval
	vector<TVector2> fxpts, bxpts;
	fxpts.reserve( nproj*(nproj-1)/2 );
	bxpts.reserve( nproj*(nproj-1)/2 );
	TVector2 fctr, bctr;
	for( Rvec_t::iterator it1 = selected.begin();
	     it1 != selected.end(); ++it1 ) {
	  for( Rvec_t::iterator it2 = it1+1; it2 != selected.end(); ++it2 ) {
	    fxpts.push_back( (*it1)->Intersect(*it2, 0.0) );
	    bxpts.push_back( (*it1)->Intersect(*it2, zback) );
	    fctr += fxpts.back();
	    bctr += bxpts.back();
#ifdef VERBOSE
	    cout << (*it1)->GetProjection()->GetName()
		 << (*it2)->GetProjection()->GetName()
		 << " front(" << fxpts.size() << ") = ";
	    fxpts.back().Print();
	    cout << (*it1)->GetProjection()->GetName()
		 << (*it2)->GetProjection()->GetName()
		 << " back (" << bxpts.size() << ") = ";
	    bxpts.back().Print();
#endif
	  }
	}
	assert( fxpts.size() == nproj*(nproj-1)/2 );
	assert( bxpts.size() == fxpts.size() );
	Double_t norm = 1.0/static_cast<Double_t>( fxpts.size() );
	fctr *= norm;
	bctr *= norm;
#ifdef VERBOSE
	cout << "fctr = "; fctr.Print();
	cout << "bctr = "; bctr.Print();
#endif
	for( vector<TVector2>::size_type k = 0; k < fxpts.size(); ++k ) {
	  matchval += (fxpts[k]-fctr).Mod2() + (bxpts[k]-bctr).Mod2();
	}
      }

      // Map of proj combos sorted by matchval
      if( matchval < f3dMatchCut ) {
	proj_combos.insert( make_pair(matchval,selected) );
      }
    } // n_combos

    //TODO:
    // now make tracks from the best combos 
  }
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
  cout << "=========== end of event ==============" << endl;
#endif
  return 0;//fNtracks;
}

//_____________________________________________________________________________
Int_t MWDC::DefineVariables( EMode mode )
{
  // Initialize global variables and lookup table for decoder


  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list
  
  RVarDef vars[] = {
#if 0
    { "trackskipped", "Tracking skipped or truncated", "fTooBusy"       },
    { "cproctime",    "Coarse Processing Time",        "fCoarseProcTime"},
    { "fproctime",    "Fine Processing Time",          "fFineProcTime"  },
    { "estngrp",      "Estimated Number of Groups",    "fEstNGroups"    },
    { "estncall",     "Estimated Number of Calls",     "fEstNCalls"     },
    { "ngrp",         "Number of Groups",              "fNGroups"       },
    { "ncall",        "Number of recursive calls",     "fNCalls"        },
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
  if( TMath::Abs(TMath::Abs(du)-45.0) > 44.0 ||
      TMath::Abs(TMath::Abs(dv)-45.0) > 44.0 ) {
    Error( Here(here), "uangle (%6.2lf) or vangle (%6.2lf) too close "
	   "to 0 or 90 degrees. Fix database.", u_angle, v_angle );
    return fStatus = kInitError;
  }

  // Check if we can use the simplified 3D matching algorithm
  if( fProj.size() == 3 ) {
    // This assumes that the first two planes are the rotated ones
    assert( kUPlane < 2 && kVPlane < 2 );
    // The abs(angle) of the two rotated planes should be (nearly) the same
    Double_t uang =
      TMath::Abs( TVector2::Phi_mpi_pi(fProj[kUPlane]->GetAngle()) );
    if( (TMath::Abs( TVector2::Phi_mpi_pi(fProj[kVPlane]->GetAngle()) )-uang )
	<  0.5*TMath::DegToRad() ) {
      f3dSimpleMatch = true;
      Double_t tan = TMath::Tan( 0.5*TMath::Pi()-uang );
      f3dMatchvalScalefact = 2.0*(1.0/3.0 + tan*tan );
    }
  } else
    f3dSimpleMatch = false;

  // Determine the width of and add the wire planes to the projections
  // TODO: this can be multithreaded, too - I think
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
	  // AddPlane() takes care of setting the plane and layer numbers in
	  // the WirePlane objects
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
    UInt_t n = theProj->GetNlayers();
    // Require at least 3 layers per projection
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
  DBRequest request[] = {
    { "planeconfig",  &planeconfig, kString },
    { "cratemap",     cmap,         kIntM,    6 },
    { "timecut",      &time_cut,    kInt,     0, 1 },
    { "pairsonly",    &pairs_only,  kInt,     0, 1 },
    { "MCdata",       &mc_data,     kInt,     0, 1 },
    { "nopartner",    &nopartner,   kInt,     0, 1 },
    { "3d_matchcut",  &f3dMatchCut, kDouble,  0, 1 },
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

  if( verbose > 1 ) {
    //  fBench->Print();
  }
}


//_____________________________________________________________________________
void MWDC::SetDebug( Int_t level )
{
  // Set debug level of this detector, including all wire planes (subdetectors)

  THaTrackingDetector::SetDebug( level );

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->SetDebug( level );
}

//_____________________________________________________________________________
// void MWDC::EnableBenchmarks( Bool_t b )
// {
//   fDoBench = b;
// }

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
