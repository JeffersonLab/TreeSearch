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
#include "TString.h"
#include "TMath.h"
#include "THaDetMap.h"
#include "THaEvData.h"

#include <algorithm>

using namespace std;
typedef string::size_type ssiz_t;

namespace TreeSearch {

typedef vector<WirePlane*>::size_type vwsiz_t;
typedef vector<Projection*>::size_type vpsiz_t;
typedef vector<WirePlane*>::iterator vwiter_t;

const Double_t kBig = 1e38;

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* desc, THaApparatus* app )
  : THaTrackingDetector(name,desc,app)
{ 
  // Constructor

  fProj.resize( kTypeEnd, NULL );

  fRefMap = new THaDetMap;

  //FIXME: test
  fDebug = 1;


#if 0
  // Default behavior for now
  // SetBit( kHardTDCcut );
  //  ResetBit( kIgnoreNegDrift );
  SetBit( kIgnoreNegDrift );
 

  fBench = new THaBenchmark;
#endif
}

//_____________________________________________________________________________
MWDC::~MWDC()
{
  // Destructor. Delete objects & subdetectors and unregister variables
 if (fIsSetup)
   RemoveVariables();
  
#if 0 
  delete fBench;
#endif

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    delete fPlanes[iplane];
  for( vpsiz_t type = 0; type < fProj.size(); ++type )
    delete fProj[type];

  delete fRefMap;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus MWDC::Init( const TDatime& date )
{
  // Initialize MWDC. Calls standard Init(), then initializes subdetectors.

  static const char* const here = "Init";

  // Initialize ourselves. This calls our ReadDatabase().
  EStatus status = THaTrackingDetector::Init(date);
  if( status )
    return fStatus = status;

  // Initialize the plane type "subdetector" - i.e. read database
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    status = fProj[type]->Init(date);
    if( status )
      return fStatus = status;
  }

  // Sanity checks of U and V angles
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
    return kInitError;
  }

  // Iinitialize the wire planes
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    status = fPlanes[iplane]->Init(date);
    if( status )
      return fStatus = status;
  }

  // Sort planes by increasing z-position
  sort( fPlanes.begin(), fPlanes.end(), WirePlane::ZIsLess() );

  // Determine per-projection plane parameters
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Projection* theProj = fProj[type];
    UInt_t n = 0;
    Double_t width = 0.0;
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      WirePlane* thePlane = fPlanes[iplane];
      if( thePlane->GetType() == type ) {
	// Add only primary planes (i.e. the first one of a partnered pair)
	// to the list and count them
	if( !thePlane->GetProjection() ) {
	  theProj->AddPlane( thePlane );
	  thePlane->SetPlaneNum(n);
	  // Save pointer to the projection object with each plane and partner
	  thePlane->SetProjection(theProj);
	  if( thePlane->GetPartner() ) {
	    WirePlane* partner = thePlane->GetPartner();
	    partner->SetProjection(theProj);
	    partner->SetPlaneNum(n);
	  }
	  ++n;
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
	     "Wire spacing too small. Fix database.", theProj->GetName() );
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
	Warning( Here(here), "For plane type \"%s\", maxslope from database = "
		 "%lf exceeds geometric maximum = %lf. Using smaller value.",
		 theProj->GetName(), theProj->GetMaxSlope(), maxslope );
	theProj->SetMaxSlope( maxslope );
      }
    } else {
      Error( Here(here), "Error calculating geometric maxslope for plane "
	     "type \"%s\". z-range of planes too small. Fix database.",
	     theProj->GetName() );
	return fStatus = kInitError;
    }
  }

  return fStatus = kOK;
}

//_____________________________________________________________________________
void MWDC::Clear( Option_t* opt )
{
  THaTrackingDetector::Clear(opt);
  
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->Clear(opt);

  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type )
    fProj[type]->Clear(opt);

  if( fNrefchan > 0 )
    fRefTime.assign( fNrefchan, kBig );  // not strictly necessary
}

//_____________________________________________________________________________
Int_t MWDC::Decode( const THaEvData& evdata )
{
  // Decode all planes and fill hitpatterns per projection
  
  static const char* const here = "Decode";

  // Decode reference channels of the VME readout (if defined)
  UInt_t iref = 0;
  for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fRefMap->GetModule(imod);	

    // Since we expect all channels to have a hit, we can efficiently loop
    // over the defined channel range
    for( Int_t chan = d->lo; chan <= d->hi; ++chan ) {
      Int_t data;
      Int_t nhits = evdata.GetNumHits( d->crate, d->slot, chan );
      if( nhits > 0 ) {
	data = evdata.GetData( d->crate, d->slot, chan, nhits-1 );
	if( nhits > 1 ) {
	  Warning( Here(here), "%d hits on reference channel %d module %d", 
		   nhits, chan, imod );
	}
	fRefTime[iref] = d->resolution * data;
      } else {
	// TODO: At this point, one could look for backup reference channels
	// to recover the data. Left for later, if needed.
	Warning( Here(here), "No hits on reference channel %d module %d.",
		 chan, imod );
	fRefTime[iref] = kBig;
      }
      ++iref;
    } // chans
  }   // modules


  // Decode the individual planes
  //TODO: multithread this
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->Decode( evdata );

  //TODO: check for excessive plane occupancy here and abort if too full
  // (via "max occupancy" parameter or similar)

  //TODO: for this to work, we need the drift distances first!
//   for( EProjType type = kTypeBegin; type < kTypeEnd; ++type )
//     fProj[type]->FillHitpattern();

  return 0;
}

//_____________________________________________________________________________
Int_t MWDC::CoarseTrack( TClonesArray& tracks )
{
  return 0;//fNtracks;
}

//_____________________________________________________________________________
Int_t MWDC::FineTrack( TClonesArray& tracks )
{
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
void MWDC::Print(const Option_t* opt) const
{
  THaTrackingDetector::Print(opt);

  //  fBench->Print();

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->Print();

  return;
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
void MWDC::EnableBenchmarks( Bool_t b )
{
  //  fDoBench = b;
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
Int_t MWDC::ReadDatabase( const TDatime& date )
{
  // Read MWDC database

  static const char* const here = "ReadDatabase";

  fIsInit = kFALSE;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Int_t status = kInitError;
  // Default plane angles and names
  // The angles can be overridden via the database
  Double_t p_angle[]   = { -60.0, 60.0, 90.0 };
  const char* p_name[] = { "u", "v", "x" };

  vector<string> planes;
  // Putting this container on the stack may cause a stack overflow!
  vector<int>* refmap = new vector<int>;

  string planeconfig;
  DBRequest request[] = {
    { "planeconfig",  &planeconfig, kString },
    { "refmap",       refmap,       kIntV,    0, 1 },
    { 0 }
  };

  Int_t err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    goto err;
  
  planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Fix database." );
    goto err;
  }

  // Delete existing configuration if re-initializing
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    // Clear previous projection objects
    delete fProj[type];
    fProj[type] = NULL;
  }
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    delete fPlanes[iplane];
  }
  fPlanes.clear();
  fRefTime.clear();

  // Fill the reference channel detector map, if given
  if( refmap->size() > 0 ) {
    UInt_t flags = THaDetMap::kFillModel;
    if( fRefMap->Fill( *refmap, flags ) <= 0 ) {
      Error( Here(here), "Invalid reference channel map. Fix database." );
      goto err;
    }
  }
  // Consistency checks of the reference channel map
  fNrefchan = fRefMap->GetTotNumChan();
  if( fNrefchan > 0 ) {
    Int_t rmin, rmax;
    fRefMap->GetMinMaxChan( rmin, rmax );
    if( rmin < 0 || rmax != static_cast<Int_t>(fNrefchan) ) {
      Error( Here(here), "Reference channel number inconsistency. "
	     "min/max/tot = %d %d %u.  Fix database.", rmin, rmax, fNrefchan );
      goto err;
    }
    // All ok - grow and fill the refchan data array with kBigs
    fRefTime.assign( fNrefchan, kBig );
  } else if( fRefMap->GetSize() > 0 ) {
    Error( Here(here), "Total number of reference channels = 0? "
	   "Fix database." );
    goto err;
  }
  // Set default resolution of channels if not already given
  for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fRefMap->GetModule(imod);
    if( THaDetMap::IsTDC(d) && d->resolution < 0.0 )
      d->resolution = 1e-9;
  }

  // Set up the definitions of the wire directions (track projections)
  for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
    Projection* proj = new Projection( type, 
				       p_name[type], 
				       p_angle[type]*TMath::DegToRad(),
				       this 
				       );
    if( !proj || proj->IsZombie() ) {
      // Urgh. Something is very bad
      Error( Here(here), "Error creating projection %s. Call expert.", 
	     p_name[type] );
      goto err;
    }
    fProj[type] = proj;
  }

  // Set up the wire planes
  for( ssiz_t i=0; i<planes.size(); ++i ) {
    TString name(planes[i].c_str());
    if( name.IsNull() )
      continue;
    vwiter_t it = find_if( fPlanes.begin(), fPlanes.end(),
			   WirePlane::NameEquals( name ) );
    if( it != fPlanes.end() ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", 
	     name.Data() );
      goto err;
    }
    WirePlane* newplane = new WirePlane( name, name, this );
    if( !newplane || newplane->IsZombie() ) {
      // Urgh. Something is very bad
      Error( Here(here), "Error creating wire plane %s. Call expert.", 
	     name.Data() );
      goto err;
    }
    fPlanes.push_back( newplane );
  }

  if( fDebug > 0 )
    Info( Here(here), "Loaded %u planes", fPlanes.size() );

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    WirePlane* thePlane = fPlanes[iplane];
    TString name( thePlane->GetName() );
    if( name.EndsWith("p") ) {
      TString other = name.Chop();
      if( other.IsNull() )
	continue;
      vwiter_t it = find_if( fPlanes.begin(), fPlanes.end(),
			     WirePlane::NameEquals( other ) );
      if( it != fPlanes.end() ) {
	WirePlane* partner = *it;
	if( fDebug > 0 )
	  Info( Here(here), "Partnering plane %s with %s",
		thePlane->GetName(), partner->GetName() );
	partner->SetPartner( thePlane );
      }
    }
  }

  status = kOK;
  fIsInit = kTRUE;
 err:
  delete refmap;
  return status;
}

//_____________________________________________________________________________
// Int_t MWDC::End(THaRunBase* run)
// {
//     //    fBench->Print();
//   return THaTrackingDetector::End(run);

// }



//_____________________________________________________________________________

}  

ClassImp(TreeSearch::MWDC)

///////////////////////////////////////////////////////////////////////////////
