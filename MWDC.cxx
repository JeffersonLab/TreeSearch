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

#include <algorithm>

// FIXME: Decoding and pattern finding in the planes should be multi-threaded

//TODO:  Decode, sort hits
//TODO:  Fill hitpatterns in CoarseTrack

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

  fProj.resize( Projection::kTypeEnd, NULL );

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

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ )
    delete fPlanes[iplane];
  for( vpsiz_t type = 0; type < fProj.size(); type++ )
    delete fProj[type];

}

//_____________________________________________________________________________
THaAnalysisObject::EStatus MWDC::Init( const TDatime& date )
{
  // Initialize MWDC. Calls standard Init(), then initializes subdetectors.

  // Initialize ourselves. This calls our ReadDatabase().
  EStatus status = THaTrackingDetector::Init(date);
  if( status )
    return fStatus = status;

  // Next, initialize each plane in turn
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ ) {
    status = fPlanes[iplane]->Init(date);
    if( status )
      return fStatus = status;
  }

  // Sort planes by increasing z-position. The sort order is defined 
  // by WirePlane::Compare().
  sort( fPlanes.begin(), fPlanes.end(), WirePlane::ZIsLess() );

  // Determine per-projection plane parameters
  for( Projection::EProjType type = Projection::kTypeBegin; 
       type < Projection::kTypeEnd; ++type ) {
    Projection* theProj = fProj[type];
    UInt_t n = 0;
    //    Double_t max_width = 0.0;
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ ) {
      WirePlane* thePlane = fPlanes[iplane];
      if( thePlane->GetType() == type ) {
	// If we have a plane of given type, make sure we have a Projection 
	// object for it
	if( !theProj )
	  fProj[type] = theProj = new Projection(type);
	// Add only primary planes (i.e. the first one of a partnered pair)
	// to the list and count them
	if( !thePlane->GetProjection() ) {
	  theProj->GetListOfPlanes().push_back(thePlane);
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

	//TODO: determine the maximum physical width of this projection 
	//if( max_width );
      }
    }
  }

  //TODO: Make sure partner planes have consistent parameters 
  //  (wire angle, maxslope etc.)
  //	thePlane->
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ ) {
    WirePlane* thePlane = fPlanes[iplane];
    TString name( thePlane->GetName() );
    if( name.EndsWith("p") ) {
      WirePlane* partner = thePlane->GetPartner();
      if( !partner ) continue;
    }
  }

  // Initialize hitpatterns for each projection plane etc.
  for( vpsiz_t iproj = 0; iproj < fProj.size(); iproj++ ) {
    if( fProj[iproj] )
      fProj[iproj]->Init();
  }



  return fStatus = kOK;
}

//_____________________________________________________________________________
void MWDC::Clear( Option_t* opt )
{
  THaTrackingDetector::Clear(opt);
  
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ )
    fPlanes[iplane]->Clear(opt);

}

//_____________________________________________________________________________
Int_t MWDC::Decode( const THaEvData& evdata )
{
  // Decode all planes

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ )
    fPlanes[iplane]->Decode( evdata );
  
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

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ )
    fPlanes[iplane]->Print();

  return;
}


//_____________________________________________________________________________
void MWDC::SetDebug( Int_t level )
{
  // Set debug level of this detector, including all wire planes (subdetectors)

  THaTrackingDetector::SetDebug( level );

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ )
    fPlanes[iplane]->SetDebug( level );
}

//_____________________________________________________________________________
void MWDC::EnableBenchmarks( Bool_t b )
{
  //  fDoBench = b;
}

//_____________________________________________________________________________
Projection::EProjType MWDC::MakePlaneType( const char* name,
					   Double_t& wire_angle )
{
  // Convert plane type name to a wire angle based on the uangle/vangle
  // database info.
  // Only the first character of name is used and must be one of "UVXY" in
  // uppercase or lowercase.
  // Returns the wire angle in degrees (90 = X, 0 = Y, uangle =, vangle = v).
  // By normal conventions, uangle = -vangle < 0, but the two angles might
  // in principle be different.
  // The wire angle is measured with respect to the X-axis, which for a 
  // magnetic spectrometer is along the dispersive direction in the direction
  // of increasing particle momentum.
  // Slight misalignments of individual chambers need to be handled in
  // the calling routine (ReadDatabase).
  // Also sets type to the plane type index defined in the Projection class.

  
  wire_angle = 0.0;
  if( !name || !*name )
    return Projection::kUndefinedType;

  TString planeTypes("UVXY");
  Ssiz_t pos = planeTypes.Index(name,1,0,TString::kIgnoreCase);
  if( pos == kNPOS )
    return Projection::kUndefinedType;

  switch(pos) {
  case 0:
    wire_angle = fUangle;
    return Projection::kUPlane;
  case 1:
    wire_angle = fVangle;
    return Projection::kVPlane;
  case 2:
    wire_angle = 90.0;
    return Projection::kXPlane;
  case 3:
    wire_angle = 0.0;
    return Projection::kYPlane;
  default:
    break;
  }
  return Projection::kUndefinedType;
}

//_____________________________________________________________________________
static void NormalizeAngle( Double_t& angle )
{
  // Put angle into [-90,90]
  
  if( angle > 0.0 )
    angle -= 180.0*TMath::Floor(angle/180.0);
  else if( angle < 0.0 )
    angle -= 180.0*TMath::Ceil(angle/180.0);
  if( angle > 90.0 )
    angle -= 180.0;
  else if( angle < -90.0 )
    angle += 180.0;
}

//_____________________________________________________________________________
Int_t MWDC::ReadDatabase( const TDatime& date )
{
  // Read MWDC database

  static const char* const here = "ReadDatabase";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  fIsInit = kFALSE;

  string planeconfig;
  DBRequest request[] = {
    { "planeconfig",    &planeconfig, kString },
    { "uangle",         &fUangle },
    { "vangle",         &fVangle },
    { 0 }
  };

  Int_t err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return kInitError;
  
  NormalizeAngle( fUangle );
  NormalizeAngle( fVangle );

  // Sanity checks of U and V angles
  if( (fUangle < 0.0 && fVangle < 0.0) || (fUangle > 0.0 && fVangle > 0.0) ) {
    Error( Here(here), "Plane misconfiguration: uangle (%6.2lf) and "
	   "vangle (%6.2lf) have the same sign. Fix database",
	   fUangle, fVangle );
    return kInitError;
  }
  if( TMath::Abs(TMath::Abs(fUangle)-45.0) > 44.0 ||
      TMath::Abs(TMath::Abs(fVangle)-45.0) > 44.0 ) {
    Error( Here(here), "uangle (%6.2lf) and vangle (%6.2lf) must be between "
	   "1-89 degrees. Fix database.", fUangle, fVangle );
    return kInitError;
  }

  // Set up the wire planes
  vector<string> planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Fix database." );
    return kInitError;
  }
  // Delete existing configuration if re-initializing
  for( vpsiz_t iproj = 0; iproj < fProj.size(); iproj++ ) {
    // Clear previous projection parameter blocks
    delete fProj[iproj];
    fProj[iproj] = NULL;
  }
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ ) {
    delete fPlanes[iplane];
    fPlanes.clear();
  }

  for( ssiz_t i=0; i<planes.size(); i++ ) {
    TString name(planes[i].c_str());
    if( name.IsNull() )
      continue;
    vwiter_t it = find_if( fPlanes.begin(), fPlanes.end(),
			   WirePlane::NameEquals( name ) );
    if( it != fPlanes.end() ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", 
	     name.Data() );
      return kInitError;
    }
    WirePlane* newplane = new WirePlane( name, name, this );
    fPlanes.push_back( newplane );
  }

  if( fDebug > 0 )
    Info( Here(here), "Loaded %u planes", fPlanes.size() );

  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); iplane++ ) {
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
  
  

  fIsInit = kTRUE;
  return kOK;
}

//_____________________________________________________________________________
Int_t MWDC::End(THaRunBase* run)
{
    //    fBench->Print();
  return THaTrackingDetector::End(run);

}



//_____________________________________________________________________________

}  

ClassImp(TreeSearch::MWDC)

///////////////////////////////////////////////////////////////////////////////
