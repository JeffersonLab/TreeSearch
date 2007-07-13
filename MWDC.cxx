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
//#include "TreeSearch.h"
#include "TList.h"
#include "TString.h"

// FIXME: Decoding and pattern finding in the planes should be multi-threaded

using namespace std;
typedef string::size_type ssiz_t;

namespace TreeSearch {

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* description,
	    THaApparatus* apparatus ) :
  THaTrackingDetector(name,description,apparatus), fPlanes(NULL)
{ 
  // Constructor

  fPlanes = new TList;

  THaDetector::SetDebug(1);
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
  // Destructor. Delete subdetectors and unregister variables
 if (fIsSetup)
   RemoveVariables();
  
#if 0 
  delete fBench;
#endif

  fPlanes->Delete();
  delete fPlanes;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus MWDC::Init( const TDatime& date )
{
  // Initialize MWDC. Calls standard Init(), then initializes subdetectors.

  EStatus status = THaTrackingDetector::Init(date);
  if( status )
    return fStatus = status;

  TIter next(fPlanes);
  WirePlane* thePlane;
  while( (thePlane = static_cast<WirePlane*>(next())) ) {
    status = thePlane->Init(date);
    if( status )
      return fStatus = status;
  }

  // Sort planes by increasing z-position. The sort order is defined 
  // by WirePlane::Compare().
  fPlanes->Sort();

  //TODO: Maintain separate lists for planes of each orientation




  return fStatus = kOK;
}

//_____________________________________________________________________________
void MWDC::Clear( Option_t* opt )
{
  THaTrackingDetector::Clear(opt);
  
  TIter next(fPlanes);
  WirePlane* thePlane;
  while( (thePlane = static_cast<WirePlane*>(next())) )
    thePlane->Clear(opt);

}

//_____________________________________________________________________________
Int_t MWDC::Decode( const THaEvData& evdata )
{
  // Decode all planes

  TIter next(fPlanes);
  WirePlane* thePlane;
  while( (thePlane = static_cast<WirePlane*>(next())) )
    thePlane->Decode( evdata );
  
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

  fPlanes->Print();
  return;
}


//_____________________________________________________________________________
void MWDC::SetDebug( Int_t level )
{
  // Set debug level of this detector, including all wire planes (subdetectors)

  THaTrackingDetector::SetDebug( level );

  TIter next(fPlanes);
  WirePlane* thePlane;
  while( (thePlane = static_cast<WirePlane*>(next())) )
    thePlane->SetDebug( level );
}

//_____________________________________________________________________________
void MWDC::EnableBenchmarks( Bool_t b )
{
  //  fDoBench = b;
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
    { 0 }
  };

  Int_t err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return kInitError;

  // Set up the wire planes
  vector<string> planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Fix database." );
    return kInitError;
  }
  // Delete existing configuration if re-initializing
  fPlanes->Delete();
  for( ssiz_t i=0; i<planes.size(); i++ ) {
    TString name(planes[i].c_str());
    if( name.IsNull() )
      continue;
    if( fPlanes->FindObject(name) ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", 
	     name.Data() );
      return kInitError;
    }
    WirePlane* newplane = new WirePlane( name, name, this );
    fPlanes->Add( newplane );
  }

  if( fDebug > 0 )
    Info( Here(here), "Loaded %d planes", fPlanes->GetSize() );

  TIter next(fPlanes);
  WirePlane* thePlane;
  while( (thePlane = static_cast<WirePlane*>(next())) ) {
    TString name( thePlane->GetName() );
    if( name.EndsWith("p") ) {
      TString other = name.Chop();
      if( other.IsNull() )
	continue;
      WirePlane* partner = static_cast<WirePlane*>(fPlanes->FindObject(other));
      if( partner ) {
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
