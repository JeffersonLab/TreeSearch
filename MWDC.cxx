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
//#include "TreeSearch.h"
#include "TList.h"
#include "WirePlane.h"

//using namespace std;
namespace TreeSearch {

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* description,
	    THaApparatus* apparatus ) :
  THaTrackingDetector(name,description,apparatus), fPlanes(NULL)
{ 
  // Constructor

  fPlanes = new TList;

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
  return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t MWDC::Decode( const THaEvData& evdata )
{

#if 0
  // Decode all planes
  for (UInt_t i=0; i<GetNPlanes() ;i++) {
    GetPlane(i)->Decode(evdata);
  }

#endif
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

#if 0
  for( UInt_t i = 0; i<GetNPlanes(); i++ ) {
    WirePlane* thePlane = static_cast<WirePlane*>( fPlanes->At(i) );
    thePlane->SetDebug( level );
  }
#endif
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


  //FIXME: quick test
  fPlanes->Add( new WirePlane("u1", "U1 Plane", this));
  fPlanes->Add( new WirePlane("u1p","U1p Plane",this));
  fPlanes->Add( new WirePlane("u2", "U2 Plane", this));
  fPlanes->Add( new WirePlane("u3", "U3 Plane", this));
  fPlanes->Add( new WirePlane("u3p","U3p Plane",this));

  fPlanes->Add( new WirePlane("x1", "X1 Plane", this));
  fPlanes->Add( new WirePlane("x1p","X1p Plane",this));
  fPlanes->Add( new WirePlane("x2", "X2 Plane", this));
  fPlanes->Add( new WirePlane("x3", "X3 Plane", this));
  fPlanes->Add( new WirePlane("x3p","X3p Plane",this));

  fPlanes->Add( new WirePlane("v1", "V1 Plane", this));
  fPlanes->Add( new WirePlane("v1p","V1p Plane",this));
  fPlanes->Add( new WirePlane("v2", "V2 Plane", this));
  fPlanes->Add( new WirePlane("v3", "V3 Plane", this));
  fPlanes->Add( new WirePlane("v3p","V3p Plane",this));




#if 0
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // load global MWDC parameters
  static const char* const here = "ReadDatabase";
  const int LEN = 200;
  char buff[LEN];
  
  // Look for the section [<prefix>.global] in the file, e.g. [ R.global ]
  TString tag1(fPrefix);
  TString tag2(fPrefix);

  Ssiz_t pos = tag1.Index("."); 
  if( pos != kNPOS )
    tag1 = tag1(0,pos+1);
  else
    tag1.Append(".");
  tag1.Prepend("[");
  tag1.Append("global]"); 
  pos = tag2.Index("."); 
  if( pos != kNPOS )
    tag2 = tag2(0,pos+1);
  else
    tag2.Append(".");
  tag2.Prepend("[");
  tag2.Append("global.done]"); 



  TString line, tag3(tag1);
  tag1.ToLower();
  tag2.ToLower();


  bool found = false;
  while (!found && fgets (buff, LEN, file) != NULL) {
    char* buf = ::Compress(buff);  //strip blanks
    line = buf;
    delete [] buf;
    if( line.EndsWith("\n") ) line.Chop();
    line.ToLower();
     if ( tag1 == line ) 
      found = true;
  }
  if( !found ) {
    Error(Here(here), "Database entry %s not found!", tag3.Data() );
    fclose(file);
    return kInitError;
  }

  if( GetNPlanes()>0 || fIsInit ) {
    Warning(Here(here), "Database has already be read in.  Using the planes we already have." );
    Warning(Here(here), "Hope you have the same configureation in the databases." );
    //fclose(file);
    //return kInitError;
  }

  // We found the section, now read the data

  // read in some basic constants first
  else
    {
      TString planename = " ";
      TString planedescr = " ";
      Int_t nPlanes = 0;  
      while ((planename!=tag2)&&(fgets(buff, LEN, file)!=NULL)) {
	char* buf = ::Compress(buff);
	planename = buf;
	delete [] buf;
	if( planename.EndsWith("\n") ) planename.Chop();
	planename.ToLower();
	if (planename!=tag2) {
	  fgets(buff, LEN, file);
	  planedescr = buff;
	  if( planedescr.EndsWith("\n") ) planedescr.Chop();
	  planedescr.ToLower();
	  WirePlane* thePlane = 
	    new((*fPlanes)[nPlanes]) WirePlane(planename,planedescr,this);
	  thePlane->SetDebug(fDebug);

	  //	  EStatus status=thePlane->Init( date );
	  //	  if (status!=0)return fStatus = status;
	  
	  nPlanes++;
	}
      }
    }

  fclose(file);
  // final sanity check
  if (MAX_PLANES<GetNPlanes()) {
    Error(Here(here),"More planes in database than permitted by MAX_PLANES!");
    return kInitError;
  }
  
#endif

  fIsInit = true;
  return kOK;
}

//_____________________________________________________________________________
void MWDC::Clear( Option_t* opt )
{
  THaTrackingDetector::Clear(opt);
  

}

//_____________________________________________________________________________
Int_t MWDC::End(THaRunBase *run)
{
    //    fBench->Print();
    return THaTrackingDetector::End(run);
}



//_____________________________________________________________________________

}  

ClassImp(TreeSearch::MWDC)

///////////////////////////////////////////////////////////////////////////////
