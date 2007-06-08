//*-- Author :    Ole Hansen, Jefferson Lab   06-Jun-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaTsMWDC                                                                 //
//                                                                           //
// Reconstruction class for horizontal drift chambers used in BigBite.       //
// Tracking is done using a tree search algorithm (DellOrso et al.,          //
// NIM A 287, 436 (1990)), in essence a recursive template matching method.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaTsMWDC.h"
//#include "TreeSearch.h"

using namespace std;
//using namespace TreeSearch;

//_____________________________________________________________________________
THaTsMWDC::THaTsMWDC( const char* name, const char* description,
		      THaApparatus* apparatus ) :
  THaTrackingDetector(name,description,apparatus), fNtracks(0),
  fDoBench(kTRUE)
{ 
  // Constructor

#if 0
  // since the number of planes needs to be read in from
  // the database, the detector is actually constructed in the 
  // readdatabase method

  // Since TCloneArrays can resize, the size here is fairly unimportant
  fPlanes = new TClonesArray("THaTsMWDCPlane",0);

  // Default behavior for now
  // SetBit( kHardTDCcut );
  //  ResetBit( kIgnoreNegDrift );
  SetBit( kIgnoreNegDrift );
 

  fBench = new THaBenchmark;
}

//_____________________________________________________________________________
THaTsMWDC::~THaTsMWDC()
{
  // Destructor. Delete subdetectors and unregister variables
#if 0 
 if (fIsSetup)
    RemoveVariables();
  
  delete fBench;
  delete fPlanes;
#endif
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus THaTsMWDC::Init( const TDatime& date )
{
  // Initialize MWDC. Calls standard Init(), then initializes subdetectors.

#if 0
  Int_t useshowercuts;

  fReconstructBits = 0;

  if( IsZombie() ) {
    return fStatus = kInitError;
  }
  EStatus status;
  status = (THaTrackingDetector::Init( date )) ;
  if (status!=0) return fStatus = status;

  for (UInt_t i=0; i<GetNPlanes(); i++) {
    status=GetPlane(i)->Init( date );
    if (status!=0)return fStatus = status;
    fReconstructBits |= (BIT(i)*GetPlane(i)->IsInReconstruction() );
  }

  // Need to pull quite a bit from THaBBe's database:
  THaBBe* bbe = dynamic_cast<THaBBe*>(GetApparatus());
  if( bbe ) {
    fMaxGroups = (UInt_t) bbe->GetMaxGroups();
    fHardMaxGroups =  (UInt_t) bbe->GetHardMaxGroups();
    fMinimumPlanes = (UInt_t) bbe->GetMinimumPlanesInRecon();
    fShowerXExt = bbe->GetShowerXExtension();
    fShowerYExt = bbe->GetShowerYExtension();
    fTargetXExt = bbe->GetTargetXExtension();
    fTargetYExt = bbe->GetTargetYExtension();
    fTargetWindowXOffset = bbe->GetTargetWindowXOffset();
    useshowercuts = bbe->UseShowerCuts();
    fMaxCallThreshold = bbe->GetMaxCallThreshold();
  }
  else {
    fMaxGroups = 1000000;
    fHardMaxGroups = 1000000;
    fMinimumPlanes =4;
    fShowerXExt = 0.0;
    fShowerYExt = 0.0;
    fTargetXExt = 0.0;
    fTargetYExt = 0.0;
    fTargetWindowXOffset = 0.0; 
    useshowercuts = false;
    fMaxCallThreshold = 1000000;
  }

  SetUseShowerCuts(useshowercuts);

  // Generate transformation matrices
  GenerateConstructMatrices();
  cout << "Requiring at least " << fMinimumPlanes << " in reconstruction" << endl;
  cout << "Only using planes ";
  for( UInt_t i = 0; i < GetNPlanes(); i++ )
    {
      if( GetPlane(i)->IsInReconstruction() )
	cout << GetPlane(i)->GetName() << " ";
    }
  cout << "in reconstruction" << endl;

  cout << "Skipping tracking on events that will have more than " << fMaxCallThreshold 
       << " estimated calls to consider."<< endl;


  // If using shower cuts print out some useful information
  if( UseShowerCuts() )
    {      
      cout << endl << "Shower cut specifications: " << endl 
	   << "-------------------------------------------" << endl;
      cout << "Shower X Range (+/- m):\t\t"   << fShowerXExt << endl;
      cout << "Shower Y Range (+/- m):\t\t"   << fShowerYExt << endl;
      cout << "Target X Extension (+/- m):\t" << fTargetXExt << endl;
      cout << "Target Y Extension (+/- m):\t" << fTargetYExt << endl;

      cout << "Target X Offset (m):\t\t"<< fTargetWindowXOffset << endl;
      cout << "-------------------------------------------" << endl;
    }
  else
    {
      cout << endl << "*Not using cuts on shower*" << endl;
    }

  // Load up target information from THaBBe
  LoadTargetData();

#endif
  return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t THaTsMWDC::Decode( const THaEvData& evdata )
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
Int_t THaTsMWDC::CoarseTrack( TClonesArray& tracks )
{
  return 0;//fNtracks;
}

//_____________________________________________________________________________
Int_t THaTsMWDC::FineTrack( TClonesArray& tracks )
{
  return 0;//fNtracks;
}

//_____________________________________________________________________________
Int_t THaTsMWDC::DefineVariables( EMode mode )
{
  // Initialize global variables and lookup table for decoder


  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list
  
#if 0
  RVarDef vars[] = {
    { "trackskipped", "Tracking skipped or truncated", "fTooBusy"       },
    { "cproctime",    "Coarse Processing Time",        "fCoarseProcTime"},
    { "fproctime",    "Fine Processing Time",          "fFineProcTime"  },
    { "estngrp",      "Estimated Number of Groups",    "fEstNGroups"    },
    { "estncall",     "Estimated Number of Calls",     "fEstNCalls"     },
    { "ngrp",         "Number of Groups",              "fNGroups"       },
    { "ncall",        "Number of recursive calls",     "fNCalls"        },
    { 0 }
  };
  DefineVarsFromList( vars, mode );
#endif  
  return 0;
}

//_____________________________________________________________________________
void THaTsMWDC::Print(const Option_t* opt) const
{
  THaTrackingDetector::Print(opt);

  //  fBench->Print();

  return;
}


//_____________________________________________________________________________
void THaTsMWDC::SetDebug( Int_t level )
{
  // Set debug level of this detector, including all wire planes (subdetectors)

  THaTrackingDetector::SetDebug( level );

#if 0
  for( UInt_t i = 0; i<GetNPlanes(); i++ ) {
    THaTsMWDCPlane* thePlane = static_cast<THaTsMWDCPlane*>( fPlanes->At(i) );
    thePlane->SetDebug( level );
  }
#endif
}

//_____________________________________________________________________________
void THaTsMWDC::EnableBenchmarks( Bool_t b )
{
  fDoBench = b;
}

//_____________________________________________________________________________
Int_t THaTsMWDC::ReadDatabase( const TDatime& date )
{
  // Read MWDC database

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
	  THaTsMWDCPlane* thePlane = 
	    new((*fPlanes)[nPlanes]) THaTsMWDCPlane(planename,planedescr,this);
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
void THaTsMWDC::Clear( Option_t* opt )
{
  THaTrackingDetector::Clear(opt);
  

}

//_______________________________________________________________________________
Int_t THaTsMWDC::End(THaRunBase *run)
{
    //    fBench->Print();
    return THaTrackingDetector::End(run);
}



//_____________________________________________________________________________


ClassImp(THaTsMWDC)
  
////////////////////////////////////////////////////////////////////////////////
