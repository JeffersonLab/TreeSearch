//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TreeSearch::WirePlane                                                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "WirePlane.h"
#include "Hit.h"
#include "THaDetMap.h"

using std::vector;

namespace TreeSearch {

//_____________________________________________________________________________
WirePlane::WirePlane( const char* name, const char* description,
		      THaDetectorBase* parent )
  : THaSubDetector(name,description,parent), fType(kUndefinedType),
    fWireStart(0.0), fWireSpacing(0.0), fSinAngle(0.0), fCosAngle(0.0),
    fTDCRes(0.0), fDriftVel(0.0), fResolution(0.0), fMaxSlope(0.0),
    fTTDConv(NULL)
{
  // Constructor

  fHits = new TClonesArray("TreeSearch::Hit", 50);
}

//_____________________________________________________________________________
WirePlane::~WirePlane()
{
  // Destructor.

  if( fIsSetup )
    RemoveVariables();

  delete fHits;
}

//_____________________________________________________________________________
Int_t WirePlane::ReadDatabase( const TDatime& date )
{


  /* Read or determine here:
     - Plane type
     -# wires (?)
     - wire start pos
     - wire spacing
     - wire angle
     - TDC resolution
     - Drift velocity (?)
     - estimated position resolution
     - ( max slope (maybe more parameters?) - FIXME: calculate this?)
     - ( array of TDC offsets )
     - ( parameters of TTD converter (opt) )
  */

  static const char* const here = "ReadDatabase";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Int_t err = ReadGeometry( file, date, kTRUE );
  if( err )
    return err;

  Double_t wire_angle = 0.0;
  vector<int> detmap;

  DBRequest request[] = {
    { "detmap",     &detmap, kIntV },
    { "nwires",     &fNelem, kInt },
    { "wire.pos",       &fWireStart },
    { "wire.spacing",   &fWireSpacing },
    { "wire.angle",     &wire_angle },
    { "tdc.res",   &fTDCRes },
    { "drift.v",   &fDriftVel },
    //    { "resolution", &fResolution },
    //    { "maxslope",   &fMaxSlope },
    { "tdc.offsets",     &fTDCOffset, kDoubleV },
    { 0 }
  };

  err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return kInitError;

  UInt_t flags = THaDetMap::kFillModel|THaDetMap::kAutoCount;
  if( FillDetMap( detmap, flags, here ) <= 0 )
    return kInitError;

  if( fNelem <= 0 ) {
    Error( Here(here), "Invalid number of wires: %d", fNelem );
    return kInitError;
  }
  Int_t nchan = fDetMap->GetTotNumChan();
  if( nchan != fNelem ) {
    Error( Here(here), "Number of detector map channels (%d) "
	   "disagrees with number of wires (%d)", nchan, fNelem );
    return kInitError;
  }
  nchan = fTDCOffset.size();
  if( nchan != fNelem ) {
    Error( Here(here), "Number of TDC offset values (%d) "
	   "disagrees with number of wires (%d)", nchan, fNelem );
    return kInitError;
  }
  wire_angle *= TMath::DegToRad();
  fSinAngle = TMath::Sin(wire_angle);
  fCosAngle = TMath::Cos(wire_angle);

  //FIXME: determine fType

  fIsInit = true;
  return kOK;
}

//_____________________________________________________________________________
Int_t WirePlane::DefineVariables( EMode mode )
{
  // initialize global variables


  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { 0 }
  };
  return DefineVarsFromList( vars, mode );

}

//_____________________________________________________________________________
THaAnalysisObject::EStatus WirePlane::Init( const TDatime& date )
{
  // Calls its own Init(), then initializes subdetectors, then calculates
  // some local geometry data.

  EStatus status = THaSubDetector::Init( date );
  if( status ) 
    return fStatus = status;

  //FIXME:
  THaDetectorBase *parent = GetDetector();
  if( parent )
    fOrigin += parent->GetOrigin();

  return fStatus = kOK;
}



//_____________________________________________________________________________
inline
void WirePlane::Clear( Option_t* opt )
{    

}

//_____________________________________________________________________________
Int_t WirePlane::Decode( const THaEvData& evData)
{    
  // Converts the raw data into hit information
  // Assumes channels & wires  are numbered in order
  // Assumes that the first channels are for reference channels
  // TODO: Make sure the wires are numbered in order, even if the channels
  //       aren't
              
#if 0
  //cout << "Decoding plane " << GetName() << endl;

  if (!evData.IsPhysicsTrigger()) return -1;
  
  Int_t nextHit = 0;

  bool no_negative;
  if( fMWDC ) {
    // If true, ignore negativ drift times completely
    no_negative      = fMWDC->TestBit(THaTsMWDC::kIgnoreNegDrift);
  } else {
    no_negative = false;
  }

  // Loop over all detector modules containing reference lines
  // they have to be at the beginning of the detector map
  // and we need one line (=one module) per reference channel
  // otherwise only the first will be used
  Int_t i=0;
  Int_t data;
  Bool_t valid;

  while ((  i < GetNRefCh() )  &&  (i < fDetMap->GetSize()) ) {
    THaDetMap::Module * d = fDetMap->GetModule(i);	
    
    // Get number of channels with hits
    valid = true;
    Int_t chan=d->lo;
    Int_t nHits = evData.GetNumHits(d->crate, d->slot, chan);
    if (nHits<1) {
      //      cout<<"Warning: Event Number "<<dec<<evData.GetEvNum()<<"; No Hit for Reference Channel "<<i<<" of detector"<<fPrefix<<endl;
      data = 0;
      valid = false;
    } else {
      if (nHits>1) {
	//cout<<"Warning: Event Number "<<evData.GetEvNum()<<"; Multiple Hits for Reference Channel "<<i<<" of detector"<<fPrefix<<endl;
	//cout<<"Using first one"<<endl;
      }
      data = evData.GetData(d->crate, d->slot, chan, nHits-1);
    }
    // Wire numbers and channels go in the same order ... 
    THaVDCWire* wire = GetRefCh(i);
    if( !wire ) { 
      cout<<"WirePlane::Decode : Ref Channels are not initialized"<<endl;
      cout<<"Skipping event"<<endl;
      return -1;
    }
    
    // FIXME:  work in resolution for reference channel which is different from TDC
    Double_t time = -0.05 * data;
    new( (*fRefHits)[i] )  THaTsMWDCHit(this, wire, data, time, 0.0, 0.0, valid );    
    i++;
  }
  if (i!=GetNRefCh()) {
    cout<<"WirePlane::Decode : Mismatch between fNRefCh and hits on RefLines"<<endl;
    cout<<i<<" "<<GetNRefCh()<<endl;
    return -1;
  }

  //  cout << "Scanning through hits" << endl;

  while (i < fDetMap->GetSize()){
    //  cout << "Module index " << i << " of " << fDetMap->GetSize() << endl;
    THaDetMap::Module * d = fDetMap->GetModule(i);
    
    // Get number of channels with hits
    Int_t nChan = evData.GetNumChan(d->crate, d->slot);
    //cout << "Number of channels here: " << nChan << endl;

    for (Int_t chNdx = 0; chNdx < nChan; chNdx++) {
      // Use channel index to loop through channels that have hits

      Int_t chan = evData.GetNextChan(d->crate, d->slot, chNdx);
      if (chan < d->lo || chan > d->hi) 
	continue; //Not part of this detector

      // Wire numbers and channels go in the same order ... 
      Int_t wireNum  = d->first + chan - d->lo;

      THaVDCWire* wire = GetWire(wireNum-GetNRefCh());
      if( !wire || wire->GetFlag() != 0 ) { continue;}
      // Get number of hits for this channel and loop through hits
      Int_t nHits = evData.GetNumHits(d->crate, d->slot, chan);
      for (Int_t hit = 0; hit < nHits; hit++) {
	// Now get the TDC data for this hit
	Int_t data = evData.GetData(d->crate, d->slot, chan, hit);   	
	Double_t refoffset = 0.0;
	if ((d->refindex)>=0) {
	  THaTsMWDCHit* ahit = GetRefHit(d->refindex);
	  if ((!ahit)||(!(ahit->IsDataValid()))) {
	    fRefOkay = false ;
	    //	    cout<<"Warning: Event Number "<<dec<<evData.GetEvNum()<<"; No Hit for Reference Channel of detector"<<fPrefix<<endl;
	    //cout<<"Warning WirePlane::Decode : no hit on RefLine"<<endl;
	  } else {
	    refoffset = ahit->GetTime();
	  }
	}

	Double_t time = 0;

	time = fTDCRes * data - wire->GetTOffset() - refoffset;
	fTotalHits++;
	if ((!no_negative)||((time>=fMinTime)&&(time<=fMaxTime))) {
	  new( (*fHits)[nextHit++] )  THaTsMWDCHit(this, wire, data, time,  wire->GetTOffset(), refoffset );
	}
	else
	  {
	    fRejectedHits++;
	  }
      } // End hit loop
    } // End channel index loop
    i++;
  }
   
  // Sort the hits in order of increasing wire number and (for the same wire
  // number) increasing time (NOT rawtime)

  fHits->Compress();
  fHits->Sort();
  
  // cbr for debugging  fill the per event data 
  
  Int_t weirno=0;
  for (i=0;i< (Int_t) GetNHits();i++) { 
    weirno=((THaTsMWDCHit*)fHits->At(i))->GetWireNum();
    fTcounter[weirno]++;
    if (fTcounter[weirno]==1) {
      fRawT[weirno]=((THaTsMWDCHit*)fHits->At(i))->GetRawTime();
      fCorT[weirno]=((THaTsMWDCHit*)fHits->At(i))->GetTime();
    }
  }
  for (i=0; i<GetNRefHits(); i++) {
    weirno=((THaTsMWDCHit*)fRefHits->At(i))->GetWireNum();
    if ((fRefT[weirno]<-1.e34)&&(((THaTsMWDCHit*)fRefHits->At(i))->IsDataValid()))
      { 
	fRefT[weirno]=((THaTsMWDCHit*)fRefHits->At(i))->GetRawTime();
      }
  }
  
#endif
  return 0;
  
}
  
//_____________________________________________________________________________

}

ClassImp(TreeSearch::WirePlane)

///////////////////////////////////////////////////////////////////////////////

