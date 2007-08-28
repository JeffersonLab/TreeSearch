//*-- Author :    Ole Hansen, Jefferson Lab   14-Jun-2007

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TreeSearch::WirePlane                                                    //
//                                                                          //
//                                                                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "WirePlane.h"
#include "Hit.h"
#include "THaDetMap.h"
#include "THaEvData.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "MWDC.h"
#include <iostream>
#include <string>

using namespace std;

// Database uses ns for TDC offsets
static const float kTDCscale = 1e-9;

namespace TreeSearch {

//_____________________________________________________________________________
WirePlane::WirePlane( const char* name, const char* description,
		      THaDetectorBase* parent )
  : THaSubDetector(name,description,parent), fPlaneNum(-1),
    fType(kUndefinedType), fWireStart(0.0), fWireSpacing(0.0), 
    fPartner(NULL), fProjection(NULL), fMWDC(NULL), 
    fDriftVel(0.0), fResolution(0.0), fTTDConv(NULL),
    fNmiss(0), fNrej(0), fWasSorted(0), fNhitwires(0), fNnohits(0)
{
  // Constructor

  static const char* const here = "WirePlane";

  fHits = new TClonesArray("TreeSearch::Hit", 200);
  if( !fHits ) {
    Fatal( Here(here), "Allocating hit array failed. Call expert." );
    MakeZombie();
    return;
  }

  fMWDC = dynamic_cast<MWDC*>( GetDetector() );
  if( !fMWDC ) {
    Error( Here(here), "No parent detector (MWDC) defined. Call expert." );
    MakeZombie();
    return;
  }
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
  // Read database

  static const char* const here = "ReadDatabase";

  // TODO: the refmap might be better in the MWDC class since, presumably,
  // we only need to read those channels once

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Int_t err = ReadGeometry( file, date, kTRUE );
  if( err )
    return err;

  Int_t status = kInitError;

  TString plane_type;
  Int_t tdc_res = 1000;  // 1ns/ch default resolution
  // Putting this on the stack may cause a stack overflow
  vector<int>* detmap = new vector<int>;

  UInt_t flags, nrefchan;
  Int_t nchan;
  TString name;
  Int_t dmin, dmax;

  DBRequest request[] = {
    { "detmap",       detmap,        kIntV },
    { "nwires",       &fNelem,       kInt },
    { "type",         &plane_type,   kTString, 0, 1 },
    { "wire.pos",     &fWireStart },
    { "wire.spacing", &fWireSpacing, kDouble,  0, 0, -1 },
    //FIXME: remove?
    { "tdc.res",      &tdc_res,      kDouble,  0, 1, -1 },
    { "drift.v",      &fDriftVel,    kDouble,  0, 0, -1 },
    { "xp.res",       &fResolution,  kDouble,  0, 0, -1 },
    { "tdc.offsets",  &fTDCOffset,   kFloatV,  0, 0 },
    { "description",  &fTitle,       kTString, 0, 1 },
    { 0 }
  };

  err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    goto err;

  if( fNelem <= 0 ) {
    Error( Here(here), "Invalid number of wires: %d", fNelem );
    goto err;
  }

  // Fill the detector map for the data channels
  flags = ( THaDetMap::kFillModel | THaDetMap::kFillRefIndex );
  if( FillDetMap( *detmap, flags, here ) <= 0 )
    goto err;

  nchan = fDetMap->GetTotNumChan();
  if( nchan != fNelem ) {
    Error( Here(here), "Number of detector map channels (%d) "
	   "disagrees with number of wires (%d)", nchan, fNelem );
    goto err;
  }
  nchan = fTDCOffset.size();
  if( nchan > 0 && nchan != fNelem ) {
    Error( Here(here), "Number of TDC offset values (%d) "
	   "disagrees with number of wires (%d)", nchan, fNelem );
    goto err;
  }
  // Convert units of TDC offsets
  for( vector<float>::size_type i = 0; i < fTDCOffset.size(); ++i ) {
    fTDCOffset[i] *= kTDCscale;
  }

  // Check consistency of reference channels and data channels
  fDetMap->GetMinMaxChan( dmin, dmax, THaDetMap::kRefIndex );
  nrefchan = fMWDC->GetNrefchan();
  if( nrefchan > 0 ) {
    if( dmax >= static_cast<Int_t>(nrefchan) ) {
      Error( Here(here), "Reference channel out of range: requested = %d, "
	     "max = %u. Fix database.", dmax, nrefchan );
      goto err;
    }
    //TODO: check if reference channel crate/slot/resolution consistent
  } else if( dmax >= 0 ) {
    Error( Here(here), "detmap specifies refindex, but no refmap defined. "
	   "Fix database." );
    goto err;
  }

  // Fill in default TDC resolution if no resolution given in module database
  // FIXME: remove?
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fDetMap->GetModule(imod);
    if( THaDetMap::IsTDC(d) && d->resolution < 0.0 )
      d->resolution = tdc_res * kTDCscale;
  }

  // Determine the type of this plane. If the optional plane type variable is
  // not in the database, use the first character of the plane name.
  name = plane_type.IsNull() ? fName[0] : plane_type[0];
  fType = fMWDC->NameToType( name );
  if( fType == kUndefinedType ) {
    Error( Here(here), "Illegal plane type (%s). Fix database.", name.Data() );
    goto err;
  }

  status = kOK;
  fIsInit = true;
 err:
  delete detmap;
  return status;
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

  if( IsZombie() )
    return fStatus = kNotinit;

  EStatus status = THaSubDetector::Init( date );
  if( status ) 
    return fStatus = status;

  //FIXME:
//   THaDetectorBase *parent = GetDetector();
//   if( parent )
//     fOrigin += parent->GetOrigin();

  return fStatus = kOK;
}



//_____________________________________________________________________________
void WirePlane::Clear( Option_t* opt )
{    
  // Clear event-by-event data (hits)

  fHits->Clear();
  fNmiss = fNrej = fWasSorted = fNhitwires = fNnohits = 0;
}

//_____________________________________________________________________________
Int_t WirePlane::Decode( const THaEvData& evData )
{    
  // Extract this plane's hit data from the raw evData.
  //
  // This routine can handle both the old Fastbus readout and the new CAEN
  // VME pipeline TDCs. The latter require a reference channel map and
  // cross-references to reference channels in the regular detector map
  // of the plane.

  //  static const char* const here = "Decode";

  UInt_t nHits = 0;
  bool positive_only = fMWDC->TestBit(MWDC::kIgnoreNegDrift);

  // Decode data. This is done fairly efficiently by looping over only the 
  // channels with hits on each module. 
  // FIXME: If a module is shared with another plane (common here), we waste
  // time skipping hits that don't belong to us.
  // NB: certain indices below are guaranteed to be in range by construction
  // in ReadDatabase, so we can avoid unneeded checks.
  bool sorted = true, revsorted = true;
  Hit* prevHit = NULL;
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module * d = fDetMap->GetModule(imod);
    Double_t ref_offset = 
      (d->refindex >= 0) ? fMWDC->GetRefTime(d->refindex) : 0.0;

    // Get number of channels with hits and loop over them, skipping channels
    // that are not part of this module
    // FIXME: this becomes very inefficient if several modules with the 
    // same crate/slot are defined - e.g. one "module" per channel...ouch
    Int_t nchan = evData.GetNumChan( d->crate, d->slot );
    for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
      Int_t chan = evData.GetNextChan( d->crate, d->slot, ichan );
      if( chan < d->lo || chan > d->hi ) {
	++fNmiss;
	continue; //Not part of this detector
      }
      // Get the wire number. Assumes that the detector map is defined in order
      // of ascending wire numbers.
      Int_t iw = d->first + chan - d->lo;
      Double_t tdc_offset = fTDCOffset[iw];

      // Get number of hits on this channel and loop over hits
      Int_t nhits = evData.GetNumHits( d->crate, d->slot, chan );
      if( nhits > 0 ) ++fNhitwires; else ++fNnohits;
      for( Int_t hit = 0; hit < nhits; hit++ ) {
	
	// Get the TDC data for this hit
	Int_t data = evData.GetData( d->crate, d->slot, chan, hit );
	
	// Convert the TDC value to the drift time. The readout uses 
	// common-stop TDCs, so drift_time = tdc_time(drift=0)-tdc_time(data).
	Double_t time = tdc_offset+ref_offset - d->resolution*(data+0.5);
	if( !positive_only || time > 0.0 ) {
	  Hit* theHit = 
	    new( (*fHits)[nHits++] ) Hit( iw, 
					  fWireStart + iw * fWireSpacing,
					  data,
					  time,
					  fResolution,
					  this
					  );
	  // We can test the order of the hits on the fly - they should
	  // come in sorted or reverse-sorted. If they are, we can avoid
	  // trying to sort an already-sorted array - which would be expensive
	  if( sorted && prevHit && theHit->Hit::Compare(prevHit) < 0 )
	    sorted = false;
	  if( revsorted && prevHit && theHit->Hit::Compare(prevHit) > 0 )
	    revsorted = false;
	  prevHit = theHit;
	}
	else
	  ++fNrej;

      } // hits
    }   // chans
  }     // modules
   
  // If ncessary, sort the hits
  fWasSorted = sorted ? 1 : revsorted ? -1 : 0;
  if( fWasSorted == 0 ) {
//   if( !sorted ) {
//     if( revsorted ) {
//       // reverse array ... urgh
//       TClonesArray* copy = new TClonesArray( fHits->GetClass(), 
// 					     fHits->GetSize() );
//       Int_t end = fHits->GetLast();
//       for( Int_t i = 0; i < end+1; ++i )
// 	new((*copy)[end-i]) Hit( *static_cast<Hit*>((*fHits)[i]) );
//       delete fHits;
//       fHits = copy;
//     } else {
    fHits->Sort();
  }

  return nHits;
}
  
//_____________________________________________________________________________
void WirePlane::SetPartner( WirePlane* p )
{    
  // Partner this plane with plane 'p'. Partner planes are expected to
  // be located close to each other and usually to have staggered wires.

  fPartner = p;
  if( p )
    p->fPartner = this;

  return;
}

//_____________________________________________________________________________
void WirePlane::Print( Option_t* opt ) const
{    
  // Print plane info

  cout << "WirePlane:  #" << GetPlaneNum() << " "
       << GetName()   << "\t"
    //       << GetTitle()        << "\t"
       << fNelem << " wires\t"
       << "z = " << GetZ();
  if( fPartner ) {
    cout << "\t partner = " 
	 << fPartner->GetName();
    //	 << fPartner->GetTitle();
  }
  cout << endl;
}

//_____________________________________________________________________________
Int_t WirePlane::Compare( const TObject* obj ) const 
{
  // Used to sort planes in a TCollection/TList by z-position

  // Fail if comparing to some other type of object
  if( !obj || IsA() != obj->IsA() )
    return -1;

  if( obj == this )
    return 0;

  const WirePlane* other = static_cast<const WirePlane*>( obj );
  if( GetZ() < other->GetZ() ) return -1;
  if( GetZ() > other->GetZ() ) return  1;
  return 0;
}

//_____________________________________________________________________________

}

ClassImp(TreeSearch::WirePlane)

///////////////////////////////////////////////////////////////////////////////

