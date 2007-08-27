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

namespace TreeSearch {

//_____________________________________________________________________________
WirePlane::WirePlane( const char* name, const char* description,
		      THaDetectorBase* parent )
  : THaSubDetector(name,description,parent), fPlaneNum(-1),
    fType(kUndefinedType), fWireStart(0.0), fWireSpacing(0.0), 
    fPartner(NULL), fProjection(NULL), fMWDC(NULL), fTDCRes(0.0), 
    fDriftVel(0.0), fResolution(0.0), fTTDConv(NULL), fRefTime(NULL),
    fNmiss(0), fNrej(0), fWasSorted(0), fNhitwires(0), fNnohits(0)
{
  // Constructor

  static const char* const here = "WirePlane";

  fRefMap = new THaDetMap;

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
  delete [] fRefTime;
  delete fRefMap;
}

//_____________________________________________________________________________
Int_t WirePlane::ReadDatabase( const TDatime& date )
{
  // Read database

  static const char* const here = "ReadDatabase";


  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Int_t err = ReadGeometry( file, date, kTRUE );
  if( err )
    return err;

  TString plane_type;
  vector<int> detmap, refmap;
  delete [] fRefTime; fRefTime = NULL;

  DBRequest request[] = {
    { "detmap",       &detmap,       kIntV },
    { "refmap",       &refmap,       kIntV,    0, 1 },
    { "nwires",       &fNelem,       kInt },
    { "type",         &plane_type,   kTString, 0, 1 },
    { "wire.pos",     &fWireStart },
    { "wire.spacing", &fWireSpacing, kDouble,  0, 0, -1 },
    //FIXME: remove?
    { "tdc.res",      &fTDCRes,      kDouble,  0, 0, -1 },
    { "drift.v",      &fDriftVel,    kDouble,  0, 0, -1 },
    { "xp.res",       &fResolution,  kDouble,  0, 0, -1 },
    { "tdc.offsets",  &fTDCOffset,   kFloatV,  0, 1 },
    { "description",  &fTitle,       kTString, 0, 1 },
    { 0 }
  };

  err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return kInitError;

  if( fNelem <= 0 ) {
    Error( Here(here), "Invalid number of wires: %d", fNelem );
    return kInitError;
  }

  // Fill the reference channel detector map, if given
  UInt_t flags;
  if( refmap.size() > 0 ) {
    flags = THaDetMap::kFillModel;
    if( fRefMap->Fill( refmap, flags ) <= 0 ) {
      Error( Here(here), "Invalid reference channel map. Fix database." );
      fRefMap->Clear();
      return kInitError;
    }
  }

  // Fill the detector map for the data channels
  flags = ( THaDetMap::kFillModel | THaDetMap::kFillRefIndex );
  if( FillDetMap( detmap, flags, here ) <= 0 )
    return kInitError;
  
  Int_t nchan = fDetMap->GetTotNumChan();
  if( nchan != fNelem ) {
    Error( Here(here), "Number of detector map channels (%d) "
	   "disagrees with number of wires (%d)", nchan, fNelem );
    return kInitError;
  }
  nchan = fTDCOffset.size();
  if( nchan > 0 && nchan != fNelem ) {
    Error( Here(here), "Number of TDC offset values (%d) "
	   "disagrees with number of wires (%d)", nchan, fNelem );
    return kInitError;
  }

  // Check consistency of reference channels and data channels
  Int_t dmin, dmax;
  fDetMap->GetMinMaxChan( dmin, dmax, THaDetMap::kRefIndex );
  Int_t nrefchan = fRefMap->GetTotNumChan();
  if( nrefchan > 0 ) {
    if( dmin < 0 && dmax < 0 ) {
      Warning( Here(here), "Reference channels defined but not used. "
	       "Check database." );
      fRefMap->Clear();
    } else {
      Int_t rmin, rmax;
      fRefMap->GetMinMaxChan( rmin, rmax );
      if( (dmin >= 0 && dmin < rmin) || dmax > rmax ) {
	Error( Here(here), "Reference channel(s) out of range: min/max=%d/%d, "
	       "requested=%d/%d. Fix database.", rmin, rmax, dmin, dmax );
	fRefMap->Clear();
	return kInitError;
      }
      fRefTime = new Double_t[ nrefchan ];
    }
  } else if ( fRefMap->GetSize() > 0 ) {
    Error( Here(here), "Total number of reference channels = %d <= 0? "
	   "Fix database.", nrefchan );
    fRefMap->Clear();
    return kInitError;
  } else if( dmax >= 0 ) {
    Error( Here(here), "detmap specifies refindex, but no refmap defined. "
	   "Fix database." );
    return kInitError;
  }

  // Determine the type of this plane. If the optional plane type variable is
  // not in the database, use the first character of the plane name.
  string name = plane_type.IsNull() ? fName.Data() : plane_type.Data();
  name = name.substr(0,1);
  fType = fMWDC->NameToType( name.c_str() );
  if( fType == kUndefinedType ) {
    Error( Here(here), "Illegal plane type (%s). Fix database.", name.c_str());
    return kInitError;
  }

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

  static const char* const here = "Decode";

  UInt_t nHits = 0;
  bool pos_only = fMWDC->TestBit(MWDC::kIgnoreNegDrift);

  // If any reference channels are specified, we need to read them first.
  // They use a separate detector map, which has the same structure as the
  // one for the data channels.
  UInt_t nref = 0;
  if( fRefMap->GetSize() > 0 ) {
    for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
      THaDetMap::Module* d = fRefMap->GetModule(imod);	

      // Since we expect all channels to have a hit, we can efficiently loop
      // over the defined channel range
      for( Int_t chan = d->lo; chan <= d->hi; ++chan ) {
	Int_t data;
	Int_t nhits = evData.GetNumHits( d->crate, d->slot, chan );
	if( nhits > 0 ) {
	  data = evData.GetData( d->crate, d->slot, chan, nhits-1 );
	  if( nhits > 1 ) {
	    Warning( Here(here), "%d hits on reference channel %d module %d", 
		     nhits, chan, imod );
	  }
	} else {
	  Error( Here(here), "No hits on reference channel %d module %d. "
		 "Event decoding failed.", chan, imod );
	  return -1;
	}
	// FIXME: separate resolution for reference channels
	//	fRefTime[nref] = d->resolution * data;
	fRefTime[nref] = fTDCRes * data;
	++nref;
      }
    }
  }

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
    Double_t ref_offset = (d->refindex >= 0) ? fRefTime[d->refindex] : 0.0;

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
      // Get the wire number. Assumes that detector map is defined in order
      // of ascending wire numbers.
      Int_t iw = d->first + chan - d->lo;
      Double_t tdc_offset = fTDCOffset[iw];

      // Get number of hits on this channel and loop over hits
      Int_t nhits = evData.GetNumHits( d->crate, d->slot, chan );
      if( nhits > 0 ) ++fNhitwires; else ++fNnohits;
      for( Int_t hit = 0; hit < nhits; hit++ ) {
	
	// Get the TDC data for this hit
	Int_t data = evData.GetData( d->crate, d->slot, chan, hit );
	
	// Convert the TDC value to the drift time
	//FIXME: implement resolution in THaDetMap
	//Double_t time = d->resolution * (data+0.5) - tdc_offset - ref_offset;
	Double_t time = fTDCRes * (data+0.5) - tdc_offset - ref_offset;
	if( !pos_only || time > 0.0 ) {
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
   
  // Sort the hits in order of increasing wire position and, for the same wire
  // position) increasing time
  fWasSorted = sorted ? 1 : revsorted ? -1 : 0;
  if( !sorted ) {
    if( revsorted ) {
      // reverse array ...
      // TODO: check if faster than qsort on revsorted data
      // TODO: can't we just leave it revsorted?
      TClonesArray* copy = new TClonesArray( fHits->GetClass(), 
					     fHits->GetSize() );
      Int_t end = fHits->GetLast();
      for( Int_t i = 0; i < end+1; ++i )
	new((*copy)[end-i]) Hit( *static_cast<Hit*>((*fHits)[i]) );
      delete fHits;
      fHits = copy;
    } else {
      fHits->Sort();
    }
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

