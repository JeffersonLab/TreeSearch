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
#include "TClass.h"
#include "MWDC.h"
#include "TimeToDistConv.h"
#include "Projection.h"
#include <iostream>
#include <string>

using namespace std;

// Database uses ns for TDC offsets
static const float kTDCscale = 1e-9;

namespace TreeSearch {

//_____________________________________________________________________________
WirePlane::WirePlane( const char* name, const char* description,
		      THaDetectorBase* parent )
  : THaSubDetector(name,description,parent), fPlaneNum(kMaxUInt),
    fType(kUndefinedType), fWireStart(0.0), fWireSpacing(0.0), 
    fPartner(0), fProjection(0), fMWDC(0), fResolution(0.0),
    fMinTime(-kBig), fMaxTime(kBig), fTTDConv(0), fHits(0), fFitCoords(0)
#ifdef TESTCODE
  , fNmiss(0), fNrej(0), fWasSorted(0), fNhitwires(0), fNmultihit(0),
    fNmaxmul(0), fNcl(0), fNdbl(0), fClsiz(0)
#endif
{
  // Constructor

  static const char* const here = "WirePlane";
  assert( name && parent );

  fMWDC = dynamic_cast<MWDC*>( GetDetector() );
  assert( fMWDC );

  if( fMWDC->TestBit(MWDC::kMCdata) ) // Monte Carlo data mode?
    fHits = new TClonesArray("TreeSearch::MCHit", 200);
  else
    fHits = new TClonesArray("TreeSearch::Hit", 200);

  fFitCoords = new TClonesArray("TreeSearch::FitCoord", 20 );

  if( !fHits or !fFitCoords ) {
    Fatal( Here(here), "Allocating hit array in wire plane %s failed. "
	   "Call expert." );
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

  delete fFitCoords;
  delete fHits;
}

//_____________________________________________________________________________
FitCoord* WirePlane::AddFitCoord( const FitCoord& coord )
{ 
  // Add given fit coordinate data to this plane's array of fit coordinates

  return
    new( (*fFitCoords)[GetNcoords()] ) FitCoord(coord);
}

//_____________________________________________________________________________
void WirePlane::Clear( Option_t* opt )
{    
  // Clear event-by-event data (hits)

  fHits->Clear();
  fFitCoords->Clear();
#ifdef TESTCODE
  fWasSorted = 0;
  fNmiss = fNrej = fNhitwires = fNmultihit = fNmaxmul = 0;
  fNcl = fNdbl = fClsiz = 0;
#endif
}

//_____________________________________________________________________________
#ifdef TESTCODE
void WirePlane::CheckCrosstalk()
{
  // Utility function to check crosstalk statistics.
  // Counts number of wire pairs (=adjacent hits) and max "cluster" size.
  // Also, marks mulithits and calculates their time differences

  UInt_t cursiz = 1;
  fClsiz = 1;
  Int_t prev_iw = -(1<<16);
  Hit* prev_hit = 0;
  for( Int_t i = 0; i < GetNhits(); ++i ) {
    Hit* hit = (Hit*)fHits->At(i);
    Int_t iw = hit->GetWireNum();
    Int_t dw = TMath::Abs(iw - prev_iw);
    if( dw == 0 ) {
      hit->fMulti = 1;
      prev_hit->fMulti = 1;
      hit->fTdiff = hit->GetDriftTime() - prev_hit->GetDriftTime();
    } else if( dw == 1 ) {
      if( cursiz == 1 ) {
	++fNcl;
	++fNdbl;
	prev_hit->fCl = 1;
      }
      ++cursiz;
      ++fNdbl;
      hit->fCl = 1;
      if( cursiz > fClsiz )
	fClsiz = cursiz;
    } else {
      cursiz = 1;
    }
    prev_iw = iw;
    prev_hit = hit;
  }
}
#endif

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
  bool no_time_cut = fMWDC->TestBit(MWDC::kDoTimeCut) == kFALSE;
  bool mc_data     = fMWDC->TestBit(MWDC::kMCdata);

  // Decode data. This is done fairly efficiently by looping over only the 
  // channels with hits on each module. 
  // FIXME: If a module is shared with another plane (common here), we waste
  // time skipping hits that don't belong to us.
  // NB: certain indices below are guaranteed to be in range by construction
  // in ReadDatabase, so we can avoid unneeded checks.
  bool sorted = true;
  Hit* prevHit = 0;
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module * d = fDetMap->GetModule(imod);
    Double_t ref_time =
      (d->refindex >= 0) ? fMWDC->GetRefTime(d->refindex) : 0.0;

    // Get number of channels with hits and loop over them, skipping channels
    // that are not part of this module
    // FIXME: this becomes very inefficient if several modules with the 
    // same crate/slot are defined - e.g. one "module" per channel...ouch
    Int_t nchan = evData.GetNumChan( d->crate, d->slot );
    // For "reversed" detector map modules, loop backwards over the channels
    // to preserve the ordering of the hits by wire number
    Int_t ichan, incr;
    if( d->reverse ) {
      ichan = nchan-1;
      incr = -1;
    } else {
      ichan = 0;
      incr = 1;
    }
    for( ; ichan < nchan && ichan >= 0; ichan += incr ) {
      Int_t chan = evData.GetNextChan( d->crate, d->slot, ichan );
      if( chan < d->lo || chan > d->hi ) {
#ifdef TESTCODE
	++fNmiss;
#endif
	continue; //Not part of this detector
      }
      // Get the wire number. Assumes that the logical channels in the detector
      // map are defined in order of ascending wire numbers.
      Int_t iw = d->first + ( (d->reverse) ? d->hi-chan : chan-d->lo );
      Double_t tdc_offset = fTDCOffset[iw];

      // Get number of hits on this wire and loop over hits
      Int_t nhits = evData.GetNumHits( d->crate, d->slot, chan );
#ifdef TESTCODE
      if( nhits > 0 ) {
	++fNhitwires; 
	if( nhits > 1 )
	  ++fNmultihit;
	if( (UInt_t)nhits > fNmaxmul )
	  fNmaxmul = nhits;
      }
#endif
      for( Int_t hit = 0; hit < nhits; hit++ ) {
	
	// Get the TDC data for this hit
	Int_t data = evData.GetData( d->crate, d->slot, chan, hit );
	
	// Convert the TDC value to the drift time. The readout is assumed to
	// use common-stop TDCs, so t_drift = t_tdc(drift=0)-t_tdc(data).
	Double_t time = tdc_offset+ref_time - d->resolution*(data+0.5);
	if( no_time_cut || (fMinTime < time && time < fMaxTime) ) {
	  Hit* theHit;
	  if( mc_data ) {
	    theHit = new( (*fHits)[nHits++] ) 
	      MCHit( iw, 
		     fWireStart + iw * fWireSpacing,
		     data,
		     time,
		     fResolution,
		     this,
		     // TODO: fill MC info here
		     0, 0.0
		     );
	  } else
	    theHit = new( (*fHits)[nHits++] )
	      Hit( iw, 
		   fWireStart + iw * fWireSpacing,
		   data,
		   time,
		   fResolution,
		   this
		   );
	  // Preliminary calculation of drift distance. Once tracks are known,
	  // the distance can be recomputed using the track slope.
	  theHit->ConvertTimeToDist( 0.0 );

	  // We can test the ordering of the hits on the fly - they should
	  // come in sorted if the lowest logical channel corresponds to
	  // the smallest wire positiion. If they do, we can skip 
	  // quicksorting an already-sorted array ;)
	  if( sorted && prevHit && theHit->Hit::Compare(prevHit) < 0 )
	    sorted = false;
	  prevHit = theHit;
	}
#ifdef TESTCODE
	else
	  ++fNrej;
#endif
      } // hits
    }   // chans
  }     // modules
   
  // If necessary, sort the hits by wire position
  if( !sorted )
    fHits->Sort();

#ifdef TESTCODE
  fWasSorted = sorted;
  CheckCrosstalk();
#endif

  return nHits;
}
  
//_____________________________________________________________________________
Int_t WirePlane::DefineVariables( EMode mode )
{
  // initialize global variables


  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { "nhits",       "Num accepted hits",  "GetNhits()" },
    { "hit.wire",    "Hit wire number",    "fHits.TreeSearch::Hit.fWireNum" },
    { "hit.tdc",     "Hit TDC value",      "fHits.TreeSearch::Hit.fRawTDC" },
    { "hit.time",    "Hit time (s)",       "fHits.TreeSearch::Hit.fTime" },
    { "hit.dist",    "Drift distance (m)",
                                    "fHits.TreeSearch::Hit.GetDriftDist()" },
    { "ncoords",     "Num fit coords",     "GetNcoords()" },
    { "coord.rank",  "Fit rank of coord",
                                "fFitCoords.TreeSearch::FitCoord.fFitRank" },
    { "coord.time", "Drift time of hit (s)",
                          "fFitCoords.TreeSearch::FitCoord.GetDriftTime()" },
    { "coord.dist", "Drift distance of hit (m)",
                          "fFitCoords.TreeSearch::FitCoord.GetDriftDist()" },
    { "coord.pos",   "Position used in uncorrected fit (m)",
                                    "fFitCoords.TreeSearch::FitCoord.fPos" },
    { "coord.trkpos","Track pos from uncorrected fit (m)",
                               "fFitCoords.TreeSearch::FitCoord.fTrackPos" },
    { "coord.trkslope","Track slope from uncorrected fit",
                             "fFitCoords.TreeSearch::FitCoord.fTrackSlope" },
    { "coord.chi2",  "Chi2 of this coord's fit",
                               "fFitCoords.TreeSearch::FitCoord.GetChi2()" },
    { "coord.resid", "Residual of trkpos (m)",
                           "fFitCoords.TreeSearch::FitCoord.GetResidual()" },
    { "coord.trkdist", "Distance of trkpos to wire (m)",
                          "fFitCoords.TreeSearch::FitCoord.GetTrackDist()" },
#ifdef TESTCODE
    { "nmiss",       "Decoder misses",     "fNmiss" },
    { "nrej",        "Time cut nopass",    "fNrej" },
    { "sorted",      "Wires were ordered", "fWasSorted" },
    { "nwhit",       "Num wires w/hits>0", "fNhitwires" },
    { "nmulti",      "Num wires w/hits>1", "fNmultihit" },
    { "maxmul",      "Max num hits/wire",  "fNmaxmul" },
    { "ncl",         "Num clusters",       "fNcl" },
    { "ndbl",        "Num double hits ",   "fNdbl" },
    { "maxclsiz",    "Max cluster size",   "fClsiz" },
    { "hit.iscl",    "Hit has neighbor",   "fHits.TreeSearch::Hit.fCl" },
    { "hit.ismulti", "Wire has multihits", "fHits.TreeSearch::Hit.fMulti" },
    { "hit.tdiff",   "multi hits tdiff",   "fHits.TreeSearch::Hit.fTdiff" },
#endif
    { 0 }
  };
  Int_t ret = DefineVarsFromList( vars, mode );

  if( fMWDC->TestBit(MWDC::kMCdata) && ret == kOK ) {
    // Additional variables for Monte Carlo data
    RVarDef mcvars[] = {
      { "mcpos", "MC track position (m)", "fHits.TreeSearch::MCHit.fMCPos" },
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
  }
  return ret;
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
//   THaDetectorBase *parent = GetDetector();
//   if( parent )
//     fOrigin += parent->GetOrigin();

  return fStatus = kOK;
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

  TString plane_type, ttd_conv;
  // Putting these containers on the stack may cause a stack overflow
  vector<Int_t>* detmap = new vector<Int_t>;
  vector<double>* ttd_param = new vector<double>;
  // Default values for optional parameters
  fMinTime = -kBig;
  fMaxTime =  kBig;
  Int_t required = 0;
  const DBRequest request[] = {
    { "detmap",        detmap,        kIntV },
    { "nwires",        &fNelem,       kInt },
    { "type",          &plane_type,   kTString, 0, 1 },
    { "wire.pos",      &fWireStart },
    { "wire.spacing",  &fWireSpacing, kDouble,  0, 0, -1 },
    { "ttd.converter", &ttd_conv,     kTString, 0, 0, -1 },
    { "ttd.param",     ttd_param,     kDoubleV, 0, 0, -1 },
    { "xp.res",        &fResolution,  kDouble,  0, 0, -1 },
    { "tdc.offsets",   &fTDCOffset,   kFloatV,  0, 0 },
    { "description",   &fTitle,       kTString, 0, 1 },
    { "drift.min",     &fMinTime,     kDouble,  0, 1, -1 },
    { "drift.max",     &fMaxTime,     kDouble,  0, 1, -1 },
    { "required",      &required,     kInt,     0, 1 },
    { 0 }
  };

  Int_t status = kInitError;
  UInt_t flags;
  err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( !err ) {
    // Parse the detector map of the data channels
    flags = THaDetMap::kFillRefChan;
    if( FillDetMap( *detmap, flags, here ) > 0 )
      status = kOK;
  }
  delete detmap; detmap = 0;

  // Create time-to-distance converter
  if( status == kOK ) {
    if( !ttd_conv.Contains("::") )
      ttd_conv.Prepend("TreeSearch::");
    const char* s = ttd_conv.Data();
    TClass* cl = TClass::GetClass( s );
    if( !cl ) {
      Error( Here(here), "Drift time-to-distance converter \"%s\" not "
	     "available. Load library or fix database.", s?s:"" );
      status = kInitError;
      goto ttderr;
    }
    if( !cl->InheritsFrom( TreeSearch::TimeToDistConv::Class() )) {
      Error( Here(here), "Class \"%s\" is not a drift time-to-distance "
	     "converter. Fix database.", s );
      status = kInitError;
      goto ttderr;
    }
    fTTDConv = static_cast<TimeToDistConv*>( cl->New() );
    if( !fTTDConv ) {
      Error( Here(here), "Unexpected error creating drift time-to-distance "
	     "converter object \"%s\". Call expert.", s );
      status = kInitError;
      goto ttderr;
    } 
    if( fTTDConv->SetParameters( *ttd_param ) != 0 ) {
      Error( Here(here), "Error initializing drift time-to-distance converter "
	     "\"%s\". Check ttd.param in database.", s );
      status = kInitError;
    }
  }
ttderr:
  delete ttd_param; ttd_param = 0;
  if( status != kOK )
    return status;

  // Retrieve TDC resolution and model number for our crateslots
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fDetMap->GetModule(imod);
    fMWDC->LoadDAQmodel(d);
    fMWDC->LoadDAQresolution(d);
    d->MakeTDC();
    UInt_t nchan = fMWDC->GetDAQnchan(d);
    if( d->hi >= nchan ) {
      Error( Here(here), "Detector map channel out of range for module "
	     "cr/sl/lo/hi = %u/%u/%u/%u. Must be < %u. Fix database.",
	     d->crate, d->slot, d->lo, d->hi, nchan );
      return kInitError;
    }
    if( d->refchan >= static_cast<Int_t>(nchan) ) {
      Error( Here(here), "Detector map reference channel %d out of range for "
	     "module cr/sl/lo/hi = %u/%u/%u/%u. Must be < %u. Fix database.",
	     d->refchan, d->crate, d->slot, d->lo, d->hi, nchan );
      return kInitError;
    }
  }

  // Sanity checks
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
  // Convert TDC offsets and timing cuts to seconds
  for( vector<float>::size_type i = 0; i < fTDCOffset.size(); ++i ) {
    fTDCOffset[i] *= kTDCscale;
  }
  if( fMinTime > -kBig )  fMinTime *= kTDCscale;
  if( fMaxTime <  kBig )  fMaxTime *= kTDCscale;

  // Determine the type of this plane. If the optional plane type variable is
  // not given, use the first character of the plane name.
  TString name = plane_type.IsNull() ? fName[0] : plane_type[0];
  fType = fMWDC->NameToType( name );
  if( fType == kUndefinedType ) {
    TString names;
    for( EProjType i = kTypeBegin; i < kTypeEnd; ++i ) {
      names += fMWDC->fProj[i]->GetName();
      if( i+1 != kTypeEnd ) 
	names += " ";
    }
    Error( Here(here), "Unsupported plane type \"%s\". Must be one of "
	   "%s. Fix database.", name.Data(), names.Data() );
    return kInitError;
  }

  if( required )
    SetBit( kIsRequired );

  fIsInit = true;
  return kOK;
}

//_____________________________________________________________________________
void WirePlane::SetPartner( WirePlane* p )
{    
  // Partner this plane with plane 'p'. Partner planes are expected to
  // be located close to each other and usually to have staggered wires.

  if( p )
    p->fPartner = this;
  else if( fPartner ) {
    assert( this == fPartner->fPartner );
    fPartner->fPartner = 0;
  }
  fPartner = p;

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
// Int_t WirePlane::Compare( const TObject* obj ) const 
// {
//   // Used to sort planes in a TCollection/TList by z-position

//   // Fail if comparing to some other type of object
//   assert( obj && IsA() == obj->IsA() );

//   if( obj == this )
//     return 0;

//   const WirePlane* other = static_cast<const WirePlane*>( obj );
//   if( GetZ() < other->GetZ() ) return -1;
//   if( GetZ() > other->GetZ() ) return  1;
//   return 0;
// }

//_____________________________________________________________________________
Double_t WirePlane::GetMaxSlope() const
{ 
  return fProjection ? fProjection->GetMaxSlope() : kBig;
}

//_____________________________________________________________________________

}

ClassImp(TreeSearch::WirePlane)

///////////////////////////////////////////////////////////////////////////////

