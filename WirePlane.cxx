//*-- Author :    Ole Hansen, Jefferson Lab   14-Jun-2007

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TreeSearch::WirePlane                                                    //
//                                                                          //
// A 1-dimensional plane of horizontal drift chamber wires                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "WirePlane.h"
#include "WireHit.h"
#include "MWDC.h"
#include "TimeToDistConv.h"
#include "Projection.h"

#include "THaDetMap.h"
#include "TClonesArray.h"
#include "TError.h"
#include "TClass.h"
#include "TString.h"
#ifndef MCDATA
#include "THaEvData.h"
#else
#include "SimDecoder.h"
#endif

#include <iostream>
#include <stdexcept>

using namespace std;
using namespace Podd;

// Database uses ns for TDC offsets
static const float kTDCscale = 1e-9;

namespace TreeSearch {

//_____________________________________________________________________________
WirePlane::WirePlane( const char* name, const char* description,
		      THaDetectorBase* parent )
  : Plane(name,description,parent),
    fMinTime(-kBig), fMaxTime(kBig), fTTDConv(0)
  , fNmiss(0), fNrej(0), fWasSorted(0), fNhitwires(0), fNmultihit(0),
    fNmaxmul(0), fNcl(0), fNdbl(0), fClsiz(0)
{
  // Constructor

  static const char* const here = "WirePlane";

  assert( dynamic_cast<MWDC*>(fTracker) );

  try {
#ifdef MCDATA
    if( fTracker->TestBit(Tracker::kMCdata) ) // Monte Carlo data mode?
      fHits = new TClonesArray("TreeSearch::MCWireHit", 200);
    else
#endif
      fHits = new TClonesArray("TreeSearch::WireHit", 200);
  }
  catch( std::bad_alloc& ) {
    Error( Here(here), "Out of memory allocating hit array for wire plane %s. "
	   "Call expert.", name );
    MakeZombie();
    return;
  }

  // Tell common code that our frontends are TDCs with reference channels
  SetBit(kTDCReadout);
  SetBit(kHaveRefChans);
}

//_____________________________________________________________________________
WirePlane::~WirePlane()
{
  // Destructor

  if( fIsSetup )
    RemoveVariables();
}

//_____________________________________________________________________________
void WirePlane::Clear( Option_t* opt )
{
  // Clear event-by-event data (hits)

  Plane::Clear(opt);

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
  WireHit* prev_hit = 0;
  for( Int_t i = 0; i < GetNhits(); ++i ) {
    WireHit* hit = static_cast<WireHit*>(fHits->At(i));
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
Int_t WirePlane::WireDecode( const THaEvData& evData )
{
  // Extract this plane's hit data from the raw evData.
  //
  // This routine can handle both the old Fastbus readout and the new CAEN
  // VME pipeline TDCs. The latter require a reference channel map and
  // cross-references to reference channels in the regular detector map
  // of the plane.

  // const char* const here = "WirePlane::Decode";

  assert( dynamic_cast<MWDC*>(fTracker) );
  MWDC* mwdc = static_cast<MWDC*>( fTracker );

  UInt_t nHits = 0;
  bool no_time_cut = !fTracker->TestBit(MWDC::kDoTimeCut);
#ifdef MCDATA
  bool mc_data     = fTracker->TestBit(Tracker::kMCdata);
  assert( !mc_data || dynamic_cast<const SimDecoder*>(&evData) != 0 );
#endif

  // Decode data. This is done fairly efficiently by looping over only the
  // channels with hits on each module.
  // FIXME: If a module is shared with another plane (common here), we waste
  // time skipping hits that don't belong to us.
  // NB: certain indices below are guaranteed to be in range by construction
  // in ReadDatabase, so we can avoid unneeded checks.
  bool sorted = true;
  WireHit* prevHit = 0;
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module * d = fDetMap->GetModule(imod);
    Double_t ref_time =
      (d->refindex >= 0) ? mwdc->GetRefTime(d->refindex) : 0.0;

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
	  WireHit* theHit;
#ifdef MCDATA
	  if( mc_data ) {
	    theHit = new( (*fHits)[nHits++] )
	      MCWireHit( iw,
			 GetStart() + iw * GetPitch(),
			 data,
			 time,
			 fResolution,
			 this,
			 // TODO: fill MC info here
			 0, 0.0, 0.0
			 );
	  } else
#endif
	    theHit = new( (*fHits)[nHits++] )
	      WireHit( iw,
		       GetStart() + iw * GetPitch(),
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
	  if( sorted && prevHit && theHit->Compare(prevHit) < 0 )
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

  // Negative return value indicates potential problem
  if( nHits > fMaxHits )
    return -nHits;

  return nHits;
}

//_____________________________________________________________________________
Hit* WirePlane::AddHitImpl( Double_t pos )
{
  // Make a dummy hit of the correct type at the given projection coordinate
  // and add it to the hit array

  assert( IsDummy() );

  WireHit* theHit = 0;

  // Emulate parameters for dummy hits
  const Int_t iw = TMath::Nint( (pos-GetStart()) / GetPitch() );
  const Int_t data = 0;
  const Double_t time = 0.0, resolution = fResolution;

#ifdef MCDATA
  const Int_t mctrack = 1;
  const Double_t mcpos = pos, mctime = 0.0;
  bool mc_data = fTracker->TestBit(Tracker::kMCdata);
  if( mc_data )
    // Monte Carlo data
    theHit = new( (*fHits)[GetNhits()] ) MCWireHit( iw,
						    GetStart() + iw * GetPitch(),
						    data,
						    time,
						    resolution,
						    this,
						    mctrack,
						    mcpos,
						    mctime
						    );
  else
#endif
    theHit = new( (*fHits)[GetNhits()] ) WireHit( iw,
						  GetStart() + iw * GetPitch(),
						  data,
						  time,
						  resolution,
						  this
						  );
  return theHit;
}

//_____________________________________________________________________________
Int_t WirePlane::Decode( const THaEvData& evData )
{
  // Convert evData to hits

  if( IsDummy() )
    // Special "decoding" for dummy planes
    return DummyDecode( evData );

  return WireDecode( evData );
}

//_____________________________________________________________________________
Int_t WirePlane::DefineVariables( EMode mode )
{
  // initialize global variables


  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { "nhits",       "Num accepted hits",                           "GetNhits()" },
    { "ncoords",     "Num fit coords",                              "GetNcoords()" },
    // { "coord.wire", "Wire number of fitted hit",                    "fFitCoords.TreeSearch::FitCoord.GetWireNum()" },
    // { "coord.time", "Drift time of hit (s)",                        "fFitCoords.TreeSearch::FitCoord.GetDriftTime()" },
    // { "coord.dist", "Drift distance of hit (m)",                    "fFitCoords.TreeSearch::FitCoord.GetDriftDist()" },
    { "coord.pos",   "Position used in fit (wirepos +/- dist) (m)", "fFitCoords.TreeSearch::FitCoord.fPos" },
    { "coord.trkpos","Track pos from projection fit (m)",           "fFitCoords.TreeSearch::FitCoord.fTrackPos" },
    { "coord.trkslope","Track slope from projection fit",           "fFitCoords.TreeSearch::FitCoord.fTrackSlope" },
    { "coord.resid", "Residual of trkpos (m)",                      "fFitCoords.TreeSearch::FitCoord.GetResidual()" },
    { "coord.trkdist", "Distance of trkpos to wire (m)",            "fFitCoords.TreeSearch::FitCoord.GetTrackDist()" },
    { "coord.3Dpos",  "Crossing position of fitted 3D track (m)",   "fFitCoords.TreeSearch::FitCoord.f3DTrkPos" },
    { "coord.3Ddist", "Distance of 3D trkpos to wire (m)",          "fFitCoords.TreeSearch::FitCoord.Get3DTrkDist()" },
    { "coord.3Dresid","Residual of 3D trkpos (m)",                  "fFitCoords.TreeSearch::FitCoord.Get3DTrkResid()" },
    { "coord.3Dslope","Slope of fitted 3D track wrt projection",    "fFitCoords.TreeSearch::FitCoord.f3DTrkSlope" },
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
#endif
    { 0 }
  };
  Int_t ret = DefineVarsFromList( vars, mode );

#ifdef MCDATA
  if( !fTracker->TestBit(MWDC::kMCdata) && ret == kOK ) {
#endif
    // Non-Monte Carlo hit data
    RVarDef nonmcvars[] = {
      { "hit.wire",    "Hit wire number",    "fHits.TreeSearch::WireHit.fWireNum" },
      { "hit.tdc",     "Hit TDC value",      "fHits.TreeSearch::WireHit.fRawTDC" },
      { "hit.time",    "Hit time (s)",       "fHits.TreeSearch::WireHit.fTime" },
      { "hit.dist",    "Drift distance (m)", "fHits.TreeSearch::WireHit.GetDriftDist()" },
#ifdef TESTCODE
      { "hit.iscl",    "Hit has neighbor",   "fHits.TreeSearch::WireHit.fCl" },
      { "hit.ismulti", "Wire has multihits", "fHits.TreeSearch::WireHit.fMulti" },
      { "hit.tdiff",   "multi hits tdiff",   "fHits.TreeSearch::WireHit.fTdiff" },
#endif
      { 0 }
    };
    ret = DefineVarsFromList( nonmcvars, mode );
#ifdef MCDATA
  } else {
    // Monte Carlo hit data includes the truth information
    RVarDef mcvars[] = {
      { "hit.wire",    "Hit wire number",    "fHits.TreeSearch::MCWireHit.fWireNum" },
      { "hit.tdc",     "Hit TDC value",      "fHits.TreeSearch::MCWireHit.fRawTDC" },
      { "hit.time",    "Hit time (s)",       "fHits.TreeSearch::MCWireHit.fTime" },
      { "hit.dist",    "Drift distance (m)", "fHits.TreeSearch::MCWireHit.GetDriftDist()" },
#ifdef TESTCODE
      { "hit.iscl",    "Hit has neighbor",   "fHits.TreeSearch::MCWireHit.fCl" },
      { "hit.ismulti", "Wire has multihits", "fHits.TreeSearch::MCWireHit.fMulti" },
      { "hit.tdiff",   "multi hits tdiff",   "fHits.TreeSearch::MCWireHit.fTdiff" },
#endif
      { "mcpos", "MC track position (m)", "fHits.TreeSearch::MCWireHit.fMCPos" },
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
  }
#endif // MCDATA
  return ret;
}

//_____________________________________________________________________________
Int_t WirePlane::ReadDatabase( const TDatime& date )
{
  // Read database

  static const char* const here = "ReadDatabase";

  // Read the database for the base class, quit if error
  Int_t status = ReadDatabaseCommon(date);
  if( status != kOK )
    return status;

  // Now read additional info needed for this derived class
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  TString ttd_conv;
  vector<double>* ttd_param = 0;
  // Default values for optional parameters
  fMinTime = -kBig;
  fMaxTime =  kBig;
  try {
    // Putting this container on the stack may cause a stack overflow
    ttd_param = new vector<double>;

    Int_t gbl = GetDBSearchLevel(fPrefix);
    const DBRequest request[] = {
      { "nwires",        &fNelem,       kInt,     0, 0, gbl },
      { "wire.pos",      &fStart },
      { "wire.spacing",  &fPitch,       kDouble,  0, 0, gbl },
      { "ttd.converter", &ttd_conv,     kTString, 0, 0, gbl },
      { "ttd.param",     ttd_param,     kDoubleV, 0, 0, gbl },
      { "tdc.offsets",   &fTDCOffset,   kFloatV,  0, 0 },
      { "drift.min",     &fMinTime,     kDouble,  0, 1, gbl },
      { "drift.max",     &fMaxTime,     kDouble,  0, 1, gbl },
      { 0 }
    };
    status = LoadDB( file, date, request, fPrefix );
  }
  // Catch exceptions here so that we can close the file
  catch(...) {
    fclose(file);
    throw;
  }

  // Finished reading the database
  fclose(file);
  if( status != kOK )
    return status;

  // Dummy planes ignore all of the parameters that are checked below,
  // so we can return right here.
  if( IsDummy() ) {
    fIsInit = true;
    return kOK;
  }

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

  fIsInit = true;
  return kOK;
}

//_____________________________________________________________________________
void WirePlane::Print( Option_t* opt ) const
{
  // Print plane info

  Plane::Print(opt);
}

//_____________________________________________________________________________
Hit* WirePlane::FindNearestHitAndPos( Double_t x, Double_t& pmin ) const
{
  // Version of FindNearestHitandPos for horizontal drift chambers. Carries out
  // additional checks to account for left/right ambiguity.

  pmin = kBig;
  if( GetNhits() == 0 )
    return 0;

  const Int_t pos = FindHitWithLowerBound( x );

  // Decide whether the wire >= x or the first one < x are closest.
  // If the track crosses between two adjacent wires, keep both.
  WireHit *hnext = 0, *hprev = 0;
  if( pos < GetNhits() ) {
    assert( dynamic_cast<WireHit*>(GetHit(pos)) );
    hnext = static_cast<WireHit*>(GetHit(pos));
    assert( hnext->GetWirePos() >= x );
  }
  if( pos > 0 ) {
    assert( dynamic_cast<WireHit*>(GetHit(pos-1)) );
    hprev = static_cast<WireHit*>(GetHit(pos-1));
    assert( hprev->GetWirePos() < x );
    if( hnext ) {
      assert( hprev->GetWireNum() < hnext->GetWireNum() );
      if( hprev->GetWireNum() + 1 < hnext->GetWireNum() ) {
	if( x - hprev->GetWirePos() < hnext->GetWirePos() - x )
	  hnext = 0;
	else
	  hprev = 0;
      }
    }
  }
  // Of the closest wire(s) found, find the closest drift distance.
  // If there are multiple hits one a wire, test all hits - without
  // making assumptions about the order of drift distances
  Double_t dmin = kBig;
  Hit* hmin = 0;
  if( hnext ) {
    hmin = hnext;
    pmin = hnext->GetPosL();
    dmin = TMath::Abs(pmin-x);
    Int_t i = pos;
    WireHit* h;
    while( ++i < GetNhits() and
	   (h = static_cast<WireHit*>(GetHit(i)))->GetWireNum() ==
	   hnext->GetWireNum() ) {
      Double_t d = TMath::Abs(h->GetPosL()-x);
      if( d < dmin ) {
	dmin = d;
	hmin = h;
	pmin = h->GetPosL();
      }
    }
  }
  if( hprev ) {
    Double_t d = TMath::Abs(hprev->GetPosR()-x);
    if( !hmin or d < dmin ) {
      dmin = d;
      hmin = hprev;
      pmin = hprev->GetPosR();
    }
    Int_t i = pos-1;
    WireHit* h;
    while( --i >= 0 and (h = static_cast<WireHit*>(GetHit(i)))->GetWireNum()
	   == hprev->GetWireNum() ) {
      d = TMath::Abs(h->GetPosR()-x);
      if( d < dmin ) {
	dmin = d;
	hmin = h;
	pmin = h->GetPosR();
      }
    }
  }
  return hmin;
}

//_____________________________________________________________________________

}

ClassImp(TreeSearch::WirePlane)

///////////////////////////////////////////////////////////////////////////////

