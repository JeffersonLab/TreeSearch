//*-- Author :    Ole Hansen, Jefferson Lab   11-Jan-2010

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TreeSearch::GEMPlane                                                     //
//                                                                          //
// A 1-dimensional readout plane of a GEM chamber.                          //
// Use two of these objects for a 2-d readout plane.                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "GEMPlane.h"
#include "GEMHit.h"
#include "Tracker.h"
#include "Projection.h"

#include "THaDetMap.h"
#include "THaEvData.h"
#include "TClonesArray.h"
#include "TH1.h"

#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
GEMPlane::GEMPlane( const char* name, const char* description,
		    THaDetectorBase* parent )
  : Plane(name,description,parent), 
    fMapType(kOneToOne), fMaxClusterSize(0), fMinAmpl(0), fSplitFrac(0), 
    fMaxSamp(1), fAmplSigma(0), fADC(0), fADCped(0), fTimeCentroid(0),
    fDnoise(0), fNhitStrips(0), fNsigStrips(0), fNRawHits(0), fRawOcc(0),
    fOccupancy(0), fADCMap(0)
{
  // Constructor

}

//_____________________________________________________________________________
GEMPlane::~GEMPlane()
{
  // Destructor.

  // Histograms in fHist should be automatically deleted by ROOT when 
  // the output file is closed

  if( fIsSetup )
    RemoveVariables();

  delete fADC;
  delete fADCped;
  delete fTimeCentroid;
}

//_____________________________________________________________________________
Int_t GEMPlane::Begin( THaRunBase* run )
{
  // Book diagnostic histos if not already done

  Plane::Begin( run );

#ifdef TESTCODE
  if( TestBit(kDoHistos) ) {
    //TODO: check if exist
    string hname, htitle;
    hname   = fName + "_adcmap";
    htitle  = fName + " adcmap";
    fADCMap = new TH1F( hname.c_str(), htitle.c_str(), fNelem, 0, fNelem );
  }
#endif
  return 0;
}

//_____________________________________________________________________________
void GEMPlane::Clear( Option_t* /* opt */ )
{    
  // Clear event-by-event data (hits)

  assert( fADC );
  memset( fADC, 0, fNelem*sizeof(Float_t) );
  assert( fADCped );
  memset( fADCped, 0, fNelem*sizeof(Float_t) );
  assert( fTimeCentroid );
  memset( fTimeCentroid, 0, fNelem*sizeof(Float_t) );

  fNhitStrips = fNsigStrips = fNRawHits = 0;
  fRawOcc = fOccupancy = 0.0;
  // FIXME: speed up with index table
  for( vector<Vflt_t>::iterator it = fADCsamp.begin(); it != fADCsamp.end();
       ++it ) {
    (*it).clear();
  }
}

//_____________________________________________________________________________
Int_t GEMPlane::MapChannel( Int_t idx ) const
{
  // Map hardware channel number to logical strip number based on mapping
  // prescription from database

  assert( idx >= 0 and idx < fNelem );

  Int_t ret;
  switch( fMapType ) {
  case kOneToOne:
    ret = idx;
    break;
  case kReverse:
    ret = fNelem-idx-1;
    break;
  case kGassiplexAdapter1:
    assert( idx < 240 );
    if( idx == 0 )
      ret = 1;
    else if( idx == 239 )
      ret = 238;
    else if( idx % 2 ) // odd
      ret = idx + 2;
    else               // even
      ret = idx - 2;
    break;
  case kGassiplexAdapter2:
    assert( idx < 240 );
    if( idx == 1 )
      ret = 0;
    else if( idx == 238 )
      ret = 239;
    else if( idx % 2 ) // odd
      ret = idx - 2;
    else               // even
      ret = idx + 2;
    break;
  case kTable:
    // Use the mapping lookup table
    assert( fChanMap.size() == static_cast<vector<Int_t>::size_type>(fNelem) );
    ret = fChanMap[idx];
    break;
  }
  assert( ret >= 0 and ret < fNelem );
  return ret;
}

//_____________________________________________________________________________
static Float_t ChargeDep( const vector<Float_t>& amp, Float_t& centroid )
{
  // Deconvolute signal given by samples in 'amp', return approximate integral.
  // Currently analyzes exactly 3 samples.
  // From Kalyan Allada
  // NIM A326, 112 (1993)

  Float_t delta_t=25; //time interval between samples (ns)
  Float_t Tp=50; // time constant (ns)
  Float_t sig[3];
  Float_t w1, w2, w3, x, sum, r1, r2;
  
  assert( amp.size() >= 3 );

  if(amp[0]==0 || amp[1] == 0 ) {
    centroid = 0.0;    
    return 0.0;
  }
  
  // TODO: calculate centroid
  centroid = 0;
  
  // Weight factors calculated based on the response of the silicon microstrip 
  // detector:
  // v(t) = (delta_t/Tp)*exp(-delta_t/Tp) 
  // Need to update this for GEM detector response(?):
  // v(t) = A*(1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2)
  // where A is the amplitude, t0 the begin of the rise, tau1 the time
  // parameter for the rising edge and tau2 the for the falling edge. 

  x = delta_t/Tp;

  w1 = TMath::Exp(x-1)/x;
  w2 = -2*TMath::Exp(-1)/x;
  w3 = TMath::Exp(-x-1)/x;

  //deconvoluted signal samples
  sig[0] = amp[0]*w1;
  sig[1] = amp[1]*w1+amp[0]*w2;
  sig[2] = amp[2]*w1+amp[1]*w2+amp[0]*w3;

  sum = delta_t*(sig[0]+sig[1]+sig[2]);
  
 // Calculate ratios for 3 samples and check for bad signals
  r1 = Float_t(amp[0])/amp[2];
  r2 = Float_t(amp[1])/amp[2];
  if(r1>1.0 || r2>1.0 || r1 > r2){
    sum = 0.0;
    centroid =0;
    return 0.0;
  }
  
  return sum;
}

//_____________________________________________________________________________
#if 0
Float_t ChargeDep( const vector<Float_t>& amp )
{
  // Deconvolute signal given by samples in 'amp', return approximate integral.
  
  Float_t dummy;

  return ChargeDep( amp, dummy );
}
#endif

//_____________________________________________________________________________
Int_t GEMPlane::Decode( const THaEvData& evData )
{    
  // Extract this plane's hit data from the raw evData.
  //
  // This routine decodes the Gassiplex analog multiplexer readout data.
  // Finds clusters of active strips (=above software threshold) and
  // computes weighted average of position. Each such cluster makes one "Hit".

  //  static const char* const here = "Decode";
  //bool mc_data = fGEM->TestBit(GEM::kMCdata);

  typedef map<Int_t,Float_t> Hitmap_t;

  assert( fADC );
  assert( fADCped );
  assert( fTimeCentroid );

#ifdef TESTCODE
  if( TestBit(kDoHistos) )
    assert( fHitMap != 0 and fADCMap != 0 );
#endif
  assert( fPed.empty() or
	  fPed.size() == static_cast<Vflt_t::size_type>(fNelem) );
  assert( fMaxSamp == 1 or
	  fADCsamp.size() == static_cast<vector<Vflt_t>::size_type>(fNelem) );

  UInt_t nHits = 0;

  // Set up pedestal and noise corrections
  Double_t noisesum = 0.0;
  UInt_t   n_noise = 0;
  fDnoise = 0.0;

  bool do_pedestal_subtraction = !fPed.empty();
  bool do_noise_subtraction    = TestBit(kDoNoise);

  // Decode data
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module * d = fDetMap->GetModule(imod);

    // Read the active channels
    Int_t nchan = evData.GetNumChan( d->crate, d->slot );
    for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
      Int_t chan = evData.GetNextChan( d->crate, d->slot, ichan );
      if( chan < d->lo or chan > d->hi ) continue; // not part of this detector

      // Map channel number to strip number
      Int_t istrip = 
	MapChannel( d->first + ((d->reverse) ? d->hi - chan : chan - d->lo) );

      // For the APV25 analog pipeline, multiple "hits" on a decoder channel
      // correspond to time samples 25 ns apart
      Int_t nsamp = evData.GetNumHits( d->crate, d->slot, chan );
      assert( nsamp > 0 );
      ++fNhitStrips;
      nsamp = TMath::Min( nsamp, static_cast<Int_t>(fMaxSamp) );

      // Integrate the signal over time, if necessary
      Float_t adc, centroid;
      if( nsamp > 1 ) {
	Vflt_t& samples = fADCsamp[istrip];
	assert( samples.empty() );
	for( Int_t isamp = 0; isamp < nsamp; ++isamp ) {
	  if( evData.GetData(d->crate, d->slot, chan, isamp) > 0 ) {
	    break;
	  }
	  continue;  // Next ichan if all requested samples are zero
	}
	for( Int_t isamp = 0; isamp < nsamp; ++isamp ) {
	  Float_t fsamp = static_cast<Float_t>
	    ( evData.GetData(d->crate, d->slot, chan, isamp) );
	  samples.push_back( fsamp );
	}
 
	adc = ChargeDep( samples, centroid );

      } else {
	adc = static_cast<Float_t>
	  ( evData.GetData(d->crate, d->slot, chan, 0) );
	centroid = 0.0;
      }
      if( adc == 0.0 ) continue;
      ++fNsigStrips;

      // Save results for cluster finding later
      fADC[istrip] = adc;
      fTimeCentroid[istrip] = centroid;

      // Strip-by-strip pedestal subtraction
      if( do_pedestal_subtraction )
	adc -= fPed[istrip];
      if( adc < 0.0 )
	adc = 0.0;

      // Sum up ADCs that are likely not a hit
      // and are not suppressed through pedestal
      if( do_noise_subtraction and adc > 0.0 and adc < fMinAmpl ){
	noisesum += adc;
	n_noise++;
      }
      fADCped[istrip] = adc;

    }  // chans
  }    // modules

  fRawOcc = static_cast<Double_t>(fNsigStrips)/static_cast<Double_t>(fNelem);

  // Calculate average common-mode noise and subtract it from corrected 
  // ADC values, if requested
  if( do_noise_subtraction and n_noise > 0 ) {
    fDnoise = noisesum/static_cast<Double_t>(n_noise);
    assert( fDnoise >= 0.0 );
    if( fDnoise > 1000.0 )
      fDnoise = 0.0;
  }

  // Put corrected ADC data into sorted map. Fill histograms.
  Hitmap_t rawhits;
  for( Int_t i = 0; i < fNelem; i++ ) {
    if( do_noise_subtraction )
      fADCped[i] -= fDnoise;

    // Apply the ADC cut
    if( fADCped[i] >= fMinAmpl ) {
      // Sort the ADC data into a map, with the strip number as index
      ++fNRawHits;
      rawhits[i] = fADCped[i];
    } 
#ifdef TESTCODE
    if( TestBit(kDoHistos) ) {
      fHitMap->Fill(i);
      fADCMap->Fill(i, fADCped[i]);
    }
#endif
  }

  fOccupancy = static_cast<Double_t>(fNRawHits)/static_cast<Double_t>(fNelem);

  // Find and analyze clusters. Clusters of active strips are considered
  // a "Hit".
  //
  // The cluster analysis is a critical part of the GEM analysis. Various
  // things can and probably need to be done right here already: splitting 
  // oversized clusters, detecting noise hits/bogus clusters, detecting and
  // fitting overlapping clusters etc. 
  //
  // This analysis may even need to be re-done after preliminary tracking to
  // see if the clustering can be improved using candidate tracks.
  // Additionally, correlated amplitude information from a second readout
  // direction in the same readout plane could be used here. These advanced
  // procedures would require significant redesign of the code: 
  // all raw strip info will have to be saved and prcessed at a later point,
  // similar to the finding of hit pairs in like-oriented planes of the MWDC.
  //
  // For the moment, we implement a very simple algorithm: any cluster of 
  // strips larger than what a single cluster should be is assumed to be two or
  // more overlapping hits, and the cluster will be split as follows: anything
  // that looks like a local peak followed by a valley will be considered an
  // actual cluster. The parameter frac = fSplitFrac (0.0 ... 1.0) can
  // be used for some crude tuning. frac > 0.0 means that a peak is
  // only a peak if the amplitude drops below (1-frac), so
  // frac = 0.1 means: trigger on a drop below 90% etc. Likewise for the
  // following valley: the bottom is found if the amplitude rises again
  // by (1+frac), so frac = 0.1 means: trigger on a rise above 110% etc.
  //
  Double_t frac_down = 1.0 - fSplitFrac, frac_up = 1.0 + fSplitFrac;
#ifndef NDEBUG
  Hit* prevHit = 0;
#endif
  Hitmap_t::iterator next = rawhits.begin();
  while( next != rawhits.end() ) {
    Hitmap_t::iterator start = next, cur = next;
    ++next;
    assert( next == rawhits.end() or ((*next).first > (*cur).first) );
    while( next != rawhits.end() and ((*next).first - (*cur).first == 1)  ) {
      ++cur;
      ++next;
    }
    // Now the cluster candidate is between start and cur
    assert( (*cur).first >= (*start).first );
    Int_t  type = 0;
    UInt_t size = (*cur).first - (*start).first + 1;
    if( size > fMaxClusterSize ) {
      Double_t maxadc = 0.0, minadc = kBig;
      Hitmap_t::iterator it = start, maxpos = start, minpos = start;
      enum EStep { kFindMax = 1, kFindMin, kDone };
      EStep step = kFindMax;
      while( step != kDone and it != next ) {
        Double_t adc = (*it).second;
        switch( step ) {
          case kFindMax:
            // Looking for maximum
            if( adc > maxadc ) {
              maxpos = it;
              maxadc = adc;
            } else if( adc < maxadc * frac_down ) {
              assert( maxadc > 0.0 );
              step = kFindMin;
              continue;
            }
            break;
          case kFindMin:
            // Looking for minimum
            if( adc < minadc ) {
              minpos = it;
              minadc = adc;
            } else if( adc > minadc * frac_up ) {
              assert( minadc < kBig );
              step = kDone;
            }
            break;
          case kDone:
            assert( false );  // should never get here
            break;
        }
        ++it;
      }
      if( step == kDone ) {
        // Found maximum followed by minimum
        assert( minpos != start );
        assert( minpos != cur );
        assert( (*minpos).first > (*maxpos).first );
        // Split the cluster at the position of the minimum, assuming that
        // the strip with the minimum amplitude is shared between both clusters
        cur  = minpos;
        next = minpos;
        // In order not to double-count amplitude, we split the signal height
        // of that strip evenly between the two clusters. This is a very
        // crude way of doing what we really should be doing: "fitting" a peak 
        // shape and using the area and centroid of the curve
        (*minpos).second /= 2.0;
      }
      type = step;
      size = (*cur).first - (*start).first + 1;
      assert( (*cur).first >= (*start).first );
    }
    // Compute weighted position average. Again, a crude (but fast) substitute
    // for fitting the centroid of the peak.
    Double_t xsum = 0.0, adcsum = 0.0;
    for( ; start != next; ++start ) {
      Double_t pos = GetStart() + (*start).first * GetPitch();
      Double_t adc = (*start).second; 
      xsum   += pos * adc;
      adcsum += adc;
    }
    assert( adcsum > 0.0 );
    // Make a new hit
    //
    // The "type" parameter indicates the result of the custer analysis:
    // 0: clean (i.e. smaller than fMaxClusterSize, no further analysis)
    // 1: large, maximum at right edge, not split
    // 2: large, no clear minimum on the right side found, not split
    // 3: split, well-defined peak found (may still be larger than maxsize)
    //
    // The resolution (sigma) of the position measurement depends on the 
    // cluster size. In particular, if the cluster consists of only a single
    // hit, the resolution is much reduced
    Double_t resolution = fResolution;
    if( size == 1 ) {
      resolution = 0.25*GetPitch();
      // The factor of 1/2*pitch is just a guess. Since with real GEMs
      // there _should_ always be more than one strip per cluster, we must
      // assume that the other strip(s) did not fire due to inefficiency.
      // As a result, the error is bigger than it would be if only ever one
      // strip fired per hit.
//       resolution = TMath::Max( 0.5*GetStripPitch(), 2.0*fResolution );
//     } else if( size == 2 ) {
//       // Again, this is a guess, to be quantified with Monte Carlo
//       resolution = 1.2*fResolution;
    }
    Double_t pos = xsum/adcsum;
#ifndef NDEBUG
    Hit* theHit = 
#endif
      new( (*fHits)[nHits++] ) GEMHit( pos,
				       adcsum,
				       size,
				       type,
				       resolution,
				       this
				       );
#ifndef NDEBUG
    // Ensure hits are ordered by position (should be guaranteed by std::map)
    assert( prevHit == 0 or theHit->Compare(prevHit) > 0 );
    prevHit = theHit;
#endif
  }

  // Negative return value indicates potential problem
  if( nHits > fMaxHits )
    nHits = -nHits;

  return nHits;
}

//_____________________________________________________________________________
Int_t GEMPlane::DefineVariables( EMode mode )
{
  // initialize global variables

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { "nstrips",     "Num strips with hits > adc.min",   "fNRawHits" },
    { "rawocc",      "strips w/data / n_all_strips",     "fRawOcc" },
    { "occupancy",   "nstrips / n_all_strips",           "fOccupancy" },
    { "strip.adc",   "Raw strip ADC values",             "fADC" },
    { "strip.adc_p", "Pedstal sub strip ADC values",     "fADCped" },
    { "strip.time",  "Centroid of strip signal (ns)",    "fTimeCentroid" },
    { "nhits",       "Num hits (clusters of strips)",    "GetNhits()" },
    { "hit.pos",     "Hit centroid (m)",    "fHits.TreeSearch::Hit.fPos" },
    { "hit.adc",     "Hit ADC sum",         "fHits.TreeSearch::Hit.fADCsum" },
    { "hit.size",    "Num strips ",         "fHits.TreeSearch::Hit.fSize" },
    { "hit.type",    "Hit analysis result", "fHits.TreeSearch::Hit.fType" },
    { "noise",       "Noise",               "fDnoise" },
    { "ncoords",     "Num fit coords",     "GetNcoords()" },
    { "coord.rank",  "Fit rank of coord",
      "fFitCoords.TreeSearch::FitCoord.fFitRank" },
    { "coord.pos",   "Position used in fit (m)",
      "fFitCoords.TreeSearch::FitCoord.fPos" },
    { "coord.trkpos","Track pos from projection fit (m)",
      "fFitCoords.TreeSearch::FitCoord.fTrackPos" },
    { "coord.trkslope","Track slope from projection fit",
      "fFitCoords.TreeSearch::FitCoord.fTrackSlope" },
    { "coord.resid", "Residual of trkpos (m)",
      "fFitCoords.TreeSearch::FitCoord.GetResidual()" },
    { "coord.3Dpos",  "Crossing position of fitted 3D track (m)",
      "fFitCoords.TreeSearch::FitCoord.f3DTrkPos" },
    { "coord.3Dresid","Residual of 3D trkpos (m)",
      "fFitCoords.TreeSearch::FitCoord.Get3DTrkResid()" },
    { "coord.3Dslope","Slope of fitted 3D track wrt projection",
      "fFitCoords.TreeSearch::FitCoord.f3DTrkSlope" },
    { 0 }
  };
  Int_t ret = DefineVarsFromList( vars, mode );

  return ret;
}

//_____________________________________________________________________________
Int_t GEMPlane::End( THaRunBase* run )
{
  // Write diagnostic histograms to file

  Plane::End( run );

#ifdef TESTCODE
  if( TestBit(kDoHistos) ) {
    if( fADCMap ) fADCMap->Write();
  }
#endif
  return 0;
}

//_____________________________________________________________________________
Int_t GEMPlane::ReadDatabase( const TDatime& date )
{
  // Read database

  static const char* const here = "ReadDatabase";

  // Read the info needed for the base class, quit if error
  Int_t status = Plane::ReadDatabase(date);
  fIsInit = false;
  if( status != kOK )
    return status;

  // Now read additional info needed for this derived class
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Set defaults
  TString mapping;
  Int_t do_noise = 1;
  fMaxClusterSize = kMaxUInt;
  fMinAmpl   = 0.0;
  fSplitFrac = 0.0;
  fMapType   = kOneToOne;
  fMaxSamp   = 1;
  fChanMap.clear();
  fPed.clear();
  fADCsamp.clear();
  fAmplSigma = 0.36; // default, an educated guess

  if( status == kOK ) {
    try {
      const DBRequest request[] = {
	//TODO: make nstrips/pos/pitch optional if already read in base class
	{ "nstrips",        &fNelem,          kInt,     0, 0, -1 },
	{ "strip.pos",      &fStart },
	{ "strip.pitch",    &fPitch,          kDouble,  0, 0, -1 },
	{ "maxclustsiz",    &fMaxClusterSize, kUInt,    0, 1, -1 },
	{ "maxsamp",        &fMaxSamp,        kUInt,    0, 1, -1 },
	{ "adc.min",        &fMinAmpl,        kDouble,  0, 1, -1 },
	{ "split.frac",     &fSplitFrac,      kDouble,  0, 1, -1 },
	{ "mapping",        &mapping,         kTString, 0, 1, -1 },
	{ "chanmap",        &fChanMap,        kIntV,    0, 1, -1 },
	{ "pedestal",       &fPed,            kFloatV,  0, 1 },
	{ "do_noise",       &do_noise,        kInt,     0, 1, -1 },
	{ "adc.sigma",      &fAmplSigma,      kDouble,  0, 1, -1 },
	{ 0 }
      };
      status = LoadDB( file, date, request, fPrefix );
    }
    // Catch exceptions here so that we can close the file
    catch(...) {
      fclose(file);
      throw;
    }
  }

  // Finished reading the database
  fclose(file);
  if( status != kOK )
    return status;

  // Sanity checks
  if( fNelem <= 0 or fNelem > 65536 ) { // arbitray upper limit
    Error( Here(here), "Invalid number of strips: %d", fNelem );
    return kInitError;
  }

  Int_t nchan = fDetMap->GetTotNumChan();
  if( nchan != fNelem ) {
    Error( Here(here), "Number of detector map channels (%d) "
        "disagrees with number of strips (%d)", nchan, fNelem );
    return kInitError;
  }

  if( do_noise )
    SetBit( kDoNoise );

  delete fADC; fADC = 0;
  delete fADCped; fADCped = 0;
  delete fTimeCentroid; fTimeCentroid = 0;

  fADC = new Float_t[fNelem];
  fADCped = new Float_t[fNelem];
  fTimeCentroid = new Float_t[fNelem];

  TString::ECaseCompare cmp = TString::kIgnoreCase;
  if( !mapping.IsNull() ) {
    if( mapping.Length() >= 3 and
        TString("one-to-one").BeginsWith(mapping,cmp) )
      fMapType = kOneToOne;
    else if( mapping.Length() >=3 and
        TString("reverse").BeginsWith(mapping,cmp) )
      fMapType = kReverse;
    else if( mapping.Length() >= 5 and
        mapping.BeginsWith(TString("gassiplex-adapter"),cmp) ) {
      if( fNelem > 240 ) {
        Error( Here(here), "Gassiplex adapter mapping allows at most 240 "
            "strips, but %d configured. Fix database.", fNelem );
        return kInitError;
      }
      if( fNelem < 240 ) {
        Warning( Here(here), "Gassiplex adapter mapping expects 240 "
            "strips, but %d configured. Database may be misconfigured "
            "(or you know what you are doing).", fNelem );
      }
      if( mapping.BeginsWith(TString("gassiplex-adapter-2"),cmp) ) {
	fMapType = kGassiplexAdapter2;
      } else {
	fMapType = kGassiplexAdapter1;
      }
    }
    else if( TString("table").CompareTo(mapping,cmp) ) {
      if( fChanMap.empty() ) {
        Error( Here(here), "Channel mapping table requested, but no map "
            "defined. Specify chanmap in database." );
        return kInitError;
      }
      if( fChanMap.size() != static_cast<UInt_t>(fNelem) ) {
        Error( Here(here), "Number of channel map entries (%u) msut equal "
            "number of strips (%d). Fix database.", 
	       static_cast<unsigned int>(fChanMap.size()), fNelem );
        return kInitError;
      }
      // check if entries in channel map are within range
      for( vector<Int_t>::const_iterator it = fChanMap.begin();
          it != fChanMap.end(); ++it ) {
        if( (*it) < 0 or (*it) >= fNelem ) {
          Error( Here(here), "Illegal chanmap entry: %d. Must be >= 0 and "
              "< %d. Fix database.", (*it), fNelem );
          return kInitError;
        }
      }
      fMapType = kTable;
    } else {
      Error( Here(here), "Unknown channel mapping type %s. Fix database.",
          mapping.Data() );
      return kInitError;
    }

  } else
    fChanMap.clear();

  if( !fPed.empty() and fPed.size() != static_cast<UInt_t>(fNelem) ) {
    Error( Here(here), "Size of pedestal array (%u) must equal "
	   "number of strips (%d). Fix database.",
	   static_cast<unsigned int>(fPed.size()), fNelem );
    return kInitError;
  }

  // Sanity checks on fMaxSamp
  static const UInt_t max_maxsamp = 32; // arbitrary sanity limit
  if( fMaxSamp == 0 )
    fMaxSamp = 1;
  else if( fMaxSamp > max_maxsamp ) {  
    Warning( Here(here), "Illegal maximum number of samples: %u. "
	     "Adjusted to maximum allowed = %u.", fMaxSamp, max_maxsamp );
    fMaxSamp = max_maxsamp;
  }
  if( fMaxSamp > 1 ) {
    // NB: this may very well allocate O(1000) vectors of O(10) floats
    try {
      fADCsamp.resize( fNelem );
      for( vector<Vflt_t>::iterator it = fADCsamp.begin();
	   it != fADCsamp.end(); ++it ) {
	(*it).reserve( fMaxSamp );
      }
    }
    catch( bad_alloc ) {
      Error( Here(here), "Out of memory trying to allocate ADC sample "
	     "buffers. Reduce maxsamp." );
      return kInitError;
    }
  }

  if( fAmplSigma < 0.0 ) {
    Warning( Here(here), "Negative adc.sigma = %lf makes no sense. Adjusted "
	     "to positive.", fAmplSigma );
    fAmplSigma = TMath::Abs( fAmplSigma );
  }
  if( fAmplSigma < 0.001 ) {
    Warning( Here(here), "adc.sigma = %lf is extremely small. "
	     "Double-check database.", fAmplSigma );
  }
	       
  fIsInit = true;
  return kOK;
}
//_____________________________________________________________________________
void GEMPlane::Print( Option_t* ) const
{    
  // Print plane info

  //TODO: add detailed info
  cout << "GEMPlane:  #" << GetPlaneNum() << " "
    << GetName()   << "\t"
    //       << GetTitle()        << "\t"
    << fNelem << " strips\t"
    << "z = " << GetZ();
  if( fPartner ) {
    cout << "\t partner = " 
      << fPartner->GetName();
    //	 << fPartner->GetTitle();
  }
  cout << endl;
}

//_____________________________________________________________________________
Double_t GEMPlane::GetAmplSigma( Double_t ampl ) const
{ 
  // Return expected sigma of the hit amplitude distribution at the given
  // amplitude 'ampl'.
  // For starters, calculate a constant relative uncertainty
  // from a database parameter.

  return fAmplSigma * ampl;
}

//_____________________________________________________________________________

}

ClassImp(TreeSearch::GEMPlane)

///////////////////////////////////////////////////////////////////////////////

