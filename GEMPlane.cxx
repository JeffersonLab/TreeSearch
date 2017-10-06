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
#include "GEMTracker.h"
#include "Projection.h"

#include "THaDetMap.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TError.h"
#include "TString.h"

#ifndef MCDATA
#include "THaEvData.h"
#else
#include "SimDecoder.h"
#endif

#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>

using namespace std;
using namespace Podd;

#define ALL(c) (c).begin(), (c).end()

namespace TreeSearch {

struct StripData_t {
  Float_t adcraw;
  Float_t adc;
  Float_t time;
  Bool_t  pass;
  StripData_t() {}
  StripData_t( Float_t _raw, Float_t _adc, Float_t _time, Bool_t _pass )
    : adcraw(_raw), adc(_adc), time(_time), pass(_pass) {}
};

// Sanity limit on number of channels (strips)
static const Int_t kMaxNChan = 20000;

//_____________________________________________________________________________
GEMPlane::GEMPlane( const char* name, const char* description,
		    THaDetectorBase* parent )
  : Plane(name,description,parent),
    fMapType(kOneToOne), fMaxClusterSize(0), fMinAmpl(0), fSplitFrac(0),
    fMaxSamp(1), fTimeCutCV(0), fTimeCutHW(0),
    fAmplSigma(0), fADCraw(0), fADC(0), fHitTime(0), fADCcor(0),
    fGoodHit(0), fDnoise(0), fNrawStrips(0), fNhitStrips(0), fHitOcc(0),
    fOccupancy(0), fADCMap(0), fHitStripTime(0), fADCsamp_sig(0), fADCsampVsStripTime_bkgd(0)
{
  // Constructor

  static const char* const here = "GEMPlane";

  assert( dynamic_cast<GEMTracker*>(fTracker) );

  try {
#ifdef MCDATA
    if( fTracker->TestBit(Tracker::kMCdata) ) // Monte Carlo data mode?
      fHits = new TClonesArray("TreeSearch::MCGEMHit", 200);
    else
#endif
      fHits = new TClonesArray("TreeSearch::GEMHit", 200);
  }
  catch( std::bad_alloc ) {
    Error( Here(here), "Out of memory allocating hit array for readout "
	   "plane %s. Call expert.", name );
    MakeZombie();
    return;
  }
}

//_____________________________________________________________________________
GEMPlane::~GEMPlane()
{
  // Destructor.

  // Histograms in fHist should be automatically deleted by ROOT when
  // the output file is closed

  if( fIsSetup )
    RemoveVariables();

  // fHits deleted in base class
  delete fGoodHit;
  delete fADCcor;
  delete fHitTime;
  delete fADC;
  delete fADCraw;
}

//_____________________________________________________________________________
Int_t GEMPlane::Begin( THaRunBase* run )
{
  // Book diagnostic histos if not already done

  Plane::Begin( run );

#ifdef TESTCODE
  if( TestBit(kDoHistos) ) {
    //TODO: check if exist
    string hname(fPrefix);
    hname.append("adcmap");
    fADCMap = new TH1F( hname.c_str(), hname.c_str(), fNelem, 0, fNelem );
    string hname_2D(fPrefix);
    hname_2D.append("HitVsStripTime");
    fHitStripTime = new TH2F( hname_2D.c_str(), hname_2D.c_str(), 350, -250, +100, 350, -250, +100);
    hname = fPrefix;
    hname.append("ADCsamp_sig");
    fADCsamp_sig = new TH1F( hname.c_str(), hname.c_str(), fMaxSamp, 0, fMaxSamp);
    hname_2D = fPrefix;
    hname_2D.append("ADCsampVsStripTime_bkgd");
    fADCsampVsStripTime_bkgd = new TH2F( hname_2D.c_str(), hname_2D.c_str(), 14, -250, +100, fMaxSamp, 0, fMaxSamp);
   
  }
#endif
  return 0;
}

//_____________________________________________________________________________
void GEMPlane::Clear( Option_t* opt )
{
  // Clear event-by-event data (hits)

  Plane::Clear(opt);

  if( !IsDummy() ) {
    assert( fADCraw and fADC and fHitTime and fADCcor and fGoodHit );
    memset( fADCraw, 0, fNelem*sizeof(Float_t) );
    memset( fADC, 0, fNelem*sizeof(Float_t) );
    memset( fHitTime, 0, fNelem*sizeof(Float_t) );
    memset( fADCcor, 0, fNelem*sizeof(Float_t) );
    memset( fGoodHit, 0, fNelem*sizeof(Byte_t) );
    fSigStrips.clear();
    fStripsSeen.assign( fNelem, false );
  }

  fNhitStrips = fNrawStrips = 0;
  fHitOcc = fOccupancy = fDnoise = 0.0;
}

//_____________________________________________________________________________
Int_t GEMPlane::MapChannel( Int_t idx ) const
{
  // Map hardware channel number to logical strip number based on mapping
  // prescription from database

  assert( idx >= 0 and idx < fNelem );

  Int_t ret = 0;
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
    assert( fChanMap.size() == static_cast<Vint_t::size_type>(fNelem) );
    ret = fChanMap[idx];
    break;
  }
  assert( ret >= 0 and ret < fNelem );
  return ret;
}

//_____________________________________________________________________________
static StripData_t ChargeDep( const vector<Float_t>& amp )
{
  // Deconvolute signal given by samples in 'amp', return approximate integral.
  // Currently analyzes exactly 3 samples.
  // From Kalyan Allada
  // NIM A326, 112 (1993)

  //FIXME: from database, proper value for Tp
  const Float_t delta_t = 25.0; // time interval between samples (ns)
  const Float_t Tp      = 50.0; // RC filter time constant (ns)

  assert( amp.size() >= 3 );
  Float_t adcraw = delta_t*(amp[0]+amp[1]+amp[2]);
  // Weight factors calculated based on the response of the silicon microstrip
  // detector:
  // v(t) = (delta_t/Tp)*exp(-delta_t/Tp)
  // Need to update this for GEM detector response(?):
  // v(t) = A*(1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2)
  // where A is the amplitude, t0 the begin of the rise, tau1 the time
  // parameter for the rising edge and tau2 the for the falling edge.

  Float_t x = delta_t/Tp;

  Float_t w1 = TMath::Exp(x-1)/x;
  Float_t w2 = -2*TMath::Exp(-1)/x;
  Float_t w3 = TMath::Exp(-x-1)/x;

  // Deconvoluted signal samples, assuming measurements of zero before the
  // leading edge
  Float_t sig[3] = { amp[0]*w1,
		     amp[1]*w1+amp[0]*w2,
		     amp[2]*w1+amp[1]*w2+amp[0]*w3 };

  Float_t adc    = delta_t*(sig[0]+sig[1]+sig[2]);
  Float_t time   = 0;     // TODO
  if(amp.size()>=6){
    //search max: temporary!!!
    int i_sampmax = 0;
    double max = 0;
    for(int i = 0; i<6; i++){
      if(amp[i]>max){
	max = amp[i];
	i_sampmax = i;
      }
    }
    time = 25.0*(i_sampmax+0.5);
  }
  Bool_t pass;
  // Calculate ratios for 3 samples and check for bad signals
  if( amp[2] > 0 ) {
    Float_t r1 = amp[0]/amp[2];
    Float_t r2 = amp[1]/amp[2];
    pass = (r1 < 1.0 and r2 < 1.0 and r1 < r2);
  } else
    pass = false;
  //printf("adcraw = %1.2f, sig[0, 1, 2] = %1.2f, %1.2f, %1.2f, adc = %1.2f \n", adcraw, sig[0], sig[1], sig[2], adc);
  return StripData_t(adcraw,adc,time,pass);
}

//_____________________________________________________________________________
void GEMPlane::AddStrip( Int_t istrip )
{
  // Record a hit on the given strip number in internal arrays.
  // Utility function used by Decode.

  Float_t adc = fADCcor[istrip];
  if( adc > 0 )
    ++fNhitStrips;
  if( fGoodHit[istrip] and adc >= fMinAmpl) {
    if(fTimeCutHW==0 || (fTimeCutCV-fTimeCutHW<=fHitTime[istrip] && fHitTime[istrip]<=fTimeCutCV+fTimeCutHW) )
      fSigStrips.push_back(istrip);
  }
#ifdef TESTCODE
  if( TestBit(kDoHistos) ) {
    fHitMap->Fill(istrip);
    fADCMap->Fill(istrip, adc);
  }
#endif
}

//_____________________________________________________________________________
Int_t GEMPlane::GEMDecode( const THaEvData& evData )
{
  // Extract this plane's hit data from the raw evData.
  //
  // This routine decodes the front-end readout data.
  // Finds clusters of active strips (=above software threshold) and
  // computes weighted average of position. Each such cluster makes one "Hit".

  const char* const here = "GEMPlane::Decode";

#ifdef MCDATA
  bool mc_data = fTracker->TestBit(Tracker::kMCdata);
  assert( !mc_data || dynamic_cast<const SimDecoder*>(&evData) != 0 );
#endif
  assert( fADCraw and fADC and fADCcor and fHitTime );

#ifdef TESTCODE
  if( TestBit(kDoHistos) )
    assert( fHitMap != 0 and fADCMap != 0 and fHitStripTime != 0 and fADCsamp_sig != 0 and fADCsampVsStripTime_bkgd != 0);
#endif
  assert( fPed.empty() or
	  fPed.size() == static_cast<Vflt_t::size_type>(fNelem) );
  assert( fSigStrips.empty() );
  assert( fStripsSeen.size() == static_cast<Vbool_t::size_type>(fNelem) );

  UInt_t nHits = 0;

  // Set up pedestal and noise corrections
  Double_t noisesum = 0.0;
  UInt_t   n_noise = 0;

  bool do_pedestal_subtraction = !fPed.empty();
  bool do_noise_subtraction    = TestBit(kDoNoise);

#ifdef MCDATA
  const SimDecoder* simdata = 0;
  if( mc_data ) {
    assert( dynamic_cast<const SimDecoder*>(&evData) );
    simdata = static_cast<const SimDecoder*>(&evData);
  }
#endif
  Vflt_t samples;
  if( fMaxSamp > 1 )
    samples.reserve(fMaxSamp);

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
      // Test for duplicate istrip, if found, warn and skip
      assert( istrip >= 0 and istrip < fNelem );
      if( fStripsSeen[istrip] ) {
	const char* inp_source = "DAQ";
#ifdef MCDATA
	if( mc_data )
	  inp_source = "digitization";
#endif
	Warning( Here(here), "Duplicate strip number %d in plane %s, event %d. "
		 "Ignorning it. Fix your %s.",
		 istrip, GetName(), evData.GetEvNum(), inp_source );
	continue;
      }
      fStripsSeen[istrip] = true;

      // For the APV25 analog pipeline, multiple "hits" on a decoder channel
      // correspond to time samples 25 ns apart
      Int_t nsamp = evData.GetNumHits( d->crate, d->slot, chan );
      // printf(" crate = %d, slot = %d, chan = %d \n", d->crate, d->slot, chan );
      // for( Int_t isamp = 0; isamp < nsamp; ++isamp ) {
      // 	Float_t fsamp = static_cast<Float_t>
      // 	  ( evData.GetData(d->crate, d->slot, chan, isamp) );
      // 	printf( " %1.0f ", fsamp );
      // }
      // printf("\n");
      assert( nsamp > 0 );
      ++fNrawStrips;
      nsamp = TMath::Min( nsamp, static_cast<Int_t>(fMaxSamp) );

      // Integrate the signal over time and analyze pulse shape
      StripData_t stripdata;
      if( nsamp > 1 ) {
	samples.clear();
	for( Int_t isamp = 0; isamp < nsamp; ++isamp ) {
	  Float_t fsamp = static_cast<Float_t>
	    ( evData.GetData(d->crate, d->slot, chan, isamp) );
	  samples.push_back( fsamp );
	  //printf( " %1.0f ", fsamp );
	}
	//printf("\n");
	// Analyze the pulse shape
	stripdata = ChargeDep(samples);
      }
      else {
	stripdata.adcraw = stripdata.adc =
	  static_cast<Float_t>( evData.GetData(d->crate, d->slot, chan, 0) );
	stripdata.time = 0;
	stripdata.pass = true;
      }
      
      // Skip null data
      if( stripdata.adcraw == 0 )
	continue;
      //printf(" crate %d, slot %d, chan %d \n", d->crate, d->slot, ichan);
      
      // Save results for cluster finding later
      fADCraw[istrip]  = stripdata.adcraw;
      fADC[istrip]     = stripdata.adc;
      fHitTime[istrip] = stripdata.time;

      // Strip-by-strip pedestal subtraction
      Float_t adc = stripdata.adc;
      if( do_pedestal_subtraction )
	adc -= fPed[istrip];

      fADCcor[istrip] = adc;
      fGoodHit[istrip] = not TestBit(kCheckPulseShape) or stripdata.pass;

      if( do_noise_subtraction ) {
	// Sum up ADCs that are likely not a hit
	if( adc < fMinAmpl ) {
	  noisesum += adc;
	  n_noise++;
	}
      }
      else {
	// If no noise subtraction is done, then we can finish up with this
	// strip number right here. Otherwise we need a second iteration below
	AddStrip( istrip );
      }

#ifdef MCDATA
      // If doing MC data, save the truth information for each strip
      if( mc_data ) {
	fMCHitList.push_back(istrip);
  	fMCHitInfo[istrip] = simdata->GetMCHitInfo(d->crate,d->slot,chan);
	//For the time being, adding the dtrip MC time.
	//fHitTime[istrip] = fMCHitInfo[istrip].fMCTime;
	
	if(fMCHitInfo[istrip].fMCTrack > 0){
	  for(int i_ = 0; i_<samples.size(); i_++)
	    fADCsamp_sig->Fill(i_, samples.at(i_));
	}else{
	  for(int i_ = 0; i_<samples.size(); i_++)
	    fADCsampVsStripTime_bkgd->Fill(fHitTime[istrip], i_, samples.at(i_));
	}
      }
#endif
    }  // chans
  }    // modules

  // Calculate average common-mode noise and subtract it from corrected
  // ADC values, if requested
  if( do_noise_subtraction ) {
    if ( n_noise > 0 ) {
      fDnoise = noisesum/n_noise;
      assert( fDnoise < fMinAmpl );
    }
    // Save strip numbers of corrected ADC data above threshold. Fill histograms.
    assert( fSigStrips.empty() );
    for( Int_t i = 0; i < fNelem; i++ ) {
      fADCcor[i] -= fDnoise;
      AddStrip( i );
    }
  }

  fHitOcc    = static_cast<Double_t>(fNhitStrips) / fNelem;
  fOccupancy = static_cast<Double_t>(GetNsigStrips()) / fNelem;

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

  // The active strip numbers must be sorted for the clustering algorithm
  sort( ALL(fSigStrips) );

  Double_t frac_down = 1.0 - fSplitFrac, frac_up = 1.0 + fSplitFrac;
#ifndef NDEBUG
  GEMHit* prevHit = 0;
#endif
  typedef Vint_t::iterator viter_t;
  Vint_t splits;  // Strips with ampl split between 2 clusters
  viter_t next = fSigStrips.begin();
  while( next != fSigStrips.end() ) {
    viter_t start = next, cur = next;
    ++next;
    assert( next == fSigStrips.end() or *next > *cur );
    while( next != fSigStrips.end() and (*next - *cur) == 1  ) {
      ++cur;
      ++next;
    }
    // Now the cluster candidate is between start and cur
    assert( *cur >= *start );
    // The "type" parameter indicates the result of the cluster analysis:
    // 0: clean (i.e. smaller than fMaxClusterSize, no further analysis)
    // 1: large, maximum at right edge, not split
    // 2: large, no clear minimum on the right side found, not split
    // 3: split, well-defined peak found (may still be larger than maxsize)
    Int_t  type = 0;
    UInt_t size = *cur - *start + 1;
    if( size > fMaxClusterSize ) {
      Double_t maxadc = 0.0, minadc = kBig;
      viter_t it = start, maxpos = start, minpos = start;
      enum EStep { kFindMax = 1, kFindMin, kDone };
      EStep step = kFindMax;
      while( step != kDone and it != next ) {
        Double_t adc = fADCcor[*it];
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
        assert( *minpos > *maxpos );
        // Split the cluster at the position of the minimum, assuming that
        // the strip with the minimum amplitude is shared between both clusters
        cur  = minpos;
        next = minpos;
        // In order not to double-count amplitude, we split the signal height
        // of that strip evenly between the two clusters. This is a very
        // crude way of doing what we really should be doing: "fitting" a peak
        // shape and using the area and centroid of the curve
	fADCcor[*minpos] /= 2.0;
	splits.push_back(*minpos);
      }
      type = step;
      size = *cur - *start + 1;
      assert( *cur >= *start );
    }
    assert( size > 0 );
    // Compute weighted position average. Again, a crude (but fast) substitute
    // for fitting the centroid of the peak.
    //printf("cluster number %d \n", *cur);
    Double_t xsum = 0.0, adcsum = 0.0, time = kBig;
    double tsum = 0.0;// time calculation by weighted sum
#ifdef MCDATA
    Double_t mcpos = 0.0, mctime = kBig;
    Int_t mctrack = 0, num_bg = 0;
#endif

#ifdef TESTCODE
    if( TestBit(kDoHistos) ){
      StripsTime.clear();
    }
#endif
    for( ; start != next; ++start ) {
      Int_t istrip = *start;
      Double_t pos = GetStart() + istrip * GetPitch();
      Double_t adc = fADCcor[istrip];
      xsum   += pos * adc;
      adcsum += adc;
      //time  = TMath::Min( time, fHitTime[istrip] ); //same method as the MC time.
      tsum   += fHitTime[istrip] * adc; // time calculation by weighted sum
#ifdef MCDATA
      // If doing MC data, analyze the strip truth information
      if( mc_data ) {
	MCHitInfo& mc = fMCHitInfo[istrip];
	// This may be smaller than the actual total number of background hits
	// contributing to the entire cluster, but counting them would involve
	// lists of secondary particle numbers ... overkill for now
	num_bg = TMath::Max( num_bg, mc.fContam );
	// All primary particle hits in the cluster are from the same track
	assert( mctrack == 0 || mc.fMCTrack == 0 || mctrack == mc.fMCTrack );
	if( mctrack == 0 ) {
	  if( mc.fMCTrack > 0 ) {
	    // If the cluster contains a signal hit, save its info and be done
	    //printf("istrip = %d, pos = %1.3f, ADC = %1.3f\n", istrip, pos, adc);
	    //printf("mc.fMCTtime = %1.3f\n", mc.fMCTime);
	    mctrack = mc.fMCTrack;
	    mcpos   = mc.fMCPos;
	    mctime  = mc.fMCTime;
#ifdef TESTCODE
	    if( TestBit(kDoHistos) ){
	      StripsTime.push_back(mc.fMCTime);
	    }
#endif
	  }
	  else {
	    // If background hits only, compute position average
	    mcpos  += mc.fMCPos;
	    mctime  = TMath::Min( mctime, mc.fMCTime );
#ifdef TESTCODE
	    if( TestBit(kDoHistos) ){
	      StripsTime.push_back(mc.fMCTime);
	    }
#endif
	    // if(mc.fMCTime>=50.0){
	    //   printf("istrip = %d, mc.fMCTime =  %f\n", istrip, mc.fMCTime);
	    //   printf("mctime = %1.3f, mc.fMCTime = %1.3f\n", mctime, mc.fMCTime);
	    // }
	    // if(-50.0<= mc.fMCTime && mc.fMCTime<=25.0){
	    //   mctime = -500.0;
	    // // if(mc.fMCPos==0.0){
	    // //   printf("istrip = %d, mc.fMCTime =  %f\n", istrip, mc.fMCTime);
	    // //   //printf("mctime = %1.3f, mc.fMCTime = %1.3f\n", mctime, mc.fMCTime);
	    // //   //printf("mcpos = %1.3f, mc.fMCpos = %1.3f\n", mcpos, mc.fMCPos);
	    // }
	  }
	}
      }
#endif // MCDATA
    }
    
#ifdef TESTCODE
    if( TestBit(kDoHistos) ) {
      for(int i_ = 0; i_<StripsTime.size(); i_++){
	fHitStripTime->Fill(mctime, StripsTime.at(i_));
      }
      StripsTime.clear();
    }
#endif
    
    
    assert( adcsum > 0.0 );
    Double_t pos = xsum/adcsum;
    time = tsum/adcsum;
    
#ifdef MCDATA
    if( mc_data && mctrack == 0 ) {
      mcpos /= static_cast<Double_t>(size);
    }
#endif
    // The resolution (sigma) of the position measurement depends on the
    // cluster size. In particular, if the cluster consists of only a single
    // hit, the resolution is much reduced
    Double_t resolution = fResolution;
    if( size == 1 ) {
      resolution = TMath::Max( 0.25*GetPitch(), fResolution );
      // The factor of 1/2*pitch is just a guess. Since with real GEMs
      // there _should_ always be more than one strip per cluster, we must
      // assume that the other strip(s) did not fire due to inefficiency.
      // As a result, the error is bigger than it would be if only ever one
      // strip fired per hit.
//       resolution = TMath::Max( 0.5*GetPitch(), 2.0*fResolution );
//     } else if( size == 2 ) {
//       // Again, this is a guess, to be quantified with Monte Carlo
//       resolution = 1.2*fResolution;
    }
    
    //if(fTimeCutHW==0 || (fTimeCutCV-fTimeCutHW<=time && time<=fTimeCutCV+fTimeCutHW) ){
    // Make a new hit
#ifndef NDEBUG
    GEMHit* theHit = 0;
#endif
#ifdef MCDATA
    if( !mc_data ) {
#endif
#ifndef NDEBUG
      theHit =
#endif
	new( (*fHits)[nHits++] ) GEMHit( pos,
					 time,
					 adcsum,
					 size,
					 type,
					 resolution,
					 this
					 );

#ifdef MCDATA
    } else {
      // Monte Carlo data
#ifndef NDEBUG
      theHit =
#endif
	new( (*fHits)[nHits++] ) MCGEMHit( pos,
					   time,
					   adcsum,
					   size,
					   type,
					   resolution,
					   this,
					   mctrack,
					   mcpos,
					   mctime,
					   num_bg
					   );
      //printf("hit pos = %1.3f, time = %1.3f, size = %d\n", pos, mctime, size);
    }
#endif // MCDATA
#ifndef NDEBUG
    // Ensure hits are ordered by position (should be guaranteed by std::map)
    assert( prevHit == 0 or theHit->Compare(prevHit) > 0 );
    prevHit = theHit;
#endif
    //}//end timing cut 
  }

  // Undo amplitude splitting, if any, so fADCcor contains correct ADC values
  for( viter_t it = splits.begin(); it != splits.end(); ++it ) {
    fADCcor[*it] *= 2.0;
  }

  // Negative return value indicates potential problem
  if( nHits > fMaxHits )
    nHits = -nHits;
  
  return nHits;
}

//_____________________________________________________________________________
Hit* GEMPlane::AddHitImpl( Double_t pos )
{
  // Make a dummy hit of the correct type at the given projection coordinate
  // and add it to the hit array

  assert( IsDummy() );

  GEMHit* theHit = 0;

  // Emulate parameters for dummy hits
  const UInt_t size = 1, type = 0;
  const Double_t adcsum = 10.*fMinAmpl, resolution = fResolution, time = 0.0;
  
#ifdef MCDATA
  const Int_t mctrack = 1, num_bg = 0;
  const Double_t mcpos = pos, mctime = 0.0;
  bool mc_data = fTracker->TestBit(Tracker::kMCdata);
  if( mc_data )
    // Monte Carlo data
    theHit = new( (*fHits)[GetNhits()] ) MCGEMHit( pos,
						   time,
						   adcsum,
						   size,
						   type,
						   resolution,
						   this,
						   mctrack,
						   mcpos,
						   mctime,
						   num_bg
						   );
    else
#endif
      theHit = new( (*fHits)[GetNhits()] ) GEMHit( pos,
						   time,
						   adcsum,
						   size,
						   type,
						   resolution,
						   this
						   );
  return theHit;
}

//_____________________________________________________________________________
Int_t GEMPlane::Decode( const THaEvData& evData )
{
  // Convert evData to hits
  if( IsDummy() )
    // Special "decoding" for dummy planes
    return DummyDecode( evData );

  return GEMDecode( evData );
}

//_____________________________________________________________________________
Int_t GEMPlane::DefineVariables( EMode mode )
{
  // initialize global variables

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { "nrawstrips",     "nstrips with decoder data",        "fNrawStrips" },
    { "nhitstrips",     "nstrips > 0",                      "fNhitStrips" },
    { "nstrips",        "Num strips with hits > adc.min",   "GetNsigStrips()" },
    { "hitocc",         "strips > 0 / n_all_strips",        "fHitOcc" },
    { "occupancy",      "nstrips / n_all_strips",           "fOccupancy" },
    { "strip.adcraw",   "Raw strip ADC sum",                "fADCraw" },
    { "strip.adc",      "Deconvoluted strip ADC sum",       "fADC" },
    { "strip.adc_c",    "Pedestal-sub strip ADC sum",       "fADCcor" },
    { "strip.time",     "Leading time of strip signal (ns)","fHitTime" },
    { "strip.good",     "Good pulse shape on strip",        "fGoodHit" },
    { "nhits",          "Num hits (clusters of strips)",    "GetNhits()" },
    { "noise",          "Noise level (avg below adc.min)",  "fDnoise" },
    { "ncoords",        "Num fit coords",                   "GetNcoords()" },
    { "coord.pos",      "Position used in fit (m)",         "fFitCoords.TreeSearch::FitCoord.fPos" },
    { "coord.trkpos",   "Track pos from projection fit (m)","fFitCoords.TreeSearch::FitCoord.fTrackPos" },
    { "coord.trkslope", "Track slope from projection fit",  "fFitCoords.TreeSearch::FitCoord.fTrackSlope" },
    { "coord.resid",    "Residual of trkpos (m)",           "fFitCoords.TreeSearch::FitCoord.GetResidual()" },
    { "coord.3Dpos",    "Crossing position of fitted 3D track (m)", "fFitCoords.TreeSearch::FitCoord.f3DTrkPos" },
    { "coord.3Dresid",  "Residual of 3D trkpos (m)",        "fFitCoords.TreeSearch::FitCoord.Get3DTrkResid()" },
    { "coord.3Dslope",  "Slope of fitted 3D track wrt projection",  "fFitCoords.TreeSearch::FitCoord.f3DTrkSlope" },
    { 0 }
  };
  Int_t ret = DefineVarsFromList( vars, mode );

  if( ret != kOK )
    return ret;

#ifdef MCDATA
  if( !fTracker->TestBit(Tracker::kMCdata) ) {
#endif
    // Non-Monte Carlo hit data
    RVarDef nonmcvars[] = {
      { "hit.pos",  "Hit centroid (m)",      "fHits.TreeSearch::GEMHit.fPos" },
      { "hit.time", "Hit reconstructed time","fHits.TreeSearch::GEMHit.fTime" },
      { "hit.adc",  "Hit ADC sum",           "fHits.TreeSearch::GEMHit.fADCsum" },
      { "hit.size", "Num strips ",           "fHits.TreeSearch::GEMHit.fSize" },
      { "hit.type", "Hit analysis result",   "fHits.TreeSearch::GEMHit.fType" },
      { 0 }
    };
    ret = DefineVarsFromList( nonmcvars, mode );
#ifdef MCDATA
  } else {
    // Monte Carlo hit data includes the truth information
    // For safety, we make sure that all hit variables are referenced with
    // respect to the MCGEMHit class and not just GEMHit - the memory layout
    // of classes under multiple inheritance might be implemetation-dependent
    RVarDef mcvars[] = {
      { "hit.pos",   "Hit centroid (m)",      "fHits.TreeSearch::MCGEMHit.fPos" },
      { "hit.time",  "Hit reconstructed time","fHits.TreeSearch::GEMHit.fTime" },
      { "hit.adc",   "Hit ADC sum",           "fHits.TreeSearch::MCGEMHit.fADCsum" },
      { "hit.size",  "Num strips ",           "fHits.TreeSearch::MCGEMHit.fSize" },
      { "hit.type",  "Hit analysis result",   "fHits.TreeSearch::MCGEMHit.fType" },
      { "hit.mctrk", "MC track number",       "fHits.TreeSearch::MCGEMHit.fMCTrack" },
      { "hit.mcpos", "MC track position (m)", "fHits.TreeSearch::MCGEMHit.fMCPos" },
      { "hit.mctime","MC track time (s)",     "fHits.TreeSearch::MCGEMHit.fMCTime" },
      { "hit.numbg", "MC num backgr hits",    "fHits.TreeSearch::MCGEMHit.fContam" },
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
  }
#endif
  return ret;
}

//_____________________________________________________________________________
Int_t GEMPlane::End( THaRunBase* run )
{
  // Write diagnostic histograms to file

  Plane::End( run );

#ifdef TESTCODE
  if( TestBit(kDoHistos) ) {
    assert( fADCMap );
    fADCMap->Write();
    assert( fHitStripTime );
    fHitStripTime->Write();
    assert( fADCsamp_sig );
    fADCsamp_sig->Write();
    assert( fADCsampVsStripTime_bkgd );
    fADCsampVsStripTime_bkgd->Write();
  }
#endif
  return 0;
}

//_____________________________________________________________________________
Int_t GEMPlane::ReadDatabase( const TDatime& date )
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

  // Set defaults
  TString mapping;
  Int_t do_noise = 1, check_pulse_shape = 1;
  fMaxClusterSize = kMaxUInt;
  fMinAmpl   = 0.0;
  fSplitFrac = 0.0;
  fMapType   = kOneToOne;
  fMaxSamp   = 1;
  fChanMap.clear();
  fPed.clear();
  fAmplSigma = 0.36; // default, an educated guess

  Int_t gbl = GetDBSearchLevel(fPrefix);
  try {
    const DBRequest request[] = {
      { "nstrips",        &fNelem,          kInt,     0, 0, gbl },
      { "strip.pos",      &fStart },
      { "strip.pitch",    &fPitch,          kDouble,  0, 0, gbl },
      { "maxclustsiz",    &fMaxClusterSize, kUInt,    0, 1, gbl },
      { "maxsamp",        &fMaxSamp,        kUInt,    0, 1, gbl },
      { "adc.min",        &fMinAmpl,        kDouble,  0, 1, gbl },
      { "timecut.cv",     &fTimeCutCV,      kDouble,  0, 1, gbl },
      { "timecut.hw",     &fTimeCutHW,      kDouble,  0, 1, gbl },
      { "split.frac",     &fSplitFrac,      kDouble,  0, 1, gbl },
      { "mapping",        &mapping,         kTString, 0, 1, gbl },
      { "chanmap",        &fChanMap,        kIntV,    0, 1, gbl },
      { "pedestal",       &fPed,            kFloatV,  0, 1 },
      { "do_noise",       &do_noise,        kInt,     0, 1, gbl },
      { "adc.sigma",      &fAmplSigma,      kDouble,  0, 1, gbl },
      { "check_pulse_shape",&check_pulse_shape, kInt, 0, 1, gbl },
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

  // Sanity checks
  if( fNelem <= 0 or fNelem > kMaxNChan ) {
    Error( Here(here), "Invalid number of channels: %d. Must be > 0 and < %d. "
	   "Fix database or recompile code with a larger limit.",
	   fNelem, kMaxNChan );
    return kInitError;
  }

  // Dummy planes ignore all of the parameters that are checked below,
  // so we can return right here.
  if( IsDummy() ) {
    fIsInit = true;
    return kOK;
  }

  Int_t nchan = fDetMap->GetTotNumChan();
  if( nchan != fNelem ) {
    Error( Here(here), "Number of detector map channels (%d) "
	   "disagrees with number of strips (%d)", nchan, fNelem );
    return kInitError;
  }

  SetBit( kDoNoise, do_noise );
  SetBit( kCheckPulseShape, check_pulse_shape );

  SafeDelete(fADCraw);
  SafeDelete(fADC);
  SafeDelete(fHitTime);
  SafeDelete(fADCcor);
  SafeDelete(fGoodHit);
#ifdef MCDATA
  delete [] fMCHitInfo; fMCHitInfo = 0;
#endif

  // Allocate arrays. The only reason that these are parallel C-arrays is
  // that the global variable system still doesn't support arrays/vectors
  // of structures/objects.
  // Out of memory exceptions from here are caught in Tracker.cxx.
  fADCraw = new Float_t[fNelem];
  fADC = new Float_t[fNelem];
  fHitTime = new Float_t[fNelem];
  fADCcor = new Float_t[fNelem];
  fGoodHit = new Byte_t[fNelem];
  fSigStrips.reserve(fNelem);
  fStripsSeen.resize(fNelem);

#ifdef MCDATA
  if( fTracker->TestBit(Tracker::kMCdata) ) {
    fMCHitInfo = new Podd::MCHitInfo[fNelem];
    fMCHitList.reserve(fNelem);
  }
#endif

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
	Error( Here(here), "Number of channel map entries (%u) must equal "
	       "number of strips (%d). Fix database.",
	       static_cast<unsigned int>(fChanMap.size()), fNelem );
	return kInitError;
      }
      // check if entries in channel map are within range
      for( Vint_t::const_iterator it = fChanMap.begin();
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
  if( fMaxSamp == 1 )
    ResetBit( kCheckPulseShape );

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
void GEMPlane::Print( Option_t* opt ) const
{
  // Print plane info

  Plane::Print( opt );
  //TODO: add detailed info
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

