//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 05-Mar-2010
//
#ifndef ROOT_TreeSearch_GEMHit
#define ROOT_TreeSearch_GEMHit

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::GEMHit                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hit.h"
#include "GEMPlane.h"
#ifdef MCDATA
#include "SimDecoder.h"  // for MC data support
#include "TSBSSimDecoder.h"
#endif

namespace TreeSearch {

  class GEMHit : public Hit {

  public:
    GEMHit() {}
  GEMHit( Int_t module, Double_t pos, Double_t adc_sum, Double_t adc_max, Double_t peaktime, 
	  UInt_t num_strips, Int_t type,
	  Double_t res, GEMPlane* pl ) :
    Hit(pos, res, static_cast<Plane*>(pl), module),
      fADCsum(adc_sum), fADCmax(adc_max), fPeaktime(peaktime),
      fSize(num_strips), fType(type) {}
    virtual ~GEMHit() {}

    virtual void Print( Option_t* opt="" ) const;

    Double_t GetADCsum()     const { return fADCsum; }
    Double_t GetADCmax()     const { return fADCmax; }
   
    Double_t GetPeaktime()   const { return fPeaktime; }
    UInt_t   GetSize()       const { return fSize; }
    Int_t    GetType()       const { return fType; }
    Int_t    GetModule()     const { return fModule;}

  protected:
    Int_t    fModule;      // Module ID of GEM hit in plane
    Double_t fADCsum;      // Sum of ADC values of active strips
    Double_t fADCmax;      // Max of cluster pulse from fit
    
    Double_t fPeaktime;    // peak time of cluster pulse from fit
    UInt_t   fSize;        // Number of active strips
    Int_t    fType;        // Result code of cluster analysis

    ClassDef(GEMHit,1)     // Hit on an ADC-based readout plane, e.g. GEMs
  };

#ifdef MCDATA
  //___________________________________________________________________________
  // Monte Carlo hit class. Same as a hit plus the MC truth info.

  class MCGEMHit : public GEMHit, public TSBSMCHitInfo {

  public:
    MCGEMHit() {}
  MCGEMHit( Int_t module, Double_t pos, Double_t adc_sum, Double_t adc_max,  Double_t prim_ratio, Double_t p_over_total, Double_t b_over_total,  Double_t pos_sigma, Double_t time_sigma, Double_t adc_sigma, Int_t nb_bg, Double_t peaktime, UInt_t num_strips, Int_t stripl, Int_t striph, Int_t type,
	      Double_t res, GEMPlane* pl, Int_t mctrk, Double_t mcpos, Double_t mccharge,
	      Double_t mctime, Int_t num_bg_strips )
    : GEMHit(module, pos, adc_sum, adc_max, peaktime, num_strips, type, res, pl),
      TSBSMCHitInfo(mctrk, mcpos, mctime, mccharge, num_bg_strips), 
      fprim_ratio(prim_ratio), fp_over_total(p_over_total), fb_over_total(b_over_total), fpos_sigma(pos_sigma), ftime_sigma(time_sigma), fadc_sigma(adc_sigma), fnb_bg(nb_bg), fstripl(stripl), fstriph(striph) {}
    
    virtual ~MCGEMHit() {}

    virtual void Print( Option_t* opt="" ) const;
    Double_t Getprim_ratio() const { return fprim_ratio;}
    Double_t Getp_over_total() const { return fp_over_total;}
    Double_t Getb_over_total() const { return fb_over_total;}
    Double_t Getpos_sigma() const { return fpos_sigma;}
    Double_t Gettime_sigma() const { return ftime_sigma;}
    Double_t Getadc_sigma() const { return fadc_sigma;}
    Int_t    Getnb_bg() const { return fnb_bg;}
    Int_t    Getstripl()       const { return fstripl;}
    Int_t    Getstriph()       const { return fstriph;}
    // Int_t Compare(const TObject *obj) const;
    
  protected:
    Double_t fprim_ratio;   // percentage of adcs of primary hit that this cluster includes 
    Double_t fp_over_total; // ratio of adc from primary hit to total adc of this cluster
    Double_t fb_over_total; // ratio of adc from background file hits to total adc of this cluster
    Double_t fpos_sigma;
    Double_t ftime_sigma;
    Double_t fadc_sigma;
    Int_t    fnb_bg;
    // Double_t
    Int_t fstripl;
    Int_t fstriph;
    ClassDef(MCGEMHit,2)    // Monte Carlo hit in ADC-based readout plane
  };
#endif // MCDATA

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
