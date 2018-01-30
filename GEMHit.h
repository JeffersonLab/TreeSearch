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
    Hit(pos, res, static_cast<Plane*>(pl)), fModule(module),
      fADCsum(adc_sum), fADCmax(adc_max), fPeaktime(peaktime),
      fSize(num_strips), fType(type) {}
    virtual ~GEMHit() {}

    virtual void Print( Option_t* opt="" ) const;

    Double_t GetADCsum()     const { return fADCsum; }
    Double_t GetADCmax()     const { return fADCmax; }
    UInt_t   GetSize()       const { return fSize; }
    Int_t    GetType()       const { return fType; }

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
  MCGEMHit( Int_t module, Double_t pos, Double_t adc_sum, Double_t adc_max, Double_t peaktime, UInt_t num_strips, Int_t type,
	      Double_t res, GEMPlane* pl, Int_t mctrk, Double_t mcpos, Double_t mccharge,
	      Double_t mctime, Int_t num_bg_strips )
    : GEMHit(module, pos, adc_sum, adc_max, peaktime, num_strips, type, res, pl),
      TSBSMCHitInfo(mctrk, mcpos, mctime, mccharge, num_bg_strips) {}
    virtual ~MCGEMHit() {}

    virtual void Print( Option_t* opt="" ) const;
    // Int_t Compare(const TObject *obj) const;
    
  protected:
    
    ClassDef(MCGEMHit,2)   // Monte Carlo hit in ADC-based readout plane
  };
#endif // MCDATA

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
