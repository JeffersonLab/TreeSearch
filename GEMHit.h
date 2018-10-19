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
#endif

namespace TreeSearch {

  class GEMHit : public Hit {

  public:
    GEMHit() : fADCsum(0), fSize(0), fType(0) {}
    GEMHit( Double_t pos, Double_t adc_sum, UInt_t num_strips, Int_t type,
	    Double_t res, GEMPlane* pl ) :
      Hit(pos, res, static_cast<Plane*>(pl)), fADCsum(adc_sum),
      fSize(num_strips), fType(type) {}
    virtual ~GEMHit() {}

    virtual void Print( Option_t* opt="" ) const;

    Double_t GetADCsum()     const { return fADCsum; }
    UInt_t   GetSize()       const { return fSize; }
    Int_t    GetType()       const { return fType; }

  protected:
    Double_t fADCsum;      // Sum of ADC values of active strips
    UInt_t   fSize;        // Number of active strips
    Int_t    fType;        // Result code of cluster analysis

    ClassDef(GEMHit,1)     // Hit on an ADC-based readout plane, e.g. GEMs
  };

#ifdef MCDATA
  //___________________________________________________________________________
  // Monte Carlo hit class. Same as a hit plus the MC truth info.

  class MCGEMHit : public GEMHit, public Podd::MCHitInfo {

  public:
    MCGEMHit() {}
    MCGEMHit( Double_t pos, Double_t adc_sum, UInt_t num_strips, Int_t type,
	      Double_t res, GEMPlane* pl, Int_t mctrk, Double_t mcpos,
	      Double_t mctime, Int_t num_bg_strips )
      : GEMHit(pos, adc_sum, num_strips, type, res, pl),
	MCHitInfo(mctrk, mcpos, mctime, num_bg_strips) {}
    virtual ~MCGEMHit() {}

    virtual void Print( Option_t* opt="" ) const;

    ClassDef(MCGEMHit,2)   // Monte Carlo hit in ADC-based readout plane
  };
#endif // MCDATA

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
