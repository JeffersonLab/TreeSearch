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

namespace TreeSearch {

  class GEMHit : public Hit {

  public:
    GEMHit() {}
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

  //___________________________________________________________________________
  // Monte Carlo hit class. Same as a hit plus the MC truth info.

  class MCGEMHit : public GEMHit, public MCHitInfo {

  public:
    MCGEMHit() {}
    MCGEMHit( Double_t pos, Double_t adc_sum, UInt_t num_strips, Int_t type,
	      Double_t res, GEMPlane* pl, MCTrack* mctrk, Double_t mcpos )
      : GEMHit(pos, adc_sum, num_strips, type, res, pl),
	MCHitInfo(mctrk, mcpos) {}
    virtual ~MCGEMHit() {}

    virtual void Print( Option_t* opt="" ) const;

    ClassDef(MCGEMHit,1)   // Monte Carlo hit in ADC-based readout plane
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
