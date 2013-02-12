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
  //TODO: Skeleton version - to be fleshed out

  class MCTrack;
  class MCGEMHit : public GEMHit {

  public:
    MCGEMHit() : fMCTrack(0) {}
    MCGEMHit( Int_t wnum, Double_t pos, Int_t tdc, Double_t time, Double_t res,
	      GEMPlane* pl, MCTrack* mctrk, Double_t mcpos )
      : GEMHit(wnum, pos, tdc, time, res, pl), fMCTrack(mctrk), fMCPos(mcpos) {}
    virtual ~MCGEMHit() {}

    virtual void Print( Option_t* opt="" ) const;

    MCTrack* GetMCTrack() const { return fMCTrack; }
    Double_t GetMCPos()   const { return fMCPos; }

  protected:
    MCTrack* fMCTrack;     // MC track generating this hit (0=noise hit)
    Double_t fMCPos;       // Exact MC track crossing position (m)

    ClassDef(MCGEMHit,1)   // Monte Carlo hit in ADC-based readout plane
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
