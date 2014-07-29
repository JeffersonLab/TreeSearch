//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 11-Jan-2010
//
#ifndef ROOT_TreeSearch_GEMPlane
#define ROOT_TreeSearch_GEMPlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::GEMPlane                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Plane.h"
#include "TBits.h"

namespace TreeSearch {

  class GEMPlane : public Plane {
  public:
    GEMPlane( const char* name, const char* description = "",
	      THaDetectorBase* parent = 0 );
    // For ROOT RTTI
    GEMPlane() : fADCraw(0), fADC(0), fHitTime(0), fADCcor(0), fGoodHit(0) {}
    virtual ~GEMPlane();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual void    Print( Option_t* opt="" ) const;

#ifdef TESTCODE
    virtual Int_t   Begin( THaRunBase* r=0 );
    virtual Int_t   End( THaRunBase* r=0 );
#endif

    Double_t        GetAmplSigma( Double_t ampl ) const;
    Double_t        GetHitOcc()      const { return fHitOcc; }
    Double_t        GetOccupancy()   const { return fOccupancy; }

  protected:

    // Hardware channel mapping
    enum EChanMapType { kOneToOne, kReverse, kGassiplexAdapter1,
			kGassiplexAdapter2, kTable };

    EChanMapType  fMapType;     // Type of hardware channel mapping to use
    Vint_t        fChanMap;     // [fNelem] Optional hardware channel mapping
    Int_t         MapChannel( Int_t idx ) const;

    // Parameters, calibration, flags
    UInt_t        fMaxClusterSize;// Maximum size of a clean cluster of strips
    Double_t      fMinAmpl;     // ADC threshold for strips to be active
    Double_t      fSplitFrac;   // Percentage of amplitude swing necessary
                                // for detecting cluster maximum/minimum
    UInt_t        fMaxSamp;     // Maximum # ADC samples per channel

    Vflt_t        fPed;         // [fNelem] Per-channel pedestal values
    TBits         fBadChan;     // Bad channel map
    Double_t      fAmplSigma;   // Sigma of hit amplitude distribution

    // Event data, hits etc.
    Float_t*      fADCraw;      // [fNelem] Integral of raw ADC samples
    Float_t*      fADC;         // [fNelem] Integral of deconvoluted ADC samples
    Float_t*      fHitTime;     // [fNelem] Leading-edge time of deconv signal (ns)
    Float_t*      fADCcor;      // [fNelem] fADC corrected for pedestal & noise
    Byte_t*       fGoodHit;     // [fNelem] Strip data passed pulse shape test
    Double_t      fDnoise;      // Event-by-event noise (avg below fMinAmpl)
    Vint_t        fSigStrips;   // Strip numbers with signal (adccor > minampl)

    UInt_t        fNrawStrips;  // Statistics: strips with any data
    UInt_t        fNhitStrips;  // Statistics: strips > 0
    Double_t      fHitOcc;      // Statistics: hit occupancy fNhitStrips/fNelem
    Double_t      fOccupancy;   // Statistics: occupancy GetNsigStrips/fNelem

    // Optional diagnostics for TESTCODE, keep for binary compatibility
    TH1*          fADCMap;      // Histogram of strip numbers weighted by ADC

    Int_t         GetNsigStrips() const { return fSigStrips.size(); }
    void          AddHit( Int_t istrip );

    // Podd interface
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );

    // Bits for GEMPlanes
    enum {
      kDoNoise         = BIT(16), // Correct data for common-mode noise
      kCheckPulseShape = BIT(17)  // Reject malformed ADC pulse shapes
    };

    ClassDef(GEMPlane,0)  // ADC-based readout plane coordinate direction
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
