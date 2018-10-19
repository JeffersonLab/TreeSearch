//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 14-Jun-2007
//
#ifndef ROOT_TreeSearch_WirePlane
#define ROOT_TreeSearch_WirePlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::WirePlane                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Types.h"
#include "Plane.h"

namespace TreeSearch {

  class TimeToDistConv;

  class WirePlane : public Plane {

  public:
    WirePlane( const char* name, const char* description = "",
	       THaDetectorBase* parent = 0 );
    WirePlane()
      : fMinTime(0), fMaxTime(0), fTTDConv(0),
        fNmiss(0), fNrej(0), fWasSorted(0), fNhitwires(0), fNmultihit(0),
        fNmaxmul(0), fNcl(0), fNdbl(0), fClsiz(0) {} // For ROOT RTTI
    virtual ~WirePlane();

    virtual void     Clear( Option_t* opt="" );
    virtual Int_t    Decode( const THaEvData& );
    virtual void     Print( Option_t* opt="" ) const;

    virtual Double_t GetMaxLRdist() const { return GetPitch(); }
    virtual Hit*     FindNearestHitAndPos( Double_t x, Double_t& pos ) const;


    TimeToDistConv*  GetTTDConv()   const { return fTTDConv; }

#ifdef TESTCODE
    void             CheckCrosstalk();
#endif

  protected:

    // Parameters, calibration, flags

    Double_t        fMinTime;   // Minimum drift time for a hit (s)
    Double_t        fMaxTime;   // Maximum drift time for a hit (s)

    TimeToDistConv* fTTDConv;   // Drift time->distance converter
    Vflt_t          fTDCOffset; // [fNelem] TDC offsets for each wire

    // Only needed for TESTCODE, but kept for binary compatibility
    UInt_t          fNmiss;     // Statistics: Decoder channel misses
    UInt_t          fNrej;      // Statistics: Rejected hits
    Int_t           fWasSorted; // Statistics: hits were sorted (0/1)
    UInt_t          fNhitwires; // Statistics: wires with one or more hits
    UInt_t          fNmultihit; // Statistics: wires with multiple hits
    UInt_t          fNmaxmul;   // Statistics: largest num hits on any wire
    UInt_t          fNcl;       // Statistics: number of hit "clusters"
    UInt_t          fNdbl;      // Statistics: num wires with neighboring hits
    UInt_t          fClsiz;     // Statistics: max cluster size

    // Support functions for dummy planes
    virtual Hit*  AddHitImpl( Double_t x );
    virtual Int_t WireDecode( const THaEvData& );

    // Podd interface
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );

    ClassDef(WirePlane,0)      // One MWDC wire plane
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
