#ifndef ROOT_TreeSearch_WirePlane
#define ROOT_TreeSearch_WirePlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::WirePlane                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaSubDetector.h"
#include "TClonesArray.h"
#include "MWDC.h"
#include <vector>
#include <string>
#include <functional>

using std::vector;
using std::string;

namespace TreeSearch {

  class TimeToDistConv;
  class Projection;
  extern const Double_t kBig;

  class WirePlane : public THaSubDetector {

  public:
    WirePlane( const char* name, const char* description = "", 
	       THaDetectorBase* parent = NULL );
    virtual ~WirePlane();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );
    virtual void    Print( Option_t* opt="" ) const;

    Int_t           Compare ( const TObject* obj ) const;
    Bool_t          IsSortable () const { return kTRUE; }

    EProjType       GetType() const { return fType; }
    Double_t        GetZ() const { return fOrigin.Z(); }
    Projection*     GetProjection() const { return fProjection; }
    WirePlane*      GetPartner() const { return fPartner; }
    void            SetPartner( WirePlane* p );

    Double_t        GetResolution() const { return fResolution; }
    Double_t        GetMaxSlope() const 
    { return fProjection ? fProjection->GetMaxSlope() : kBig; }
    Double_t        GetWireStart() const { return fWireStart; }
    Double_t        GetWireSpacing() const { return fWireSpacing; }

    TimeToDistConv* GetTTDConv() const { return fTTDConv; }

    TSeqCollection* GetHits() const { return fHits; }
    Int_t           GetNhits() const { return fHits->GetLast()+1; }
    UInt_t          GetPlaneNum() const { return fPlaneNum; }

    void            SetPlaneNum( UInt_t n ) { fPlaneNum = n; }
    void            SetProjection( Projection* p ) { fProjection = p; }

    // Helper functors for STL algorithms...
    struct ZIsLess
      : public std::binary_function< WirePlane*, WirePlane*, bool >
    {
      bool operator() ( const WirePlane* a, const WirePlane* b ) const
      { return ( a && b && a->GetZ() < b->GetZ() ); }
    };

    class NameEquals : public std::unary_function< WirePlane*, bool > {
    public:
      NameEquals( const char* s ) : name(s?s:"") {}
      bool operator() ( const WirePlane* wp ) const
      { return ( wp && name == wp->GetName() ); }
    private:
      string name;
    };

  protected:

    // Geometry, configuration

    UInt_t      fPlaneNum;     // Ordinal of this plane within its projection
    EProjType   fType;         // Plane type (x,y,u,v)
    Double_t    fWireStart;    // Position of 1st wire (along wire coord) (m)
    Double_t    fWireSpacing;  // Wire spacing (assumed constant) (m)
    WirePlane*  fPartner;      //! Partner plane (usually with staggered wires)
    Projection* fProjection;   //! The projection that we belong to
    MWDC*       fMWDC;         //! Our parent detector

    // Parameters, calibration, flags

    //FIXME: per daq experts, the tdc res is per module not per plane
    Double_t    fTDCRes;       // TDC Resolution ( s / channel)
    Double_t    fDriftVel;     // Drift velocity in the wire plane (m/s)
    Double_t    fResolution;   // Drift distance resolution (sigma) (m)

    TimeToDistConv* fTTDConv;   // Drift time->distance converter
    vector<float>   fTDCOffset; // [fNelem] TDC offsets for each wire
    THaDetMap*      fRefMap;    // Map of reference channels, if any

    // Event data, hits etc.

    TClonesArray*   fHits;      // Hit data
    Double_t*       fRefTime;   // [fRefMap->GetSize()] reference channel data
    UInt_t          fNmiss;     // Statistics: Decoder channel misses
    UInt_t          fNrej;      // Statistics: Rejected hits
    Int_t           fWasSorted; // Statistics: hits were sorted fwd/rev (1/-1)
    UInt_t          fNhitwires; // Statistics: number of _wires_ hit
    UInt_t          fNnohits;   // Statistics: channels without hits

    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );

    ClassDef(WirePlane,0)    // One MWDC wire plane
  };
}

///////////////////////////////////////////////////////////////////////////////

#endif
