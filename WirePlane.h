#ifndef ROOT_TreeSearch_WirePlane
#define ROOT_TreeSearch_WirePlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::WirePlane                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaSubDetector.h"
#include "TClonesArray.h"
#include "Types.h"
#include <vector>
#include <string>
#include <cassert>
#include <functional>

using std::vector;
using std::string;

namespace TreeSearch {

  class TimeToDistConv;
  class Projection;
  class MWDC;
  class Hit;
  class FitCoord;
  extern const Double_t kBig;

  class WirePlane : public THaSubDetector {

  public:
    WirePlane( const char* name, const char* description = "", 
	       THaDetectorBase* parent = 0 );
    virtual ~WirePlane();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );
    virtual void    Print( Option_t* opt="" ) const;

//     virtual Int_t   Compare ( const TObject* obj ) const;
//     virtual Bool_t  IsSortable () const { return kTRUE; }

    FitCoord*       AddFitCoord( const FitCoord& coord );
    EProjType       GetType() const { return fType; }
    Double_t        GetZ() const { return fOrigin.Z(); }
    Projection*     GetProjection() const { return fProjection; }
    WirePlane*      GetPartner() const { return fPartner; }
    void            SetPartner( WirePlane* p );

    Double_t        GetResolution() const { return fResolution; }
    Double_t        GetMaxSlope() const; 
    Double_t        GetWireStart() const { return fWireStart; }
    Double_t        GetWireSpacing() const { return fWireSpacing; }
    virtual Double_t GetMaxLRdist() const { return GetWireSpacing(); }

    TimeToDistConv* GetTTDConv() const { return fTTDConv; }

    TSeqCollection* GetHits()     const { return fHits; }
    Int_t           GetNhits()    const { return fHits->GetLast()+1; }
    TSeqCollection* GetCoords()   const { return fFitCoords; }
    Int_t           GetNcoords()  const { return fFitCoords->GetLast()+1; }
    UInt_t          GetPlaneNum() const { return fPlaneNum; }
    UInt_t          GetLayerNum() const { return fLayerNum; }
    Bool_t          IsRequired()  const { return false; } //TODO: use user bit

    void            SetPlaneNum( UInt_t n ) { fPlaneNum = n; }
    void            SetLayerNum( UInt_t n ) { fLayerNum = n; }
    void            SetProjection( Projection* p ) { fProjection = p; }
#ifdef TESTCODE
    void            CheckCrosstalk();
#endif

    // Helper functors for STL algorithms...
    struct ZIsLess
      : public std::binary_function< WirePlane*, WirePlane*, bool >
    {
      bool operator() ( const WirePlane* a, const WirePlane* b ) const
      { assert(a&&b); return ( a->GetZ() < b->GetZ() ); }
    };

    class NameEquals : public std::unary_function< WirePlane*, bool > {
    public:
      NameEquals( const char* s ) : name(s?s:"") {}
      bool operator() ( const WirePlane* wp ) const
      { assert(wp); return ( name == wp->GetName() ); }
    private:
      string name;
    };

  protected:

    // Geometry, configuration

    UInt_t      fPlaneNum;     // Ordinal of this plane within its projection
    UInt_t      fLayerNum;     // Layer number of this plane within its proj
    EProjType   fType;         // Plane type (x,y,u,v)
    Double_t    fWireStart;    // Position of 1st wire (along wire coord) (m)
    Double_t    fWireSpacing;  // Wire spacing (assumed constant) (m)
    WirePlane*  fPartner;      //! Partner plane (usually with staggered wires)
    Projection* fProjection;   //! The projection that we belong to
    MWDC*       fMWDC;         //! Our parent detector

    // Parameters, calibration, flags

    Double_t    fResolution;   // Drift distance resolution (sigma) (m)
    Double_t    fMinTime;      // Minimum drift time for a hit (s)
    Double_t    fMaxTime;      // Maximum drift time for a hit (s)

    TimeToDistConv* fTTDConv;   // Drift time->distance converter
    vector<float>   fTDCOffset; // [fNelem] TDC offsets for each wire

    // Event data, hits etc.

    TClonesArray*   fHits;      // Hit data
    TClonesArray*   fFitCoords; // Hit coordinates used by good fits in roads

#ifdef TESTCODE
    UInt_t          fNmiss;     // Statistics: Decoder channel misses
    UInt_t          fNrej;      // Statistics: Rejected hits
    Int_t           fWasSorted; // Statistics: hits were sorted (0/1)
    UInt_t          fNhitwires; // Statistics: wires with one or more hits
    UInt_t          fNmultihit; // Statistics: wires with multiple hits
    UInt_t          fNmaxmul;   // Statistics: largest num hits on any wire
    UInt_t          fNcl;       // Statistics: number of hit "clusters"
    UInt_t          fNdbl;      // Statistics: num wires with neighboring hits
    UInt_t          fClsiz;     // Statistics: max cluster size
#endif

    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );

    ClassDef(WirePlane,0)    // One MWDC wire plane
  };
}

///////////////////////////////////////////////////////////////////////////////

#endif
