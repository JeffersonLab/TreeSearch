#ifndef ROOT_TreeSearch_WirePlane
#define ROOT_TreeSearch_WirePlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::WirePlane                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaSubDetector.h"
#include "Projection.h"
#include <vector>
#include <string>
#include <functional>

class TSeqCollection;

using std::vector;

namespace TreeSearch {

  class TimeToDistConv;
  typedef Projection::EProjType EProjType;

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
    Double_t        GetMaxSlope() const { return fMaxSlope; }

    TimeToDistConv* GetTTDConv() const { return fTTDConv; }

    TSeqCollection* GetHits() const { return fHits; }
    Int_t           GetNhits() const;
    UInt_t          GetPlaneNum() const { return fPlaneNum; }

    void            SetPlaneNum( UInt_t n ) { fPlaneNum = n; }
    void            SetProjection( Projection* p ) { fProjection = p; }

    // Helper functors for STL algorithms...
    //    static bool     ComparePlaneZ( WirePlane* a, WirePlane* b )
    //    { return ( a && a->Compare(b) < 0 ); }
    struct ComparePlaneZ
      : public std::binary_function< WirePlane*, WirePlane*, bool >
    {
      bool operator() ( const WirePlane* a, const WirePlane* b ) const
      { return ( a && a->Compare(b) < 0 ); }
    };

    class NameIs : public std::unary_function< WirePlane*, bool > {
    public:
      NameIs( const char* s ) : name(s) {}
      bool operator() ( const WirePlane* wp ) const
      { return ( wp && name == wp->GetName() ); }
    private:
      std::string name;
    };

  protected:

    // Geometry, configuration

    UInt_t     fPlaneNum;     // Number of this plane (0..Nplanes)
    EProjType  fType;         // Plane type (x,y,u,v)
    Double_t   fWireStart;    // Position of 1st wire (along wire coord) (m)
    Double_t   fWireSpacing;  // Wire spacing (assumed constant) (m)
    Double_t   fSinAngle;     // Sine of angle between dispersive direction (x)
                              //  and direction of decreasing wire number
    Double_t   fCosAngle;     // Cosine of wire angle
    WirePlane* fPartner;      //! Partner plane with staggered wires
    Projection* fProjection;  //! Parameters of this plane type

    // Parameters, calibration, flags

    Double_t   fTDCRes;       // TDC Resolution ( s / channel)
    Double_t   fDriftVel;     // Drift velocity in the wire plane (m/s)
    Double_t   fResolution;   // Drift distance resolution (sigma) (m)
    Double_t   fMaxSlope;     // Maximum physical slope of track (0=perp)

    TimeToDistConv* fTTDConv;   // Drift time->distance converter
    vector<double>  fTDCOffset; // [fNelem] TDC offsets for each wire

    // Event data, hits etc.

    TSeqCollection* fHits;      // Hit data

    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );

    ClassDef(WirePlane,0)    // One MWDC wire plane
  };
}

///////////////////////////////////////////////////////////////////////////////

#endif
