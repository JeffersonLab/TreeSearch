#ifndef ROOT_TreeSearch_WirePlane
#define ROOT_TreeSearch_WirePlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::WirePlane                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaSubDetector.h"
#include "TClonesArray.h"
#include <vector>

using std::vector;

class THaVDCTimeToDistanceConv;

namespace TreeSearch {

  class WirePlane : public THaSubDetector {

    enum EType { kUndefinedType = 0, kXPlane, kYPlane, kUPlane, kVPlane };

  public:
    WirePlane( const char* name, const char* description = "", 
	       THaDetectorBase* parent = NULL );
    virtual ~WirePlane();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );

    EType   GetType() const { return fType; }
  

  protected:

    // Geometry, configuration

    EType     fType;         // Plane type (x,y,u,v)
    Double_t  fWireStart;    // Position of 1st wire (along wire coord) (m)
    Double_t  fWireSpacing;  // Wire spacing (assumed constant) (m)
    Double_t  fSinAngle;     // Sine of angle between dispersive direction (x)
                             //  and direction of decreasing wire number (rad)
    Double_t  fCosAngle;     // Cosine of wire angle (rad)

    // Parameters, calibration, flags

    Double_t  fTDCRes;       // TDC Resolution ( s / channel)
    Double_t  fDriftVel;     // Drift velocity in the wire plane (m/s)
    Double_t  fResolution;   // Drift distance resolution (sigma) (m)
    Double_t  fMaxSlope;     // Max physical track slope in this plane (FIXME?)

    THaVDCTimeToDistanceConv* fTTDConv; // Drift time->distance converter
    vector<double>  fTDCOffset;    // [fNelem] TDC offsets for each wire

    // Event data, hits etc.

    TClonesArray*   fHits;    // Hit data

    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );

    ClassDef(WirePlane,0)    // One MWDC wire plane
  };
}

///////////////////////////////////////////////////////////////////////////////

#endif
