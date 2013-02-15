//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 11-Jan-2010
//
#ifndef ROOT_TreeSearch_Plane
#define ROOT_TreeSearch_Plane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Plane                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Types.h"
#include "THaSubDetector.h"
#include "TClonesArray.h"
#include "TVector2.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <cassert>
#include <functional>

using std::vector;
using std::string;

class TH1;

namespace TreeSearch {

  class Projection;
  class Tracker;
  class Hit;
  class FitCoord;
  extern const Double_t kBig;

  class Plane : public THaSubDetector {

  public:
    Plane( const char* name, const char* description = "", 
	   THaDetectorBase* parent = 0 );
    // For ROOT RTTI
    Plane() : fPartner(0), fProjection(0), fHits(0), fFitCoords(0),
	      fHitMap(0) {}
    virtual ~Plane();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& ) = 0;
    virtual void    Print( Option_t* opt="" ) const;

    virtual Int_t   Begin( THaRunBase* r=0 );
    virtual Int_t   End( THaRunBase* r=0 );

    virtual Double_t GetMaxLRdist() const { return 0; }

//     virtual Int_t   Compare ( const TObject* obj ) const;
//     virtual Bool_t  IsSortable () const { return kTRUE; }

    FitCoord*       AddFitCoord( const FitCoord& coord );
    Bool_t          Contains( Double_t x, Double_t y ) const;
    Bool_t          Contains( const TVector2& point ) const;
    void            EnableCalibration( Bool_t enable = true );
    EProjType       GetType()        const { return fType; }
    Double_t        GetZ()           const { return fOrigin.Z(); }
    Projection*     GetProjection()  const { return fProjection; }
    Plane*          GetPartner()     const { return fPartner; }
    void            SetPartner( Plane* p );
    const char*     GetPartnerName() const { return fPartnerName.Data(); }
    Tracker*        GetTracker()     const { return fTracker; }
    Double_t        GetResolution()  const { return fResolution; }
    Double_t        GetMaxSlope()    const; 
    Double_t        GetStart()       const { return fStart+fCoordOffset; }
    Double_t        GetPitch()       const { return fPitch; }

    TSeqCollection* GetHits()        const { return fHits; }
    Int_t           GetNhits()       const { return fHits->GetLast()+1; }
    TSeqCollection* GetCoords()      const { return fFitCoords; }
    Int_t           GetNcoords()     const { return fFitCoords->GetLast()+1; }
    UInt_t          GetPlaneNum()    const { return fPlaneNum; }
    UInt_t          GetDefinedNum()  const { return fDefinedNum; }
    Bool_t          IsCalibrating()  const { return TestBit(kCalibrating); }
    Bool_t          IsRequired()     const;

    void            SetPlaneNum( UInt_t n )   { fPlaneNum = n; }
    void            SetDefinedNum( UInt_t n ) { fDefinedNum = n; }
    void            SetProjection( Projection* p )
    { fProjection = p; UpdateOffset(); }
    void            SetRequired( Bool_t enable = true );
    void            UpdateOffset();

    // Helper functors for STL algorithms...
    struct ZIsLess
      : public std::binary_function< Plane*, Plane*, bool >
    {
      bool operator() ( const Plane* a, const Plane* b ) const
      { 
	assert( a && b );
	assert( a->GetDefinedNum() != b->GetDefinedNum() );
	if( a->GetZ() < b->GetZ() ) return true;
	if( b->GetZ() < a->GetZ() ) return false;
	return( a->GetDefinedNum() < b->GetDefinedNum() );
      }
    };

    class NameEquals : public std::unary_function< Plane*, bool > {
    public:
      NameEquals( const char* s ) : name(s?s:"") {}
      bool operator() ( const Plane* pl ) const
      { assert(pl); return ( name == pl->GetName() ); }
    private:
      string name;
    };

  protected:

    // Geometry, configuration
    UInt_t        fPlaneNum;    // Ordinal of this plane within its projection
    UInt_t        fDefinedNum;  // Ordinal of this plane in planeconfig string
    EProjType     fType;        // Plane type (e.g. x,y,u,v)
    Double_t      fStart;       // Position of 1st sensor (e.g. wire/strip) (m)
    Double_t      fPitch;       // Sensor spacing (assumed constant) (m)
    Double_t      fCoordOffset; // Coord offset wrt Tracker due to fOrigin
    TString       fPartnerName; // Name of partner plane in same 2-D readout
    Plane*        fPartner;     //! Partner plane
    Projection*   fProjection;  //! The projection that we belong to
    Tracker*      fTracker;     //! Our parent detector

    // Parameters, calibration, flags
    Double_t      fResolution;  // Sensor position resolution (sigma) (m)
    UInt_t        fMaxHits;     // Maximum # hits before flagging decode error

    // Event data, hits etc.
    TClonesArray* fHits;        // Cluster data (groups of hits)
    TClonesArray* fFitCoords;   // Cluster coordinates used by good road fits

    Int_t ReadDatabaseCommon( const TDatime& date );

    // Podd interface
    virtual Int_t ReadDatabase( const TDatime& date ) = 0;

    // Only for TESTCODE, but kept for binary compatibility
    TH1*          fHitMap;     // Histogram of active sensor numbers

    // Bits for Planes
    enum {
      kIsRequired  = BIT(14), // Tracks must have a hit in this plane
      kCalibrating = BIT(15), // Plane in calibration mode (implies !required)
      kTDCReadout  = BIT(21), // Readout is TDC-based
      kHaveRefChans= BIT(22)  // Hardware modules have reference channels
#ifdef TESTCODE
     ,kDoHistos    = BIT(23)  // Generate diagnostic histograms
#endif
    };

  private:
#ifndef NDEBUG
    Tracker* GetTrackerSafe()  const;
#endif
    
    ClassDef(Plane,0)  // One Tracker plane coordinate direction
  };

  //___________________________________________________________________________
  inline
  Bool_t Plane::Contains( Double_t x, Double_t y ) const
  {
    // Check if the given point is within the active area of this wire plane.
    // Coordinates are relative to the Tracker origin. Time-critical, may be
    // be called O(1e5) per event

    //TODO: allow for (small) rotation due to misalignment?

    return ( TMath::Abs( x-fOrigin.X() ) < fSize[0] and
	     TMath::Abs( y-fOrigin.Y() ) < fSize[1] );
  }

  //___________________________________________________________________________
  inline
  Bool_t Plane::Contains( const TVector2& point ) const
  {
    // Same as Contains(x,y), but using a TVector2 as input

    return Contains( point.X(), point.Y() );
  }

  //___________________________________________________________________________
  inline
  void Plane::EnableCalibration( Bool_t enable )
  {
    // Enable/disable calibration mode flag

    SetBit( kCalibrating, enable );
  }
  
  //___________________________________________________________________________
  inline
  Bool_t Plane::IsRequired() const
  { 
    return ( not IsCalibrating() and TestBit(kIsRequired) );
  }

  //___________________________________________________________________________
  inline
  void Plane::SetRequired( Bool_t enable )
  {
    // Enable/disable calibration mode flag

    SetBit( kIsRequired, enable );
  }
  
///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
