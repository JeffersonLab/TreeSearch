//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 27-Jun-2007
//
#ifndef ROOT_TreeSearch_Hit
#define ROOT_TreeSearch_Hit

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Plane.h"
#include "TBits.h"
#include "Node.h"   // for NodeDescriptor
#include "Helper.h" // for NumberOfSetBits
#include <utility>
#include <set>
#include <cassert>
#include <functional>
#include <iostream>

class TSeqCollection;
class TIterator;

namespace TreeSearch {

  class Road;
  extern const Double_t kBig;

  class Hit : public TObject {

  public:
    Hit() : fPlane(0) {}
  Hit( Double_t pos, Double_t res, Plane* pl, Int_t moduleID )
    : fPos(pos), fResolution(res), fPlane(pl), fModuleID(moduleID) { assert(fPlane); }
    // Default copy and assignment are fine
    //    Hit( const Hit& );
    //    Hit& operator=( const Hit& );
    virtual ~Hit() {}

    virtual Int_t Compare( const TObject* obj ) const;
    virtual Int_t Compare( const Hit* rhs, Double_t maxdist ) const;
    virtual Bool_t IsSortable () const { return kTRUE; }
    virtual void Print( Option_t* opt="" ) const;

    // Abstraction of the GetPosL/R interface in WireHits, currently
    // only used in Roads::CollectCoordinates
    virtual UInt_t   GetNumPos()         const;
    virtual Double_t GetPosI( UInt_t i ) const;

    Double_t GetPos()        const { return fPos; }
    Double_t GetResolution() const { return fResolution; }
    Double_t GetZ()          const { return fPlane->GetZ(); }

    Plane*   GetPlane()      const { return fPlane; }
    Int_t    GetModuleID()   const { return fModuleID; } 
    UInt_t   GetPlaneNum()   const { return fPlane->GetPlaneNum(); }
    UInt_t   GetAltPlaneNum() const { return fPlane->GetAltPlaneNum(); }

    // Functor for ordering hits in sets
    struct PosIsLess : public std::binary_function< Hit*, Hit*, bool >
    {
      bool operator() ( const Hit* a, const Hit* b ) const
      {
	assert( a && b );
	assert( a->GetPlaneNum() != kMaxUInt && b->GetPlaneNum() != kMaxUInt );
	if( a->GetPlane()->GetType() != b->GetPlane()->GetType() ) {
	  return (a->GetPlane()->GetType() < b->GetPlane()->GetType());
	}
	return ( a->GetPlaneNum() != b->GetPlaneNum() ) ?
	  (a->GetPlaneNum() < b->GetPlaneNum()) : (a->GetPos() < b->GetPos());
      }
    };

  protected:
    Double_t fPos;         // Hit position along plane coordinate axis (m)
    Double_t fResolution;  // Resolution of fPos (sigma, m)
    Plane*   fPlane;       //! Pointer to the plane obj where this hit occurred
    Int_t    fModuleID;

    ClassDef(Hit,1)        // Generic tracker plane hit
  };

  //___________________________________________________________________________
  // Utility class for iterating over one or two collections of hits.
  // Used for generating hit patterns. If two collections are given, they
  // are assumed to contain hits from adjacent planes with parallel
  // (and usually staggered) wires, and hit pairs are returned for hits
  // whose positions are within 'maxdist' of each other.

  typedef std::pair<TObject*,TObject*> ObjPair_t;

  class HitPairIter {

  public:
    HitPairIter( const TSeqCollection* collA, const TSeqCollection* collB,
		 Double_t maxdist );
    HitPairIter( const HitPairIter& rhs );
    HitPairIter& operator=( const HitPairIter& rhs );
    virtual ~HitPairIter();

    const TSeqCollection* GetCollection( Int_t n=0 ) const
    { return (n==0) ? fCollA : fCollB; }
    void Reset();

    // Iterator functions.
    HitPairIter& Next();
    HitPairIter& operator++() { return Next(); }
    const HitPairIter operator++(int) {
      HitPairIter clone(*this);  Next();  return clone;
    }
    // Current value.
    ObjPair_t  operator()() const  { return fCurrent; }
    ObjPair_t& operator* ()        { return fCurrent; }
    // Comparisons
    bool operator==( const HitPairIter& rhs ) const {
      return( fCollA == rhs.fCollA && fCollB == rhs.fCollB &&
	      fCurrent == rhs.fCurrent );
    }
    bool operator!=( const HitPairIter& rhs ) const { return !(*this==rhs); }
    operator bool() const
    { return (fCurrent.first != 0 || fCurrent.second != 0); }
    bool operator!() const { return !((bool)*this); }

  private:
    const TSeqCollection* fCollA;
    const TSeqCollection* fCollB;
    TIterator* fIterA;
    TIterator* fIterB;
    TIterator* fSaveIter;
    Hit* fSaveHit;
    Double_t fMaxDist;
    Bool_t fStarted;
    Bool_t fScanning;
    ObjPair_t fCurrent;
    ObjPair_t fNext;

    HitPairIter();

    ClassDef(HitPairIter,0)  // Iterator over two lists of hits
  };

  //___________________________________________________________________________
  // Coordinate information derived from fitting hits in a wire plane.
  // A given raw Hit may be associated with any number of FitCoord objects.

  class FitCoord : public TObject {

  public:
    FitCoord( Hit* hit, Road* road, Double_t pos,
	      Double_t trkpos2d, Double_t trkslope2d,
	      Double_t trkpos, Double_t trkslope )
      : fHit(hit), fRoad(road), fPos(pos), fTrackPos(trkpos2d),
	fTrackSlope(trkslope2d), f3DTrkPos(trkpos), f3DTrkSlope(trkslope) {}
    FitCoord() : fHit(0), fRoad(0) {} // For ROOT RTTI
    virtual ~FitCoord() {}

    Double_t  GetChi2()       const;
    Hit*      GetHit()        const { return fHit; }
    Road*     GetRoad()       const { return fRoad; }
    Double_t  GetPos()        const { return fPos; }
    // Double_t  GetDriftTime()  const { return fHit ? fHit->GetDriftTime():kBig;}
    // Double_t  GetDriftDist()  const { return fHit ? fHit->GetDriftDist():kBig;}
    Double_t  GetTrackPos()   const { return fTrackPos; }
    Double_t  GetTrackSlope() const { return fTrackSlope; }
    Double_t  GetTrackDist()  const
    { return fHit ? fTrackPos-fHit->GetPos() : kBig; }
    Double_t  GetResidual()   const { return fHit ? fTrackPos-fPos : kBig; }
    Double_t  Get3DTrkPos()   const { return f3DTrkPos; }
    Double_t  Get3DTrkSlope() const { return f3DTrkSlope; }
    Double_t  Get3DTrkDist()  const
    { return fHit ? f3DTrkPos-fHit->GetPos() : kBig; }
    Double_t  Get3DTrkResid() const { return fHit ? f3DTrkPos-fPos : kBig; }
    // Int_t     GetWireNum()    const { return fHit ? fHit->GetWireNum() : -1; }

  private:
    Hit*      fHit;        // Decoded raw hit data
    Road*     fRoad;       // Road that created this fit
    Double_t  fPos;        // Uncorrected hit position (posL/R) used in fit (m)
    Double_t  fTrackPos;   // Track crossing position from projection fit (m)
    Double_t  fTrackSlope; // Track slope from projection fit
    Double_t  f3DTrkPos;   // Crossing position of fitted 3D track (m)
    Double_t  f3DTrkSlope; // Slope of fitted 3D track
    //TODO: add more members (for correction data etc)

    ClassDef(FitCoord,2) // Coordinate information from road fit
  };

  //___________________________________________________________________________
  // Utility structure for storing sets of hits along with NodeDescriptors

  typedef std::set<Hit*,Hit::PosIsLess> Hset_t;
  struct HitSet {
    Hset_t  hits;          // Hits associated with a pattern
    UInt_t  plane_pattern; // Bit pattern of plane numbers occupied by hits
    UInt_t  nplanes;       // number of active planes
    mutable UInt_t used;   // Pattern has been assigned to a road

    HitSet() : plane_pattern(0), nplanes(0), used(0) {}
    virtual ~HitSet() {}
    void          CalculatePlanePattern();
    static Bool_t CheckMatch( const Hset_t& hits, const TBits* bits );
    Bool_t        CheckMatch( const TBits* bits ) const;
    static UInt_t GetMatchValue( const Hset_t& hits );
    static UInt_t GetAltMatchValue( const Hset_t& hits );
    Bool_t        IsSimilarTo( const HitSet& tryset, Int_t maxdist=0 ) const;

    ClassDef(HitSet, 0)  // A set of hits associated with a pattern
  };

  typedef std::pair<NodeDescriptor,HitSet> Node_t;

  //___________________________________________________________________________
  inline
  Int_t Hit::Compare( const TObject* obj ) const
  {
    // Used for sorting hits in a TSeqCollection (e.g. TList, TClonesArray).
    // A hit is "less than" another hit if its position is smaller.
    // Returns -1 if this is smaller than rhs, 0 if equal, +1 if greater.

    assert( dynamic_cast<const Hit*>(obj) );
    const Hit* rhs = static_cast<const Hit*>(obj);
    assert( fPlane == rhs->fPlane );

    if( fPos  < rhs->fPos )  return  -1;
    if( fPos  > rhs->fPos )  return  1;
    return 0;
  }

  //___________________________________________________________________________
  inline
  Int_t Hit::Compare( const Hit* rhs, Double_t maxdist ) const
  {
    // Determine if two hits are within maxdist of each other.
    // Returns -1 if this<rhs, 0 if overlap, +1 if this>rhs.
    if( fPos+maxdist < rhs->fPos )
      // this hit is "smaller than" (to the left of) rhs
      return -1;
    if( rhs->fPos+maxdist < fPos )
      // this hit is "larger than" (to the right of) rhs
      return +1;
    // The hits overlap within the maxdist tolerance
    return 0;
  }

  //___________________________________________________________________________
  inline
  UInt_t HitSet::GetMatchValue( const Hset_t& hits )
  {
    // Return plane occupancy pattern of given hitset

    UInt_t curpat = 0;
    for( Hset_t::const_iterator it = hits.begin(); it != hits.end(); ++it )
      curpat |= 1U << (*it)->GetPlaneNum();

    return curpat;
  }

  //___________________________________________________________________________
  inline
  Bool_t HitSet::CheckMatch( const Hset_t& hits, const TBits* bits )
  {
    // Check if the plane occupancy pattern of the given hits is marked as
    // allowed in the given bitfield

    return bits->TestBitNumber( GetMatchValue(hits) );
  }

  //___________________________________________________________________________
  inline
  Bool_t HitSet::CheckMatch( const TBits* bits ) const
  {
    // Check if the plane occupancy pattern of the hits in this hitset is
    // marked as allowed in the given bitfield

    assert( plane_pattern || hits.empty() );
    return bits->TestBitNumber(plane_pattern);
  }

  //___________________________________________________________________________
  inline
  void HitSet::CalculatePlanePattern()
  {
    // Set plane_pattern and nplanes according to the hits stored in this HitSet

    plane_pattern = GetMatchValue(hits);
    nplanes = NumberOfSetBits(plane_pattern);
  }

  //___________________________________________________________________________
#ifdef VERBOSE
  inline
  void PrintHits( const Hset_t& hits )
  {
    //  cout << hits.size() << " hits" << endl;
    for( Hset_t::reverse_iterator it = hits.rbegin(); it != hits.rend();
	 ++it ) {
      std::cout << " ";
      (*it)->Print();
    }
  }
#endif

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
