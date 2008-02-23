#ifndef ROOT_TreeSearch_Hit
#define ROOT_TreeSearch_Hit

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TimeToDistConv.h"
#include "WirePlane.h"
#include <utility>
#include <set>
#include <cassert>
#include <functional>

class TSeqCollection;
class TIterator;

namespace TreeSearch {

  class Road;
  extern const Double_t kBig;

  class Hit : public TObject {

  public:
    Hit() : fWirePlane(0) {}
    Hit( Int_t wnum, Double_t pos, Int_t tdc, Double_t time, Double_t res,
	 WirePlane* wp ) :
      fWireNum(wnum), fRawTDC(tdc), fTime(time), fPos(pos), fPosL(pos), 
      fPosR(pos), fResolution(res), fWirePlane(wp)
#ifdef TESTCODE
      , fCl(0), fMulti(0), fTdiff(0.0)
#endif
    { assert(fWirePlane); }
    // Default copy and assignment are fine
    //    Hit( const Hit& );
    //    Hit& operator=( const Hit& );
    virtual ~Hit() {}

    virtual Int_t Compare( const TObject* obj ) const;
    Int_t  Compare( const Hit* rhs, Double_t maxdist ) const;
    Bool_t IsSortable () const { return kTRUE; }
    virtual void Print( Option_t* opt="" ) const;
 
    Double_t ConvertTimeToDist( Double_t slope );

    Int_t    GetWireNum()    const { return fWireNum; }
    Double_t GetWirePos()    const { return fPos; }
    Double_t GetZ()          const;
    Double_t GetRawTDC()     const { return fRawTDC; }
    Double_t GetDriftTime()  const { return fTime; }
    Double_t GetDriftDist()  const { return fPosR-fPos; }
    Double_t GetPosL()       const { return fPosL; }
    Double_t GetPosR()       const { return fPosR; }
    Double_t GetResolution() const { return fResolution; }

    WirePlane* GetWirePlane() const { return fWirePlane; }

    // Functors for ordering hits in sets
    struct WireNumLess : public std::binary_function< Hit*, Hit*, bool >
    {
      bool operator() ( const Hit* a, const Hit* b ) const
      { 
	assert( a && b );
	const WirePlane* A = a->GetWirePlane(), *B = b->GetWirePlane();
	assert( A && B );
	if( A->GetPlaneNum() < B->GetPlaneNum() ) return true;
	if( A->GetPlaneNum() > B->GetPlaneNum() ) return false;
	if( a->GetWireNum()  < b->GetWireNum()  ) return true;
	if( a->GetWireNum()  > b->GetWireNum()  ) return false;
	return ( a->GetDriftTime() < b->GetDriftTime() );
      }
    };

    struct WireDistLess : public std::binary_function< Hit*, Hit*, bool >
    {
      WireDistLess( Int_t maxdist ) : fMaxDist(maxdist) {}
      bool operator() ( const Hit* a, const Hit* b ) const
      { 
	assert( a && b );
	const WirePlane* A = a->GetWirePlane(), *B = b->GetWirePlane();
	assert( A && B );
	if( A->GetPlaneNum() < B->GetPlaneNum() ) return true;
	if( A->GetPlaneNum() > B->GetPlaneNum() ) return false;
	return ( a->GetWireNum() + fMaxDist < b->GetWireNum() );
      }
      Int_t GetMaxDist() const { return fMaxDist; }
    private:
      Int_t fMaxDist;      // Max allowed distance between hits in a cluster
    };

  protected:
    Int_t    fWireNum;     // Wire number
    Int_t    fRawTDC;      // Raw TDC value (channels)
    Double_t fTime;        // Hit time corrected for TDC offset (s)
    Double_t fPos;         // Wire position along plane axis (m)
    Double_t fPosL;        // fPos - raw drift distance (m)
    Double_t fPosR;        // fPos + raw drift distance (m)
    Double_t fResolution;  // Resolution of fPosR/fPosL (m)

    WirePlane* fWirePlane; //! Pointer to parent wire plane

#ifdef TESTCODE
    Int_t    fCl;          // Neighboring wire also fired
    Int_t    fMulti;       // Additional hits present on same wire
    Double_t fTdiff;       // Time difference to previous multihit
    friend void WirePlane::CheckCrosstalk();
#endif

    ClassDef(Hit,1)        // Horizontal drift chamber hit
  };


  //___________________________________________________________________________
  // Monte Carlo hit class. Same as a hit plus the MC truth info.
  //TODO: Skeleton version - to be fleshed out

  class MCTrack;
  class MCHit : public Hit {

  public:
    MCHit() : fMCTrack(0) {}
    MCHit( Int_t wnum, Double_t pos, Int_t tdc, Double_t time, Double_t res,
	   WirePlane* wp, MCTrack* mctrk, Double_t mcpos )
      : Hit(wnum, pos, tdc, time, res, wp), fMCTrack(mctrk), fMCPos(mcpos) {}
    virtual ~MCHit() {}

    virtual void Print( Option_t* opt="" ) const;

    MCTrack* GetMCTrack() const { return fMCTrack; }
    Double_t GetMCPos()   const { return fMCPos; }

  protected:
    MCTrack* fMCTrack;     // MC track generating this hit (0=noise hit)
    Double_t fMCPos;       // Exact MC track crossing position (m)

    ClassDef(MCHit,1)      // Monte Carlo hit in horizontal drift chamber
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
  // A given raw Hit may be associated with zero or more FitCoord objects.

  class FitCoord : public TObject {

  public:
    FitCoord( Hit* hit, Road* road, UInt_t ifit,
	      Double_t pos, Double_t trackpos, Double_t trackslope ) 
      : fHit(hit), fRoad(road), fFitRank(ifit), fPos(pos), fTrackPos(trackpos),
	fTrackSlope(trackslope)
    { assert(fHit&&fRoad); }
    FitCoord() {} // For internal ROOT use
    virtual ~FitCoord() {}

    Double_t  GetChi2()       const;
    Hit*      GetHit()        const { return fHit; }
    Road*     GetRoad()       const { return fRoad; }
    UInt_t    GetRank()       const { return fFitRank; }
    Double_t  GetPos()        const { return fPos; }
    Double_t  GetDriftDist()  const { return fHit->GetDriftDist(); }
    Double_t  GetTrackPos()   const { return fTrackPos; }
    Double_t  GetTrackSlope() const { return fTrackSlope; }
    Double_t  GetTrackDist()  const { return fTrackPos-fHit->GetWirePos(); }
    Double_t  GetResidual()   const { return fTrackPos-fPos; }

  private:
    Hit*      fHit;      // Decoded raw hit data
    Road*     fRoad;     // Road that created this fit
    UInt_t    fFitRank;  // Fit rank by chi2 within the road (0 = best)
    Double_t  fPos;      // Uncorrected hit position (posL/R) used in fit (m)
    Double_t  fTrackPos; // Track crossing position from original track fit (m)
    Double_t  fTrackSlope; // Track slope from original track fit (m)
    //TODO: add more members (for correction data etc)

    ClassDef(FitCoord,1) // Coordinate information from road fit
  };

  //___________________________________________________________________________
  // Utility structure for storing sets of hits along with NodeDescriptors

  typedef std::set<Hit*,Hit::WireNumLess> Hset_t;
  struct HitSet {
    Hset_t  hits;  // Hits associated with a pattern
    UInt_t  used;  // 0=not in any road, 1=some hits used, 2=all hits used
    HitSet() : used(0) {}
    virtual ~HitSet() {}
    ClassDef(HitSet, 0)  // A set of hits associated with a pattern
  };

  //___________________________________________________________________________
  inline
  Int_t Hit::Compare( const TObject* obj ) const 
  {
    // Used for sorting hits in a TSeqCollection (e.g. TList, TClonesArray).
    // A hit is "less than" another hit if its position is smaller.
    // Returns -1 if this is smaller than rhs, 0 if equal, +1 if greater.
    // Also, for hits on the same wire, the first hit on the wire (the one with
    // the smallest time) is "less than" one with a higher time.  If the hits
    // are sorted according to this scheme, they will be in order of increasing
    // wire number and, for each wire, will be in the order in which they hit
    // the wire

    // For dissimilar objects, compare addresses for lack of anything better
//     if( !obj || IsA() != obj->IsA() )
//       return (this<obj) ? -1 : (this==obj) ? 0 : 1;

    if( fPos  < static_cast<const Hit*>(obj)->fPos )  return -1;
    if( fPos  > static_cast<const Hit*>(obj)->fPos )  return  1;
    if( fTime < static_cast<const Hit*>(obj)->fTime ) return -1;
    if( fTime > static_cast<const Hit*>(obj)->fTime ) return  1;
    return 0;
  }

  //___________________________________________________________________________
  inline
  Int_t Hit::Compare( const Hit* rhs, Double_t maxdist ) const {
    // Determine if two hits are within maxdist of each other.
    // Returns -1 if this<rhs, 0 if overlap, +1 if this>rhs.
    // Overlap is necessary but NOT sufficient for two hits to be a pair;
    // additional comparisons of the L/R positions are necessary for that.
    if( fPosR+maxdist < rhs->fPosL )
      // this hit is "smaller than" (to the left of) rhs
      return -1;
    if( rhs->fPosR+maxdist < fPosL )
      // this hit is "larger than" (to the right of) rhs
      return +1;
    // The hits overlap within the maxdist tolerance
    return 0;
  }

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
