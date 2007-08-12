#ifndef ROOT_TreeSearch_Hit
#define ROOT_TreeSearch_Hit

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "WirePlane.h"
#include <utility>

class TSeqCollection;
class TIterator;

namespace TreeSearch {

  class MCTrack;

  extern const Double_t kBig;

  class Hit : public TObject {

  public:
    Hit() : fWirePlane(NULL) {}
    Hit( Int_t wnum, Double_t pos, Int_t tdc, Double_t time, Double_t res,
	 const WirePlane* wp )
      : fWireNum(wnum), fRawTDC(tdc), fTime(time), fPos(pos), fPosL(pos), 
	fPosR(pos), fResolution(res), fTrackPos(kBig), fWirePlane(wp) {}
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
    Double_t GetZ() const { return fWirePlane ? fWirePlane->GetZ() : kBig; }
    Double_t GetRawTDC()     const { return fRawTDC; }
    Double_t GetDriftTime()  const { return fTime; }
    Double_t GetDriftDist()  const { return fPosR-fPos; }
    Double_t GetPosL()       const { return fPosL; }
    Double_t GetPosR()       const { return fPosR; }
    Double_t GetTrackPos()   const { return fTrackPos; }
    Double_t GetTrackDist()  const { return fTrackPos-fPos; }
    Double_t GetResolution() const { return fResolution; }

  protected:
    Int_t    fWireNum;     // Wire number
    Int_t    fRawTDC;      // Raw TDC value (channels)
    Double_t fTime;        // Hit time corrected for TDC offset (s)
    Double_t fPos;         // Wire position along plane axis (m)
    Double_t fPosL;        // fPos - raw drift distance (m)
    Double_t fPosR;        // fPos + raw drift distance (m)
    Double_t fResolution;  // Resolution of fPosR/fPosL (m)
    Double_t fTrackPos;    // Track crossing position from track fit (m)

    const WirePlane* fWirePlane; //! Pointer to parent wire plane

    ClassDef(Hit,1)        // Horizontal drift chamber hit
  };


  //___________________________________________________________________________
  // Monte Carlo hit class. Same as a hit plus the MC truth info.

  class MCHit : public Hit {

  public:
    MCHit() : fMCTrack(NULL) {}
    MCHit( Int_t wnum, Double_t pos, Int_t tdc, Double_t time, Double_t res,
	   const WirePlane* wp, MCTrack* mctrk, Double_t mcpos )
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

    ObjPair_t& Next();
    ObjPair_t& operator() () { return Next(); }
    const ObjPair_t& Current() const { return fCurrent; }
    const TSeqCollection* GetCollection( Int_t n=0 ) const 
    { return (n==0) ? fCollA : fCollB; }
    void Reset();

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
