#ifndef ROOT_TreeSearch_Hit
#define ROOT_TreeSearch_Hit

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include <utility>

class TSeqCollection;
class TIterator;

namespace TreeSearch {

  class WirePlane;

  class Hit : public TObject {

  public:
    Hit() {}
    Hit( const Hit& );
    Hit& operator=( const Hit& );
    virtual ~Hit() {}

    virtual Int_t Compare( const TObject* obj ) const;
    Int_t Compare( const Hit* rhs, Double_t maxdist = 0.0 ) const;
    Bool_t IsSortable () const { return kTRUE; }

  protected:

    Double_t fPos;
    Double_t fPosL;
    Double_t fPosR;

    WirePlane* fWirePlane;


#if 0
    virtual Double_t ConvertTimeToDist(Double_t slope);
  
    // Get and Set Functions
    THaVDCWire* GetWire() const { return fWire; }
    Int_t    GetWireNum() const { return fWire ? fWire->GetNum() : -1; }
    Int_t    GetRawTime() const { return fRawTime; }
    Double_t GetTime()    const { return fTime; }
    Double_t GetDist()    const { return fDist; }
    //Position of hit wire
    Double_t GetPos()     const { return fWire ? fWire->GetPos() : 1e38; } 
    Double_t GetdDist()   const { return fdDist; }

    void     SetWire(THaVDCWire * wire) { fWire = wire; }
    void     SetRawTime(Int_t time)     { fRawTime = time; }
    void     SetTime(Double_t time)     { fTime = time; }
    void     SetDist(Double_t dist)     { fDist = dist; }
    void     SetdDist(Double_t ddist)   { fdDist = ddist; }
    void     SetFitDist(Double_t dist)  { ftrDist = dist; }

  
  protected:
    static const Double_t kBig;  //!
  
    THaVDCWire* fWire;     // Wire on which the hit occurred
    Int_t       fRawTime;  // TDC value (channels)
    Double_t    fTime;     // Time corrected for time offset of wire (s)
    Double_t    fDist;     // (Perpendicular) Drift Distance
    Double_t    fdDist;    // uncertainty in fDist (for chi2 calc)
    Double_t    ftrDist;   // (Perpendicular) distance from the track
#endif
  
    ClassDef(Hit,1)        // MWDC Hit class
  };

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
    ObjPair_t& Current()     { return fCurrent; }
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

    ClassDef(HitPairIter,0)  // Iterator over collections of hits
  };

//_____________________________________________________________________________
  inline
  Int_t Hit::Compare( const Hit* rhs, Double_t maxdist ) const {
    // Compare if two hits are within maxdist of each other
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
