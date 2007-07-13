#ifndef ROOT_TreeSearch_Hit
#define ROOT_TreeSearch_Hit

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

namespace TreeSearch {

  class Hit : public TObject {

  public:
    Hit() {}
    Hit( const Hit& );
    Hit& operator=( const Hit& );
    virtual ~Hit() {}

#if 0
    virtual Double_t ConvertTimeToDist(Double_t slope);
    Int_t  Compare ( const TObject* obj ) const;
    Bool_t IsSortable () const { return kTRUE; }
  
    // Get and Set Functions
    THaVDCWire* GetWire() const { return fWire; }
    Int_t    GetWireNum() const { return fWire->GetNum(); }
    Int_t    GetRawTime() const { return fRawTime; }
    Double_t GetTime()    const { return fTime; }
    Double_t GetDist()    const { return fDist; }
    Double_t GetPos()     const { return fWire->GetPos(); } //Position of hit wire
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
}

///////////////////////////////////////////////////////////////////////////////

#endif
