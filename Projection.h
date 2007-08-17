#ifndef ROOT_TreeSearch_Projection
#define ROOT_TreeSearch_Projection

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>

using std::vector;

namespace TreeSearch {

  class Hitpattern;
  class PatternTree;
  class WirePlane;

  class Projection {
  public:
    enum EProjType { kUndefinedType = -1, kUPlane, kTypeBegin = kUPlane,
		     kVPlane, kXPlane, kYPlane, kTypeEnd };

    Projection( EProjType type );
    Projection( const Projection& orig );
    Projection& operator=( const Projection& rhs );
    virtual ~Projection();

    void  Reset();
    Int_t Init( UInt_t depth = 0 );

    EProjType   GetType() const { return fType; }
    //FIXME: needed?
    vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }
    Double_t    GetMaxSlope() const { return fMaxSlope; }
    Double_t    GetWidth() const { return fWidth; }
    Hitpattern* GetHitpattern() const { return fHitpattern; }

    void        SetMaxSlope( Double_t m ) { fMaxSlope = m; }
    void        SetPatternTree( PatternTree* pt ) { fPatternTree = pt; }

  private:
    EProjType           fType;
    vector<WirePlane*>  fPlanes;
    Double_t            fMaxSlope;
    Double_t            fWidth;
    Hitpattern*         fHitpattern;
    PatternTree*        fPatternTree;

    ClassDef(Projection,1)  // Track projection plane
  };

  // Wondrous C++...
  inline
  Projection::EProjType& operator++( Projection::EProjType& e )
  { switch(e) { 
    case Projection::kUndefinedType: return e=Projection::kUPlane;
    case Projection::kUPlane: return e=Projection::kVPlane;
    case Projection::kVPlane: return e=Projection::kXPlane;
    case Projection::kXPlane: return e=Projection::kYPlane;
    case Projection::kYPlane: return e=Projection::kTypeEnd;
    case Projection::kTypeEnd: default: return e=Projection::kUndefinedType;
    }
  }
  inline
  const Projection::EProjType operator++( Projection::EProjType e, int )
  { Projection::EProjType r(e); ++e; return r; }

}

#endif
