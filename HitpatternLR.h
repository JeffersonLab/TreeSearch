//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 11-Feb-2013
//
#ifndef ROOT_TreeSearch_HitpatternLR
#define ROOT_TreeSearch_HitpatternLR

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::HitpatternLR                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hitpattern.h"

namespace TreeSearch {

  class HitpatternLR : public Hitpattern {

  public:
    explicit HitpatternLR( const PatternTree& pt );  // preferred constructor
    HitpatternLR( UInt_t nlevels, UInt_t nplanes, Double_t width );
    explicit HitpatternLR( const Hitpattern& orig );
    Hitpattern& operator=( const Hitpattern& rhs );
    virtual ~HitpatternLR();

    virtual Int_t Fill( const std::vector<TreeSearch::Plane*>& planes );
    virtual Int_t ScanHits( Plane* A, Plane* B = 0 );

    ClassDef(HitpatternLR,0)  // Hitpattern filled by L/R-ambiguous wire hits
  };

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch


#endif
