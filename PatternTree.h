#ifndef ROOT_TreeSearch_PatternTree
#define ROOT_TreeSearch_PatternTree

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

namespace TreeSearch {

  class Link;

  class PatternTree {
  public:
    PatternTree( UInt_t nlevels, UInt_t nplanes, Link* root );
    virtual ~PatternTree();

    Link* GetRoot() const { return fRoot; }


  private:
    UInt_t     fNlevels;     // Number of levels of the tree (0-nlevels-1)
    UInt_t     fNplanes;     // Number of hitpattern planes

    Link*      fRoot;        // Root node 

    ClassDef(PatternTree,0)   // Precomputed template database
  };


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
