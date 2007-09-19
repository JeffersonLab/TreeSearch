#ifndef ROOT_TreeSearch_PatternTree
#define ROOT_TreeSearch_PatternTree

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

namespace TreeSearch {

  class ListNode;

  class PatternTree {
  public:
    PatternTree( UInt_t depth, UInt_t nplanes, ListNode* root );
    virtual ~PatternTree();

    ListNode* GetRoot() const { return fRoot; }


  private:
    UInt_t         fDepth;       // Depth of tree (levels = 0 - depth-1)
    UInt_t         fNplanes;     // Number of hitpattern planes

    ListNode*      fRoot;        // Root node 

    ClassDef(PatternTree,0)   // Precomputed template database
  };


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
