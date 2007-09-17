#ifndef ROOT_TreeSearch_PatternTree
#define ROOT_TreeSearch_PatternTree

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

namespace TreeSearch {

  struct ListNode;

  class TreeNode {
  public:

    //  protected:
    UInt_t      fMinDepth;    // Minimum valid depth for this pattern
    ListNode*   fChild;       // Linked list of child patterns
    union {
      UShort_t* fBitP;        // Bit pattern of this node (unique)
      UShort_t  fBit[2];      // Bit pattern array (NB: size may be >2!!)
    };
  };

  struct ListNode {
    Int_t       fOp;          // Operation to be applied to pattern
    TreeNode*   fTreeNode;    // Bit pattern treenode
    ListNode*   fNext;        // Next list element
  };

  class PatternTree {
  public:
    PatternTree();
    virtual ~PatternTree();

    Int_t Build();

    void Make();

  private:
    UInt_t      fNplanes;     // Number of hitpattern planes
    UInt_t      fDepth;       // Depth of tree (levels = 0 - depth-1)
    Double_t    fMaxSlope;    // Max allowed slope, in normalized units (0-1)
    Double_t*   fZ;           // z positions of planes, normalized (0-1)

    ListNode*   root;         // Root node 

    ClassDef(PatternTree,0)   // Precomputed template database
  };

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
