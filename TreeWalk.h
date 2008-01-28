#ifndef TreeSearch__TreeWalk
#define TreeSearch__TreeWalk

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::TreeWalk                                                      //
//                                                                           //
// Generic function object to traverse a PatternTree.                        //
//                                                                           //
// This implements an internal iterator pattern [E. Gamma et al.,            //
// "Design Patterns", Addison-Wesley, 1995] that applies generic operation   //
// to each tree element.  The operation itself represents a simplified form  //
// of a visitor pattern [ibid.], where the simplification lies in the fact   //
// that the tree contains only a single class of nodes instead of many.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "NodeVisitor.h"

namespace TreeSearch {

  class TreeWalk {
  private:
    UInt_t   fNlevels;  // Number of levels in tree
  public:
    TreeWalk( UInt_t nlevels = 0 ) : fNlevels(nlevels) {}
    virtual ~TreeWalk() {}
    void SetNlevels( UInt_t n ) { fNlevels = n; }
    
    NodeVisitor::ETreeOp 
    operator() ( Link* link, NodeVisitor& op, Pattern* parent = 0,
		 UInt_t depth = 0, UInt_t shift = 0,
		 Bool_t mirrored = false ) const;
    ClassDef(TreeWalk, 0)  // Generic traversal function for a PatternTree
  };

} // end namespace TreeSearch

#endif
