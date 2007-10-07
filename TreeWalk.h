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

#include "Pattern.h"

namespace TreeSearch {

  // Structure to describe a tree node (an actual pattern, including shifts
  // and mirroring) 
  struct NodeDescriptor {
    Link*    link;     // Linked-list node pointing to a base pattern
    Pattern* parent;   // Parent node
    UShort_t shift;    // Shift of the base pattern to its actual position
    UChar_t  depth;    // Current recursion depth
    Bool_t   mirrored; // Pattern is mirrored

    NodeDescriptor( Link* ln, Pattern* p, UShort_t shft, Bool_t mir, 
		    UChar_t dep )
      : link(ln), parent(p), shift(shft), depth(dep), mirrored(mir) {}
  };    

  // "Operation" must be a function object whose operator() takes a
  // NodeDescriptor or const NodeDescriptor& as argument
  class TreeWalk {
  private:
    UInt_t   fNlevels;  // Number of levels in tree
  public:
    TreeWalk( UInt_t nlevels = 0 ) : fNlevels(nlevels) {}
    virtual ~TreeWalk() {}
    void SetNlevels( UInt_t n ) { fNlevels = n; }
    
    // Supported return codes from Operation function object
    enum ETreeOp { kError, kRecurse, kRecurseUncond, kSkipChildNodes };

    template<typename Operation>
    ETreeOp operator() ( Link* link, Operation& op, Pattern* parent = 0,
			 UInt_t depth = 0, UInt_t shift = 0,
			 Bool_t mirrored = false ) const;
    ClassDef(TreeWalk, 0)  // Generic traversal function for a PatternTree
  };

  //___________________________________________________________________________
  template<typename Operation>
  inline TreeWalk::ETreeOp
  TreeWalk::operator()( Link* link, Operation& action, Pattern* parent, 
			UInt_t depth, UInt_t shift, Bool_t mirrored ) const
  {
    // Traverse the tree and call function object "action" for each link. 
    // The return value from action determines the behavior:
    //  kError: error, return immediately
    //  kRecurseToMaxdepth: process child nodes until reaching maxdepth
    //  kRecurse: process child nodes (regardless of depth)
    //  fSkipChildNodes: ignore child nodes

    if( !link ) return TreeWalk::kError;
    ETreeOp ret = action(NodeDescriptor(link, parent, shift, mirrored, depth));
    if( ret == TreeWalk::kRecurseUncond or
	( ret == TreeWalk::kRecurse and depth+1 < fNlevels ) ) {
      Pattern* pat = link->GetPattern();
      Link* ln = pat->GetChild();
      while( ln ) {
	// Set up parameters of child pattern based on current position in the
	// tree. The mirroring flag for the child pattern is the pattern's
	// mirroring flag xor the mirror state of the parent (so that 
	// mirrored+mirrored = unmirrored)
	ret = (*this)( ln, action, pat, depth+1, (shift << 1) + ln->Shift(),
		       mirrored xor ln->Mirrored() );
	if( ret == TreeWalk::kError ) return ret;
	// Continue along the linked list of child nodes
	ln = ln->Next();
      }
    }
    return ret;
  }

  /////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
