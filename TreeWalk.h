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
    UShort_t shift;    // Shift of the base pattern to its actual position
    UChar_t  depth;    // Current recursion depth
    Bool_t   mirrored; // Pattern is mirrored

    NodeDescriptor( Link* ln, UShort_t shft, Bool_t mir, UChar_t dep )
      : link(ln), shift(shft), depth(dep), mirrored(mir) {}
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
    
    template<typename Operation>
    Int_t operator() ( Link* link, Operation& op, UInt_t depth = 0,
		       UInt_t shift = 0, Bool_t mirrored = false ) const;

    ClassDef(TreeWalk, 0)  // Generic traversal function for a PatternTree
  };

  //___________________________________________________________________________
  // Base class for "Visitor" classes to the tree nodes. Left empty instead of 
  // using a pure virtual operator() because the template iterator TreeWalk 
  // performs better than a non-template version using inheritance from 
  // NodeVisitor. Inheriting from this base class is still useful for visitors
  // to gain access to certain classes (like PatternGenerator) as friends.
  class NodeVisitor {
    // virtual Int_t operator() ( const NodeDescriptor& nd ) = 0;
  };

  //___________________________________________________________________________
  template<typename Operation>
  inline
  Int_t TreeWalk::operator()( Link* link, Operation& action, UInt_t depth,
			      UInt_t shift, Bool_t mirrored ) const
  {
    // Traverse the tree and call function object "action" for each link. 
    // The return value from action determines the behavior:
    //  <0: error, return immediately
    //   0: process child nodes recursively
    //  >0: ignore child nodes

    if( depth >= fNlevels )
      return 0;

    Int_t ret = action( NodeDescriptor(link, shift, mirrored, depth) );
    if( ret == 0 ) {
      Link* ln = link->GetPattern()->GetChild();
      while( ln ) {
	// Set up parameters of child pattern based on current position in the
	// tree. The mirroring flag for the child pattern is the pattern's
	// mirroring flag xor the mirror state of the parent (so that 
	// mirrored+mirrored = unmirrored)
	ret = (*this)( ln, action, depth+1, (shift << 1) + ln->Shift(),
		       mirrored xor ln->Mirrored() );
	if( ret ) return ret;
	// Continue along the linked list of child nodes
	ln = ln->Next();
      }
    } else if( ret > 0 )
      ret = 0;
    return ret;
  }

  /////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
