#ifndef TreeSearch__NodeVisitor
#define TreeSearch__NodeVisitor

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::NodeVisitor                                                   //
//                                                                           //
// Base class for "Visitors" to the patterm tree nodes                       //
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

    // Comparison operators for inserting into ordered STL containers.
    // Order nodes by ascending start bin.
    // NB: pattern bit[0] = 0, so nothing needs to be added to "shift"
    bool operator<( const NodeDescriptor& rhs ) const {
      return ( shift < rhs.shift );
    }
    bool operator==( const NodeDescriptor& rhs ) const {
      return ( shift == rhs.shift );
    }
  };    

  class NodeVisitor {
  public:
    enum ETreeOp { kRecurse, kRecurseUncond, kSkipChildNodes, kError };

    virtual ETreeOp operator() ( const NodeDescriptor& nd ) = 0;
    virtual ~NodeVisitor() {}
  protected:
    // Functions to allow derived classes to modify Links and Patterns.
    // (needed since C++ friend access doesn't extend to derived classes)
    static void SetLinkPattern( Link* link, Pattern* pattern );
    static void SetLinkType( Link* link, Int_t type );
    static void SetLinkNext( Link* link, Link* next );
    static void SetPatternChild( Pattern* pat, Link* link );
  };

} // end namespace TreeSearch

#endif
