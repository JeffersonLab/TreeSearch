#ifndef TreeSearch__NodeVisitor
#define TreeSearch__NodeVisitor

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::NodeVisitor                                                   //
//                                                                           //
// Base class for "Visitors" to the patterm tree nodes                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TreeWalk.h"
#include "Pattern.h"

namespace TreeSearch {

  class NodeVisitor {
  public:
    virtual TreeWalk::ETreeOp operator() ( const NodeDescriptor& nd ) = 0;
    virtual ~NodeVisitor() {}
  protected:
    // Functions to allow derived classes to modify Links and Patterns.
    // (needed since C++ friend access doesn't extend to derived classes)
    static void SetLinkPattern( Link* link, Pattern* pattern ) {
      link->fPattern = pattern;
    }
    static void SetLinkType( Link* link, Int_t type ) {
      link->fOp = type;
    }
    static void SetLinkNext( Link* link, Link* next ) {
      link->fNext = next;
    }
    static void SetPatternChild( Pattern* pat, Link* link ) {
      pat->fChild = link;
      // Explicitly set pointer is managed externally
      pat->fDelChld = false;
    }
  };

} // end namespace TreeSearch

#endif
