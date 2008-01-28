///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::NodeVisitor abstract base class                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "NodeVisitor.h"

using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
void NodeVisitor::SetLinkPattern( Link* link, Pattern* pattern ) {
  link->fPattern = pattern;
}

void NodeVisitor::SetLinkType( Link* link, Int_t type ) {
  link->fOp = type;
}

void NodeVisitor::SetLinkNext( Link* link, Link* next ) {
  link->fNext = next;
}

void NodeVisitor::SetPatternChild( Pattern* pat, Link* link ) {
  pat->fChild = link;
  // Explicitly set pointer is managed externally
  pat->fDelChld = false;
}
 
///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
