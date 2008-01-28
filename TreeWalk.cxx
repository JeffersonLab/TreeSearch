///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::TreeWalk                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TreeWalk.h"

ClassImp(TreeSearch::TreeWalk)

namespace TreeSearch {

//_____________________________________________________________________________
NodeVisitor::ETreeOp
TreeWalk::operator()( Link* link, NodeVisitor& action, Pattern* parent, 
		      UInt_t depth, UInt_t shift, Bool_t mirrored ) const
{
  // Traverse the tree and call function object "action" for each link. 
  // The return value from action determines the behavior:
  //  kRecurse: process child nodes until reaching maxdepth
  //  kRecurseUncond: process child nodes (regardless of depth)
  //  fSkipChildNodes: ignore child nodes
  //  kError: error, return immediately

  if( !link ) return NodeVisitor::kError;
  NodeVisitor::ETreeOp ret = 
    action(NodeDescriptor(link, parent, shift, mirrored, depth));
  if( ret == NodeVisitor::kRecurseUncond or 
      ( ret == NodeVisitor::kRecurse and depth+1 < fNlevels ) ) {
    Pattern* pat = link->GetPattern();
    Link* ln = pat->GetChild();
    while( ln ) {
      // Set up parameters of child pattern based on current position in the
      // tree. The mirroring flag for the child pattern is the pattern's
      // mirroring flag xor the mirror state of the parent (so that 
      // mirrored+mirrored = unmirrored). The shift corresponds either
      // to the pattern's left or right edge for unmirrored or mirrored
      // patterns, respectively.
      Bool_t new_mir = mirrored xor ln->Mirrored();
      ret = (*this)( ln, action, pat, depth+1, 
		     (shift << 1) + (new_mir xor ln->Shift()), new_mir );
      if( ret == NodeVisitor::kError ) return ret;
      // Continue along the linked list of child nodes
      ln = ln->Next();
    }
  }
  return ret;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
