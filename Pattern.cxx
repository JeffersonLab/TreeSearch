///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Pattern                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Pattern.h"

namespace TreeSearch {

//_____________________________________________________________________________
Pattern::Pattern( UInt_t size )
  : fChild(0), fNbits(size), fMinDepth(-1), fMaxDepth(0), fRefIndex(-1)
{
  // Constructor. Puts bit array on heap. Minimum size is 1.

  if( fNbits == 0 )
    ++fNbits;
  fBits = new UShort_t[fNbits];
  memset( fBits, 0, fNbits*sizeof(UShort_t) );
}

//_____________________________________________________________________________
// Pattern::Pattern( UShort_t* loc )
//: fBits(loc), fChild(0), fNbits(0), fMinDepth(0), fMaxDepth(0), fRefIndex(-1)
// {
//   // Constructor. Puts bit array at loc, a location managed by the caller.
// }

//_____________________________________________________________________________
Pattern::Pattern( const Pattern& orig )
  : fChild(0), fNbits(orig.fNbits), fMinDepth(-1), fMaxDepth(0), fRefIndex(-1)
{
  // Copy constructor. Copy has same bits, no children

  fBits = new UShort_t[fNbits];
  memcpy( fBits, orig.fBits, fNbits*sizeof(UShort_t) );
}

//_____________________________________________________________________________
Pattern& Pattern::operator=( const Pattern& rhs )
{
  // Assignment. Copies only this bits, not the pointers to the children.

  if( this != &rhs ) {
    fChild = 0;
    fNbits = rhs.fNbits;
    fMinDepth = -1;
    fMaxDepth = 0;
    fRefIndex = -1;
    delete fBits;
    fBits = new UShort_t[fNbits];
    memcpy( fBits, rhs.fBits, fNbits*sizeof(UShort_t) );
  }
  return *this;
}

//_____________________________________________________________________________
ListNode* Pattern::AddChild( Pattern* child, Int_t type )
{
  assert(child);
  return (fChild = new ListNode( child, fChild, type ));
}


//_____________________________________________________________________________
Pattern* Pattern::FindChild( const Pattern* pat, Int_t type )
{
  // Search for the given pattern with the given type in the list of child
  // nodes of this pattern

  ListNode* child = fChild;
  while( child ) {
    Pattern* rhs = child->GetPattern();
    assert(rhs);
    if( child->Type() == type && *pat == *rhs )
      return rhs;
    child = child->Next();
  }
  return 0;
}

//_____________________________________________________________________________
Pattern::~Pattern()
{
  // Destroys a Pattern. Deletes all child listnodes (but not the treenodes
  // they point to)

  while( fChild ) {
    ListNode* child = fChild;
    fChild = fChild->Next();
    delete child;
  }
  delete [] fBits;
}


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
