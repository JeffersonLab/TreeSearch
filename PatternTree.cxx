///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternTree.h"
#include "TMath.h"
#include <iostream>

ClassImp(TreeSearch::PatternTree)

using namespace std;
namespace TreeSearch {

// Iterator over all suitable child patterns of a parent pattern.
// Child pattern bits are either 2*bit or 2*bit+1 of the parent bits,
// yielding 2^nbits (=2^nplanes) different combinations.
// The bits of suitable patterns must monotonically increase.
// bit[0] is always zero (otherwise the pattern could be shifted).
// Note that type() is an important part of the result.
// type & 1 indicates a pattern shifted by one to the right
// type & 2 indicates a mirrored pattern 
// To recover the pattern, mirror first, then shift left, as approriate.
class ChildIter {
private:
  const Pattern fParent;
  Pattern   fChild;
  Int_t     fCount;
  Int_t     fType;
public:
  ChildIter( const Pattern& parent ) 
    : fParent(parent), fChild(parent.GetNbits()), 
      fCount(1<<(parent.GetNbits())), fType(0)  { ++(*this); }
  ChildIter& operator++() { 
    if( fCount > 0 ) {
      while( fCount-- ) {
	Int_t maxbit = 0;
	Int_t minbit = 1;
	UInt_t nbits = fChild.GetNbits();
	for( UInt_t ibit = nbits; ibit--; ) {
	  Int_t bit = fParent[ibit] << 1;
	  if( fCount & (1<<ibit) )
	    ++bit;
	  fChild[ibit] = bit;
	  if( bit < minbit ) minbit = bit;
	  if( bit > maxbit ) maxbit = bit;
	}
	Int_t width = fChild.GetWidth();
	if( maxbit-minbit > TMath::Abs(width) )
	  continue;
	if( minbit == 0 )
	  fType = 0;
	else {
	  fType = 1;
	  for( UInt_t ibit = nbits; ibit; )
	    --fChild[--ibit];
	}
	if( width < 0 ) {
	  fType += 2;
	  width = -width;
	  for( UInt_t ibit = nbits; ibit--; )
	    fChild[ibit] = width-fChild[ibit];
	}
	break;
      }
    } else
      fCount = -1;
    return *this;
  }
  const ChildIter operator++(int) { 
    ChildIter clone(*this);
    ++(*this); 
    return clone;
  }
  Pattern& operator*() { return fChild; }
  operator bool() const { return (fCount >= 0); }
  bool operator!() const { return !((bool)*this); }
  Int_t type() const { return fType; }
};

//_____________________________________________________________________________
Pattern::Pattern( UInt_t size )
  : fChild(0), fNbits(size), fMinDepth(-1), fMaxDepth(0), fRefIndex(-1)
{
  // Constructor. Puts bit array on heap.

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
    Pattern* rhs = child->fPattern;
    assert(rhs);
    if( child->fOp == type && *pat == *rhs )
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

//_____________________________________________________________________________
PatternTree::PatternTree( UInt_t depth, const vector<double>& zpos )
  : fDepth(depth), fNplanes(zpos.size()), fMaxSlope(kBig), fZ(zpos), fRoot(0)
{
  // Constructor. 
  // The number of planes is defined by the size of the zpos vector.

}

//_____________________________________________________________________________
PatternTree::~PatternTree()
{
  // Destructor

}

//_____________________________________________________________________________
void PatternTree::DeleteBuildTree()
{
  for( vector<ListNode*>::iterator it = fHashTable.begin();
       it != fHashTable.end(); ++it ) {
    ListNode* listnode = *it;
    while( listnode ) {
      delete listnode->fPattern;
      ListNode* prev_node = listnode;
      listnode = listnode->Next();
      delete prev_node;
    }
  }
}

//_____________________________________________________________________________
void PatternTree::AddHash( Pattern* pat )
{
  // Add given pattern to the hash table (used during database generation)

  UInt_t hashsize = fHashTable.size();
  assert(pat);
  assert(hashsize);
  Int_t hash = pat->Hash()%hashsize;
  ListNode* prev_node = fHashTable[hash];
  fHashTable[hash] = new ListNode( pat, prev_node, 0 );
}

//_____________________________________________________________________________
Int_t PatternTree::Init()
{
  // Initialize the pattern tree

  // TODO: try to read from disk first. If not successful, then do this:
  Pattern* root = new Pattern( fNplanes );
  fHashTable.clear();
  fHashTable.resize( 1<<(fDepth-1), 0 );
  AddHash( root );

  MakeChildNodes( root, 1 );

  // TODO: write the tree to disk and read it back

  DeleteBuildTree();

  return 0;
}

//_____________________________________________________________________________
bool PatternTree::TestSlope( const Pattern& pat, UInt_t depth )
{
  UInt_t width = pat.GetWidth();
  return ( width < 2 ||
	   TMath::Abs((double)(width-1) / (double)(1<<depth)) <= fMaxSlope );
}

//_____________________________________________________________________________
bool PatternTree::LineCheck( const Pattern& pat )
{
  // Check if the gievn bit pattern is consistent with a straight line.
  // The intersection plane positions are given by fZ[].
  // Assumes a normalized pattern, for which pat[0] is always zero.

  assert(fNplanes);
  Double_t xL   = pat[fNplanes-1];
  Double_t xRm1 = xL;               // xR-1
  Double_t zL   = fZ[fNplanes-1];
  Double_t zR   = zL;

  for( Int_t i = fNplanes-2; i > 0; --i ) {
    // Compare the intersection point with the i-th plane of the left edge 
    // of the band, (xL-x0) * z[i]/zL, to the left edge of the bin, pat[i]-x0. 
    // If the difference is larger than one bin width (=1), the bin is
    // outside of the allowed band.
    // Multiply with zL (to avoid division) and recall x0 = 0.
    Double_t dL = xL*fZ[i] - pat[i]*zL;
    if( TMath::Abs(dL) >= zL )
      return false;
    // Likewise for the right edge
    Double_t dR = xRm1*fZ[i] - pat[i]*zR;
    if( TMath::Abs(dR) >= zR )
      return false;

    if( i > 1 ) {
      // If dL>0, the right edge of the bin is inside the band, so set a
      // new right-side limit.
      if( dL > 0 ) {
	xRm1 = pat[i];
	zR   = fZ[i];
      }
      // Likewise for the left-side limit
      if( dR < 0 ) {
	xL = pat[i];
	zL = fZ[i];
      }
    }
  } // planes
  return true;
}

//_____________________________________________________________________________
Pattern* PatternTree::Find( const Pattern& pat )
{
  // Search for the given pattern in the hash table (used during build)

  UInt_t hashsize = fHashTable.size();
  assert(hashsize);
  Int_t hash = pat.Hash()%hashsize;
  ListNode* listnode = fHashTable[hash];
  while( listnode ) {
    Pattern* rhs = listnode->fPattern;
    if( pat == *rhs )
      return rhs;
    listnode = listnode->Next();
  }
  return 0;
}

//_____________________________________________________________________________
void PatternTree::MakeChildNodes( Pattern* parent, UInt_t depth )
{
  // Generate child nodes for the given parent pattern

  // Requesting child nodes for the parent at this depth implies that the 
  // parent is being used at the level above
  if( depth > 0 )
    parent->UsedAtDepth( depth-1 );

  // Base case of the recursion: no child nodes beyond fDepth-1
  if( depth >= fDepth )
    return;

  // Iterate over child patterns of the parent
  ChildIter it( *parent );
  while( it ) {
    Pattern& child = *it;
    bool insert = true;

    // Pattern already exists?
    Pattern* node = Find( child );
    if( node ) {
      // If this pattern exists, but has not been used at this depth before,
      // mark it for this depth and generate child nodes for it
      bool is_lower  = depth < node->fMinDepth;
      bool is_higher = depth > node->fMaxDepth;
      // If the pattern has only been tested at a higher depth, we need 
      // to redo the slope test since the absolute slope increases with
      // decreasing depth
      if( is_lower )
	is_lower = is_lower && TestSlope( *node, depth );

      if( is_higher || is_lower )
	MakeChildNodes( node, depth+1 );

    } else {
      // If the pattern is new, check it for consistency with maxslope
      // and the straight line condition. If good, add it to the database
      // and generate its child nodes.
      if( TestSlope(child, depth) && LineCheck(child) ) {
	node = new Pattern( child );
	// This pattern is guaranteed to be a added as a child node at
	// this depth or below, either here through the recursive call
	// or below as a child node of the current parent. Therefore, we
	// can add it to the hashtable here.
	AddHash( node );
	MakeChildNodes( node, depth+1 );
	
      } else
	insert = false;
    }
    // If this candidate child pattern was found suitable, add it as a
    // child node to the parent - provided, it isn't already there.
    // (How could it get there?!?)
    if( insert && !parent->FindChild(node, it.type()) )
      parent->AddChild(node, it.type());

    ++it;
  }
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
