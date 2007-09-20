///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternGenerator                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternGenerator.h"
#include "PatternTree.h"
#include "TMath.h"

ClassImp(TreeSearch::PatternGenerator)

namespace TreeSearch {

  //typedef PatternGenerator::ChildIter chldIter_t;

//_____________________________________________________________________________
inline
PatternGenerator::ChildIter& PatternGenerator::ChildIter::operator++()
{ 
  // Iterator over all suitable child patterns of a given parent pattern.
  // Child pattern bits are either 2*bit or 2*bit+1 of the parent bits,
  // yielding 2^nbits (=2^nplanes) different combinations.
  // The bits of suitable patterns must monotonically increase.
  // bit[0] is always zero (otherwise the pattern could be shifted).
  // Note that type() is an important part of the result.
  // type & 1 indicates a pattern shifted by one to the right
  // type & 2 indicates a mirrored pattern 
  // To recover the pattern, mirror first, then shift left, as approriate.

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

//_____________________________________________________________________________
PatternGenerator::PatternGenerator() : fDepth(0), fNplanes(0), fMaxSlope(0)
{
  // Constructor. 

}

//_____________________________________________________________________________
PatternGenerator::~PatternGenerator()
{
  // Destructor

  DoTree( kDelete );

}

//_____________________________________________________________________________
Int_t PatternGenerator::DoTree( EOperation op )
{
  // Execute given operation on all unique Pattern nodes of the build tree.
  // Internal utility function.

  Int_t result = 0;
  if( op == kBytesRequired ) {
    Int_t npatt = DoTree( kCountPatterns );
    result += sizeof(Pattern) * npatt;
    result += sizeof(UShort_t) * npatt * fNplanes;
    result += sizeof(ListNode) * DoTree( kCountChildNodes );
  } else {
    for( vector<ListNode*>::iterator it = fHashTable.begin();
	 it != fHashTable.end(); ++it ) {
      // Each hashnode points to a unique pattern by construction of the table
      ListNode* hashnode = *it;
      Int_t count = 0;
      while( hashnode ) {
	ListNode* cur_node = hashnode;
	hashnode = hashnode->Next();
	assert( cur_node->GetPattern() );
	switch( op ) {
	case kDelete:
	  delete cur_node->GetPattern();
	  delete cur_node;
	  break;
	case kResetRefIndex:
	  cur_node->GetPattern()->fRefIndex = -1;
	  break;
	case kCountPatterns:
	  ++result;
	  break;
	case kCountChildNodes:
	  {
	    Pattern* pat = cur_node->GetPattern();
	    assert(pat);
	    ListNode* ln = pat->fChild;
	    while( ln ) {
	      ++result;
	      ln = ln->Next();
	    }
	  }
	  break;
	case kMaxHashDepth:
	  ++count;
	  break;
	default:
	  break;
	}
      } // while hashnode
      if( op == kMaxHashDepth && count > result )
	result = count;
    } // for hashtable elements
    if( op == kDelete )
      fHashTable.clear();
  }
  return result;
}

//_____________________________________________________________________________
void PatternGenerator::Print( Option_t* opt ) const
{
  // Print information about the tree, depending on option

}

//_____________________________________________________________________________
void PatternGenerator::AddHash( Pattern* pat )
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
PatternTree* PatternGenerator::Generate( UInt_t depth, Double_t width,
					 const vector<double>& zpos,
					 Double_t maxslope )
{
  // Generate a new pattern tree for the given parameters. Returns a pointer
  // to the generated tree, or zero if error

  // Clear out previous build, if any
  DoTree( kDelete );

  // Set parameters for the new build.
  //TODO: normalize zpos and maxslope
  fDepth    = depth;
  fZ        = zpos;
  fNplanes  = fZ.size();
  fMaxSlope = maxslope;
  // TODO: check for unreasonable input

  // 2^(depth-1) * 2^(nplanes-2) is the upper limit for the number of patterns.
  // The following table size should give decent speed. Anything larger would
  // require a cleverer hash function.
  fHashTable.resize( 1<<(fDepth-1), 0 );

  // Start with the trivial all-zero root node at depth 0. 
  Pattern* root = new Pattern( fNplanes );
  AddHash( root );

  // Generate the tree recursively
  MakeChildNodes( root, 1 );

  
  return 0;
}

//_____________________________________________________________________________
bool PatternGenerator::TestSlope( const Pattern& pat, UInt_t depth )
{
  UInt_t width = pat.GetWidth();
  return ( width < 2 ||
	   TMath::Abs((double)(width-1) / (double)(1<<depth)) <= fMaxSlope );
}

//_____________________________________________________________________________
bool PatternGenerator::LineCheck( const Pattern& pat )
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
      // If dL>0, the right edge of the bin is inside the band,
      // so set a new right-side limit.
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
Pattern* PatternGenerator::Find( const Pattern& pat )
{
  // Search for the given pattern in the current database

  UInt_t hashsize = fHashTable.size();
  assert(hashsize);
  Int_t hash = pat.Hash()%hashsize;
  ListNode* listnode = fHashTable[hash];
  while( listnode ) {
    Pattern* rhs = listnode->GetPattern();
    if( pat == *rhs )
      return rhs;
    listnode = listnode->Next();
  }
  return 0;
}

//_____________________________________________________________________________
void PatternGenerator::MakeChildNodes( Pattern* parent, UInt_t depth )
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
