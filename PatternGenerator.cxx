///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternGenerator                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternGenerator.h"
#include "PatternTree.h"
#include "TMath.h"
#include <iostream>
#include <fstream>

//FIXME: TEST
#include <ctime>

using namespace std;

ClassImp(TreeSearch::PatternGenerator)

namespace TreeSearch {

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
PatternGenerator::PatternGenerator() : fNlevels(0), fNplanes(0), fMaxSlope(0)
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
void PatternGenerator::DoTree( EOperation op )
{
  // Execute given operation on all unique Pattern nodes of the build tree.
  // Internal utility function.

  for( vector<Link*>::iterator it = fHashTable.begin();
       it != fHashTable.end(); ++it ) {
    Link* hashnode = *it;
    while( hashnode ) {
      Link* cur_node = hashnode;
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
      }
    }
  }
  if( op == kDelete )
    fHashTable.clear();
}

//_____________________________________________________________________________
void PatternGenerator::GetTreeStatistics( Statistics_t& stats ) const
{
  // Collect statistics on the build tree. Used by Print().

  memset( &stats, 0, sizeof(Statistics_t) );

  for( vector<Link*>::const_iterator it = fHashTable.begin();
       it != fHashTable.end(); ++it ) {
    // Each hashnode points to a unique pattern by construction of the table
    Link* hashnode = *it;
    UInt_t hash_length = 0;
    while( hashnode ) {
      // Count patterns
      stats.nPatterns++;
      // Count child nodes and length of child list
      Pattern* pat = hashnode->GetPattern();
      assert(pat);
      Link* ln = pat->fChild;
      UInt_t list_length = 0;
      while( ln ) {
	stats.nLinks++;
	list_length++;
	ln = ln->Next();
      }
      if( list_length > stats.MaxChildListLength )
	stats.MaxChildListLength = list_length;
      // Count collision list depth
      hash_length++;

      hashnode = hashnode->Next();

    } // while hashnode

    if( hash_length > stats.MaxHashDepth )
      stats.MaxHashDepth = hash_length;

  } // hashtable elements

  stats.nBytes = 
    stats.nPatterns * sizeof(Pattern)
    + stats.nPatterns * fNplanes * sizeof(UShort_t)
    + stats.nLinks * sizeof(Link);
  stats.nHashBytes = fHashTable.size() * sizeof(Link*)
    + stats.nPatterns * sizeof(Link);
}

//_____________________________________________________________________________
void PatternGenerator::Print( Option_t* opt, ostream& os ) const
{
  // Print information about the tree, depending on option

  // Dump all nodes
  if( *opt == 'D' ) {
    for( vector<Link*>::const_iterator it = fHashTable.begin();
	 it != fHashTable.end(); ++it ) {
      Link* hashnode = *it;
      while( hashnode ) {
	Pattern* pat = hashnode->GetPattern();
	pat->Print( true, os );
	hashnode = hashnode->Next();
      }
    }
    return;
  }

  // Basic info
  os << "tree: nlevels = " << fNlevels
     << ", nplanes = " << fNplanes
     << ", zpos = ";
  for( UInt_t i=0; i<fZ.size(); i++ ) {
    os << fZ[i];
    if( i+1 != fZ.size() )
      os << ",";
  }
  os << endl;

  Statistics_t stats;
  GetTreeStatistics( stats );
  os << "patterns = " << stats.nPatterns
     << ", links = "   << stats.nLinks
     << ", bytes = " << stats.nBytes
     << endl;
  os << "maxlinklen = " << stats.MaxChildListLength
     << ", maxhash = " << stats.MaxHashDepth
     << ", hashbytes = " << stats.nHashBytes
     << endl;
 
 //TODO: add more features

}

//_____________________________________________________________________________
void PatternGenerator::AddHash( Pattern* pat )
{
  // Add given pattern to the hash table

  assert(pat);
  UInt_t hashsize = fHashTable.size();
  if( hashsize == 0 ) {
    // Set the size of the hash table.
    // 2^(nlevels-1)*2^(nplanes-2) is the upper limit for the number of
    // patterns, so a size of 2^(nlevels-1) will give 2^(nplanes-2) collisions
    // per entry (i.e. 2, 4, 8), with which we can live. Anything better would
    // require a cleverer hash function.
    fHashTable.resize( 1<<(fNlevels-1), 0 );
    hashsize = fHashTable.size();
  }
  Int_t hash = pat->Hash()%hashsize;
  fHashTable[hash] = new Link( pat, fHashTable[hash], 0 );
}

//_____________________________________________________________________________
PatternTree* PatternGenerator::Generate( UInt_t maxdepth, Double_t width,
					 const vector<double>& zpos,
					 Double_t maxslope )
{
  // Generate a new pattern tree for the given parameters. Returns a pointer
  // to the generated tree, or zero if error

  // Clear out previous build, if any
  DoTree( kDelete );

  // Set parameters for the new build.
  //TODO: normalize zpos and maxslope
  fNlevels  = maxdepth+1;
  fZ        = zpos;
  fNplanes  = fZ.size();
  fMaxSlope = maxslope;
  // TODO: check for unreasonable input


  // FIXME: test
  clock_t start = clock();
  double cpu_secs;

  // Start with the trivial all-zero root node at depth 0. 
  Pattern* root = new Pattern( fNplanes );
  AddHash( root );

  // Generate the tree recursively
  MakeChildNodes( root, 1 );
  
  // FIXME: TEST TEST
  cpu_secs = ((double)(clock()-start))/CLOCKS_PER_SEC;

  Print();
  cout << "time = " << cpu_secs << " s" << endl;
  // Dump the entire database for inspection
  ofstream outfile("nodes.txt");
  if( outfile ) {
    Print("D",outfile);
    outfile.close();
  }

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
  // Assumes identical bin sizes and positions in each plane.

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
  Link* link = fHashTable[hash];
  while( link ) {
    Pattern* rhs = link->GetPattern();
    if( pat == *rhs )
      return rhs;
    link = link->Next();
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

  // Base case of the recursion: no child nodes beyond fNlevels-1
  if( depth >= fNlevels )
    return;

  // If not already done, generate the child patterns of the parent
  if( !parent->fChild ) {
    ChildIter it( *parent );
    while( it ) {
      Pattern& child = *it;

      // Pattern already exists?
      Pattern* node = Find( child );
      if( node ) {
	// If the pattern has only been tested at a higher depth, we need to
	// redo the slope test since the slope is larger now at lower depth
	if( depth >= node->fMinDepth || TestSlope(*node, depth)) {
	  // Only add a reference to the existing pattern
	  parent->AddChild( node, it.type() );
	}
      } else if( TestSlope(child, depth) && LineCheck(child) ) {
	// If the pattern is new, check it for consistency with maxslope and 
	// the straight line condition.
	node = new Pattern( child );
	AddHash( node );
	parent->AddChild( node, it.type() );
      }
      ++it;
    }
  }

  // Recursively generate child nodes down the tree
  Link* ln = parent->GetChild();
  while( ln ) {
    Pattern* node = ln->GetPattern();
    // We only need to go deeper if either this pattern does not have children
    // yet OR (important!), children were previously generated only from a 
    // deeper location in the tree and so this pattern's subtree needs to be
    // extended deeper down now.
    if( !node->fChild || node->fMinDepth > depth )
      MakeChildNodes( node, depth+1 );
    ln = ln->Next();
  }
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
