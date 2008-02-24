///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternGenerator                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternGenerator.h"
#include "PatternTree.h"
#include "TMath.h"
#include "TString.h"
#include "TError.h"
#include <iostream>
#include <stdexcept>
#include <ctime>
#include <cassert>
#if ROOT_VERSION_CODE < ROOT_VERSION(5,8,0)
#include <cstdlib>   // for atof()
#endif

using namespace std;

ClassImp(TreeSearch::PatternGenerator)

// Private definitions, used internally
namespace {

using namespace TreeSearch;

// Iterator over child patterns of a given parent patten
class ChildIter {
private:
  const Pattern fParent;  // copy of parent pattern
  Pattern   fChild;       // current child pattern
  Int_t     fCount;       // trial iterations left to do
  Int_t     fType;        // current pattern type (normal/shifted/mirrored)
public:
  ChildIter( const Pattern& parent ) 
    : fParent(parent), fChild(parent), fType(0) { reset(); }
  ChildIter&      operator++();
  const ChildIter operator++(int) { 
    ChildIter clone(*this);
    ++(*this); 
    return clone;
  }
  Pattern& operator*()            { return fChild; }
           operator bool()  const { return (fCount >= 0); }
  Int_t    type()           const { return fType; }
  void     reset() { 
    fCount = 1<<fParent.GetNbits();
    ++(*this);
  }
};

//_____________________________________________________________________________
inline
ChildIter& ChildIter::operator++()
{ 
  // Iterator over all suitable child patterns of a given parent pattern
  // that occur when the bin resolution is doubled.
  // Child pattern bits are either 2*bit or 2*bit+1 of the parent bits,
  // yielding 2^nbits (=2^nplanes) different combinations.
  // The bits of suitable patterns must monotonically increase.
  // bit[0] is always zero (otherwise the pattern could be shifted).
  // Note that type() is an important part of the result.
  // type & 1 indicates a pattern shifted by one to the right
  // type & 2 indicates a mirrored pattern 
  // To recover the actual pattern, mirror first, then shift, as appropriate.
  // With the self-referential tree structure used here, mirrored patterns
  // only ever occur as children of the root of the tree. Simultaneously
  // mirrored and shifted patterns never occur.
  // Hence, type = 0, 1, and, very rarely, 2.

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

///////////////////////////////////////////////////////////////////////////////

} //end namespace

namespace TreeSearch {

//_____________________________________________________________________________
PatternGenerator::PatternGenerator()
  : fNlevels(0), fNplanes(0), fMaxSlope(0)
{
  // Constructor

}

//_____________________________________________________________________________
PatternGenerator::~PatternGenerator()
{
  // Destructor

  DeleteTree();

}

//_____________________________________________________________________________
void PatternGenerator::DeleteTree()
{
  // Delete the build tree along with its hash table.
  // Internal utility function.

  for( vector<HashNode>::iterator it = fHashTable.begin();
       it != fHashTable.end(); ++it ) {
    delete (*it).GetPattern();
  }
  fHashTable.clear();
  ClearStatistics();
}

//_____________________________________________________________________________
void PatternGenerator::CalcStatistics()
{
  // Collect statistics on the build tree. This is best done separately here
  // because some things (averages, memory requirements) can only be 
  // calculated once the tree is complete.

  ClearStatistics();

  for( vector<HashNode>::const_iterator it = fHashTable.begin();
       it != fHashTable.end(); ++it ) {
    // Each hashnode points to a unique pattern by construction of the table
    Pattern* pat = (*it).GetPattern();
    if( pat ) {
      // Count patterns
      fStats.nPatterns++;
      // Count child nodes and length of child list
      Link* ln = pat->fChild;
      UInt_t list_length = 0;
      while( ln ) {
	fStats.nLinks++;
	list_length++;
	ln = ln->Next();
      }
      if( list_length > fStats.MaxChildListLength )
	fStats.MaxChildListLength = list_length;
    } // if pattern
  } // hashtable elements

  // Count the root node's link, too
  fStats.nLinks++;

  fStats.nBytes = 
    fStats.nPatterns * sizeof(Pattern)
    + fStats.nPatterns * fNplanes * sizeof(UShort_t)
    + fStats.nLinks * sizeof(Link);
  fStats.nHashBytes = fHashTable.size() * sizeof(HashNode);
}

//_____________________________________________________________________________
void PatternGenerator::ClearStatistics()
{
  memset( &fStats, 0, sizeof(Statistics_t) );
}

//_____________________________________________________________________________
void PatternGenerator::Print( Option_t* opt, ostream& os ) const
{
  // Print information about the tree, depending on option

  // Dump all stored patterns, using Pattern::print()
  if( *opt == 'D' ) {
    for( vector<HashNode>::const_iterator it = fHashTable.begin();
	 it != fHashTable.end(); ++it ) {
      Pattern* pat = (*it).GetPattern();
      if( pat )
	pat->Print( true, os );
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

  os << "patterns = " << fStats.nPatterns
     << ", links = "   << fStats.nLinks
     << ", bytes = " << fStats.nBytes
     << endl;
  os << "maxlinklen = " << fStats.MaxChildListLength
     << ", hashsize = " << fHashTable.size()
     << ", hashbytes = " << fStats.nHashBytes
     << endl;
  os << "time = " << fStats.BuildTime << " s" << endl;
}

//_____________________________________________________________________________
PatternTree* PatternGenerator::Generate( TreeParam_t parameters )
{
  // Generate a new pattern tree for the given parameters. Returns a pointer
  // to the generated tree, or zero if error

  static const char* const here = "PatternGenerator::Generate";

  // Set parameters for the new build.
  if( parameters.Normalize() != 0 )
    return 0;

  DeleteTree();

  fNlevels  = parameters.maxdepth()+1;
  fZ        = parameters.zpos();
  fNplanes  = fZ.size();
  fMaxSlope = parameters.maxslope();

  // Benchmark the build
  clock_t cpu_ticks = clock();

  // Start with the trivial all-zero root node at depth 0. 
  Pattern* root = new Pattern( fNplanes );
  HashNode* hroot = AddHash( root );

  // Generate the tree recursively
  MakeChildNodes( hroot, 1 );
  
  // Calculate tree statistics (number of patterns, links etc.)
  cpu_ticks = clock() - cpu_ticks;
  CalcStatistics();
  fStats.BuildTime = static_cast<Double_t>(cpu_ticks)/CLOCKS_PER_SEC;

  //FIXME: TEST
  // Print tree statistics
  Print();

  // Create a PatternTree object from the build tree
  PatternTree* tree = 0;
  try {
    tree = new PatternTree( parameters, fStats.nPatterns, fStats.nLinks );
  }
  catch( bad_alloc ) {
    delete tree;
    tree = 0;
  }
  if( tree ) {
    // Copy all base patterns of the build tree to the PatternTree
    Link root_link(root,0,0);
    TreeWalk walk(fNlevels);
    PatternTree::CopyPattern copy(tree);
    cout << "Filling tree..." << flush;
    NodeVisitor::ETreeOp ret = walk( &root_link, copy );
    cout << "done" << endl;
    if( ret == kError ) {
      ::Error( here, "Unable to create PatternTree structure" );
      delete tree;
      return 0;
    }
    //FIXME: TEST
    // print the copied tree to test file
//     ofstream outf( "pt.txt", ios::out|ios::trunc );
//     PrintPattern print(outf);
//     fTreeWalk( tree->GetRoot(), print );
//     cout << "Pattern count = " << print.GetCount() << endl;
  }

  return tree;
}

//_____________________________________________________________________________
PatternTree* PatternGenerator::Generate( UInt_t maxdepth, 
					 Double_t detector_width, 
					 const char* zpos,
					 Double_t maxslope )
{
  // Convenience function for use at the command line and in scripts

  vector<Double_t> zvec;
  TString zlist( zpos );
  TString tok;
  Ssiz_t pos = 0;
  while( zlist.Tokenize(tok, pos, ",") ) {
    if( tok.IsNull() )
      continue;
#if ROOT_VERSION_CODE < ROOT_VERSION(5,8,0)
    Double_t val = atof(tok.Data());
#else
    Double_t val = tok.Atof();
#endif
    zvec.push_back( val );
  }
  
  return Generate( TreeParam_t( maxdepth, detector_width, maxslope, zvec ));
}

//_____________________________________________________________________________
UInt_t PatternGenerator::Hash( const Pattern& pat ) const
{
  // Calculate unique a hash value for the given pattern. The maximum value is
  // 2^(fNlevels-1 + fNplanes-2) - 1. With the hash computed here, collisions 
  // never occur for patterns that pass LineTest(). 

  UInt_t hash = pat[fNplanes-1] * (1U << fNplanes-2);
  const Double_t x = pat[fNplanes-1];
  const Double_t z = fZ[fNplanes-1];

  // In each intermediary plane, patterns can have at most two valid bit
  // positions, yielding a unique signature for each valid pattern
  // (assuming equal bin width for all planes).
  for( Int_t i = fNplanes-2; i > 0; --i ) {
    Double_t d = pat[i]*z - x*fZ[i];
    if( d > 0 )
      hash += (1U << (i-1));
    // The case d==0, where a bin's edge falls _exactly_ on the boundary of
    // the allowed region, is ambiguous and may be affected by floating point
    // rounding behavior. The test here relies on LineTest() to reject at least
    // patterns with an active bin to the left of a d==0 bin.
  }
  return hash;
}

//_____________________________________________________________________________
PatternGenerator::HashNode* PatternGenerator::AddHash( Pattern* pat )
{
  // Add given pattern to the hash table

  // This implements a perfect hash table (no collisions), the fastest way to
  // do the pattern lookup during build. 
  assert(pat);
  if( fHashTable.empty() )
    fHashTable.resize( (1U << (fNlevels-1)) * (1U << (fNplanes-2)) );

  HashNode& h = fHashTable[ Hash(*pat) ];
  assert(h.fPattern == 0); // Overwriting exisiting entries should never happen
  h.fPattern = pat;

  return &h;
}

//_____________________________________________________________________________
PatternGenerator::HashNode* PatternGenerator::Find( const Pattern& pat )
{
  // Search for the given pattern in the current database

  HashNode& h = fHashTable[ Hash(pat) ];
  if( h.fPattern ) {
    if( pat == *h.fPattern )
      return &h;
    // A hash collision for valid patterns should never happen
    assert( LineTest(pat) == false );
  }
  return 0;
}

//_____________________________________________________________________________
inline
bool PatternGenerator::SlopeTest( const Pattern& pat, UInt_t depth ) const
{
  UInt_t width = pat.GetWidth();
  return ( width < 2 or
	   TMath::Abs((double)(width-1) / (double)(1<<depth)) <= fMaxSlope );
}

//_____________________________________________________________________________
bool PatternGenerator::LineTest( const Pattern& pat ) const
{
  // Check if the gievn bit pattern is consistent with a straight line.
  // The intersection plane positions are given by fZ[]. Assumes fZ[0]=0.
  // The other z values must increase strictly monotonically, fZ{i] > fZ[i-1].
  // In the parent class, the z-values are normalized so that 
  // fZ[fNplanes-1] = 1.0, but that is not required.
  // Assumes a normalized pattern, for which pat[0] is always zero.
  // Assumes identical bin sizes and bin boundaries in each plane.
  // If any of these assumption are relaxed, the algorithm requires additional
  // tests and might need to save some data in temporary variables.
  //
  // The algorithm works by computing left and right boundary lines within
  // which a bin must lie. The boundaries are the connections between the
  // top and bottom reference points, (SL, SR) and (BL, BR), respectively.
  // Starting from the plane one below the top (fZ[fNplanes-2]), bins in
  // successively lower planes are tested to lie in the allowed region. If they
  // do, the boundaries are narrowed if necessary to ensure that subsequent
  // bins are consistent with all bins above. If the test succeeds in all
  // planes, the pattern is consistent with a line.
  //
  // Algorithm originally developed by Brandon Belew, SULI summer intern,
  // Jefferson Lab, 2007.

  // FIXME: for certain z-values, the following can be _very_ sensitive 
  // to the floating point rounding behavior! (Can this be fixed?)

  struct Point {
    Double_t x, z;
  };

  Point SL = { pat[fNplanes-1], fZ[fNplanes-1] };
  Point SR = { SL.x + 1, SL.z };
  // Since fZ[0] = 0, we only need the x-coordinate of BL and BR
  Double_t xBL = 0.0;
  Double_t xBR = 1.0;
  // mL and mR are the current slopes of the left and right boundaries, resp.
  Double_t mL  = SL.x/SL.z;
  Double_t mR  = mL;
  
  for( Int_t j = fNplanes-2; j > 0; --j ) {
    // xL and xR are the edges of the active bin in the current plane
    // The bin width is assumed to be unity in all planes
    Double_t xL = pat[j];
    Double_t xR = xL+1;
    // z is the current plane's z position
    Double_t z  = fZ[j];
    // jL and jR define the left and right boundary of the allowed x-region
    // at the current plane
    Double_t jL = z*mL + xBL;
    Double_t jR = z*mR + xBR;

    // If the current bin is outside of the boundaries, we're done. 
    // The test fails.
    if( xL >= jR or xR <= jL )
      return false;

    // If this was the last plane before the bottom, all tests succeeded,
    // and we are done with the loop.
    if( j == 1 )
      break;

    // Recalculate the right and left boundaries for the next hit, given that
    // it passes through this bin. This is only necessary, of course, if the
    // present bin "cuts" into the allowed region (either its left or right
    // edge are between the boundaries).
    assert( not( xR < jR and xL > jL )); //both edges must never be inside
    if( xR < jR ) {
      assert((SL.z-z)>1e-3);
      mR = (SL.x - xR) / (SL.z - z);
      Double_t new_xBR = xR - z*mR;
      if( new_xBR < xBR )
	xBR = new_xBR;
      else {
	assert(z>1e-3);
	mR = (xR - xBR) / z;
      }
      SR.x = xR;
      SR.z = z;
    }
    // By construction (though not totally obvious), we always have 
    // jR <= jL+1. Also, we are guaranteed xR = xL+1 (with bin width = 1).
    // Thus, if the above test, xR < jR, is true, then also xL < jL,
    // and the following can always be skipped. Hence the "else" here:
    // (This is important because the previous block reassigns SR!)
    else if( xL > jL ) {
      assert((SR.z-z)>1e-3);
      mL = (SR.x - xL) / (SR.z - z);
      Double_t new_xBL = xL - z*mL;
      if( new_xBL > xBL )
	xBL = new_xBL;
      else {
	assert(z>1e-3);
	mL = (xL - xBL) / z;
      }
      SL.x = xL;
      SL.z = z;
    }
  } // planes

  return true;
}

//_____________________________________________________________________________
void PatternGenerator::MakeChildNodes( HashNode* pnode, UInt_t depth )
{
  // Generate child nodes for the given parent pattern

  // Requesting child nodes for the parent at this depth implies that the 
  // parent is being used at the level above
  if( depth > 0 )
    pnode->UsedAtDepth( depth-1 );

  // Base case of the recursion: no child nodes beyond fNlevels-1
  if( depth >= fNlevels )
    return;

  // If not already done, generate the child patterns of this parent
  Pattern* parent = pnode->GetPattern();
  assert(parent);
  if( !parent->fChild ) {
    ChildIter it( *parent );
    while( it ) {
      Pattern& child = *it;

      // Pattern already exists?
      HashNode* node = Find( child );
      if( node ) {
	Pattern* pat = node->GetPattern();
	assert(pat);
	// If the pattern has only been tested at a higher depth, we need to
	// redo the slope test since the slope is larger now at lower depth
	if( depth >= node->fMinDepth or SlopeTest(*pat, depth)) {
	  // Only add a reference to the existing pattern
	  parent->AddChild( pat, it.type() );
	}
      } else if( SlopeTest(child, depth) and LineTest(child) ) {
	// If the pattern is new, check it for consistency with maxslope and 
	// the straight line condition.
	Pattern* pat = new Pattern( child );
	AddHash( pat );
	parent->AddChild( pat, it.type() );
      }
      ++it;
    }
  }

  // Recursively generate child nodes down the tree
  Link* ln = parent->GetChild();
  while( ln ) {
    Pattern* pat = ln->GetPattern();
    // This Find() is necessary because we don't want to store the state 
    // parameter fMinDepth with the pattern itself - if we did so, we'd
    // get bigger patterns but minimally faster code here
    HashNode* node = Find( *pat );
    assert(node);
    // We only need to go deeper if either this pattern does not have children
    // yet OR (important!) children were previously generated from a deeper
    // location in the tree and so this pattern's subtree needs to be extended
    // deeper down now.
    if( !pat->fChild or node->fMinDepth > depth )
      MakeChildNodes( node, depth+1 );
    ln = ln->Next();
  }
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
