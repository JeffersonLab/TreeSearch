///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternTree.h"
#include "Pattern.h"
#include "TError.h"
#include <iostream>
#include <stdexcept>

using namespace std;

ClassImp(TreeSearch::PatternTree)

namespace TreeSearch {

//_____________________________________________________________________________
PatternTree::PatternTree( const TreeParam_t& parameters, UInt_t nPatterns,
			  UInt_t nLinks )
try
  : fParameters(parameters), fParamOK(false), fNpat(0), fNlnk(0), fNbit(0), 
    fID(0)
{
  // Constructor. 

  if( fParameters.Normalize() != 0 )
    return;

  fPatterns.resize( nPatterns );
  fLinks.resize( nLinks );
  fBits.resize( nPatterns * fParameters.zpos.size() );

  fParamOK = true;
}
catch ( bad_alloc ) {
  ::Error( "PatternTree::PatternTree", "Out of memory trying to create "
	   "%u patterns, %u links", nPatterns, nLinks );
}

//_____________________________________________________________________________
PatternTree::~PatternTree()
{
  // Destructor

}

//_____________________________________________________________________________
Int_t PatternTree::Read( const char* filename )
{
  // Read tree from file

  try {

  }
  catch ( bad_alloc ) {
    ::Error( "PatternTree::Read", "Out of memory trying to read "
	     "%u patterns, %u links", 0,0 );
    // TODO: clean up
  }

  return 0;
}


//_____________________________________________________________________________
Int_t TreeParam_t::Normalize()
{
  // Test, then normalize parameter values to width = 1 and z_max = 1.
  // Returns 0 on success, != 0 on error. On success only, the input data
  // are irreversibly modified. zpos.front(), zpos.back(), and width will 
  // have trivial values (0, 1, 1, respectively) that algorithms can rely on.

  static const char* const here = "TreeParam_t::Normalize";

  if( maxdepth >= 16 ) {
    ::Error( here, "Illegal maxdepth = %u. Must be < 16", maxdepth );
    return -1;
  }
  if( width < 1e-2 ) {
    ::Error( here, "Illegal detector width %lf. Must be >= 0.01.", 
	     width );
    return -2;
  }
  if( zpos.size() < 3 || zpos.size() > 16 ) {
    ::Error( here, "Illegal number of planes = %u. Must be between 3-16.",
	     zpos.size() );
    return -3;
  }
  if( maxslope < 0.0 )
    maxslope = -maxslope;

  // Check zpos array for sorting and minimum spacing
  for( vector<Double_t>::size_type i = 1; i < zpos.size(); ++i ) {
    if( zpos[i] < zpos[i-1] + 1e-3 ) {
      ::Error( here, "Array of z-positions not sorted or planes not "
	       "spaced by at least 0.001 at index = %u.", i );
      return -4;
    }
  }

  // Normalize the z-positions to the interval [0,1]
  Double_t zsize = zpos.back() - zpos.front();
  for( vector<Double_t>::size_type i = 1; i+1 < zpos.size(); ++i ) {
    zpos[i] -= zpos.front();
    zpos[i] /= zsize;
  }
  zpos.back()  = 1.0;
  zpos.front() = 0.0;

  // Scale maxslope to the new aspect ratio
  maxslope *= zsize / width;

  width = 1.0;

  return 0;
}

//_____________________________________________________________________________
void PatternTree::CopyPattern::AddChild( Pattern* node, Pattern* child, 
					 Int_t type )
{
  // Add child to node's child pattern list

  Link* ln = node->GetChild();
  while( ln ) {
    if( !ln->fPattern ) {
      ln->fPattern = child;
      ln->fOp = type;
      break;
    }
    assert( ln->GetPattern() != child or ln->Type() != type );
    ln = ln->Next();
    assert(ln);  // An empty slot must exist somewhere in the list
  }
}

//_____________________________________________________________________________
Int_t PatternTree::CopyPattern::operator() ( const NodeDescriptor& nd )
try {
  // Add pattern to the PatternTree fTree

  Pattern* copied_parent = 0;
  if( nd.parent ) {
    map<Pattern*,Int_t>::iterator ip = fMap.find(nd.parent);
    assert( ip != fMap.end() );
    copied_parent = &(fTree->fPatterns.at( ip->second ));
  }
  Pattern* node = nd.link->GetPattern();
  map<Pattern*,Int_t>::iterator idx = fMap.find(node);
  if( idx == fMap.end() ) {
    // New pattern: add pattern to pattern and bits arrays, and add its 
    // child links to links array
    map< Pattern*,Int_t>::size_type np = fMap.size();
    fMap[node] = np;
    assert( fTree->fNpat == np );
    Pattern* cur_pat = &(fTree->fPatterns.at(fTree->fNpat));
    // Copy the pattern to the array element via its assignment operator,
    // which sets the fChild pointer to zero and copies the bits
    *cur_pat = *node;
    // Per 2003 C++ Standard TC, std::vector is guaranteed to store its 
    // elements in contiguous memory, so this is safe:
    UShort_t* bitloc = &(fTree->fBits.at(fTree->fNbit));
    // Tell the pattern to store its bits at address bitloc. This copies the
    // bits into the fBits vector.
    cur_pat->SetBitloc( bitloc );
    // Link this pattern to its parent
    if( copied_parent ) {
      AddChild( copied_parent, cur_pat, nd.link->Type() );
    } else {
      // If there is no parent, this had better be the root node
      assert( fTree->fNlnk == 0 );
      fTree->fLinks.at(fTree->fNlnk++) = Link(cur_pat, 0, nd.link->Type());
    }
    // Create the child node links, but with empty pattern pointers
    Int_t nchild = node->GetNchildren();
    PatternTree::vlsz_t lpos = fTree->fNlnk;
    // Set the child pointer of the copied pattern to the first child node
    if( nchild > 0 ) {
      cur_pat->fChild = &(fTree->fLinks.at(lpos));
      cur_pat->fDelChld = false;
      // Link the child node list (not really needed, but we don't want the
      // fNext pointers to dangle)
      for( Int_t i = nchild-2; i >= 0; --i )
	fTree->fLinks.at(lpos+i).fNext = &(fTree->fLinks.at(lpos+i+1));
    }
    // Update cursors
    fTree->fNpat++;
    fTree->fNlnk += nchild;
    fTree->fNbit += cur_pat->GetNbits();
    // Proceed with this pattern's child nodes
    return 0;
  } 
  else {
    // Existing pattern: add link to the parent pattern's child node list,
    // pointing to the referenced pattern
    Pattern* ref_node = &(fTree->fPatterns.at( idx->second ));
    AddChild( copied_parent, ref_node, nd.link->Type() );
    // Skip this pattern's child nodes since they are already in the tree
    return 1;
  }
}
catch ( out_of_range ) {
  ::Error( "TreeSearch::CopyPattern", "Array index out of range at %u %u %u "
	   "(internal logic error). Tree not copied. Call expert.",
	   fTree->fNpat, fTree->fNlnk, fTree->fNbit );
  return -1;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
