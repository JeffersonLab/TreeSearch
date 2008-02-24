///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternTree.h"
#include "Pattern.h"
//#include "TreeFile.h"
#include "TError.h"
#include <iostream>
#include <stdexcept>

using namespace std;

ClassImp(TreeSearch::PatternTree)

namespace TreeSearch {

//_____________________________________________________________________________
PatternTree::PatternTree( const TreeParam_t& param, UInt_t nPatterns,
			  UInt_t nLinks )
try
  : fParameters(param), fParamOK(false), fNpat(0), fNlnk(0), fNbit(0)
{
  // Constructor. 

  if( fParameters.Normalize() != 0 )
    return;

  fPatterns.resize( nPatterns );
  fLinks.resize( nLinks );
  fBits.resize( nPatterns * fParameters.zpos().size() );

  fParamOK = true;
}
catch ( bad_alloc ) {
  ::Error( "PatternTree::PatternTree", "Out of memory trying to create "
	   "%u patterns, %u links. Tree not created.", nPatterns, nLinks );
  throw;
}

//_____________________________________________________________________________
PatternTree::~PatternTree()
{
  // Destructor

}

//_____________________________________________________________________________
PatternTree* PatternTree::Read( const char* filename, const TreeParam_t& tp )
{
  // Read tree from file

  try {

    // TODO: implement
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
  // are irreversibly modified: zpos.front() and zpos.back() will be 0 and 1,
  // respectively, and maxslope is scaled to the normalized aspect ratio.
  // width is not changed

  static const char* const here = "TreeParam_t::Normalize";

  if( fNormalized )  return 0;

  if( fMaxdepth >= 16 ) {
    ::Error( here, "Illegal maxdepth = %u. Must be < 16", fMaxdepth );
    return -1;
  }
  if( fWidth < 1e-2 ) {
    ::Error( here, "Illegal detector width %lf. Must be >= 0.01.", 
	     fWidth );
    return -2;
  }
  if( fZpos.size() < 3 || fZpos.size() > 16 ) {
    ::Error( here, "Illegal number of planes = %u. Must be between 3-16.",
	     fZpos.size() );
    return -3;
  }
  if( fMaxslope < 0.0 )
    fMaxslope = -fMaxslope;

  // Check fZpos array for sorting and minimum spacing
  for( vector<Double_t>::size_type i = 1; i < fZpos.size(); ++i ) {
    if( fZpos[i] < fZpos[i-1] + 1e-3 ) {
      ::Error( here, "Array of z-positions not sorted or planes not "
	       "spaced by at least 0.001 at index = %u.", i );
      return -4;
    }
  }

  // Normalize the z-positions to the interval [0,1]
  Double_t zsize = fZpos.back() - fZpos.front();
  for( vector<Double_t>::size_type i = 1; i+1 < fZpos.size(); ++i ) {
    fZpos[i] -= fZpos.front();
    fZpos[i] /= zsize;
  }
  fZpos.back()  = 1.0;
  fZpos.front() = 0.0;

  // Scale maxslope to the new aspect ratio
  fMaxslope *= zsize / fWidth;

  // NB: the width is not changed because we need it later
  fNormalized = true;

  return 0;
}

//_____________________________________________________________________________
void PatternTree::Print( Option_t* opt, ostream& os )
{
  // Print information about the tree, depending on option

  TreeWalk walk( GetNlevels() );
  // Print ASCII pictures of ALL actual patterns in the tree ("P") or
  // dump n-tuples of all actual patterns, one per line ("L")
  if( *opt == 'P' or *opt == 'L' ) {
    PrintPattern print(os, (*opt == 'L'));
    walk( GetRoot(), print );
    return;
  }

  // Count all actual patterns
  if( *opt == 'C' ) {
    CountPattern count;
    walk( GetRoot(), count );
    os << "Total pattern count = " << count.GetCount() << endl;
    return;
  }

  // Basic info
//   os << "tree: nlevels = " << fNlevels
//      << ", nplanes = " << fNplanes
//      << ", zpos = ";
//   for( UInt_t i=0; i<fZ.size(); i++ ) {
//     os << fZ[i];
//     if( i+1 != fZ.size() )
//       os << ",";
//   }
//   os << endl;

//   os << "patterns = " << fStats.nPatterns
//      << ", links = "   << fStats.nLinks
//      << ", bytes = " << fStats.nBytes
//      << endl;
//   os << "maxlinklen = " << fStats.MaxChildListLength
//      << ", hashsize = " << fHashTable.size()
//      << ", hashbytes = " << fStats.nHashBytes
//      << endl;
//   os << "time = " << fStats.BuildTime << " s" << endl;
}

//_____________________________________________________________________________
Int_t PatternTree::Write( const char* filename )
{
  // Write tree to binary file

  // TODO: write header

  size_t index_size = sizeof(Int_t);
  vpsz_t npatt = fPatterns.size();
  if( npatt < (1U<<8) )
    index_size = 1;
  else if( npatt < (1U<<16) )
    index_size = 2;
  WritePattern write(filename,index_size);
  TreeWalk walk( GetNlevels() );
  Int_t ret = walk( GetRoot(), write );
  if( ret != kError )
    ret = 0;
  return ret;
}

//_____________________________________________________________________________
void 
PatternTree::CopyPattern::AddChild( Pattern* node, Pattern* child, Int_t type )
{
  // Add child to node's child pattern list

  Link* ln = node->GetChild();
  while( ln ) {
    if( !ln->GetPattern() ) {
      SetLinkPattern( ln, child );
      SetLinkType( ln, type );
      break;
    }
    assert( ln->GetPattern() != child or ln->Type() != type );
    ln = ln->Next();
    assert(ln);  // An empty slot must exist somewhere in the list
  }
}

//_____________________________________________________________________________
NodeVisitor::ETreeOp
PatternTree::CopyPattern::operator() ( const NodeDescriptor& nd )
try {
  // Add pattern to the PatternTree fTree

  map<Pattern*,Int_t>::iterator idx;
  Pattern* copied_parent = 0;
  if( nd.parent ) {
    idx = fMap.find(nd.parent);
    assert( idx != fMap.end() );
    copied_parent = &(fTree->fPatterns.at( idx->second ));
  }
  Pattern* node = nd.link->GetPattern();
  idx = fMap.find(node);
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
    // bits into the fTree->fBits vector.
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
    // Set the child pointer of the copied pattern to the first child node
    if( nchild > 0 ) {
      PatternTree::vlsz_t lpos = fTree->fNlnk;
      SetPatternChild( cur_pat, &(fTree->fLinks.at(lpos)) );
      // Link the child node list (not really needed, but we don't want the
      // fNext pointers to dangle)
      for( Int_t i = nchild-2; i >= 0; --i )
	SetLinkNext( &(fTree->fLinks.at(lpos+i)), 
		     &(fTree->fLinks.at(lpos+i+1)) );
    }
    // Update cursors
    fTree->fNpat++;
    fTree->fNlnk += nchild;
    fTree->fNbit += cur_pat->GetNbits();
    // Proceed with this pattern's child nodes
    return kRecurseUncond;
  } 
  else {
    // Existing pattern: add link to the parent pattern's child node list,
    // pointing to the referenced pattern
    Pattern* ref_node = &(fTree->fPatterns.at( idx->second ));
    AddChild( copied_parent, ref_node, nd.link->Type() );
    // Skip this pattern's child nodes since they are already in the tree
    return kSkipChildNodes;
  }
}
catch ( out_of_range ) {
  ::Error( "TreeSearch::CopyPattern", "Array index out of range at %u %u %u "
	   "(internal logic error). Tree not copied. Call expert.",
	   fTree->fNpat, fTree->fNlnk, fTree->fNbit );
  return kError;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
