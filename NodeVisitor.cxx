///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::NodeVisitor classes                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "NodeVisitor.h"
#include "TError.h"
#include <iostream>

using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
WritePattern::WritePattern( const char* filename, size_t index_size )
  : os(0), fIdxSiz(index_size)
{
  static const char* here = "PatternGenerator::WritePattern";

  if( filename && *filename ) {
    os = new ofstream(filename, ios::out|ios::binary|ios::trunc);
    if( !os || !(*os) ) {
      ::Error( here, "Error opening treefile %s", filename );
      if( os )
	delete os;
      os = 0;
    }
  } else {
    ::Error( here, "Invalid file name" );
  }
  if( (fIdxSiz & (fIdxSiz-1)) != 0 ) {
    ::Warning( here, "Invalid index_size = %d. Must be a power of 2", fIdxSiz);
    fIdxSiz = sizeof(Int_t);
  }    
}

//_____________________________________________________________________________
template< typename T>
static inline
void swapped_binary_write( ostream& os, const T& data, size_t n = 1, 
			   size_t start = 0 )
{
  // Write "n" elements of "data" to "os" in binary big-endian (MSB) format.
  // "start" indicates an optional _byte_ skip count from the beginning of
  // the big-endian format data - designed to write only the non-trivial
  // part of scalar data for which the upper bytes are known to be always zero.
  size_t size = sizeof(T);
  const char* bytes = reinterpret_cast<const char*>( &data );
#ifdef R__BYTESWAP
  size_t k = size-1;
  for( size_t i = start; i < n * size; ++i )
    os.put( bytes[(i&~k)+(k-i&k)] );
#else
  os.write( bytes+start, n*size-start );
#endif
}

//_____________________________________________________________________________
TreeWalk::ETreeOp WritePattern::operator() ( const NodeDescriptor& nd )
{
  // Write the single pattern referenced by "nd" to the binary output stream.
  // In essence, this implements the serialization method for cyclical graphs
  // described in
  // http://www.parashift.com/c++-faq-lite/serialization.html#faq-36.11

  if( !os )
    return TreeWalk::kError;
  Pattern* node = nd.link->GetPattern();
  map<Pattern*,Int_t>::iterator idx = fMap.find(node);
  if( idx == fMap.end() ) {
    map<Pattern*,Int_t>::size_type n = fMap.size();
    fMap[node] = n;
    // Header for new pattern: link type + 128 (=128-130)
    os->put( nd.link->Type() | 0x80 );
    if( os->fail() ) return TreeWalk::kError;
    // Pattern data. NB: fBits[0] is always 0, so we can skip it
    swapped_binary_write( *os, node->GetBits()[1], node->GetNbits()-1 );
    // Child node count
    UShort_t nchild = node->GetNchildren();
    swapped_binary_write( *os, nchild );
    if( os->fail() ) return TreeWalk::kError;
    // Write child nodes regardless of depth
    return TreeWalk::kRecurseUncond;
  } else {
    // Reference pattern header: the plain link type (=0-2)
    os->put( nd.link->Type() );
    // Reference index
    swapped_binary_write( *os, idx->second, 1, sizeof(Int_t)-fIdxSiz );
    if( os->fail() ) return TreeWalk::kError;
    // Don't write child nodes
    return TreeWalk::kSkipChildNodes;
  }
}

//_____________________________________________________________________________
TreeWalk::ETreeOp PrintPattern::operator() ( const NodeDescriptor& nd )
{
  // Print actual (shifted & mirrored) pattern described by "nd".

  ++fCount;
  if( fDump )
    os << setw(2) << nd.depth;

  Pattern* node = nd.link->GetPattern();
  for( UInt_t i = 0; i < node->GetNbits(); i++ ) {
    UInt_t v = node->GetBits()[i];
    if( nd.mirrored )  // Mirrored pattern
      v = node->GetWidth() - v;
    v += nd.shift;

    // In dump mode, write out one pattern n-tuple per line
    if( fDump )
      os << " " << setw(5) << v;

    // Otherwise draw a pretty ASCII picture of the pattern
    else {
      UInt_t op = (nd.mirrored ? 2 : 0) + nd.link->Shift();
      os << static_cast<UInt_t>(nd.depth) << "-" << op;
      for( UInt_t k = 0; k < nd.depth; ++k )
	os << " ";
      os << " |";
      for( UInt_t k = 0; k < v; ++k )
	os << ".";
      os << "O";
      for( UInt_t k = (1<<nd.depth)-1; k > v; --k )
	os << ".";
      os << "|" << endl;
    }
  }
  os << endl;

  // Process child nodes through maxdepth
  return TreeWalk::kRecurse;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
