///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Pattern                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Pattern.h"

#include <iostream>
using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
Pattern::Pattern( UInt_t size )
  : fChild(0), fNbits(size)
{
  // Constructor. Puts bit array on heap. Minimum size is 1.

  assert( fNbits <= 16 );
  if( fNbits == 0 )
    ++fNbits;
  fBits = new UShort_t[fNbits];
  memset( fBits, 0, fNbits*sizeof(UShort_t) );
}

//_____________________________________________________________________________
// Pattern::Pattern( UShort_t* loc )
//: fBits(loc), fChild(0), fNbits(0), fMinDepth(0), fRefIndex(-1)
// {
//   // Constructor. Puts bit array at loc, a location managed by the caller.
// }

//_____________________________________________________________________________
Pattern::Pattern( const Pattern& orig )
  : fChild(0), fNbits(orig.fNbits)
{
  // Copy constructor. Copy has same bits, no children

  fBits = new UShort_t[fNbits];
  memcpy( fBits, orig.fBits, fNbits*sizeof(UShort_t) );
}

//_____________________________________________________________________________
const Pattern& Pattern::operator=( const Pattern& rhs )
{
  // Assignment. Copies only this bits, not the pointers to the children.

  if( this != &rhs ) {
    fChild = 0;
    fNbits = rhs.fNbits;
    delete fBits;
    fBits = new UShort_t[fNbits];
    memcpy( fBits, rhs.fBits, fNbits*sizeof(UShort_t) );
  }
  return *this;
}

//_____________________________________________________________________________
Pattern::~Pattern()
{
  // Destroys a Pattern. Deletes all child links (but not the Patterns
  // they point to)

  while( fChild ) {
    Link* child = fChild;
    fChild = fChild->Next();
    delete child;
  }
  delete [] fBits;
}

//_____________________________________________________________________________
void Pattern::Print( bool print_links, ostream& os, bool end_line ) const 
{
  // Print this pattern and, if requested, its child patterns

  for( UInt_t i=0; i<fNbits; i++ ) {
    os << fBits[i] << " ";
  }
  //  os << "=" << fMinDepth;
  if( print_links ) {
    os << ":  ";
    Link* ln = fChild;
    while( ln ) {
      // Print each link
      os << "(" << ln->Type() << ") ";
      ln->GetPattern()->Print( false, os, false );
      ln = ln->Next();
      if( ln )
	os << ", ";
    }
  }
  if( end_line )
    os << endl;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
