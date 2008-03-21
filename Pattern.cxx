//*-- Author :    Ole Hansen, Jefferson Lab   19-Sep-2007

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
  : fChild(0), fNbits(size), fDelBits(true), fDelChld(false)
{
  // Constructor. If size > 0, create array for the bits on the heap.

  assert( fNbits <= 16 );
  if( fNbits == 0 ) {
    fDelBits = false;
    fBits = 0;
  } else {
    fBits = new UShort_t[fNbits];
    memset( fBits, 0, fNbits*sizeof(UShort_t) );
  }
}

//_____________________________________________________________________________
Pattern::Pattern( const Pattern& orig )
  : fChild(0), fNbits(orig.fNbits), fDelBits(orig.fDelBits), fDelChld(false)
{
  // Copy constructor. Copy has same bits, no children

  if( fDelBits ) {
    fBits = new UShort_t[fNbits];
    memcpy( fBits, orig.fBits, fNbits*sizeof(UShort_t) );
  } else {
    fBits = orig.fBits;
  }
}

//_____________________________________________________________________________
const Pattern& Pattern::operator=( const Pattern& rhs )
{
  // Assignment. Copies only this bits, not the pointers to the children.

  if( this != &rhs ) {
    fChild = 0;
    fNbits = rhs.fNbits;
    fDelBits = rhs.fDelBits;
    fDelChld = false;
    if( fDelBits ) {
      delete fBits;
      if( fNbits ) {
	fBits = new UShort_t[fNbits];
	memcpy( fBits, rhs.fBits, fNbits*sizeof(UShort_t) );
      } else {
	fBits = 0;
      }
    } else {
      fBits = rhs.fBits;
    }
  }
  return *this;
}

//_____________________________________________________________________________
Pattern::~Pattern()
{
  // Destroys a Pattern. Deletes all child links (but not the Patterns
  // they point to)

  if( fDelChld ) {
    while( fChild ) {
      Link* child = fChild;
      fChild = fChild->Next();
      delete child;
    }
  }
  if( fDelBits )
    delete [] fBits;
}

//_____________________________________________________________________________
Int_t Pattern::GetNchildren() const
{
  // Returns number of child nodes of the pattern

  Int_t n = 0;
  Link* ln = fChild;
  while( ln ) {
    ++n;
    ln = ln->Next();
  }
  return n;
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

//_____________________________________________________________________________
void Pattern::SetBitloc( UShort_t* bitloc )
{
  // Tell the pattern to store its bits at the externally-managed location
  // bitloc. The current fBits[fNbits] will be copied to bitloc. 
  // If bitloc is NULL, revert to internally-managed bits.

  UShort_t* oldbits = fBits;
  if( bitloc ) {
    fBits = bitloc;
    if( oldbits ) {
      memcpy( fBits, oldbits, fNbits*sizeof(UShort_t) );
      if( fDelBits )
	delete [] oldbits;
    }
    fDelBits = false;
  } else if( !fDelBits ) {
    fBits = new UShort_t[fNbits];
    if( oldbits )
      memcpy( fBits, oldbits, fNbits*sizeof(UShort_t) );
    fDelBits = true;
  }
}


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
