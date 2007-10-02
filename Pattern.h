#ifndef ROOT_TreeSearch_Pattern
#define ROOT_TreeSearch_Pattern

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Pattern                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <cassert>
#include <cstring>
#include <iostream>

namespace TreeSearch {

  class Pattern;

  class Link {
  private:
    Pattern*   fPattern;    // Bit pattern treenode
    Link*      fNext;       // Next list element
    Int_t      fOp;         // Operation to be applied to pattern (bits 0-1)
  public:
    Link( Pattern* pat, Link* next, Int_t op )
      : fPattern(pat), fNext(next), fOp(op) { assert( op >= 0 && op < 4 ); }
    Pattern*   GetPattern() const { return fPattern; }
    Link*      Next()       const { return fNext; }
    Int_t      Type()       const { return fOp; }
    Bool_t     Mirrored()   const { return ((fOp & 2) != 0); }
    UInt_t     Shift()      const { return (fOp & 1); }
  }; // end class Link

  class Pattern {
    friend class PatternGenerator;
  private:
    UShort_t*  fBits;        // [fNbits] Bit pattern array
    Link*      fChild;       // Linked list of child patterns
    UShort_t   fNbits;       // Bit count (size of fBits array, <=16)
    UShort_t   fMinDepth;    // Minimum valid depth for this pattern (<=16)
    Int_t      fRefIndex;    // Reference index for serializing the tree

    Link*      AddChild( Pattern* child, Int_t type ) {
      assert(child);
      return (fChild = new Link( child, fChild, type ));
    }
    void       UsedAtDepth( UInt_t depth ) {
      if( depth < fMinDepth ) fMinDepth = depth;
    }

  public:
    Pattern( UInt_t size );
    Pattern( const Pattern& orig );
    Pattern& operator=( const Pattern& rhs );
    //    Pattern( UShort_t* bitloc );
    ~Pattern();
    bool operator==( const Pattern& rhs ) const {
      assert( fNbits == rhs.fNbits );
      return (0 == memcmp( fBits, rhs.fBits, fNbits*sizeof(UShort_t)));
    }
    bool operator!=( const Pattern& rhs ) const { 
      return !(*this == rhs );
    }
    UShort_t&  operator[](UInt_t i) { 
      assert(i<fNbits);
      return fBits[i];
    }
    UShort_t   operator[](UInt_t i) const { 
      assert(i<fNbits);
      return fBits[i];
    }
    Link*      GetChild()    const  { return fChild; }
    UShort_t*  GetBits()            { return fBits; }
    //    UInt_t     GetMinDepth() const  { return fMinDepth; }
    UInt_t     GetWidth()    const  { return fBits[fNbits-1]-fBits[0]; }
    UInt_t     GetNbits()    const  { return fNbits; }
    Int_t      GetRefIndex() const  { return fRefIndex; }
    Int_t      Hash()        const  { return GetWidth(); }
    void       SetRefIndex( Int_t i ) { fRefIndex = i; }
    void       Print( bool print_links = true, std::ostream& os = std::cout,
		      bool end_line = true ) const;
  }; // end class Pattern

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
