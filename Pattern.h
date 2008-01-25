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
    friend class NodeVisitor;
  private:
    Pattern*   fPattern;    // Bit pattern treenode
    Link*      fNext;       // Next list element
    Int_t      fOp;         // Operation to be applied to pattern (bits 0-1)
  public:
    Link( Pattern* pat = 0, Link* next = 0, Int_t op = 0 )
      : fPattern(pat), fNext(next), fOp(op) { assert( op >= 0 && op < 4 ); }
    Pattern*   GetPattern() const { return fPattern; }
    Link*      Next()       const { return fNext; }
    Int_t      Type()       const { return fOp; }
    Bool_t     Mirrored()   const { return ((fOp & 2) != 0); }
    UInt_t     Shift()      const { return (fOp & 1); }
  }; // end class Link

  class Pattern {
    friend class PatternGenerator;
    friend class NodeVisitor;
  private:
    UShort_t*  fBits;        // [fNbits] Bit pattern array
    Link*      fChild;       // Linked list of child patterns
    UShort_t   fNbits;       // Bit count (size of fBits array, <=16)
    Bool_t     fDelBits;     // If true, delete fBits in d'tor
    Bool_t     fDelChld;     // If true, delete fChild list in d'tor

    Link*      AddChild( Pattern* child, Int_t type ) {
      assert(child);
      fDelChld = true;
      return (fChild = new Link( child, fChild, type ));
    }
  public:
    explicit Pattern( UInt_t size = 0 );
    Pattern( const Pattern& orig );
    const Pattern& operator=( const Pattern& rhs );
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

    Link*      GetChild()     const { return fChild; }
    UShort_t*  GetBits()            { return fBits; }
    UInt_t     GetWidth()     const { return fBits[fNbits-1]-fBits[0]; }
    UInt_t     GetNbits()     const { return fNbits; }
    Int_t      GetNchildren() const;
    void       Print( bool print_links = true, std::ostream& os = std::cout,
		      bool end_line = true ) const;
    void       SetBitloc( UShort_t* bitloc );

  }; // end class Pattern

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
