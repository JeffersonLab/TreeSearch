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

  class Link;

  class Pattern {
    friend class PatternGenerator;
  private:
    UShort_t*  fBits;        // [fNbits] Bit pattern array
    Link*      fChild;       // Linked list of child patterns
    UChar_t    fNbits;       // Bit count
    UChar_t    fMinDepth;    // Minimum valid depth for this pattern
    UChar_t    fMaxDepth;    // Maximum valid depth (used during building)
    Int_t      fRefIndex;    // Reference index for serializing the tree

    Link*      AddChild( Pattern* child, Int_t type );
    Pattern*   FindChild( const Pattern* pat, Int_t type );
    void       UsedAtDepth( UInt_t depth );

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
    Link*      GetChild()    const  { return fChild; }
    UShort_t*  GetBits()            { return fBits; }
    UInt_t     GetMinDepth() const  { return fMinDepth; }
    UInt_t     GetWidth()    const  { return fBits[fNbits-1]-fBits[0]; }
    UInt_t     GetNbits()    const  { return fNbits; }
    Int_t      Hash()        const  { return GetWidth(); }
    UShort_t&  operator[](UInt_t i) { 
      assert(i<fNbits);
      return fBits[i];
    }
    UShort_t   operator[](UInt_t i) const { 
      assert(i<fNbits);
      return fBits[i];
    }
    void print( bool print_links = true, std::ostream& os = std::cout,
		bool end_line = true ) const;
  }; // end class Pattern

  class Link {
  private:
    Pattern*   fPattern;    // Bit pattern treenode
    Link*      fNext;       // Next list element
    Int_t      fOp;         // Operation to be applied to pattern
  public:
    Link( Pattern* pat, Link* next, Int_t op )
      : fPattern(pat), fNext(next), fOp(op) {}
    Pattern*   GetPattern() const { return fPattern; }
    Link*      Next() const { return fNext; }
    Int_t      Type() const { return fOp; }
    void print( std::ostream& os = std::cout, bool end_line = true ) const;
  }; // end class Link

  //___________________________________________________________________________
  inline
  void Pattern::UsedAtDepth( UInt_t depth )
  {
    // Mark pattern as having been used at given depth

    if( depth < fMinDepth ) fMinDepth = depth;
    if( depth > fMaxDepth ) fMaxDepth = depth;
  }


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
