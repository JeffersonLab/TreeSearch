#ifndef ROOT_TreeSearch_Hitpattern
#define ROOT_TreeSearch_Hitpattern

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TBits.h"
#include <cstring>

namespace TreeSearch {

  // Modified version of TBits that clears the bits without deleting the
  // reserved memory
  class Bits : public TBits {
  public:
    Bits( UInt_t nbits = 8 ) : TBits(nbits) {}
    Bits( const TBits& orig ) : TBits(orig) {}
    Bits& operator=( const TBits& rhs ) {
      TBits::operator=(rhs); return *this;
    }
    virtual ~Bits() {}
    void ResetBitRange( UInt_t lo, UInt_t hi );
    void SetBitRange( UInt_t lo, UInt_t hi );
    void FastClear() { memset(fAllBits,0,fNbytes); }
    ClassDef(Bits,1)  // Bit container with additional methods over TBits
  };

  class PatternTree;
  class WirePlane;

  class Hitpattern {

  public:
    //    Hitpattern( const PatternTree& pt );
    Hitpattern( UInt_t depth, UInt_t nplanes, Double_t width );
    virtual ~Hitpattern();

    Bits*  GetRow( UInt_t i ) const { return i<fNplanes ? fPattern[i] : 0; }
    Bits** GetPattern() const { return fPattern; }
    Double_t GetWidth() const { return fWidth; }

    Int_t  SetPoints( WirePlane* hitsA, WirePlane* hitsB,
		      UInt_t plane );

    void   Clear( Option_t* opt="" );
    void   Print( Option_t* opt="" ) const;

//FIXME: add Draw() (=event display)

  protected:

    UInt_t   fDepth;    // Depth of the pattern tree
    UInt_t   fNplanes;  // Number of wire planes contained in the pattern
    Double_t fWidth;    // Physical width of region encompassing hitpattern (m)
    Bits**   fPattern;  // [fNPlanes] pattern at all fDepth resolutions

    ClassDef(Hitpattern,0)  // Wire chamber hitpattern at multiple resolutions
  };

}  // end namespace TreeSearch

//_____________________________________________________________________________
inline
void TreeSearch::Bits::SetBitRange( UInt_t lo, UInt_t hi )
{
  // Set range of bits from lo to hi to value.

  if( hi<lo ) return;
  SetBitNumber( hi ); // expand array if necessary
  if( lo==hi ) return;
  UChar_t mask  = ~((1U<<(lo&7))-1);
  UChar_t mask2 = (1U<<(hi&7))-1;
  lo >>= 3;
  hi >>= 3;
  if( lo < hi ) {
    fAllBits[hi] |= mask2;
    memset( fAllBits+lo+1, 0xFF, hi-lo-1 );
  } else {
    mask &= mask2;
  }
  fAllBits[lo] |= mask;
}

//_____________________________________________________________________________
inline
void TreeSearch::Bits::ResetBitRange( UInt_t lo, UInt_t hi )
{
  // Set range of bits from lo to hi to value.

  if( hi<lo ) return;
  UChar_t mask = ~((1U<<(lo&7))-1);
  lo >>= 3;
  if( lo >= fNbytes ) return;
  UChar_t mask2 = (1U<<((hi&7)+1))-1;
  hi >>= 3;
  if( hi >= fNbytes ) {
    hi = fNbytes-1;
    mask2 = 0xFF;
  }
  if( lo < hi ) {
    fAllBits[hi] &= (0xFF ^ mask2);
    memset( fAllBits+lo+1, 0, hi-lo-1 );
  } else {
    mask &= mask2;
  }
  fAllBits[lo] &= (0xFF ^ mask);
}

///////////////////////////////////////////////////////////////////////////////

#endif
