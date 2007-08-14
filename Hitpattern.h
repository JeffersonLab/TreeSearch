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
    void SetBitRange( UInt_t lo, UInt_t hi, Bool_t value = kTRUE );
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
void TreeSearch::Bits::SetBitRange( UInt_t lo, UInt_t hi, Bool_t value )
{
  // Set range of bits from lo to hi to value.

  if( hi<lo ) return;
  SetBitNumber( hi, value ); // expand array if necessary
  if( lo==hi ) return;
  UInt_t  loc_lo = lo/8;
  UChar_t bit_lo = lo%8;
  UInt_t  loc_hi = hi/8;
  UChar_t bit_hi = hi%8;
  UChar_t mask_lo = ~((1<<bit_lo)-1);
  if( loc_lo < loc_hi ) {
    UChar_t mask = (1U<<bit_hi)-1;
    if( value ) {
      fAllBits[loc_hi] |= mask;
      mask = 0xFF;
    } else {
      fAllBits[loc_hi] &= (0xFF ^ mask);
      mask = 0;
    }
    memset( fAllBits+loc_lo+1, mask, loc_hi-loc_lo-1 );
  } else {
    mask_lo &= (1U<<bit_hi)-1;
  }
  if( value )
    fAllBits[loc_lo] |= mask_lo;
  else
    fAllBits[loc_lo] &= (0xFF ^ mask_lo);
}

///////////////////////////////////////////////////////////////////////////////

#endif
