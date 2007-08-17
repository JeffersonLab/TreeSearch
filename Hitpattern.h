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
    Hitpattern( const Hitpattern& orig );
    Hitpattern& operator=( const Hitpattern& rhs );
    virtual ~Hitpattern();

    Bits*  GetRow( UInt_t i ) const { return i<fNplanes ? fPattern[i] : 0; }

    Double_t GetWidth() const { return fWidth; }
    UInt_t   GetDepth() const { return fDepth; }
    UInt_t   GetNplanes() const { return fNplanes; }

    void     SetPositionRange( Double_t start, Double_t end, UInt_t plane )
    { if(plane<fNplanes && start<=end) uSetPositionRange(start,end,plane); }
    Int_t    ScanHits( WirePlane* hitsA, WirePlane* hitsB );
    Bool_t   TestPosition( Double_t pos, UInt_t plane, UInt_t depth=32 ) const;
    Bool_t   TestBin( UInt_t bin, UInt_t plane, UInt_t depth=32 ) const;

    void     Clear( Option_t* opt="" );
    void     Print( Option_t* opt="" ) const;

//FIXME: add Draw() (=event display)

  protected:

    UInt_t   fDepth;    // Depth of the pattern tree (e.g. 10: max 2^9 bits)
    UInt_t   fNplanes;  // Number of wire planes contained in the pattern
    Double_t fWidth;    // Physical width of region encompassing hitpattern (m)
    Double_t fScale;    // 1/(bin resolution) = 2^(fDepth-1)/fWidth
    Double_t fOffset;   // Offset of zero hit position wrt zero det coord (m)
    Bits**   fPattern;  // [fNPlanes] pattern at all fDepth resolutions

    // More efficient "unchecked" versions of these routines for internal use
    void   uSetPositionRange( Double_t start, Double_t end, UInt_t plane );
    void   uSetPosition( Double_t pos, Double_t res, UInt_t plane )
    { uSetPositionRange( pos-res, pos+res, plane ); }

    ClassDef(Hitpattern,0)  // Wire chamber hitpattern at multiple resolutions
  };

}  // end namespace TreeSearch

//_____________________________________________________________________________
inline
Bool_t TreeSearch::Hitpattern::TestBin( UInt_t bin, UInt_t plane, 
					UInt_t depth ) const
{
  // Test if point is set at the given depth of the pattern tree in the
  // given plane.

  if( depth >= fDepth )
    depth = fDepth-1;
  UInt_t offset = 1U<<depth;
  if ( plane >= fNplanes || bin >= offset )
    return kFALSE;
  return fPattern[plane]->TestBitNumber( bin + offset );
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

///////////////////////////////////////////////////////////////////////////////

#endif
