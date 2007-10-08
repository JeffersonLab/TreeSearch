#ifndef ROOT_TreeSearch_Hitpattern
#define ROOT_TreeSearch_Hitpattern

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TBits.h"
#include "TMath.h"
#include "TreeWalk.h"
#include "Pattern.h"
#include <cstring>
#include <cassert>

namespace TreeSearch {

  // Extended version of TBits that supports setting and resetting bit ranges
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
    ClassDef(Bits,1)  // TBits container with range setting methods
  };

  class PatternTree;
  class WirePlane;

  class Hitpattern {

  public:
    //    Hitpattern( const PatternTree& pt );
    Hitpattern( UInt_t nlevels, UInt_t nplanes, Double_t width );
    Hitpattern( const Hitpattern& orig );
    Hitpattern& operator=( const Hitpattern& rhs );
    virtual ~Hitpattern();

    Bits*    GetRow( UInt_t i ) const { return i<fNplanes ? fPattern[i] : 0; }

    Double_t GetWidth()   const { return (1U<<(fNlevels-1))/fScale; }
    UInt_t   GetNlevels() const { return fNlevels; }
    UInt_t   GetNplanes() const { return fNplanes; }
    Double_t GetOffset()  const { return fOffset; }

    Bool_t   IsError()    const { return (fNplanes == 0); }

    void     SetPositionRange( Double_t start, Double_t end, UInt_t plane );
    void     SetPosition( Double_t pos, Double_t res, UInt_t plane )
    { SetPositionRange( pos-res, pos+res, plane ); }
    Int_t    ScanHits( WirePlane* A, WirePlane* B );
    Bool_t   TestPosition( Double_t pos, UInt_t plane, UInt_t depth ) const;
    Bool_t   TestBin( UInt_t bin, UInt_t plane, UInt_t depth ) const;
    UInt_t   ContainsPattern( const NodeDescriptor& nd ) const;

    void     Clear( Option_t* opt="" );
    void     Print( Option_t* opt="" ) const;

    void     SetOffset( Double_t off ) { fOffset = off; }

//FIXME: add Draw() (=event display)

  protected:

    UInt_t   fNlevels;  // Number of levels in the pattern tree 
    UInt_t   fNplanes;  // Number of wire planes contained in the pattern
    Double_t fScale;    // 1/(bin resolution) = 2^(fNlevels-1)/width (1/m)
    Double_t fOffset;   // Offset of zero hit position wrt zero det coord (m)
    Bits**   fPattern;  // [fNplanes] pattern at all fNlevels resolutions


    ClassDef(Hitpattern,0)  // Wire chamber hitpattern at multiple resolutions
  };

  //___________________________________________________________________________
  inline
  void Hitpattern::Clear( Option_t* opt )
  {
    // Clear the hitpattern

    for( UInt_t i=fNplanes; i; )
      fPattern[--i]->FastClear();
  }

  //___________________________________________________________________________
  inline
  Bool_t TreeSearch::Hitpattern::TestBin( UInt_t bin, UInt_t plane, 
					  UInt_t depth ) const
  {
    // Test if point is set at the given depth and plane.

    assert( depth < fNlevels && plane < fNplanes );
    UInt_t offset = 1U<<depth;
    assert( bin < offset );
    return fPattern[plane]->TestBitNumber( bin + offset );
  }

  //___________________________________________________________________________
  inline
  Bool_t Hitpattern::TestPosition( Double_t pos, UInt_t plane, 
				   UInt_t depth ) const
  {
    // Test if position 'pos' (in m) is marked in the hit pattern.
    // The pattern will be tested at the given depth.

    assert( depth < fNlevels && plane < fNplanes );
    Int_t bin = TMath::FloorNint( fScale*pos );
    if( bin < 0 || bin >= 1<<(fNlevels-1) )
      return kFALSE;
    return 
      fPattern[plane]->TestBitNumber( (bin>>(fNlevels-depth-1))+(1U<<depth) );
  }

  //___________________________________________________________________________
  inline
  UInt_t Hitpattern::ContainsPattern( const NodeDescriptor& nd ) const
  {
    // Check if the hitpattern contains the pattern specified by the
    // NodeDescriptor. Returns the count of planes where the pattern's bit
    // was found set in the hitpattern.
    // Used to compare with the patterns stored in the PatternTree class.
      
    assert( nd.depth < fNlevels );
    // The offset of the hitpattern bits at this depth
    UInt_t offs = 1U<<nd.depth;
    UInt_t n_found = 0;
    // The start bit number of the tree pattern we are comparing to
    UInt_t startpos = offs + nd.shift;
    Pattern* pat = nd.link->GetPattern();
    assert(pat);
    // Pointer to the last element of the pattern's bit array + 1
    UShort_t* theBin = pat->GetBits() + fNplanes;
    // Check if the pattern's bits are set in the hitpattern, plane by plane
    if( nd.mirrored ) {
      assert( startpos < (offs<<1) );
      for( UInt_t i=fNplanes; i; ) {
	if( fPattern[--i]->TestBitNumber(startpos - *--theBin) )
	  ++n_found;
      }
    } else {
      assert( startpos + pat->GetWidth() < (offs<<1) );
      for( UInt_t i=fNplanes; i; ) {
	if( fPattern[--i]->TestBitNumber(startpos + *--theBin) )
	  ++n_found;
      }
    }
    return n_found;
  }


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch


#endif
