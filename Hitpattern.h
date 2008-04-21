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
#include <vector>
#include <set>
#include <utility>

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
  class Hit;

  class Hitpattern {

  public:
    explicit Hitpattern( const PatternTree& pt );  // preferred constructor
    Hitpattern( UInt_t nlevels, UInt_t nplanes, Double_t width );
    Hitpattern( const Hitpattern& orig );
    Hitpattern& operator=( const Hitpattern& rhs );
    virtual ~Hitpattern();

    std::pair<UInt_t,UInt_t> ContainsPattern( const NodeDescriptor& nd ) const;

    const std::vector<TreeSearch::Hit*>&  GetHits( UInt_t plane,
						   UInt_t bin ) const {
      // Get array of hits that set the given bin in the given plane
      return fHits[ MakeIdx(plane,bin) ];
    }
    UInt_t   GetNbins()   const { return 1U<<(fNlevels-1); }
    UInt_t   GetNlevels() const { return fNlevels; }
    UInt_t   GetNplanes() const { return fNplanes; }
    Double_t GetOffset()  const { return fOffset; }
    Double_t GetWidth()   const { return GetNbins()/fScale; }
    Double_t GetBinWidth() const { return fBinWidth; }  // meters per bin
    Double_t GetBinScale() const { return fScale; }     // bins per meter

    Bool_t   IsError()    const { return (fNplanes == 0); }

    void     SetPositionRange( Double_t start, Double_t end, UInt_t plane,
			       Hit* hit );
    void     SetPosition( Double_t pos, Double_t res, UInt_t plane,
			  Hit* hit )
    { SetPositionRange( pos-res, pos+res, plane, hit ); }
    Int_t    ScanHits( WirePlane* A, WirePlane* B );
    //    Bool_t   TestPosition( Double_t pos, UInt_t plane, UInt_t depth ) const;
    //    Bool_t   TestBin( UInt_t bin, UInt_t plane, UInt_t depth ) const;

    void     Clear( Option_t* opt="" );
    void     Print( Option_t* opt="" ) const;

    void     SetOffset( Double_t off ) { fOffset = off; }

#ifdef TESTCODE
    // Number of bins set at the highest resolution
    UInt_t   GetBinsSet() const;
    // Number of hits recorded
    UInt_t   GetNhits()   const { return (UInt_t)fHitList.size(); }
    // Maximum number of hits recorded per bin
    UInt_t   GetMaxhitBin() const { return fMaxhitBin; }
#endif

//TODO: add Draw() (=event display)

  protected:

    UInt_t   fNlevels;  // Number of levels in the pattern tree 
    UInt_t   fNplanes;  // Number of wire planes contained in the pattern
    Double_t fScale;    // 1/(bin resolution) = 2^(fNlevels-1)/width (1/m)
    Double_t fBinWidth; // 1/fScale (meters per bin)
    Double_t fOffset;   // Offset of zero hit position wrt zero det coord (m)
    Bits**   fPattern;  // [fNplanes] pattern at all fNlevels resolutions

    // Storage for saving pointers to the hits that set each active bin at
    // max level in each plane. Since each plane has the same number of
    // levels, each plane/bin combination can be represented with a 
    // single index (see MakeIdx below).
    std::vector< std::vector<Hit*> > fHits;
    // The plane/bin indices set in fHits. Used for fast clearing
    std::vector<UInt_t> fHitList;
    
    UInt_t MakeIdx( UInt_t plane, UInt_t bin ) const {
      // Return index into fHits corresponding to the given plane and bin
      assert( plane<fNplanes && bin<GetNbins() );
      UInt_t idx = (plane<<(fNlevels-1)) + bin;
      assert( idx < fHits.size());
      return idx;
    }

    void AddHit( UInt_t plane, UInt_t bin, Hit* hit );

#ifdef TESTCODE
    UInt_t  fMaxhitBin;  // Maximum depth of hit array per bin
#endif

  private:
    void Init( Double_t width );

    ClassDef(Hitpattern,0)  // Wire chamber hitpattern at multiple resolutions
  };

  //___________________________________________________________________________
//   inline
//   Bool_t TreeSearch::Hitpattern::TestBin( UInt_t bin, UInt_t plane, 
// 					  UInt_t depth ) const
//   {
//     // Test if point is set at the given depth and plane.

//     assert( depth < fNlevels && plane < fNplanes );
//     UInt_t offset = 1U<<depth;
//     assert( bin < offset );
//     return fPattern[plane]->TestBitNumber( bin + offset );
//   }

  //___________________________________________________________________________
//   inline
//   Bool_t Hitpattern::TestPosition( Double_t pos, UInt_t plane, 
// 				   UInt_t depth ) const
//   {
//     // Test if position 'pos' (in m) is marked in the hit pattern.
//     // The pattern will be tested at the given depth.

//     assert( depth < fNlevels && plane < fNplanes );
//     UInt_t bin = TMath::FloorNint( fScale*pos );
//     if( bin < 0 || bin >= GetNbins() )
//       return kFALSE;
//     return 
//       fPattern[plane]->TestBitNumber( (bin>>(fNlevels-depth-1))+(1U<<depth) );
//   }

  //___________________________________________________________________________
  inline std::pair<UInt_t,UInt_t> 
  Hitpattern::ContainsPattern( const NodeDescriptor& nd ) const
  {
    // Check if the hitpattern contains the pattern specified by the
    // NodeDescriptor. Returns the plane occupancy bitpattern and the 
    // count of planes where the pattern's bit was found set in the hitpattern.
    //
    // Used to compare with the patterns stored in the PatternTree class.
      
    assert( nd.depth < fNlevels );
    // The offset of the hitpattern bits at this depth
    UInt_t offs = 1U<<nd.depth;
    UInt_t matchval = 0, nmatch = 0;
    // The start bit number of the tree pattern we are comparing to
    UInt_t startpos = offs + nd.shift;
    Pattern* pat = nd.link->GetPattern();
    // Pointer to the last element of the pattern's bit array + 1
    UShort_t* bitnum = pat->GetBits() + fNplanes;
    // Check if the pattern's bits are set in the hitpattern, plane by plane
    if( nd.mirrored ) {
      assert( startpos < (offs<<1) );
      for( UInt_t i=fNplanes; i; ) {
	if( fPattern[--i]->TestBitNumber(startpos - *--bitnum) ) {
	  matchval |= (1U<<i);
	  ++nmatch;
	}
      }
    } else {
      assert( startpos + pat->GetWidth() < (offs<<1) );
      for( UInt_t i=fNplanes; i; ) {
	if( fPattern[--i]->TestBitNumber(startpos + *--bitnum) ) {
	  matchval |= (1U<<i);
	  ++nmatch;
	}
      }
    }
    return std::make_pair(matchval,nmatch);
  }


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch


#endif
