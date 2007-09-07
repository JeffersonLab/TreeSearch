///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hitpattern.h"
#include "Hit.h"
#include "WirePlane.h"
#include "TError.h"
#include "TMath.h"

#include <iostream>

using std::cout;
using std::endl;

ClassImp(TreeSearch::Hitpattern)

namespace TreeSearch {

//_____________________________________________________________________________
Hitpattern::Hitpattern( UInt_t depth, UInt_t nplanes, Double_t width )
  : fDepth(depth), fNplanes(0), fWidth(width), fScale(0.0), 
    fOffset(0.5*width), fPattern(NULL)
{
  // Constructor

  static const char* const here = "TreeSearch::Hitpattern";
  
  if( nplanes == 0 || nplanes > 100 || fDepth == 0 || fDepth > 20 ) {
    ::Error( here, "Illegal number of planes or depth: %d %d.\n"
	     "Both > 0, nplanes <= 100, depth <= 20.", nplanes, depth );
  } else if( fWidth < 1e-2 ) { // Negative or very small width?
    ::Error( here, "Illegal detector width %lf. Must be >= +1cm.", fWidth );
  } else {
    fNplanes = nplanes;
    fPattern = new Bits*[fNplanes];
    UInt_t nbins = 1U<<fDepth;
    for( UInt_t i=0; i<fNplanes; i++ )
      fPattern[i] = new Bits( nbins );
    fScale = static_cast<Double_t>(nbins>>1) / fWidth;
  }
}

//_____________________________________________________________________________
Hitpattern::Hitpattern( const Hitpattern& orig ) 
  : fDepth(orig.fDepth), fNplanes(orig.fNplanes), fWidth(orig.fWidth),
    fScale(orig.fScale), fOffset(orig.fOffset), fPattern(NULL)
{
  // Copy ctor

  if( fNplanes > 0 ) {
    fPattern = new Bits*[fNplanes];
    for( UInt_t i=fNplanes; i--; )
      fPattern[i] = new Bits(*orig.fPattern[i]);
  }
}

//_____________________________________________________________________________
Hitpattern& Hitpattern::operator=( const Hitpattern& rhs )
{
  // Assignment

  if( this != &rhs ) {
    fDepth   = rhs.fDepth;
    fNplanes = rhs.fNplanes;
    fWidth   = rhs.fWidth;
    fScale   = rhs.fScale;
    fOffset  = rhs.fOffset;
    delete fPattern; fPattern = NULL;
    if( fNplanes > 0 ) {
      fPattern = new Bits*[fNplanes];
      for( UInt_t i=fNplanes; i--; )
	fPattern[i] = new Bits(*rhs.fPattern[i]);
    }
  }
  return *this;
}

//_____________________________________________________________________________
Hitpattern::~Hitpattern()
{
  // Destructor

  for( UInt_t i=fNplanes; i; )
    delete fPattern[--i];
  delete [] fPattern;
}

//_____________________________________________________________________________
void Hitpattern::uSetPositionRange( Double_t start, Double_t end,
				    UInt_t plane )
{
  // Set pattern bins corresponding to the exact physical positions
  // between start and end (in m) in the given plane. 
  // Positions may range from 0.0 to fWidth.

  Int_t hi = TMath::FloorNint( fScale*end );
  if( hi < 0 ) return;
  Int_t lo = TMath::FloorNint( fScale*start );
  // At the deepest tree level, there are 2^(fDepth-1) bins.
  Int_t nbins = 1<<(fDepth-1);
  if( lo >= nbins ) return;
  if( lo < 0 )
    lo = 0;
  if( hi >= nbins )
    hi = nbins-1;

  // Loop through the tree levels, starting at the highest resolution.
  // In practice, we usually have hi-lo <= 1 even at the highest resolution.
  while (true) {
    fPattern[plane]->SetBitRange( lo+nbins, hi+nbins );
    nbins >>= 1;
    if( nbins == 0 ) break;
    lo >>= 1;
    hi >>= 1;
  }
}

//_____________________________________________________________________________
Int_t Hitpattern::ScanHits( WirePlane* A, WirePlane* B )
{
  // Set the points at all depths of the hit pattern that correspond to
  // the hits in plane A. The plane number is extracted from A.
  // The second, optional plane B represents an optional partner
  // plane of A, i.e. a nearby plane with (usually) staggered wires. 
  // If B is present, test for pairs of hits in A and B.
  // Naturally, the two planes and their wires must be parallel.
  //
  // Returns number of hits processed

  static const Double_t inv2sqrt2 = 1.0/(2.0*TMath::Sqrt(2.0));

  if( !A ) return 0;
  UInt_t plane = A->GetPlaneNum();
  if( plane >= fNplanes ) {
    ::Warning( "TreeSearch::Hitpattern::ScanHits", 
	       "Plane number %d out of range (max = %d). Call expert.", 
	       plane, fNplanes-1 );
    return 0;
  }
  Double_t dz = B ? B->GetZ() - A->GetZ() : 0.0;
  Double_t maxdist = A->GetMaxSlope() * dz;
  Double_t maxdist2 = 0.5*maxdist;
  bool do_single_hits =
    ( A->GetDetector()->TestBit(MWDC::kPairsOnly) == kFALSE || B == NULL );

  Int_t nhits = 0;

  HitPairIter it( A->GetHits(), B ? B->GetHits() : NULL, maxdist );
  while( it ) {
    nhits++;
    const Hit* hitA = static_cast<const Hit*>((*it).first);
    const Hit* hitB = static_cast<const Hit*>((*it).second);
    bool found = false;
    if( hitA && hitB ) {
      // Combined resolution of the two hits (assuming their individual
      // resolutions are similar):
      Double_t res = inv2sqrt2*(hitA->GetResolution()+hitB->GetResolution());
      for( int i=4; i--; ) {
	Double_t posA = (i&2 ? hitA->GetPosL() : hitA->GetPosR()) + fOffset;
	Double_t posB = (i&1 ? hitB->GetPosL() : hitB->GetPosR()) + fOffset;
	if( TMath::Abs( posA-posB ) < maxdist ) {
	  found = true;
	  uSetPosition( 0.5*(posA+posB), res, plane );
	}
      }
    }
    if( !found && do_single_hits ) {
      // Here, we have either an unpaired hit or a pair whose positions do
      // not match within maxdist (probably rare unless maxdist is very small).
      // Either way, we set the hits individually by projecting their left
      // and right positions onto the reference plane (either plane A, if no B,
      // or the midplane between A and B).
      for( int i=2; i--; ) {
	const Hit* hit = i ? hitA : hitB;
	if( hit ) {
	  // If projecting onto the midplane (maxdist != 0), the effective
	  // resolution is larger by 0.5*maxdist so that we set all the bits
	  // the hit position can project to. If maxdist is big enough,
	  // this will probably cause some ghosts, which is unavoidable to
	  // maintain high tracking efficiency. The full hit resolution
	  // is recovered when the hit positions are fit within each road.
	  //
	  // The bigger issue with unpaired hits is that the LR-ambiguity
	  // is not resolved, so two entries have to be made into the pattern.
	  Double_t res = hit->GetResolution() + maxdist2;
	  uSetPosition( hit->GetPosL() + fOffset, res, plane );
	  uSetPosition( hit->GetPosR() + fOffset, res, plane );
	}
      }
    }
    ++it;
  }
  return nhits;
}


//_____________________________________________________________________________
Bool_t Hitpattern::TestPosition( Double_t pos, UInt_t plane, 
				 UInt_t depth ) const
{
  // Test if position 'pos' (in m) is marked in the hit pattern.
  // The pattern will be tested at the given depth (default: deepest level).

  if( plane >= fNplanes )
    return kFALSE;
  Int_t bin = TMath::FloorNint( fScale*pos );
  if( bin < 0 || bin >= 1<<(fDepth-1) )
    return kFALSE;
  if( depth >= fDepth )
    depth = fDepth-1;
  return fPattern[plane]->TestBitNumber( (bin>>(fDepth-depth-1))+(1<<depth) );
}

//_____________________________________________________________________________
void Hitpattern::Print( Option_t* opt ) const
{
  // Print basic info about hitpattern.

  //TODO
}


//_____________________________________________________________________________
void Bits::ResetBitRange( UInt_t lo, UInt_t hi )
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
void Bits::SetBitRange( UInt_t lo, UInt_t hi )
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

} // end namespace TreeSearch

