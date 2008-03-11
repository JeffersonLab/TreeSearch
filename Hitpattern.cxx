///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hitpattern.h"
#include "Hit.h"
#include "WirePlane.h"
#include "PatternTree.h"
#include "TError.h"
#include "TMath.h"
#include "MWDC.h"

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <set>

using namespace std;

typedef std::vector<TreeSearch::Hit*>::size_type vsiz_t;

ClassImp(TreeSearch::Hitpattern)

namespace TreeSearch {

//_____________________________________________________________________________
Hitpattern::Hitpattern( const PatternTree& pt )
  : fNlevels(pt.GetNlevels()), fNplanes(pt.GetNplanes()), fScale(0), 
    fOffset(0.5*pt.GetWidth()), fPattern(0)
#ifdef TESTCODE
  , fMaxhitBin(0)
#endif
{
  // Construct Hitpattern using paramaters of pattern tree

  Init( pt.GetWidth() );
}

//_____________________________________________________________________________
Hitpattern::Hitpattern( UInt_t nlevels, UInt_t nplanes, Double_t width )
  : fNlevels(nlevels), fNplanes(0), fScale(0), fOffset(0.5*width), fPattern(0)
#ifdef TESTCODE
  , fMaxhitBin(0)
#endif
{
  // Constructor

  static const char* const here = "TreeSearch::Hitpattern";
  
  if( nplanes == 0 || nplanes > 16 || fNlevels == 0 || fNlevels > 16 ) {
    ::Error( here, "Illegal number of planes or tree levels: %d %d.\n"
	     "Both must be 0 < n <= 16.", nplanes, nlevels );
  } else if( width < 1e-2 ) { // Negative or very small width?
    ::Error( here, "Illegal detector width %lf. Must be >= +1cm.", width );
  } else {
    fNplanes = nplanes;
    Init( width );
  }
}

//_____________________________________________________________________________
void Hitpattern::Init( Double_t width )
{
  // Allocate memory for new Hitpattern object.
  // Internal utility function called by constructors.

  assert( width >= 1e-2 );
  fScale = GetNbins() / width;
  fBinWidth = 1.0/fScale;

  try {
    fPattern = new Bits*[fNplanes];
    UInt_t nbins2 = 2*GetNbins();  // 2*number of bins at deepest level
    for( UInt_t i=0; i<fNplanes; i++ )
      fPattern[i] = new Bits( nbins2 );
    fHits.resize( fNplanes*GetNbins() );
  }
  catch ( std::bad_alloc ) {
    ::Error( "Hitpattern::Hitpattern", "Out of memory trying to construct "
	     "new Hitpattern object. Call expert." );
    throw;
  }
}

//_____________________________________________________________________________
Hitpattern::Hitpattern( const Hitpattern& orig ) 
try
  : fNlevels(orig.fNlevels), fNplanes(orig.fNplanes),
    fScale(orig.fScale), fBinWidth(orig.fBinWidth), fOffset(orig.fOffset),
    fPattern(0), fHits(orig.fHits), fHitList(orig.fHitList)
#ifdef TESTCODE
  , fMaxhitBin(orig.fMaxhitBin)
#endif
{
  // Copy ctor

  if( fNplanes > 0 ) {
    fPattern = new Bits*[fNplanes];
    for( UInt_t i=fNplanes; i--; )
      fPattern[i] = new Bits(*orig.fPattern[i]);
  }
  assert( fHits.size() == fNplanes*GetNbins() );
}
catch ( std::bad_alloc ) {
  ::Error( "Hitpattern::Hitpattern", "Out of memory trying to copy Hitpattern "
	   "object. Call expert." );
  throw;
}

//_____________________________________________________________________________
Hitpattern& Hitpattern::operator=( const Hitpattern& rhs )
{
  // Assignment

  if( this != &rhs ) {
    fNlevels = rhs.fNlevels;
    fNplanes = rhs.fNplanes;
    fScale   = rhs.fScale;
    fBinWidth= rhs.fBinWidth;
    fOffset  = rhs.fOffset;
    delete fPattern; fPattern = 0;
    if( fNplanes > 0 ) {
      fPattern = new Bits*[fNplanes];
      for( UInt_t i=fNplanes; i--; )
	fPattern[i] = new Bits(*rhs.fPattern[i]);
    }
    fHits = rhs.fHits;
    assert( fHits.size() == fNplanes*GetNbins() );
    fHitList = rhs.fHitList;
#ifdef TESTCODE
    fMaxhitBin = rhs.fMaxhitBin;
#endif
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
inline
void Hitpattern::AddHit( UInt_t plane, UInt_t bin, Hit* hit ) {
  // Add hit for given bin in plane to the hit arrays
  assert(hit);
  UInt_t idx = MakeIdx( plane, bin );
  fHits[idx].push_back( hit );
  fHitList.push_back( idx );
#ifdef TESTCODE
  if( fMaxhitBin < (UInt_t)fHits[idx].size() )
    fMaxhitBin = fHits[idx].size();
#endif
}

//_____________________________________________________________________________
void Hitpattern::Clear( Option_t* )
{
  // Clear the hitpattern

  for( UInt_t i=fNplanes; i; )
    fPattern[--i]->FastClear();

  // For speed, clear only arrays that are actually filled
  for( vector<UInt_t>::size_type i = fHitList.size(); i; ) {
    UInt_t idx = fHitList[--i];
    assert( idx < fHits.size());
    fHits[idx].clear();
  }
  fHitList.clear();

#ifdef TESTCODE
  fMaxhitBin = 0;
#endif
}

//_____________________________________________________________________________
#ifdef TESTCODE
UInt_t Hitpattern::GetBinsSet() const
{
  // Return number of bins set at the highest resolution
  UInt_t n = 0, nbins = GetNbins();
  for( UInt_t i=fNplanes; i; ) {
    n += fPattern[--i]->CountBits( nbins );
  }
  return n;
}
#endif

//_____________________________________________________________________________
void Hitpattern::SetPositionRange( Double_t start, Double_t end,
				   UInt_t plane, Hit* hitA, Hit* hitB )
{
  // Set pattern bins corresponding to the exact physical positions
  // between start and end (in m) in the given plane. 
  // Positions may range from 0.0 to width.
  // Associate these bins with hitA and (optional) hitB.

  assert( plane<fNplanes && start<=end );
  Int_t hi = TMath::FloorNint( fScale*end );
  if( hi < 0 ) return;
  Int_t lo = TMath::FloorNint( fScale*start );
  // At the deepest tree level, there are 2^(fNlevels-1) bins.
  Int_t nbins = GetNbins();
  if( lo >= nbins ) return;
  if( lo < 0 )
    lo = 0;
  if( hi >= nbins )
    hi = nbins-1;

  // Save the hit pointer(s) in the hit array so that we can efficiently
  // retrieve later the hit(s) that caused the bits to be set.
  // NB: This is quite expensive
  for( Int_t i = lo; i <= hi; ++i ) {
    AddHit( plane, i, hitA );
    if( hitB )
      AddHit( plane, i, hitB );
  }

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
  // NB: a "plane" in the hitpattern is a "layer" in the projection
  // (= single plane or plane pair, as appropriate) 
  UInt_t plane = A->GetLayerNum();
  assert( plane < fNplanes );
  Double_t dz = B ? B->GetZ() - A->GetZ() : 0.0;
  Double_t maxdist = A->GetMaxSlope() * dz;
  Double_t maxdist2 = 0.5*maxdist;
  bool do_single_hits =
    ( A->GetMWDC()->TestBit(MWDC::kPairsOnly) == kFALSE || B == 0 );

  Int_t nhits = 0;

  HitPairIter it( A->GetHits(), B ? B->GetHits() : 0, maxdist );
  while( it ) {
    nhits++;
    Hit* hitA = static_cast<Hit*>((*it).first);
    Hit* hitB = static_cast<Hit*>((*it).second);
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
	  SetPosition( 0.5*(posA+posB), res, plane, hitA, hitB );
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
	Hit* hit = i ? hitA : hitB;
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
	  SetPosition( hit->GetPosL() + fOffset, res, plane, hit );
	  SetPosition( hit->GetPosR() + fOffset, res, plane, hit );
	}
      }
    }
    ++it;
  }
  return nhits;
}


//_____________________________________________________________________________
void Hitpattern::Print( Option_t* ) const
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
  // Set range of bits from lo to hi (inclusive, i.e. [lo,hi])

  if( hi<lo ) return;
  SetBitNumber( hi ); // expands array if necessary
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

