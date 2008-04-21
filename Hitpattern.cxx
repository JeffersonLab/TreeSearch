//*-- Author :    Ole Hansen, Jefferson Lab   13-Aug-2007

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
  // Even though this may record the same idx multiple times, it performs
  // better than anything more fancy (like e.g. a set)
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
  for( vector<UInt_t>::iterator it = fHitList.begin(); it != fHitList.end();
       ++it ) {
    UInt_t idx = *it;
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
				   UInt_t plane, Hit* hit )
{
  // Set pattern bins corresponding to the exact physical positions
  // between start and end (in m) in the given plane. 
  // Positions may range from 0.0 to width.
  // Associate these bins with given hit

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
  for( Int_t i = lo; i <= hi; ++i )
    AddHit( plane, i, hit );

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

  if( !A ) return 0;
  UInt_t planeA = A->GetPlaneNum();
  assert( planeA < fNplanes );
  Double_t maxdist = 0.0;
  UInt_t planeB = fNplanes;
  if( B ) {
    maxdist = A->GetMaxSlope() * (B->GetZ() - A->GetZ());
    planeB = B->GetPlaneNum();
    assert( planeB < fNplanes );
  }
  bool do_single_hits =
    ( A->GetMWDC()->TestBit(MWDC::kPairsOnly) == kFALSE || B == 0 );

  Int_t nhits = 0;

  HitPairIter it( A->GetHits(), B ? B->GetHits() : 0, maxdist );
  while( it ) {
    // At least one hit found
    nhits++;
    Hit* hitA = static_cast<Hit*>((*it).first);
    Hit* hitB = static_cast<Hit*>((*it).second);
    assert( hitA || hitB );
    if( hitA && hitB ) {
      // A pair of hits registered in partner planes. One or more combinations
      // of hit positions may be within maxdist of each other.
      // Determine which combinations these are and set their bins.
      int set = 0, bitA, bitB;
      Double_t posA, posB;
      for( int i = 4; i; ) { --i;
	if( i&2 ) { posA = hitA->GetPosL(); bitA = 8; }
	else      { posA = hitA->GetPosR(); bitA = 4; }
	if( i&1 ) { posB = hitB->GetPosL(); bitB = 2; }
	else      { posB = hitB->GetPosR(); bitB = 1; }
	if( TMath::Abs( posA-posB ) <= maxdist ) {
	  if( (bitA & set) == 0 ) {
	    SetPosition( posA+fOffset, hitA->GetResolution(), planeA, hitA );
	    // Prevent duplicate entries for zero-drift hits
	    if( hitA->GetDriftDist() == 0 )
	      set |= 12;
	    else
	      set |= bitA;
	  }
	  if( (bitB & set) == 0 ) {
	    SetPosition( posB+fOffset, hitB->GetResolution(), planeB, hitB );
	    if( hitB->GetDriftDist() == 0 )
	      set |= 3;
	    else
	      set |= bitB;
	  }
	}
      }
    }
    else if( do_single_hits ) {
      // Unpaired hit in only one plane
      if( hitA ) {
	SetPosition( hitA->GetPosL()+fOffset, hitA->GetResolution(),
		     planeA, hitA );
	SetPosition( hitA->GetPosR()+fOffset, hitA->GetResolution(),
		     planeA, hitA );
      } else {
	SetPosition( hitB->GetPosL()+fOffset, hitB->GetResolution(),
		     planeB, hitB );
	SetPosition( hitB->GetPosR()+fOffset, hitB->GetResolution(),
		     planeB, hitB );
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

