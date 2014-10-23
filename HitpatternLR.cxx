//*-- Author :    Ole Hansen, Jefferson Lab   11-Feb-2013

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::HitpatternLR                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "HitpatternLR.h"
#include "WireHit.h"
#include "Plane.h"
#include "TMath.h"
#include "MWDC.h" // for MWDC bits test

#include <iostream>
#include <stdexcept>
#include <algorithm>

using namespace std;

ClassImp(TreeSearch::HitpatternLR)

namespace TreeSearch {

//_____________________________________________________________________________
HitpatternLR::HitpatternLR( const PatternTree& pt )
  : Hitpattern(pt)
{
  // Construct Hitpattern using paramaters of pattern tree

}

//_____________________________________________________________________________
HitpatternLR::HitpatternLR( UInt_t nlevels, UInt_t nplanes, Double_t width )
  : Hitpattern(nlevels, nplanes, width)
{
  // Constructor

}

//_____________________________________________________________________________
HitpatternLR::HitpatternLR( const Hitpattern& orig )
  : Hitpattern(orig)
{
  // Copy ctor

}

//_____________________________________________________________________________
Hitpattern& HitpatternLR::operator=( const Hitpattern& rhs )
{
  // Assignment

  return Hitpattern::operator=( rhs );
}

//_____________________________________________________________________________
HitpatternLR::~HitpatternLR()
{
  // Destructor

}

//_____________________________________________________________________________
Int_t HitpatternLR::Fill( const vector<Plane*>& planes )
{
  // Fill this hitpattern from hits in the given planes, all assumed to belong
  // to the same projecion. Returns the total number of hits processed.

#ifndef NDEBUG
  Projection* proj = 0;
#endif

  Int_t ntot = 0;
  for( vector<Plane*>::const_iterator it = planes.begin();
       it != planes.end(); ++it ) {
    Plane* plane = *it;
    assert( plane );
    if( plane->IsDummy() )
      SETBIT(fDummyPlanePattern, plane->GetPlaneNum());
    // Hitpatterns interprets hits as wire hits with L/R ambiguity
    Plane* partner = plane->GetPartner();
    // If a plane has a partner (usually with staggered wires), scan them
    // together to resolve some of the L/R ambiguity of the hit positions
    ntot += ScanHits( plane, partner );
    // If the partner plane was just scanned, don't scan it again
    if( partner ) {
      ++it;
      assert( it != planes.end() );
      assert( *it == partner );
    }
#ifndef NDEBUG
    // Bugcheck for projection mismatch
    if( proj ) {
      assert( plane->GetProjection() == proj );
      if( partner )
	assert( partner->GetProjection() == proj );
    } else {
      proj = plane->GetProjection();
      assert( proj );
    }
#endif
  }

  return ntot;
}

//_____________________________________________________________________________
Int_t HitpatternLR::ScanHits( Plane* A, Plane* B )
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
  assert( A->GetTracker() );
  bool do_single_hits =
    ( A->GetTracker()->TestBit(MWDC::kPairsOnly) == kFALSE || B == 0 );

  Int_t nhits = 0;

  HitPairIter it( A->GetHits(), B ? B->GetHits() : 0, maxdist );
  while( it ) {
    // At least one hit found
    nhits++;
    assert( !(*it).first  || dynamic_cast<WireHit*>((*it).first)  );
    assert( !(*it).second || dynamic_cast<WireHit*>((*it).second) );
    WireHit* hitA = static_cast<WireHit*>((*it).first);
    WireHit* hitB = static_cast<WireHit*>((*it).second);
    assert( hitA || hitB );
    // Don't record the pseudo-hits in dummy planes
    WireHit* recA = A->IsDummy() ? 0 : hitA;
    WireHit* recB = B->IsDummy() ? 0 : hitB;
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
	    // Cover +/- 2 sigma position range
	    SetPosition( posA+fOffset, 2.*hitA->GetResolution(), planeA, recA );
	    // Prevent duplicate entries for zero-drift hits
	    if( hitA->GetDriftDist() == 0 )
	      set |= 12;
	    else
	      set |= bitA;
	  }
	  if( (bitB & set) == 0 ) {
	    SetPosition( posB+fOffset, 2.*hitB->GetResolution(), planeB, recB );
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
	SetPosition( hitA->GetPosL()+fOffset, 2.*hitA->GetResolution(),
		     planeA, recA );
	SetPosition( hitA->GetPosR()+fOffset, 2.*hitA->GetResolution(),
		     planeA, recA );
      } else {
	SetPosition( hitB->GetPosL()+fOffset, 2.*hitB->GetResolution(),
		     planeB, recB );
	SetPosition( hitB->GetPosR()+fOffset, 2.*hitB->GetResolution(),
		     planeB, recB );
      }
    }
    ++it;
  }
  return nhits;
}


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

