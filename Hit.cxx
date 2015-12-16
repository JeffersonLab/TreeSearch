//*-- Author :    Ole Hansen, Jefferson Lab   27-Jun-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hit.h"
#include "Road.h"
#include "TSeqCollection.h"

#include <iostream>

using std::cout;
using std::endl;
using std::make_pair;

ClassImp(TreeSearch::Hit)
ClassImp(TreeSearch::HitPairIter)
ClassImp(TreeSearch::HitSet)

namespace TreeSearch {

//_____________________________________________________________________________
void Hit::Print( Option_t* opt ) const
{
  // Print hit info

  cout << "Hit: plane=" << (fPlane ? fPlane->GetName() : "??")
       << " pos="       << GetPos()
       << " z="         << GetZ()
       << " res="       << GetResolution();
  if( *opt != 'C' )
    cout << endl;
}

//_____________________________________________________________________________
UInt_t Hit::GetNumPos() const
{
  // Number of positions available from this hit. For generic hits, this is 1.
  // Wire chamber hits may override this to provide 2 positions due to L/R
  // ambiguity.

  return 1;
}

//_____________________________________________________________________________
#ifdef NDEBUG
Double_t Hit::GetPosI( UInt_t ) const
#else
Double_t Hit::GetPosI( UInt_t i ) const
#endif
{
  // Get the i-th position. For generic hits, i must be 0, and the function
  // always returns GetPos()

  assert( i<1 );
  return GetPos();
}

//_____________________________________________________________________________
Double_t FitCoord::GetChi2() const
{
  // Return chi2 of the fit that used this coordinate

  return fRoad ? fRoad->GetChi2() : kBig;
}

//_____________________________________________________________________________
HitPairIter::HitPairIter( const TSeqCollection* collA,
			  const TSeqCollection* collB,
			  Double_t maxdist )
  : fCollA(collA), fCollB(collB), fIterA(0), fIterB(0),
    fSaveIter(0), fSaveHit(0), fMaxDist(maxdist), fStarted(kFALSE),
    fScanning(kFALSE)
{
  // Constructor

  if( !fIterA && fCollA )
    fIterA = fCollA->MakeIterator();
  if( !fIterB && fCollB ) {
    fIterB = fCollB->MakeIterator();
    if( !fSaveIter )
      fSaveIter = fCollB->MakeIterator();
  }
  // Initialize our state so we point to the first item
  Next();
}

//_____________________________________________________________________________
HitPairIter::HitPairIter( const HitPairIter& rhs )
  : fCollA(rhs.fCollA), fCollB(rhs.fCollB), fIterA(0), fIterB(0),
    fSaveIter(0), fSaveHit(rhs.fSaveHit),
    fMaxDist(rhs.fMaxDist), fStarted(rhs.fStarted), fScanning(rhs.fScanning),
    fCurrent(rhs.fCurrent), fNext(rhs.fNext)
{
  // Copy ctor

  if( fCollA ) {
    fIterA = fCollA->MakeIterator();
    *fIterA = *rhs.fIterA;
  }
  if( fCollB ) {
    fIterB = fCollB->MakeIterator();
    *fIterB = *rhs.fIterB;
    fSaveIter = fCollB->MakeIterator();
    *fSaveIter = *rhs.fSaveIter;
  }
}

//_____________________________________________________________________________
HitPairIter& HitPairIter::operator=( const HitPairIter& rhs )
{
  // Assignment operator

  if( this != &rhs ) {
    fCollA     = rhs.fCollA;
    fCollB     = rhs.fCollB;
    fMaxDist   = rhs.fMaxDist;
    fStarted   = rhs.fStarted;
    fScanning  = rhs.fScanning;
    fCurrent   = rhs.fCurrent;
    fNext      = rhs.fNext;
    delete fIterA;
    delete fIterB;
    delete fSaveIter;
    if( fCollA ) {
      fIterA = fCollA->MakeIterator();
      *fIterA = *rhs.fIterA;
    }
    if( fCollB ) {
      fIterB = fCollB->MakeIterator();
      *fIterB = *rhs.fIterB;
      fSaveIter = fCollB->MakeIterator();
      *fSaveIter = *rhs.fSaveIter;
    }
  }
  return *this;
}

//_____________________________________________________________________________
HitPairIter::~HitPairIter()
{
  // Destructor

  delete fIterA;
  delete fIterB;
  delete fSaveIter;
}


//_____________________________________________________________________________
void HitPairIter::Reset()
{
  // Reset the iterator to the start

  fStarted = fScanning = kFALSE;
  if( fIterA )
    fIterA->Reset();
  if( fIterB )
    fIterB->Reset();
  // Our initial state is to point to the first object
  Next();
}


//_____________________________________________________________________________
HitPairIter& HitPairIter::Next()
{
  // Return next pair of hits along the wire plane. If a hit in either
  // plane is unpaired (no matching hit on the other plane within maxdist)
  // then only that hit is set in the returned pair object. If both
  // hits returned are zero, then there are no more hits in either plane.

  if( !fStarted ) {
    fNext = make_pair( fIterA ? fIterA->Next() : 0,
		       fIterB ? fIterB->Next() : 0 );
    fStarted = kTRUE;
  }

  fCurrent = fNext;
  Hit* hitA = static_cast<Hit*>( fCurrent.first );
  Hit* hitB = static_cast<Hit*>( fCurrent.second );

  if( hitA && hitB ) {
    switch( hitA->Compare(hitB,fMaxDist) ) {
    case -1: // A<B
      fNext.first  = fIterA->Next();
      fCurrent.second = 0;
      break;
    case  1: // A>B
      fNext.second = fIterB->Next();
      fCurrent.first = 0;
      break;
    default: // A==B
      {
	// Found a pair
	Hit* nextB = static_cast<Hit*>( fIterB->Next() );
	if( !nextB || hitA->Compare(nextB,fMaxDist) < 0 ) {
	  if( fScanning ) {
	    // End of a scan of plane B with fixed hitA
	    fScanning = kFALSE;
	    // Return B to the point where we started the scan
	    *fIterB = *fSaveIter;
	    hitB = fSaveHit;
	    // Advance to the next hit in A
	    hitA = static_cast<Hit*>( fIterA->Next() );
	    // Advance B until either B >= A or B == nextB (the hit in B that
	    // ended the prior scan), whichever comes first.
	    // The Bs for which saveB <= B < nextB have been paired with the
	    // prior A in the prior scan and so can't be considered unpaired,
	    // but they might pair with the new A, still.
	    if( hitA ) {
	      while( hitB != nextB && hitB->Compare(hitA,fMaxDist) < 0 )
		hitB = static_cast<Hit*>( fIterB->Next() );
	    } else {
	      // Of course, if there are no more hits in A, we only have to
	      // scan the rest of the Bs.
	      hitB = nextB;
	    }
	    fNext = make_pair( hitA, hitB );

	  } else {
	    // This is the normal case: nextB > hitA (usually true for
	    // small maxdist). So hitA/hitB are a pair, and
	    // the next hits to consider are the next ones in each plane.
	    fNext = make_pair( fIterA->Next(), nextB );
	  }
	} else {
	  // A==B and A==nextB, so more than one B matches this A.
	  // Start scanning mode where we keep hitA fixed and walk along
	  // B as long as B==A.
	  if( !fScanning ) {
	    fScanning = kTRUE;
	    // Save the starting point of the scan. We have to save both
	    // the iterator and the object because ROOT's TIterator
	    // is rather dumb (no Current() method).
	    *fSaveIter = *fIterB;
	    fSaveHit = hitB;
	  }
	  // nextB != 0 and hitA == nextB guaranteed here, so the
	  // next iteration will give hitA==hitB again.
	  fNext.second = nextB;
	}
	break;
      }
    }
  } else if( hitA ) {
    fNext.first  = fIterA->Next();
  } else if( hitB ) {
    fNext.second = fIterB->Next();
  }

  return *this;
}

//_____________________________________________________________________________
UInt_t HitSet::GetAltMatchValue( const Hset_t& hits )
{
  // Return plane occupancy pattern of given hitset using
  // Plane::GetAltPlaneNum(), i.e. the plane index including dummy planes.
  // Currently only used for asserting that TreeSearch found the correct
  // plane pattern.

  UInt_t curpat = 0;
  for( Hset_t::const_iterator it = hits.begin(); it != hits.end(); ++it )
    curpat |= 1U << (*it)->GetPlane()->GetAltPlaneNum();

  return curpat;
}

//_____________________________________________________________________________
Bool_t HitSet::IsSimilarTo( const HitSet& tryset, Int_t /* maxdist */ ) const
{
  // If maxdist == 0:
  // Similar to STL includes() algorithm, but allows tryset to have additional
  // hits in a given wire plane if there is at least one included hit in that
  // plane.
  //
  // Example: the following matches, despite the extra hit in "try":
  //   this:  30/   32/40/50/51
  //   try:   --/31 32/40/50/51
  //
  // Standard includes() implies intersection == set2.
  // This algorithm tests planepattern(intersection) == planepattern(set2)
  //
  // If maxdist > 0:
  // Same as above, but consider hits "equal" not only if they are identical
  // but also if their wire numbers are at most maxdist apart.
  //
  // Example: the following two patters are "similar" with maxdist = 1:
  //   this:  30/32/40/50/51
  //   try:   31/32/40/50/51
  //
  // This mode can be used to build "clusters" of patterns.

  assert( tryset.plane_pattern );

  Hset_t::const_iterator ihits = hits.begin();
  Hset_t::const_iterator ehits = hits.end();
  Hset_t::const_iterator itry  = tryset.hits.begin();
  Hset_t::const_iterator etry  = tryset.hits.end();
  //  Hset_t::key_compare    comp  = hits.key_comp();
  Hit::PosIsLess comp;

  UInt_t intersection_pattern = 0;

  while( ihits != ehits and itry != etry ) {
    if( comp(*itry, *ihits) )
      ++itry;
    else if( comp(*ihits, *itry) )
      ++ihits;
    else {
      intersection_pattern |= 1U << (*itry)->GetPlaneNum();
      ++ihits;
      ++itry;
    }
  }
  return tryset.plane_pattern == intersection_pattern;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

