///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hit.h"
#include "TimeToDistConv.h"
#include "TSeqCollection.h"

#include <iostream>

using std::cout;
using std::endl;
using std::make_pair;

ClassImp(TreeSearch::Hit)
ClassImp(TreeSearch::MCHit)

namespace TreeSearch {

//_____________________________________________________________________________
Int_t Hit::Compare( const TObject* obj ) const 
{
  // Used for sorting hits in a TSeqCollection (e.g. TList, TClonesArray).
  // A hit is "less than" another hit if its position is smaller.

  if( !obj || IsA() != obj->IsA() )
    return -1;

  const Hit* rhs = static_cast<const Hit*>( obj );
 
  if( fPos < rhs->fPos ) return -1;
  if( fPos > rhs->fPos ) return  1;
  return 0;
}

#if 0
  // Also, for hits on the same wire, the first hit on the wire (the one with
  // the smallest time) is "less than" one with a higher time.  If the hits
  // are sorted according to this scheme, they will be in order of increasing
  // wire number and, for each wire, will be in the order in which they hit
  // the wire
  Int_t myWireNum = fWire->GetNum();
  Int_t hitWireNum = hit->GetWire()->GetNum();
  // Compare wire numbers
  if (myWireNum < hitWireNum) return -1;
  if (myWireNum > hitWireNum) return  1;
  if (myWireNum == hitWireNum) {
    // If wire numbers are the same, compare times
    Double_t hitTime = hit->GetTime();
    if (fTime < hitTime) return -1;
    if (fTime > hitTime) return  1;
  }
  return 0;
}
#endif


//_____________________________________________________________________________
Double_t Hit::ConvertTimeToDist( Double_t slope )
{
  // Convert drift time to drift distance. 'slope' is the approximate
  // slope of the track.
  // Updates the internal variables fPosL, fPosR and fResolution.
  // Must be called before doing analysis of drift chamber hits.

  TimeToDistConv* ttd;
  Double_t dist;
  if( fWirePlane && (ttd = fWirePlane->GetTTDConv()) ) {
    dist = ttd->ConvertTimeToDist( fTime, slope );
  } else {
    dist = 0.0;
  }
  if( dist > 0.0 ) {
    fPosL = fPos-dist;
    fPosR = fPos+dist;
  } else {
    fPosL = fPosR = fPos;
  }
  return dist;
}

//_____________________________________________________________________________
void Hit::Print( Option_t* opt ) const
{
  // Print hit info

  cout << "Hit: wnum=" << GetWireNum()
       << " wpos=" << GetWirePos()
       << " z=" << GetZ()
       << " res=" << GetResolution()
       << " time=" << GetDriftTime()
       << " drift=" << GetDriftDist()
       << " trk="  << GetTrackDist();
  if( *opt != 'C' )
    cout << endl;
}

//_____________________________________________________________________________
void MCHit::Print( Option_t* opt ) const
{
  // Print hit info

  Hit::Print("C");
  cout << " MCpos=" << GetMCPos()
       << endl;
}

//_____________________________________________________________________________

}

ClassImp(TreeSearch::HitPairIter)

namespace TreeSearch {

//_____________________________________________________________________________
HitPairIter::HitPairIter( const TSeqCollection* collA,
			  const TSeqCollection* collB,
			  Double_t maxdist ) 
  : fCollA(collA), fCollB(collB), fIterA(NULL), fIterB(NULL),
    fSaveIter(NULL), fSaveHit(NULL), fMaxDist(maxdist), fStarted(kFALSE),
    fScanning(kFALSE)
{
  // Constructor
}

//_____________________________________________________________________________
HitPairIter::HitPairIter( const HitPairIter& rhs )
  : fCollA(rhs.fCollA), fCollB(rhs.fCollB), fSaveHit(rhs.fSaveHit),
    fMaxDist(fMaxDist), fStarted(rhs.fStarted), fScanning(rhs.fScanning),
    fCurrent(rhs.fCurrent), fNext(rhs.fNext)
{
  // Copy ctor

  if( fCollA ) {
    fIterA = fCollA->MakeIterator();
    *fIterA = *rhs.fIterA;
  }
  if( fCollB ) {
    fIterB = fCollA->MakeIterator();
    *fIterB = *rhs.fIterB;
    fSaveIter = fCollA->MakeIterator();
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
      fIterB = fCollA->MakeIterator();
      *fIterB = *rhs.fIterB;
      fSaveIter = fCollA->MakeIterator();
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
}


//_____________________________________________________________________________
ObjPair_t& HitPairIter::Next()
{
  // Return next pair of hits along the wire plane. If a hit in either
  // plane is unpaired (no matching hit on the other plane within maxdist)
  // then only that hit is set in the returned pair object. If both
  // hits returned are zero, then there are no more hits in either plane.

  if( !fIterA && fCollA )
    fIterA = fCollA->MakeIterator();
  if( !fIterB && fCollB ) {
    fIterB = fCollB->MakeIterator();
    if( !fSaveIter )
      fSaveIter = fCollB->MakeIterator();
  }
  if( !fStarted ) {
    fNext = make_pair( fIterA ? fIterA->Next() : NULL, 
		       fIterB ? fIterB->Next() : NULL );
    fStarted = kTRUE;
  }

  fCurrent = fNext;
  Hit* hitA = static_cast<Hit*>( fCurrent.first );
  Hit* hitB = static_cast<Hit*>( fCurrent.second );

  if( hitA && hitB ) {
    switch( hitA->Compare(hitB,fMaxDist) ) {
    case -1: // A<B 
      fNext.first  = fIterA->Next();
      fCurrent.second = NULL;
      break;
    case  1: // A>B
      fNext.second = fIterB->Next();
      fCurrent.first = NULL;
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
	  // nextB != NULL and hitA == nextB guaranteed here, so the
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

  return fCurrent;
}

} // end namespace TreeSearch

///////////////////////////////////////////////////////////////////////////////
