///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hit.h"
#include "TSeqCollection.h"

#include <iostream>

using std::cout;
using std::endl;
using std::make_pair;

ClassImp(TreeSearch::Hit)

namespace TreeSearch {

//_____________________________________________________________________________
Int_t Hit::Compare( const TObject* obj ) const 
{
  // Used for sorting hits in a TSeqCollection (e.g. TList, TClonesArray).
  // A hit is "less than" another hit if it occurred on a lower wire number.
  // Also, for hits on the same wire, the first hit on the wire (the one with
  // the smallest time) is "less than" one with a higher time.  If the hits
  // are sorted according to this scheme, they will be in order of increasing
  // wire number and, for each wire, will be in the order in which they hit
  // the wire

  if( !obj || IsA() != obj->IsA() )
    return -1;

  const Hit* hit = static_cast<const Hit*>( obj );
 
  if( fPos < hit->fPos ) return -1;
  if( fPos > hit->fPos ) return  1;
  return 0;
}

#if 0
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


//_____________________________________________________________________________
Double_t THaVDCHit::ConvertTimeToDist(Double_t slope)
{
  // Converts TDC time to drift distance
  // Takes the (estimated) slope of the track as an argument
  
  THaVDCTimeToDistConv* ttdConv = (fWire) ? fWire->GetTTDConv() : NULL;
  
  if (ttdConv) {
    // If a time to distance algorithm exists, use it to convert the TDC time 
    // to the drift distance
    fDist = ttdConv->ConvertTimeToDist(fTime, slope, &fdDist);
    return fDist;
  }
  
  Error("ConvertTimeToDist()", "No Time to dist algorithm available");
  return 0.0;

}
#endif


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

  //FIXME: deal with fIterA, fIterB, fSaveIter

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
    //FIXME: deal with fIterA fIterB fSaveIter
    fIterA = fIterB = NULL;
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

  fStarted = kFALSE;
  if( fIterA )
    fIterA->Reset();
  if( fIterB )
    fIterB->Reset();
}


//_____________________________________________________________________________
ObjPair_t& HitPairIter::Next()
{
  // Return next pair of hits along the wire plane. If a hit in either
  // plane is unpaired (no matching hit on the other plane within maxdist
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
	    // ended the prior scan, whichever comes first.
	    // The Bs for which saveB <= B < nextB have been paired with the 
	    // prior A in the prior scan and so can't be considered unpaired,
	    // but they might pair with the new A, too.

	    if( hitA ) {
	      // NB: hitB < nextB guaranteed here, and if nextB is the end 
	      // of the list (i.e. nextB == NULL) this loop will end there.
	      while( hitB != nextB && hitB->Compare(hitA,fMaxDist) < 0 )
		hitB = static_cast<Hit*>( fIterB->Next() );
	    } else {
	      // Of course, if there are no more hits in A, we only have to 
	      // look at the rest of the Bs.
	      hitB = nextB;
	    }
	    fNext = make_pair( hitA, hitB );
	    // We are back to the same situation as at the start of the scan.
	    return Next();

	  } else {
	    // This is the normal case: nextB > hitA (always true for 
	    // maxdist = 0). So hitA/hitB are a pair, and
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
