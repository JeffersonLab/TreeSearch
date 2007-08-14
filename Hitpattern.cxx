///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hitpattern.h"
#include "Hit.h"
#include "TSeqCollection.h"
#include "MWDC.h"

#include <iostream>

using std::cout;
using std::endl;

ClassImp(TreeSearch::Hitpattern)

namespace TreeSearch {

//_____________________________________________________________________________
Hitpattern::Hitpattern( UInt_t depth, UInt_t nplanes, Double_t width )
  : fDepth(depth), fNplanes(nplanes), fWidth(width)
{
  // Constructor

  if( fNplanes > 0 ) {
    fPattern = new Bits*[fNplanes];
    for( UInt_t i=0; i<fNplanes; i++ )
      fPattern[i] = new Bits( 1U<<fDepth );
  } else {
    fPattern = NULL;
    fNplanes = 0;
  }
}

//_____________________________________________________________________________
Hitpattern::~Hitpattern()
{
  // Destructor

  for( UInt_t i=0; i<fNplanes; i++ )
    delete fPattern[i];
  delete [] fPattern;
}

//_____________________________________________________________________________
void Hitpattern::Clear( Option_t* opt )
{
  // Clear the hitpattern

  for( UInt_t i=0; i<fNplanes; i++ )
    fPattern[i]->FastClear();
}

//_____________________________________________________________________________
Int_t Hitpattern::SetPoints( WirePlane* A, WirePlane* B, UInt_t plane )
{
  // Set all bits at all depths of the hit pattern corresponding to
  // the hits in 'A' in plane number 'plane'.
  // The second, optional plane 'B' represents an optional partner
  // plane of 'A', usually a nearby plane with staggered wires. 
  // Naturally, the two planes and their wires must be parallel.
  //
  // Returns number of hits processed

  static const ObjPair_t null(0,0);
  
  Clear();

  if( !A )
    return 0;
  MWDC* mwdc = dynamic_cast<MWDC*>( A->GetDetector() );
  if( !mwdc )
    return 0;
  Double_t dz = B ? B->GetZ() - A->GetZ() : 0.0;
  Double_t maxslope = mwdc->GetMaxSlope();
  Double_t maxdist = maxslope * dz;

  Int_t nhits = 0;

  HitPairIter next( A->GetHits(), B ? B->GetHits() : NULL, maxdist );
  ObjPair_t& hitpair = next();
  while( hitpair != null ) {
    nhits++;
    
    hitpair = next();
  }


  return nhits;
}

//_____________________________________________________________________________
void Hitpattern::Print( Option_t* opt ) const
{
  // Print basic info about hitpattern.

}



///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

