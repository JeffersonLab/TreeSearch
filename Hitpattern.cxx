///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hitpattern.h"
#include "Hit.h"
#include "TSeqCollection.h"

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

  if( fNplanes > 0 && fDepth > 0 ) {
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
  delete fPattern;
}

//_____________________________________________________________________________
Int_t Hitpattern::SetPoints( TSeqCollection* hitsA, TSeqCollection* hitsB, 
			     UInt_t plane )
{
  // Set all bits at all depths of the hit pattern corresponding to
  // the hits in 'hitsA' in plane number 'plane'.
  // The second collection, 'hitsB', represents hits in a optional partner
  // plane of the 'hitsA' plane, usually a nearby second plane with staggered
  // wires. Naturally, the two planes and their wires must be parallel.

  return 0;


}

//_____________________________________________________________________________
void Hitpattern::Clear( Option_t* opt )
{
  // Clear the hitpattern

  for( UInt_t i=0; i<fNplanes; i++ )
    fPattern[i]->Clear();
}

//_____________________________________________________________________________
void Hitpattern::Print( Option_t* opt ) const
{
  // Print basic info about hitpattern.

}



///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

