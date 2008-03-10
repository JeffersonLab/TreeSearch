#ifndef ROOT_TreeSearch_BigBite
#define ROOT_TreeSearch_BigBite

//////////////////////////////////////////////////////////////////////////
//
// TreeSearch::BigBite
//
// Bare-bones BigBite spectrometer class containing the TreeSearch::MWDC.
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h"

namespace TreeSearch {

  class BigBite : public THaSpectrometer {
  
  public:
    BigBite( const char* name, const char* description );
    virtual ~BigBite();

    virtual Int_t   FindVertices( TClonesArray& tracks );
    virtual Int_t   TrackCalc();
  
  protected:
    
    ClassDef(BigBite,0) // BigBite spectrometer
  };

} // end namespace TreeSearch

#endif

