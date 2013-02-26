#ifndef ROOT_SoLID_SolSpec
#define ROOT_SoLID_SolSpec

//////////////////////////////////////////////////////////////////////////
//
// SoLID::SolSpec
//
// SoLID spectrometer class containing GEM trackers
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h"

namespace SoLID {

  class SolSpec : public THaSpectrometer {
  
  public:
    SolSpec( const char* name, const char* description );
    virtual ~SolSpec();

    virtual Int_t   FindVertices( TClonesArray& tracks );
    virtual Int_t   TrackCalc();
  
  protected:
    
    ClassDef(SolSpec,0) // SoLID spectrometer
  };

} // end namespace SoLID

#endif

