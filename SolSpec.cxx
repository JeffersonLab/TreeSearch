//*-- Author :    Ole Hansen  03-December-2012

//////////////////////////////////////////////////////////////////////////
//
// SoLID::SolSpec
//
// SoLID spectrometer class containing GEM trackers
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SolSpec.h"
#include "SoLIDGEMTracker.h"

using namespace std;

ClassImp(SoLID::SolSpec)

namespace SoLID {

//_____________________________________________________________________________
SolSpec::SolSpec( const char* name, const char* description )
  : THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

  //  AddDetector( new MWDC("mwdc","SoLID multi-wire drift chamber") );
}

//_____________________________________________________________________________
SolSpec::~SolSpec()
{
  // Destructor
}

//_____________________________________________________________________________
Int_t SolSpec::FindVertices( TClonesArray& /* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t SolSpec::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

//_____________________________________________________________________________

} // end namespace SoLID
