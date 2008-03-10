//*-- Author :    Ole Hansen  10-March-2008

//////////////////////////////////////////////////////////////////////////
//
// TreeSearch::BigBite
//
// Bare-bones BigBite spectrometer class containing the TreeSearch::MWDC.
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "BigBite.h"
#include "MWDC.h"

using namespace std;

ClassImp(TreeSearch::BigBite)

namespace TreeSearch {

//_____________________________________________________________________________
BigBite::BigBite( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

  AddDetector( new MWDC("mwdc","BigBite multi-wire drift chamber") );
}

//_____________________________________________________________________________
BigBite::~BigBite()
{
  // Destructor
}

//_____________________________________________________________________________
  Int_t BigBite::FindVertices( TClonesArray& /* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t BigBite::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

//_____________________________________________________________________________

} // end namespace TreeSearch
