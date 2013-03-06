//*-- Author :    Ole Hansen  03-December-2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// SoLID::SolSpec                                                       //
//                                                                      //
// SoLID spectrometer class containing GEM trackers                     //
// Used for testing.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "SolSpec.h"
#include "SoLIDGEMTracker.h"

#include <string>
#include <sstream>
#include <exception>

using namespace std;

ClassImp(SoLID::SolSpec)

namespace SoLID {

//_____________________________________________________________________________
SolSpec::SolSpec( const char* name, const char* description,
		  UInt_t nsectors )
  : THaSpectrometer( name, description )
{
  // Constructor. Define a GEM tracker detector for each of the 'nsectors'
  // sectors

  static const char* const here = "SolSpec";

  if( nsectors > 64 ) {
    Error( Here(here), "Number of sectors = %u too large. Must be <= 64. "
	   "Creating SoLID spectrometer failed.", nsectors );
    MakeZombie();
    throw range_error("SolSpec: number of sectors too large");
  }
  if( nsectors == 0 ) {
    Warning( Here(here), "Creating SoLID spectrometer with zero sectors" );
  }

  for( UInt_t i = 0; i < nsectors; ++i ) {
    stringstream sn, sd;
    sn << "tracker." << i;
    sd << "SoLID tracker in sector " << i;
    Int_t ret = AddDetector( new GEMTracker( sn.str().c_str(), 
					     sd.str().c_str()) );
    if( ret != 0 ) {
      stringstream s;
      s << "Error adding GEM tracker \"" << sn.str() << "\" "
	<< "to SoLID spectrometer \"" << GetName() << "\"";
      MakeZombie();
      throw logic_error(s.str());
    }
  }
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
