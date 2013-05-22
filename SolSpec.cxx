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

#include "TList.h"
#include "THaGlobals.h"
#include "THaTextvars.h"

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

  for( UInt_t i = 1; i <= nsectors; ++i ) {
    stringstream sn, sd;
    sn << "tracker." << i;
    sd << "SoLID tracker in sector " << i;
    GEMTracker* theTracker = new GEMTracker( sn.str().c_str(), 
					     sd.str().c_str() );
    theTracker->SetSectorNumber(i-1);
    Int_t ret = AddDetector( theTracker );
    if( ret != 0 ) {
      stringstream s;
      s << "Error adding GEM tracker \"" << sn.str() << "\" "
	<< "to SoLID spectrometer \"" << GetName() << "\"";
      MakeZombie();
      throw logic_error(s.str());
    }
  }

  fTrackCoord = new TClonesArray( "SoLID::SolTrackCoords", kInitTrackMultiplicity );

  // For now, don't require run database
  // We might need it later to read things like the field setting
  fProperties &= ~kNeedsRunDB;
}

//_____________________________________________________________________________
SolSpec::~SolSpec()
{
  // Destructor

  delete fTrackCoord; fTrackCoord = 0;

  DefineVariables( kDelete );
}

//_____________________________________________________________________________
void SolSpec::Clear( Option_t* opt )
{
  // Clear event-by-event data

  THaSpectrometer::Clear(opt);
  fTrackCoord->Clear();
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus SolSpec::Init( const TDatime& run_time )
{
  // Extra initialization for the SoLID spectrometer.
  //
  // To simplify database access, automatically define the following
  // gHaTextvars macros:
  // "allsectors": all sector numbers
  // "plane1": all supported projection types suffixed by "1"
  // "plane2"..."plane5": dto. with suffixes "2"..."5"

  Int_t nsect = fDetectors->GetSize();
  if( gHaTextvars != 0 ) {
    if( nsect > 0 ) {
      stringstream s;
      for( Int_t i = 1; i <= nsect; ++i ) {
	s << i;
	if( i+1 <= nsect )
	  s << ",";
      }
      assert( s && !s.str().empty() );
      gHaTextvars->Set( "allsectors", s.str() );
    }

    for( Int_t i = 1; i <= 5; ++i ) {
      stringstream s, p;
      Int_t ntypes = TreeSearch::kTypeEnd-TreeSearch::kTypeBegin;
      assert( ntypes > 0 );
      for( Int_t j = 0; j < ntypes; ++j ) {
	s << TreeSearch::kProjParam[j].name << i;
	if( j != ntypes-1 )
	  s << ",";
      }
      p << "plane" << i;
      assert( s && !s.str().empty() );
      assert( p && !p.str().empty() );
      gHaTextvars->Set( p.str(), s.str() );
    }
  }

  // Proceed with normal spectrometer initialization
  // EStatus ret = THaSpectrometer::Init( run_time );
  // if( ret != kOK )
  //   return ret;
  return THaSpectrometer::Init( run_time );

  // TODO: set up text variables like "plane1" etc. for output 
  // definitions & cuts, using actual plane names defined
  // (quite a job)

  // return kOK;
}

//_____________________________________________________________________________
Int_t SolSpec::DefineVariables( EMode mode )
{
  // Define/delete standard variables for a spectrometer (tracks etc.)
  // Can be overridden or extended by derived (actual) apparatuses

  if( mode == kDefine && fIsSetup ) return kOK;

  // Define standard spectrometer variables (track array)
  if( mode == kDefine )
    THaSpectrometer::DefineVariables(mode);

  RVarDef vars[] = {
    { "tr.sect",      "Sector number",
      "fTrackCoord.SoLID::SolTrackCoords.fSector" },
    { "tr.r",         "Transv dist (m)",
      "fTrackCoord.SoLID::SolTrackCoords.fRtransv" },
    { "tr.theta",     "Polar angle (rad)",
      "fTrackCoord.SoLID::SolTrackCoords.fTheta" },
    { "tr.phi",       "Azimuth (rad)",
      "fTrackCoord.SoLID::SolTrackCoords.fPhi" },
    { "tr.phi_rot",   "Azimuth in sector (rad)",
      "fTrackCoord.SoLID::SolTrackCoords.fPhiRot" },
    { "tr.thdir",     "Direction polar angle (rad)",
      "fTrackCoord.SoLID::SolTrackCoords.fThetaDir" },
    { "tr.phdir",     "Direction azimuth (rad)",
      "fTrackCoord.SoLID::SolTrackCoords.fPhiDir" },
    { "tr.phdir_rot", "Direction azimuth in sector (rad)",
      "fTrackCoord.SoLID::SolTrackCoords.fPhiDirRot" },
    { 0 }
  };

  return DefineVarsFromList( vars, mode );
}

//_____________________________________________________________________________
Int_t SolSpec::ReadRunDatabase( const TDatime& date )
{
  // Dummy run database reader to override the THaSpectrometer function.

  // This is a bit of a kludge. SolSpec actually shouldn't inherit from
  // THaSpectrometer, but from a different spectrometer base class that
  // doesn't assume small-angle focusing optics.

  return THaApparatus::ReadRunDatabase( date );
}

//_____________________________________________________________________________
#ifdef NDEBUG
Int_t SolSpec::FindVertices( TClonesArray& /* tracks */)
#else
Int_t SolSpec::FindVertices( TClonesArray& tracks )
#endif
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  assert( tracks.GetLast() == fTrackCoord->GetLast() );

  return 0;
}

//_____________________________________________________________________________
Int_t SolSpec::TrackCalc()
{
  // Additional track calculations

  // TODO

  return 0;
}

//_____________________________________________________________________________

} // end namespace SoLID
