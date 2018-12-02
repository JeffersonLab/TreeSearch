//*-- Author :    Ole Hansen  03-December-2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// SBS::SBSSpec                                                         //
//                                                                      //
// SBS spectrometer class containing GEM trackers                       //
// Used for testing.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "SBSSpec.h"
#include "SBSGEMTracker.h"

#include "TList.h"
#include "THaGlobals.h"
#include "THaTextvars.h"
#include "THaTrack.h"

#include <string>
#include <sstream>
#include <exception>

using namespace std;

ClassImp(SBS::SBSSpec)

namespace SBS {

//_____________________________________________________________________________
SBSSpec::SBSSpec( const char* name, const char* description,
		  UInt_t nsectors )
  : THaSpectrometer( name, description )
{
  // Constructor. Define a GEM tracker detector for each of the 'nsectors'
  // sectors

  static const char* const here = "SBSSpec";

  if( nsectors > 64 ) {
    Error( Here(here), "Number of sectors = %u too large. Must be <= 64. "
	   "Creating SBS spectrometer failed.", nsectors );
    MakeZombie();
    throw range_error("SBSSpec: number of sectors too large");
  }
  if( nsectors == 0 ) {
    Warning( Here(here), "Creating SBS spectrometer with zero sectors" );
  }

  for( UInt_t i = 1; i <= nsectors; ++i ) {
    stringstream sn, sd;
    sn << "tracker." << i;
    sd << "SBS tracker in sector " << i;
    SBS::GEMTracker* theTracker = new SBS::GEMTracker( sn.str().c_str(),
							   sd.str().c_str() );
    theTracker->SetSectorNumber(i-1);
    Int_t ret = AddDetector( theTracker );
    if( ret != 0 ) {
      stringstream s;
      s << "Error adding GEM tracker \"" << sn.str() << "\" "
	<< "to SBS spectrometer \"" << GetName() << "\"";
      MakeZombie();
      throw logic_error(s.str());
    }
  }

  fSBSTrackInfo = new TClonesArray( "SBS::SBSTrackInfo", kInitTrackMultiplicity );

  // For now, don't require run database
  // We might need it later to read things like the field setting
  fProperties &= ~kNeedsRunDB;
}

//_____________________________________________________________________________
SBSSpec::~SBSSpec()
{
  // Destructor

  delete fSBSTrackInfo; fSBSTrackInfo = 0;

  DefineVariables( kDelete );
}

//_____________________________________________________________________________
void SBSSpec::Clear( Option_t* opt )
{
  // Clear event-by-event data

  THaSpectrometer::Clear(opt);
  fSBSTrackInfo->Clear();
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus SBSSpec::Init( const TDatime& run_time )
{
  // Extra initialization for the SBS spectrometer.
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
Int_t SBSSpec::DefineVariables( EMode mode )
{
  // Define/delete standard spectrometer variables (tracks etc.) plus
  // SBS-specific additional track variables (cylindrical coords etc.)

  if( mode == kDefine && fIsSetup ) return kOK;

  // Define standard spectrometer variables (track array)
  if( mode == kDefine )
    THaSpectrometer::DefineVariables(mode);

  RVarDef vars[] = {
    { "tr.sect",        "Sector number",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fSector" },
    { "tr.x_t",         "x transport (m)",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fX_T" },
    { "tr.y_t",         "y transport (m)",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fY_T" },
    { "tr.xdir_t",      "dx/dz transport",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fXdir_T" },
    { "tr.ydir_t",      "dy/dz transport",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fYdir_T" },
    { "tr.3dmatchbits", "bits patterns of 3D ampcorr matched 3d hits",
      "fSBSTrackInfo.SBS::SBSTrackInfo.f3dMatchBits" },
    { "tr.n3dmatch",    "number of 3D ampcorr matched 3d hits",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fN3dMatch" },
#ifdef MCDATA
    { "tr.mchitbits", "Planepattern of true MC hits in this track",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fMCHitBits" },
    { "tr.nmchits",   "Number of MC hits in this track",
      "fSBSTrackInfo.SBS::SBSTrackInfo.fNMCHits" },
#endif
    { 0 }
  };

  return DefineVarsFromList( vars, mode );
}

//_____________________________________________________________________________
Int_t SBSSpec::ReadRunDatabase( const TDatime& date )
{
  // Dummy run database reader to override the THaSpectrometer function.

  // This is a bit of a kludge. SBSSpec actually shouldn't inherit from
  // THaSpectrometer, but from a different spectrometer base class that
  // doesn't assume small-angle focusing optics.

  return THaApparatus::ReadRunDatabase( date );
}

//_____________________________________________________________________________
Int_t SBSSpec::FindVertices( TClonesArray& tracks )
{
  // Reconstruct target coordinates for all tracks found.

  assert( tracks.GetLast() == fSBSTrackInfo->GetLast() );

  for( Int_t i = 0; i < GetNTracks(); ++i ) {
    // VERTEX RECONSTRUCTION FOR STRAIGHT TRACKS
#ifdef NDEBUG  // this is asking for Heisenbugs....
    THaTrack* theTrack = static_cast<THaTrack*>( fTracks->UncheckedAt(i) );
    // The track coordinates (fX etc., called "focal plane coordinates"
    // in THaTrack) are in the lab frame, i.e. wrt to this spectrometer's
    // origin. We do need the tracker origin for its z-position, however ...
    Double_t z = theTrack->GetCreator()->GetOrigin().Z();
#else
    THaTrack* theTrack = dynamic_cast<THaTrack*>( fTracks->At(i) );
    assert( theTrack );
    THaTrackingDetector* tracker = theTrack->GetCreator();
    assert( tracker );
    Double_t z = tracker->GetOrigin().Z();
#endif
    // Now, define the vertex as the point of closest approach to the beam.
    // For the time being, take the beam to be exactly at (x,y) = (0,0)
    // At the closest approach of two skewed lines, the connecting segment
    // is simultaneously perpendicular to both.
    // Call the lines g(t) = g0 + g1*t, h(s) = h0 + h1*s
    // where g1 and h1 are unit vectors. Let h(s) be the beam.
    // Then at the closest approach,
    // tc = - (g0-h0)*(g1-h1(g1*h1))/(1-(g1*h1)^2)
    // sc =   (g0-h0)*(h1-g1(g1*h1))/(1-(g1*h1)^2)
    TVector3 g0( theTrack->GetX(), theTrack->GetY(), z );
    TVector3 g1( theTrack->GetTheta(), theTrack->GetPhi(), 1.0 );
    g1 = g1.Unit();
    TVector3 h0( 0, 0, 0  );
    TVector3 h1( 0, 0, 1. );
    Double_t gh = g1*h1;
    Double_t denom = 1.-gh*gh;
    if( TMath::Abs(denom) > 1e-6 ) {
      TVector3 D0 = g0-h0;
      Double_t tc = -D0*(g1-h1*gh)/denom;
      Double_t sc =  D0*(h1-g1*gh)/denom;
      TVector3 vertex = g0 + g1*tc;
      theTrack->SetVertex( vertex );
      theTrack->SetVertexError( vertex - h0-h1*sc );
    }
    // With straight tracks, we only know the momentum vector's direction
    theTrack->SetPvect(g1);
  }
  return 0;
}

//_____________________________________________________________________________
Int_t SBSSpec::TrackCalc()
{
  // Additional track calculations


  return 0;
}

//_____________________________________________________________________________

} // end namespace SBS
