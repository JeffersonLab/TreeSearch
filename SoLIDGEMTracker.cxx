//*-- Author :    Ole Hansen, Jefferson Lab   06-Feb-2012

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SoLID::GEMTracker                                                         //
//                                                                           //
// This is class GEMTracker in namespace SoLID. It inherits from class       //
// GEMTracker in namespace TreeSearch. Although confusing, this is fine.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SoLIDGEMTracker.h"
#include "SoLIDGEMPlane.h"
#include "SolSpec.h"
#include "TClonesArray.h"

using namespace std;

namespace SoLID {

typedef vector<Plane*>::size_type vrsiz_t;

//_____________________________________________________________________________
GEMTracker::GEMTracker( const char* name, const char* desc, THaApparatus* app )
  : TreeSearch::GEMTracker(name,desc,app), fSector(-1), fPhi(0)
{
  // Constructor
}

//_____________________________________________________________________________
GEMTracker::~GEMTracker()
{
  // Destructor
}

//_____________________________________________________________________________
const char* GEMTracker::GetDBFileName() const
{
  // Return database file name prefix. For SoLID trackers, this is the detector
  // name with any trailing sector number removed. In this way, all trackers
  // use the same database file.

  // fDBPrefix is set in MakePrefix()
  return fDBPrefix.Data();
}

//_____________________________________________________________________________
void GEMTracker::MakePrefix()
{
  // Set up name prefixes for global variables and database.
  // Global variables and database keys get the standard prefix,
  // e.g. "solid.tracker.3."
  // The database file name gets a different, shorter prefix, allowing
  // the trackers to share a single database file. For the example above,
  // "solid.tracker."

  TreeSearch::GEMTracker::MakePrefix();

  TString prefix( GetPrefix() );
  assert( prefix.EndsWith(".") );
  Int_t ndot = prefix.CountChar('.');
  if( ndot > 1 ) {
    prefix.Chop();
    if( !prefix.EndsWith(".") )
      prefix.Remove( prefix.Last('.')+1 );
    else
      Warning( Here("MakePrefix"), "Double dot in detector prefix = "
	       "\"%s\"?", GetPrefix() );
  }
  fDBPrefix = prefix;
}

//_____________________________________________________________________________
Plane* GEMTracker::MakePlane( const char* name, const char* description,
			      THaDetectorBase* parent ) const
{
  // Create an object of the plane class used by this implementation

  return new SoLID::GEMPlane( name, description, parent );
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus GEMTracker::PartnerPlanes()
{
  // Before doing PartnerPlanes() for GEM trackers, do additional geometry
  // calculations. This is done here because the planes need to have read
  // their database and be sorted by z.

  // TODO: protect against repeated calls? (Don't keep subtracting vectors)

  // Take the origin of the first plane as the origin of this Tracker and
  // make all the plane coordinates relative to the Tracker origin
  if( !fPlanes.empty() ) {
    // Bypass the const& - presumably we know what we're doing...
    TVector3& org = const_cast<TVector3&>( fPlanes.front()->GetOrigin() );
    fOrigin = org;
    org.SetXYZ( 0.0, 0.0, 0.0 );
    for( vrsiz_t iplane = 1; iplane < fPlanes.size(); ++iplane ) {
      Plane* thePlane = fPlanes[iplane];
      assert( thePlane->GetZ() >= fOrigin.Z() ); // else not sorted
      TVector3& other_org = const_cast<TVector3&>( thePlane->GetOrigin() );
      other_org -= fOrigin;
    }
  }
  // fOrigin needs to be given in the global (lab) reference frame
  fOrigin *= fRotation;

  // Now continue with standard PartnerPlanes() of the base class
  return TreeSearch::GEMTracker::PartnerPlanes();
}

//_____________________________________________________________________________
Int_t GEMTracker::ReadGeometry( FILE* file, const TDatime& date,
				Bool_t /* required */ )
{
  // Read basic geometry for a GEM tracker sector
  //
  // The only geometry parameter needed for the a SoLID Tracker (which
  // represents the collection of GEM trackers in a sector) is 'phi', the
  // central phi angle of the sector. The origin of the tracker coordinate
  // system will be determined automatically in PartnerPlanes() using the
  // positions of the GEM planes that are defined for this tracker.
  //
  // 'phi' is always required. The 'required' argument is ignored.

  //  static const char* const here = "ReadGeometry";

  DBRequest request[] = {
    { "phi",  &fPhi,  kDouble, 0, 0, 0, "central phi of sector [deg]" },
    { 0 }
  };
  Int_t err = LoadDB( file, date, request, fPrefix );
  if( err )
    return err;

  // Keep phi in rad and normalize to -pi..pi
  fPhi = TVector2::Phi_mpi_pi( fPhi*TMath::DegToRad() );

  // fOrigin will be set later PartnerPlanes (after Plane init)

  // Define this detector's rotation wrt the global coordinate system.
  // The rotation rotates the axes and not the vectors, hence it is a rotation
  // about z by _positive_ phi.
  // Any vector in this Tracker's reference frame can be transformed to the
  // global system simply by multiplication with the rotation matrix:
  // v_global = fRotation * v

  fRotation.SetToIdentity();
  fRotation.RotateZ( fPhi );

  return kOK;
}

//_____________________________________________________________________________
Int_t GEMTracker::NewTrackCalc( THaTrack*, const TVector3& pos,
				const TVector3& dir )
{
  // For every new track, convert track coordinates to the cylindrical system
  // appropriate for SoLID. This is a temporary solution; the right way to do
  // this is to have the SoLID spectrometer use its own track class.

  assert( dynamic_cast<SolSpec*>(GetApparatus()) );
  SolSpec* solid = static_cast<SolSpec*>(GetApparatus());

  TClonesArray* trackInfo = solid->GetTrackInfo();
  Int_t i = trackInfo->GetLast()+1;

#ifdef MCDATA
  SolTrackInfo* trkinfo =
#endif
  new( (*trackInfo)[i] ) SolTrackInfo( fSector, pos, dir, fPhi );

#ifdef MCDATA
  assert( static_cast<vector<UInt_t>::size_type>(i) < fMCHitBits.size() );
  trkinfo->SetMCHitBits( fMCHitBits[i] );
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace SoLID

ClassImp(SoLID::GEMTracker)

