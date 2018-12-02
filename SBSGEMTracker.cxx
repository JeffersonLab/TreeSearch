//*-- Author :    Ole Hansen, Jefferson Lab   06-Feb-2012

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBS::GEMTracker                                                         //
//                                                                           //
// This is class GEMTracker in namespace SBS. It inherits from class       //
// GEMTracker in namespace TreeSearch. Although confusing, this is fine.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSGEMTracker.h"
#include "SBSGEMPlane.h"
#include "SBSSpec.h"
#include "Helper.h"
#include "TClonesArray.h"
#include "TMatrixD.h"

using namespace std;

namespace SBS {

typedef vector<Plane*>::size_type vrsiz_t;

//_____________________________________________________________________________
GEMTracker::GEMTracker( const char* name, const char* desc, THaApparatus* app )
  : TreeSearch::GEMTracker(name,desc,app), fSector(-1), fXOffset(0)
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
  // Return database file name prefix. For SBS trackers, this is the detector
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

  return new SBS::GEMPlane( name, description, parent );
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
    org.SetXYZ( 0.0, 0.0, 0.0 );  // update first plane
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
				Bool_t required )
{
  // Read basic geometry for a SBS GEM tracker sector
  //
  // The geometry parameters needed for a SBS Tracker 
  // (which represents the collection of GEM trackers in a sector) are: 
  // X_offset, the shift in X_transport of the plane wrt the central ray of the spectrometer, 
  // and two angles: one of horizontal rotation (thetaH), which translates the SBS angle,
  // one of vertical rotation (thetaV), translating the "bending" of the SBS wrt the xOz plane.
  //
  // Those will allow to make the rotation to calculate the 
  
  
  DBRequest request[] = {
    {"dmag",        &fDMag,         kDouble, 0, 1},
    {"xoff",        &fXOffset,      kDouble, 0, 1},
    {"thetaH",      &fThetaH,       kDouble, 0, 1},
    {"thetaV",      &fThetaV,       kDouble, 0, 1},
    { 0 }
  };
  Int_t err = LoadDB( file, date, request, fPrefix );
  if( err )
    return err;
  
  fOrigin.SetX(fXOffset);
  
  Double_t torad = atan(1) / 45.0;
  fThetaH *= torad;
  fThetaV *= torad;
  
  //Setting up the rotations to calcuate the plane origin.
  Double_t arr_roty0[9] = {cos(fThetaH),  0, sin(fThetaH),
			  0,             1,            0,
			  -sin(fThetaH), 0, cos(fThetaH)};
  Double_t arr_rotx1[9] = {1,            0,            0,
			   0, cos(fThetaV), -sin(fThetaV),
			   0, sin(fThetaV),  cos(fThetaV)};
  Double_t arr_rotz2[9] = {0, -1,  0,
  			   1,  0,  0,
  			   0,  0,  1};
  
  TMatrixD Roty0(3,3,arr_roty0);// rotation along hall pivot (y): spectrometer theta
  TMatrixD Rotx1(3,3,arr_rotx1);// rotation along x': spectrometer bending
  TMatrixD Rotz2(3,3,arr_rotz2);// rotation along z": box rotation
  TMatrixD Rotzx(3,3,arr_rotz2);
  TMatrixD fRotMat_LB(3,3, arr_roty0);
  
  Rotzx.Mult(Rotz2, Rotx1);
  fRotMat_LB.Mult(Rotzx, Roty0);// Box to Lab transformation
  
  TMatrixD fRotMat_BL = fRotMat_LB;
  fRotMat_BL.Invert();// Lab to Box transformation
  // fOrigin will be set later PartnerPlanes (after Plane init)
  
  Double_t r_temp[3] = {fOrigin.X(), fOrigin.Y(), fOrigin.Z()};
  TMatrixD m_temp(3, 1, r_temp); 
  
  TMatrixD m_res(3, 1, r_temp);
  m_res.Mult(fRotMat_BL, m_temp);
  
  fOrigin.SetX(m_res(0, 0)-fDMag*sin(fThetaH)*1.0e3);
  fOrigin.SetY(m_res(1, 0));
  fOrigin.SetY(m_res(2, 0)+fDMag*cos(fThetaH)*1.0e3);
  
  return kOK;
}

//_____________________________________________________________________________
Int_t GEMTracker::NewTrackCalc( Int_t idx, THaTrack*, const TVector3& pos,
				const TVector3& dir, const FitRes_t& )
{
  // For every new track, convert track coordinates to the system
  // appropriate for SBS. This is a temporary solution; the right way to do
  // this is to have the SBS spectrometer use its own track class.

  assert( dynamic_cast<SBSSpec*>(GetApparatus()) );
  SBSSpec* sbs = static_cast<SBSSpec*>(GetApparatus());

  TClonesArray* trackInfo = sbs->GetTrackInfo();

#ifdef MCDATA
  SBSTrackInfo* theInfo =
#endif
  new( (*trackInfo)[idx] ) SBSTrackInfo( fSector, pos, dir);
  
  using TreeSearch::NumberOfSetBits;
  assert( static_cast<vrsiz_t>(idx) < f3dMatchBits.size() );
  theInfo->Set3DMatchBits( f3dMatchBits[idx] );
  theInfo->SetN3DMatch( NumberOfSetBits(theInfo->Get3DMatchBits()) );
  //cout << idx << " " << f3dMatchBits.size() << " " << theInfo->Get3DMatchBits() << " " << theInfo->GetN3DMatch() << endl;
#ifdef MCDATA
  using TreeSearch::NumberOfSetBits;
  assert( static_cast<vrsiz_t>(idx) < fMCHitBits.size() );
  theInfo->SetMCHitBits( fMCHitBits[idx] );
  theInfo->SetNMCHits( NumberOfSetBits(theInfo->GetMCHitBits()) );
#endif

  return 0;
}

#ifdef MCDATA
//_____________________________________________________________________________
Int_t GEMTracker::FitMCPoints( Podd::MCTrack* mctrk ) const
{
  // SBS version of FitMCPoints to a MC track. In addition to fitting
  // the track, also reconstruct the track to the target.
  // Currently assumes stright tracks.

  Int_t npt = TreeSearch::Tracker::FitMCPoints( mctrk );
  if( npt < 0 )
    return npt;

  // See SBSSpec::FindVertices. Find distance of closest approach of
  // the fitted MC track to the beam (assumed to be exactly along z)
  // At this point, the track parameters have been converted to the lab frame
  Double_t* coef = mctrk->fMCFitPar;
  TVector3 g0( coef[0], coef[2], fOrigin.Z() );
  TVector3 g1( coef[1], coef[3], 1.0  );
  g1 = g1.Unit();
  TVector3 h0( 0, 0, 0  );
  TVector3 h1( 0, 0, 1. );
  Double_t gh = g1*h1;
  Double_t denom = 1.-gh*gh;
  if( TMath::Abs(denom) > 1e-6 ) {
    TVector3 D0 = g0-h0;
    Double_t tc = -D0*(g1-h1*gh)/denom;
    //Double_t sc =  D0*(h1-g1*gh)/denom;
    TVector3 vertex = g0 + g1*tc;
    // Save the results as additional fit results with mctrk
    coef[6] = vertex.X();
    coef[7] = vertex.Y();
    coef[8] = vertex.Z();
  }
  return npt;
}
#endif // MCDATA

///////////////////////////////////////////////////////////////////////////////

} // end namespace SBS

ClassImp(SBS::GEMTracker)

