//*-- Author :    Ole Hansen, Jefferson Lab   11-Jan-2010

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::GEMTracker                                                    //
//                                                                           //
// Reconstruction class for a set of GEM tracking chambers in the HRS.       //
// Adapted from BigBite MWDC tracking code.                                  //
// Tracking is done using a tree search algorithm (DellOrso et al.,          //
// NIM A 287, 436 (1990)), in essence a recursive template matching method.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GEMTracker.h"
#include "GEMPlane.h"
#include "Projection.h"
#include "Hit.h"
#include "Road.h"
#include "Helper.h"

#include "THaDetMap.h"
#include "THaEvData.h"
#include "THaTrack.h"

#include "TString.h"
#include "TMath.h"
#include "THashTable.h"
#include "TVector2.h"
#include "TBits.h"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <string>
#include <stdexcept>

using namespace std;
typedef string::size_type ssiz_t;

namespace TreeSearch {

typedef vector<Plane*>::size_type vrsiz_t;
typedef vector<Plane*>::iterator  vriter_t;
typedef vector<Projection*>::size_type vpsiz_t;
typedef vector<Projection*>::iterator  vpiter_t;
typedef vector<vector<Int_t> >::iterator vviter_t;

#define ALL(c) (c).begin(), (c).end()

//_____________________________________________________________________________
GEMTracker::GEMTracker( const char* name, const char* desc, THaApparatus* app )
  : Tracker(name,desc,app),
    fMaxCorrMismatches(0), fMaxCorrNsigma(1.0)
{ 
  // Constructor - nothing special for standard GEM trackers
}

//_____________________________________________________________________________
GEMTracker::~GEMTracker()
{
  // Destructor - nothing special for standard GEM trackers
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus GEMTracker::Init( const TDatime& date )
{
  // Initialize GEMTracker. Calls standard Init(), then initializes subdetectors.

  static const char* const here = "Init";

  // The "cratemap" is only needed during Init of GEM and Plane
  assert( fCrateMap == 0 );
  fCrateMap = new THashTable(100);
  fCrateMap->SetOwner();

  // Initialize ourselves. This calls our ReadDatabase() and DefineVariables()
  EStatus status = THaTrackingDetector::Init(date);

  // Initialize the readout planes
  if( !status ) {
    for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      status = fPlanes[iplane]->Init(date);
      if( status )
	break;
    }
  }
  delete fCrateMap; fCrateMap = 0;
  if( status )
    return fStatus = status;

  // Sort planes by increasing z-position
  sort( ALL(fPlanes), Plane::ZIsLess() );

  // Associate planes and partners
  bool all_partnered = true;
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    Plane* thePlane = fPlanes[iplane];
    TString other( thePlane->GetPartnerName() );
    if( other.IsNull() ) {
      thePlane->SetPartner( 0 );
      all_partnered = false;
      continue;
    }
    vriter_t it = find_if( ALL(fPlanes), Plane::NameEquals(other) );
    if( it != fPlanes.end() ) {
      Plane* partner = *it;
      // Sanity checks
      if( thePlane == partner ) {
	Error( Here(here), "Plane %s: cannot partner a plane with itself. "
	       "Fix database.", other.Data() );
	return fStatus = kInitError;
      }
      // Partner planes must be of different types (multi-dimensional readout)
      if( thePlane->GetType() == partner->GetType() ) {
	Error( Here(here), "Partner planes %s and %s have the same type!"
	       " Fix database.", thePlane->GetName(), partner->GetName() );
	return fStatus = kInitError;
      }
      // 2D readouts must have essentially the same z-position
      if( TMath::Abs( thePlane->GetZ() - partner->GetZ() ) > 1e-3 ) {
	Error( Here(here), "Partner planes %s and %s must have the same "
	       "z-position within 1 mm. Fix database.", 
	       thePlane->GetName(), partner->GetName() );
	return fStatus = kInitError;
      }
      // Check for consistency
      TString ppname( partner->GetPartnerName() );
      if( ppname.IsNull() or ppname != thePlane->GetName() ) {
	Error( Here(here), "Inconsistent plane partnering. Partner(%s) "
	       "= %s, but partner(%s) = %s. Fix database.",
	       thePlane->GetName(), other.Data(), 
	       partner->GetName(), ppname.IsNull() ? "(none)":ppname.Data() );
	return fStatus = kInitError;
      }
      // Nothing to do if partner already set (prevents duplicate printouts)
      if( thePlane->GetPartner() == partner ) {
	assert( partner->GetPartner() == thePlane );
	continue;
      }
      // Mutually associate thePlane and partner
      if( fDebug > 0 )
	Info( Here(here), "Partnering plane %s with %s",
	      thePlane->GetName(), partner->GetName() );

      partner->SetPartner( thePlane );

    } else {
      Error( Here(here), "Partner plane %s of %s is not defined!"
	     " Fix database.", other.Data(), thePlane->GetName() );
      return fStatus = kInitError;
    }
  }

  // Set up the projections based on which plane types are defined
  assert( fProj.empty() );
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    Plane* thePlane = fPlanes[iplane];
    assert( !thePlane->GetProjection() );
    EProjType type = thePlane->GetType();
    vpiter_t it = find_if( ALL(fProj), Projection::TypeEquals(type) );
    Projection* proj;
    if( it != fProj.end() ) {
      proj = *it;
    } else {
      try {
	proj = new Projection( type,
			       kProjParam[type].name, 
			       kProjParam[type].angle*TMath::DegToRad(),
			       this 
			       );
      }
      catch( bad_alloc ) { proj = 0; }
      if( !proj or proj->IsZombie() ) {
	// Urgh. Something very bad is going on
	Error( Here(here), "Error creating projection %s. Call expert.", 
	       kProjParam[type].name );
	return fStatus = kInitError;
      }
      fProj.push_back( proj );
    }
    proj->AddPlane( thePlane );
  }

  // Sort projections by ascending EProjType
  sort( ALL(fProj), Projection::ByType() );

  // If exactly two projections are defined, switch on the 3D amplitude
  // correlation matching mode
  if( fProj.size() == 2 ) {
    SetBit( k3dCorrAmpl );
    if( fDebug > 0 )
      Info( Here(here), "2 projections defined: Enabling amplitude "
	    "correlation matching.");

    if( not all_partnered ) {
      Error( Here(here), "Amplitude correlation mode enabled, but one or "
	     "more readout planes are 1D (=no partner). Fix database." );
      return fStatus = kInitError;
    }
  }

  // Initialize the projections. This will read the database and set
  // the projections' angle and maxslope, which we need in the following
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    status = (*it)->Init(date);
    if( status )
      return fStatus = status;
  }

  // Sanity checks of U and V angles which the projections just read via Init
  vpiter_t iu = find_if( ALL(fProj), Projection::TypeEquals(kUPlane) );
  vpiter_t iv = find_if( ALL(fProj), Projection::TypeEquals(kVPlane) );
  if( iu != fProj.end() and iv != fProj.end() ) {
    // Both u and v planes are defined
    //TODO: isn't this an unnecessary limitation?
    //TODO: instead, check if any projections have (nearly) identical angles
    Double_t u_angle = (*iu)->GetAngle()*TMath::RadToDeg();
    Double_t v_angle = (*iv)->GetAngle()*TMath::RadToDeg();
    Int_t qu = TMath::FloorNint( u_angle/90.0 );
    Int_t qv = TMath::FloorNint( v_angle/90.0 );
    if( (qu&1) == (qv&1) ) {
      Error( Here(here), "Plane misconfiguration: uangle (%6.2lf) and vangle "
	     "(%6.2lf) are in equivalent quadrants. Fix database.",
	     u_angle, v_angle );
      return fStatus = kInitError;
    }
    //TODO: put this check into Projection::Init?
    //FIXME: this does not require both u and v
    Double_t du = u_angle - 90.0*qu;
    Double_t dv = v_angle - 90.0*qv;
    if( TMath::Abs(TMath::Abs(du)-45.0) > 44.0 or
	TMath::Abs(TMath::Abs(dv)-45.0) > 44.0 ) {
      Error( Here(here), "uangle (%6.2lf) or vangle (%6.2lf) too close "
	     "to 0 or 90 degrees. Fix database.", u_angle, v_angle );
      return fStatus = kInitError;
    }

    // Check if we can use the simplified 3D matching algorithm
    //FIXME: fast_3d should work for x,y,u(45) and similar, too
    //TODO: generalize to 4 projections with symmetry
    if( fProj.size() == 3 ) {
      assert( !TestBit(k3dCorrAmpl) );
      // This assumes that u and v are the first two defined projections
      assert( kUPlane < kVPlane && kVPlane < 2 );
      // The abs(angle) of the two rotated planes must be (nearly) the same
      Double_t uang = TMath::Abs( TVector2::Phi_mpi_pi( (*iu)->GetAngle() ));
      if( (TMath::Abs( TVector2::Phi_mpi_pi( (*iv)->GetAngle() ))-uang )
	  <  0.5*TMath::DegToRad() ) {
	SetBit(k3dFastMatch);
	Double_t tan = TMath::Tan( 0.5*TMath::Pi()-uang );
	// The scale factor converts the fast_3d matchvalue to the one computed
	// by the general algorithm
	f3dMatchvalScalefact = 2.0 * (1.0/3.0 + tan*tan );
	// Avoid scaling for every event
	f3dMatchCut /= f3dMatchvalScalefact;
      }
    }
  }

  // Determine width and maxslope of the projections
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    Projection* theProj = *it;
    EProjType type = theProj->GetType();
    Double_t width = 0.0;
    //TODO: can this be part of Projection::Init?
    for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      Plane* thePlane = fPlanes[iplane];
      //FIXME: loop over the projection's planes only -> do in Projection::Init
      if( thePlane->GetType() == type ) {
	// Determine the "width" of this projection plane (=width along the
	// projection coordinate).
	// The idea is that all possible hit positions in all planes of a 
	// given projection must fall within the range [-W/2,W/2]. It is normal
	// that some planes cover less than this width (e.g. because they are 
	// smaller or offset) - it is meant to be the enclosing range for all
	// planes. In this way, the tree search can use one fixed bin width.
	// The total width found here divided by the number of bins used in
	// tree search is the resolution of the roads.
	// Of course, this will become inefficient for grotesque geometries.
	Double_t s = thePlane->GetStripStart();
	Double_t d = thePlane->GetStripPitch();
	Double_t n = static_cast<Double_t>( thePlane->GetNelem() );
	Double_t lo = s - 0.5*d;
	Double_t hi = s + (n-0.5)*d;
	Double_t w = max( TMath::Abs(hi), TMath::Abs(lo) );
	if( w > width )
	  width = w;
      }
    }

    //TODO: make this work if projection does not have enough planes (>=2)
    // Set width of this projection
    width *= 2.0;
    if( width > 0.01 )
      theProj->SetWidth( width );
    else {
      Error( Here(here), "Error calculating width of projection \"%s\". "
	     "Strip pitch too small? Fix database.", theProj->GetName() );
      return fStatus = kInitError;
    }
    // maxslope is the maximum expected track slope in the projection.
    // width/depth is the maximum geometrically possible slope. It may be
    // further limited by the trigger acceptance, optics, etc.
    Double_t dz = TMath::Abs(theProj->GetZsize());
    if( dz > 0.01 ) {
      Double_t maxslope = width/dz;
      if( theProj->GetMaxSlope() < 0.01 ) {  // Consider unset
	theProj->SetMaxSlope( maxslope );
      } else if( theProj->GetMaxSlope() > maxslope ) {
	if( fDebug > 0 ) {
	  Warning( Here(here), "For projection \"%s\", maxslope from "
		   "database = %lf exceeds geometric maximum = %lf. "
		   "Using smaller value.",
		   theProj->GetName(), theProj->GetMaxSlope(), maxslope );
	}
	theProj->SetMaxSlope( maxslope );
      }
    } else {
      Error( Here(here), "Error calculating geometric maxslope for plane "
	     "type \"%s\". z-range of planes too small. Fix database.",
	     theProj->GetName() );
	return fStatus = kInitError;
    }

    // Now that the projection's list of planes, width, and maxslope is known,
    // do the level-2 initialization of the projections - load the pattern
    // database and initialize the hitpattern
    //FIXME: no longer need 2-part Init
    status = theProj->InitLevel2(date);
    if( status )
      return fStatus = status;

  }

  return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t GEMTracker::ReadDatabase( const TDatime& date )
{
  // Read GEM database

  static const char* const here = "ReadDatabase";

  fIsInit = kFALSE;
  // Delete existing configuration (in case we are re-initializing)
  DeleteContainer( fProj );
  DeleteContainer( fPlanes );
  fCalibPlanes.clear();

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Read fOrigin (detector position) and fSize. fOrigin is the position of
  // the GEM detector relative to some superior coordinate system
  // (typically the spectrometer detector stack reference frame). 
  // fOrigin will be added to all tracks generated; if fOrigin.Z() is not
  // zero, tracks will be projected into the z=0 plane.
  fOrigin.SetXYZ(0,0,0);
  Int_t err = ReadGeometry( file, date );
  if( err )
    return err;

  Int_t do_adc_cut = 1;
  f3dMatchCut = 1e-4;
  Int_t event_display = 0, disable_tracking = 0, disable_finetrack = 0;
  Int_t maxmiss = -1, maxthreads = -1, ampcorr_maxmiss = -1;
  fMaxCorrNsigma = 1.0;
  Double_t conf_level = 1e-9;
  ResetBit( k3dCorrAmpl );  // Set in Init()
  ResetBit( k3dFastMatch );
  DBRequest request[] = {
    { "do_adc_cut",        &do_adc_cut,        kInt,    0, 1 },
    { "3d_matchcut",       &f3dMatchCut,       kDouble, 0, 1 },
    { "event_display",     &event_display,     kInt,    0, 1 },
    { "disable_tracking",  &disable_tracking,  kInt,    0, 1 },
    { "disable_finetrack", &disable_finetrack, kInt,    0, 1 },
    { "3d_maxmiss",        &maxmiss,           kInt,    0, 1 },
    { "3d_chi2_conflevel", &conf_level,        kDouble, 0, 1 },
    { "3d_ampcorr_maxmiss",&ampcorr_maxmiss,   kInt,    0, 1 },
    { "3d_ampcorr_nsigma", &fMaxCorrNsigma,    kDouble, 0, 1 },
    { "maxthreads",        &maxthreads,        kInt,    0, 1 },
    { 0 }
  };

  Int_t status = kInitError;
  err = LoadDB( file, date, request, fPrefix );
  fclose(file);

  // Set analysis control flags
  SetBit( kDoADCCut,      do_adc_cut );
  SetBit( kEventDisplay,  event_display );
  SetBit( kDoCoarse,      !disable_tracking );
  SetBit( kDoFine,        !(disable_tracking or disable_finetrack) );
  bool doing_tracking = TestBit(kDoCoarse);

  vector<string> planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Set \"planeconfig\" in database." );
    return kInitError;
  }
  vector<string> calibplanes = vsplit(calibconfig);

  // Set up the readout planes
  for( ssiz_t i=0; i<planes.size(); ++i ) {
    assert( !planes[i].empty() );
    const char* name = planes[i].c_str();
    vriter_t it = find_if( ALL(fPlanes), Plane::NameEquals( name ) );
    if( it != fPlanes.end() ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", name );
      return kInitError;
    }
    Plane* newplane = 0;
    try { newplane = new Plane( name, name, this ); }
    catch( bad_alloc ) { newplane = 0; }
    if( !newplane or newplane->IsZombie() ) {
      // Urgh. Something is very bad
      Error( Here(here), "Error creating readout plane %s. Call expert.",
	     name );
      return kInitError;
    }
    fPlanes.push_back( newplane );
    newplane->SetDebug( fDebug );
    newplane->SetDefinedNum( i );

    // Set calibration mode if requested for any planes
    if( !calibplanes.empty() ) {
      vector<string>::iterator its =
	find_if( ALL(calibplanes), bind2nd(equal_to<string>(), planes[i]) );
      if( its != calibplanes.end() ) {
	Info( Here(here), "Plane %s in calibration mode", name );
	newplane->EnableCalibration();
	fCalibPlanes.push_back( newplane );
	calibplanes.erase( its );
      }
    }
  }

  // Warn if any requested calibration plane(s) do not exist
  if( !calibplanes.empty() ) {
    string s("plane");
    if( calibplanes.size() > 1 )
      s.append("s");
    for( vector<string>::size_type i = 0; i < calibplanes.size(); ++i ) {
      s.append(" ");
      s.append(calibplanes[i]);
    }
    Warning( Here(here), "Requested calibration for undefined plane %s. "
	     "Error in database?", s.c_str() );
  }

  UInt_t nplanes = fPlanes.size();
  if( nplanes < 5 ) {
    Error( Here(here), "Insufficient number of planes = %u. Need at least 5. "
	   "Fix database.", nplanes );
    return kInitError;
  }
  if( fDebug > 0 )
    Info( Here(here), "Loaded %u planes", nplanes );

  if( doing_tracking ) {
    // Convert maximum number of missing hits to ndof of fits
    fMinNdof = ( maxmiss >= 0 ) ? nplanes - 4 - maxmiss : 1;

    // Set parameter for maximum number of amplitude correlation mismatches
    fMaxCorrMismatches = ( ampcorr_maxmiss >= 0 ) ? ampcorr_maxmiss : 0;

    // Determine Chi2 confidence interval limits for the selected CL and the
    // possible degrees of freedom of the 3D track fit
    if( conf_level < 0.0 or conf_level > 1.0 ) {
      Error( Here(here), "Illegal fit confidence level = %lf. "
	     "Must be 0-1. Fix database.", conf_level );
      return kInitError;
    }
    fChisqLimits.clear();
    fChisqLimits.resize( nplanes-3, make_pair<Double_t,Double_t>(0,0) );
    for( vec_pdbl_t::size_type dof = fMinNdof; dof < fChisqLimits.size();
	 ++dof ) {
      fChisqLimits[dof].first  = TMath::ChisquareQuantile( conf_level, dof );
      fChisqLimits[dof].second = 
	TMath::ChisquareQuantile( 1.0-conf_level, dof );
    }
  }
  
  fIsInit = kTRUE;
  return kOK;
}

//_____________________________________________________________________________
Plane* GEMTracker::MakePlane( const char* name, const char* description, 
			      THaDetectorBase* parent ) const
{
  // Create an object of the plane class used by this implementation

  return new GEMPlane( name, description, parent );
}

//_____________________________________________________________________________
Projection* GEMTracker::MakeProjection( EProjType type, const char* name,
					Double_t angle, 
					THaDetectorBase* parent ) const
{
  // Create an object of the projection class used by this implementation.
  // The type of projection determines basically the method with which the
  // Hitpattern will be filled.

  return new Projection( type, name, angle, parent );
}

//_____________________________________________________________________________
UInt_t GEMTracker::GetCrateMapDBcols() const
{
  // Return number of columns for the detector crate map in the database.
  // 5 columns means that the no resolution value is needed to analyze the
  // data, as is typical for ADCs

  return 5;
}

//_____________________________________________________________________________
UInt_t GEMTracker::MatchRoadsCorrAmpl( vector<Rvec_t>& roads,
		       UInt_t /* ncombos */,
		       list< pair<Double_t,Rvec_t> >& combos_found,
		       Rset_t& unique_found )
{
  // Matching of Roads via amplitude correlation. Obviously, this requires
  // detectors that measure some sort of amplitude (e.g. energy deposited)
  // for hits in some sort of shared readout plane.
  // Currently requires exactly two projections with identical number of planes.

  vector<Rvec_t>::size_type nproj = roads.size();
  assert( nproj == 2 );
  assert( not (roads[0].empty() or roads[1].empty()) );

#ifdef VERBOSE
  cout << "amplitude correlation algo):";
#endif

  // Vector holding a combination of roads to test. One road from each
  // projection
  Rvec_t selected( nproj, 0 );

  UInt_t nfound = 0;

  const Rpvec_t& xplanes =
    roads[0].front()->GetProjection()->GetListOfPlanes();
  const Rpvec_t& yplanes =
    roads[1].front()->GetProjection()->GetListOfPlanes();

  assert( not xplanes.empty() and (xplanes.size() == yplanes.size()) );
  assert( xplanes.front()->GetType() != yplanes.front()->GetType() );

  UInt_t nplanes = xplanes.size();
  TBits xybits(nplanes), ybits(nplanes); 
  // Look at all possible combinations of x-roads and y-roads.
  // Try to match them via the ADC amplitudes of the hits in the shared
  // readout planes. Amplitudes should correlate well in the absence
  // of pileup. 
  for( Rvec_t::iterator itx = roads[0].begin(); itx != roads[0].end();
       ++itx ) {
    Road* xroad = *itx;
    const Road::Pvec_t& xpoints = xroad->GetPoints();
    UInt_t xpat = xroad->GetPlanePattern();
    for( Rvec_t::iterator ity = roads[1].begin(); ity != roads[1].end();
	 ++ity ) {
      Road* yroad = *ity;

      // xpat and ypat are bitpatterns of the plane numbers that have hits.
      // The AND of these patters is the pattern of planes where both read-
      // out directions have hits (matching or not). Check here if there
      // are enough such common active planes, else all following work can
      // be skipped.
      UInt_t ypat = yroad->GetPlanePattern();
      xybits.Set( nplanes, &xpat );
      ybits.Set( nplanes, &ypat );
      xybits &= ybits;
      assert( xybits.CountBits() <= nplanes );
      if( xybits.CountBits() + fMaxCorrMismatches < nplanes )
	continue;

      // For all points (=hits that yield the best fit) of this xroad, 
      // get the yroad point in the same plane (if available), then check
      // if their ADC amplitudes approximately match (within a hard cut)
      const Road::Pvec_t& ypoints = yroad->GetPoints();
      Road::Pvec_t::const_iterator ityp = ypoints.begin();
      UInt_t nmatches = 0;
      Double_t matchval = 0.0;
      for( Road::Pvec_t::const_iterator itxp = xpoints.begin();
	   itxp != xpoints.end(); ++itxp ) {
	const Road::Point* xp = *itxp;
	assert( xp and xp->hit );
	const Plane* xplane = xp->hit->GetPlane();
	assert( xplane );
	UInt_t xnum = xplane->GetPlaneNum();
	assert( xnum < xplanes.size() );
	// No hit in the other readout direction of this plane?
	if( not ybits.TestBitNumber(xnum) )
	  continue;
	// Move y-iterator forward until it gets to the same plane
	while( ityp != ypoints.end() and
	       (*ityp)->hit->GetPlaneNum() != xnum ) {
	  ++ityp;
	}
	// Must find a corresponding hit here, else bug in ybits test
	assert( ityp != ypoints.end() );
	if( ityp == ypoints.end() ) break;
	const Road::Point* yp = *ityp;
	assert( yp and yp->hit );
	const Plane* yplane = yp->hit->GetPlane();
	assert( yplane );
	assert( yplane->GetPartner() == xplane );

	// Get ratio of hit amplitudes
	Double_t xampl = xp->hit->GetADCsum();
	Double_t yampl = yp->hit->GetADCsum();
	assert( xampl > 0.0 and yampl > 0.0 ); // ensured in Decoder
	if( xampl < 1.0 or yampl < 1.0 )
	  continue;
	Double_t ratio = xampl/yampl;
	Double_t asym = (xampl - yampl)/(xampl + yampl);
	// Compute the cutoff. In general, the sigma of the distribution can
	// be amplitude-dependent, so we get it via a calibration function
	// that is a property of each readout plane
	// 	  Double_t xsigma = xplane->GetAmplSigma( xampl );
	// 	  Double_t ysigma = yplane->GetAmplSigma( yampl );
	// Apply overall scale factor ("number of sigmas") from the database
	// 	  Double_t cutoff = fMaxCorrNsigma / yampl *
	// 	    TMath::Sqrt( xsigma*xsigma + ysigma*ysigma*ratio*ratio );
	Double_t cutoff = fMaxCorrNsigma;
#ifdef VERBOSE
	if( fDebug > 3 ) {
	  cout << xplane->GetName() << yplane->GetName()
	       << " ampl = (" << xampl << ", " << yampl << ")"
	       << ",\tratio = " << ratio
	       << ", asym = " << asym << endl;
	}
#endif
	// Count readout planes whose x/y-hit amplitudes match
	if( TMath::Abs(asym) < cutoff ) {
	  matchval += TMath::Abs(asym);  // not really used later (yet)
	  ++nmatches;
	}
      } //xpoints
	  
#ifdef VERBOSE
      if( fDebug > 3 ) {
	cout <<   "nmatches = " << nmatches 
	     << ", matchval = " << matchval << "   ";
      }
#endif
      // If enough planes have correlated x/y hit amplitudes, save this
      // x/y road pair as a possible 3D track candidate. The 3D track fit
      // will separate the wheat from the chaff later
      if( nmatches + fMaxCorrMismatches >= nplanes ) {
	selected[0] = xroad;
	selected[1] = yroad;
	Add3dMatch( selected, matchval, combos_found, unique_found );
	++nfound;
      }
#ifdef VERBOSE
      else if( fDebug > 3 ) { cout << endl; }
#endif
    } //yroads
  }   //xroads

  return nfound;
}

//_____________________________________________________________________________
UInt_t GEMTracker::MatchRoadsImpl( vector<Rvec_t>& roads, UInt_t ncombos,
				   list<std::pair<Double_t,Rvec_t> >& combos_found,
				   Rset_t& unique_found )
{
  // Implementation of MatchRoads for GEM trackers.

  vector<Rvec_t>::size_type nproj = roads.size();

  bool fast_3d = TestBit(k3dFastMatch);
  bool correlate_amplitudes = TestBit(k3dCorrAmpl);

  // Limitation: amplitude correlation algo only implemented for 2 projections
  assert( nproj == 2 or not correlate_amplitudes );

  // Fast 3D matching and amplitude correlation are mutually exclusive
  assert( not (fast_3d and correlate_amplitudes) );

  // if( nproj == 2 and not correlate_amplitudes )
  //   return 0;

  UInt_t nfound;

  if( correlate_amplitudes ) {
    nfound = MatchRoadsCorrAmpl( roads, ncombos, combos_found, unique_found );
  } else if( fast_3d ) {
    nfound = MatchRoadsFast3D( roads, ncombos, combos_found, unique_found );
  } else {
    nfound = MatchRoadsGeneric( roads, ncombos, combos_found, unique_found );
  }

  return nfound;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

ClassImp(TreeSearch::GEMTracker)

