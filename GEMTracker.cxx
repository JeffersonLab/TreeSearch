//*-- Author :    Ole Hansen, Jefferson Lab   11-Jan-2010

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::GEMTracker                                                    //
//                                                                           //
// Reconstruction class for a set of GEM tracking chambers (GEMPlanes).      //
// Tracking is done using a tree search algorithm (DellOrso et al.,          //
// NIM A 287, 436 (1990)), in essence a recursive template matching method.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GEMTracker.h"
#include "GEMPlane.h"
#include "Projection.h"
#include "GEMHit.h"
#include "Road.h"

#include "TString.h"
#include "TMath.h"
#include "TBits.h"

#include <iostream>
#include <algorithm>

using namespace std;

namespace TreeSearch {

typedef vector<Plane*>::size_type vrsiz_t;
typedef vector<Plane*>::iterator  vriter_t;

#define ALL(c) (c).begin(), (c).end()

//_____________________________________________________________________________
GEMTracker::GEMTracker( const char* name, const char* desc, THaApparatus* app )
  : Tracker(name,desc,app), fMaxCorrMismatches(0), fMaxCorrNsigma(1.0)
{ 
  // Constructor - nothing special for standard GEM trackers

  // GEMs support matching 2 projections only (via amplitude correlations
  // in MatchRoadsCorrAmpl). Override the default (3) from Tracker base class.
  fMinReqProj = 2;
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

  // Initialize ourselves. This runs common tasks such as instantiating and
  // initializing planes and projections and calls our ReadDatabase(),
  // DefineVariables() and PartnerPlanes().
  EStatus status = Tracker::Init(date);
  if( status )
    return fStatus = status;

  // If exactly two projections are defined, switch on the 3D amplitude
  // correlation matching mode
  if( fProj.size() == 2 ) {
    assert( !TestBit(k3dFastMatch) );
    SetBit( k3dCorrAmpl );
    if( fDebug > 0 )
      Info( Here(here), "2 projections defined: Enabling amplitude "
	    "correlation matching.");

    if( not fAllPartnered ) {
      Error( Here(here), "Amplitude correlation mode enabled, but one or "
	     "more readout planes are 1D (=no partner). Fix database." );
      return fStatus = kInitError;
    }
  }

  return fStatus = kOK;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus GEMTracker::PartnerPlanes()
{
  // Associate planes and partners for GEM trackers. Partner planes are those
  // that share the same readout plane (2D readouts). Requires planes to be
  // initialized, i.e. to know their names, type and geometry.

  static const char* const here = "PartnerPlanes";

  fAllPartnered = true;
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    Plane* thePlane = fPlanes[iplane];
    assert( thePlane->IsInit() );
    TString other( thePlane->GetPartnerName() );
    if( other.IsNull() ) {
      thePlane->SetPartner( 0 );
      fAllPartnered = false;
      continue;
    }
    vriter_t it = find_if( ALL(fPlanes), Plane::NameEquals(other) );
    if( it != fPlanes.end() ) {
      Plane* partner = *it;
      assert( partner->IsInit() );
      // Sanity checks
      if( thePlane == partner ) {
	Error( Here(here), "Plane %s: cannot partner a plane with itself. "
	       "Fix database.", other.Data() );
	return kInitError;
      }
      // Partner planes must be of different types (multi-dimensional readout)
      if( thePlane->GetType() == partner->GetType() ) {
	Error( Here(here), "Partner planes %s and %s have the same type!"
	       " Fix database.", thePlane->GetName(), partner->GetName() );
	return kInitError;
      }
      // 2D readouts must have essentially the same z-position
      if( TMath::Abs( thePlane->GetZ() - partner->GetZ() ) > 1e-3 ) {
	Error( Here(here), "Partner planes %s and %s must have the same "
	       "z-position within 1 mm. Fix database.", 
	       thePlane->GetName(), partner->GetName() );
	return kInitError;
      }
      // Check for consistency
      TString ppname( partner->GetPartnerName() );
      if( ppname.IsNull() or ppname != thePlane->GetName() ) {
	Error( Here(here), "Inconsistent plane partnering. Partner(%s) "
	       "= %s, but partner(%s) = %s. Fix database.",
	       thePlane->GetName(), other.Data(), 
	       partner->GetName(), ppname.IsNull() ? "(none)":ppname.Data() );
	return kInitError;
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
      return kInitError;
    }
  }
  return kOK;
}

//_____________________________________________________________________________
Int_t GEMTracker::ReadDatabase( const TDatime& date )
{
  // Read GEM database

  // Read common database info
  Int_t err = Tracker::ReadDatabase( date );
  fIsInit = kFALSE;
  if( err != kOK )
    return err;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // These bits will be set later in Init()
  ResetBit( k3dCorrAmpl );
  ResetBit( k3dFastMatch );

  Int_t do_adc_cut = 1;
  Int_t ampcorr_maxmiss = -1;
  fMaxCorrNsigma = 1.0;
  DBRequest request[] = {
    { "3d_ampcorr_maxmiss", &ampcorr_maxmiss, kInt,    0, 1 },
    { "3d_ampcorr_nsigma",  &fMaxCorrNsigma,  kDouble, 0, 1 },
    { 0 }
  };

  err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return err;

  // Set analysis control flags
  SetBit( kDoADCCut, do_adc_cut );

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
	assert( dynamic_cast<GEMHit*>(xp->hit) );
	assert( dynamic_cast<GEMHit*>(yp->hit) );
	Double_t xampl = static_cast<GEMHit*>(xp->hit)->GetADCsum();
	Double_t yampl = static_cast<GEMHit*>(yp->hit)->GetADCsum();
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

  assert( nproj >= 2 );

  bool fast_3d = TestBit(k3dFastMatch);
  bool correlate_amplitudes = TestBit(k3dCorrAmpl);

  if( nproj == 2 and not correlate_amplitudes )
    return 0;

  // Limitation: amplitude correlation algo only implemented for 2 projections
  assert( nproj == 2 or not correlate_amplitudes ); // else bug in Init

  // Fast 3D matching and amplitude correlation are mutually exclusive
  assert( not (fast_3d and correlate_amplitudes) );

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

