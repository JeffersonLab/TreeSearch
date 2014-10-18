//*-- Author :    Ole Hansen, Jefferson Lab   06-Jun-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
// Reconstruction class for horizontal drift chambers used in BigBite.       //
// Tracking is done using a tree search algorithm (DellOrso et al.,          //
// NIM A 287, 436 (1990)), in essence a recursive template matching method.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "MWDC.h"
#include "WirePlane.h"
#include "ProjectionLR.h"
#include "WireHit.h"

#include "THaDetMap.h"
#include "THaEvData.h"

#include "TString.h"
#include "TMath.h"

#include <algorithm>

using namespace std;

namespace TreeSearch {

typedef vector<Plane*>::size_type vwsiz_t;
typedef vector<Plane*>::iterator  vwiter_t;

#define ALL(c) (c).begin(), (c).end()

//_____________________________________________________________________________
MWDC::MWDC( const char* name, const char* desc, THaApparatus* app )
  : Tracker(name,desc,app)

{
  // Constructor

  fRefMap = new THaDetMap;
}

//_____________________________________________________________________________
MWDC::~MWDC()
{
  // Destructor

  delete fRefMap;
}

//_____________________________________________________________________________
Int_t MWDC::Decode( const THaEvData& evdata )
{
  // Decode all planes and fill hitpatterns per projection

  static const char* const here = "Decode";

  // Decode reference channels of the VME readout (if any)
  for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fRefMap->GetModule(imod);
    // By construction, this map has one channel per module
    Int_t chan = d->lo;
    UInt_t nhits = evdata.GetNumHits( d->crate, d->slot, chan );
    if( nhits > 0 ) {
      Int_t data = evdata.GetData( d->crate, d->slot, chan, nhits-1 );
      if( nhits > 1 ) {
	Warning( Here(here), "%d hits on reference channel %d module %d",
		 nhits, chan, imod );
      }
      fRefTime[imod] = d->resolution * data;
    } else {
      // TODO: At this point, one could look for backup reference channels
      // to recover the data. Left for later, if needed.
      Warning( Here(here), "No hits on reference channel %d module %d.",
	       chan, imod );
      fRefTime[imod] = kBig;
    }
  }   // modules

  // Standard decoding: decode all projections and fill hitpatterns
  Tracker::Decode( evdata );

  return 0;
}

//_____________________________________________________________________________
UInt_t MWDC::MatchRoadsImpl( vector<Rvec_t>& roads, UInt_t ncombos,
			     list<std::pair<Double_t,Rvec_t> >& combos_found,
			     Rset_t& unique_found )
{
  // Implementation of MatchRoads for MWDCs. Supports the generic and fast
  // geometric matching algorithms. Fast matching is selected automatically
  // in Init if the detector configuration is appropriate.

  vector<Rvec_t>::size_type nproj = roads.size();

  if( nproj < 3 )
    return 0;

  UInt_t nfound;
  if( TestBit(k3dFastMatch) ) {
    nfound = MatchRoadsFast3D( roads, ncombos, combos_found, unique_found );
  } else {
    nfound = MatchRoadsGeneric( roads, ncombos, combos_found, unique_found );
  }
  return nfound;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus MWDC::Init( const TDatime& date )
{
  // Initialize MWDC. Calls standard Init(), then initializes subdetectors.

  static const char* const here = "Init";

  fRefMap->Reset();
  fRefTime.clear();

  // Initialize ourselves. This runs common tasks such as instantiating and
  // initializing planes and projections and calls our ReadDatabase(),
  // DefineVariables() and PartnerPlanes().
  EStatus status = Tracker::Init(date);
  if( status )
    return fStatus = status;

  // Automatically create a separate detector map for the reference channels
  // used by the data channels. By construction, each "module" in this map
  // has exactly one channel.
  // This map is decoded in our Decode() for efficiency.
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    THaDetMap* planeMap = fPlanes[iplane]->GetDetMap();
    for( Int_t imod = 0; imod < planeMap->GetSize(); ++imod ) {
      THaDetMap::Module* d = planeMap->GetModule(imod);
      if( d->refchan < 0 ) continue;
      THaDetMap::Module* refmod = fRefMap->Find( d->crate, d->slot,
						 d->refchan );
      // To allow the data channels to look up their reference channel data
      // quickly, we save the index into the refmap in d->refindex.
      if( refmod ) {
	d->refindex = refmod->first;
      }	else {
	d->refindex = fRefMap->GetSize();
	fRefMap->AddModule( d->crate, d->slot, d->refchan,
			    d->refchan, d->refindex, d->model );
	refmod = fRefMap->GetModule(d->refindex);
	if( !refmod ) {
	  Error( Here(here), "Error when adding reference channel cr/sl/ch="
		 "%u/%u/%u. Call expert.", d->crate, d->slot, d->refchan );
	  return fStatus = kInitError;
	}
	// The reference channel module's model and resolution have already
	// been determined when the plane's detector map was filled. Also,
	// refindex was checked to be a vaid channel number for this model.
	refmod->SetResolution( d->resolution );
	refmod->MakeTDC();
      }
    }
  }
  // Check if any reference channels are also defined as data channels
  for( Int_t imod = 0; imod < fRefMap->GetSize(); ++imod ) {
    THaDetMap::Module* r = fRefMap->GetModule(imod);
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      THaDetMap* planeMap = fPlanes[iplane]->GetDetMap();
      if( planeMap->Find( r->crate, r->slot, r->lo ) ) {
	Error( Here(here), "Reference channel cr/sl/ch=%u/%u/%u also defined "
	       "as data channel in plane \"%s\". Fix database.",
	       r->crate, r->slot, r->lo, fPlanes[iplane]->GetName() );
	return fStatus = kInitError;
      }
    }
  }
  // Reference map set up correctly and without conflicts
  fRefTime.assign( fRefMap->GetSize(), kBig );

#ifdef MCDATA
  // Set up handler for MC data aware of LR hits
  if( TestBit(kMCdata) ) {
    delete fMCPointUpdater;
    fMCPointUpdater = new LRPointUpdater;
  }
#endif

  return fStatus = kOK;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus MWDC::PartnerPlanes()
{
  // Associate planes and partners for wire chambers. Partner planes are
  // adjacent planes with the same wire direction and usually staggered
  // wire positions.
  // Requires planes to be initialized, i.e. to know their names,
  // type and geometry.

  static const char* const here = "PartnerPlanes";

  // Associate planes and partners (either via partner name or, if not given,
  // identical names except trailing "p")
  if( !TestBit(kNoPartner) ) {
    fAllPartnered = true;
    for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
      Plane* thePlane = fPlanes[iplane];
      assert( thePlane->IsInit() );
      TString other( thePlane->GetPartnerName() );
      if( other.IsNull() ) {
	TString name( thePlane->GetName() );
	if( name.EndsWith("p") )
	  other = name.Chop();
      }
      if( other.IsNull() )
	continue;
      vwiter_t it = find_if( ALL(fPlanes), Plane::NameEquals(other) );
      if( it != fPlanes.end() ) {
	Plane* partner = *it;
	assert( partner->IsInit() );
	// Sanity checks
	if( thePlane == partner ) {
	  Error( Here(here), "Plane %s: cannot partner a plane with itself. "
		 "Fix database.", other.Data() );
	  return kInitError;
	}
	// Partner planes must be of the same type!
	if( thePlane->GetType() != partner->GetType() ) {
	  Error( Here(here), "Partner planes %s and %s have different types!"
		 " Fix database.", thePlane->GetName(), partner->GetName() );
	  return kInitError;
	}
	// Check for consistency
	TString ppname( partner->GetPartnerName() );
	if( !ppname.IsNull() and ppname != thePlane->GetName() ) {
	  Error( Here(here), "Inconsistent plane partnering. Partner(%s) "
		 "= %s, but partner(%s) = %s. Fix database.",
		 thePlane->GetName(), other.Data(),
		 partner->GetName(), ppname.Data() );
	  return kInitError;
	}
	// Nothing to do if partner already set (prevents duplicate printouts)
	if( thePlane->GetPartner() == partner ) {
	  assert( partner->GetPartner() == thePlane );
	  continue;
	}
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
  } else {
    fAllPartnered = false;
  }

  return kOK;
}

//_____________________________________________________________________________
Int_t MWDC::ReadDatabase( const TDatime& date )
{
  // Read MWDC database

  Int_t err = Tracker::ReadDatabase( date );
  fIsInit = kFALSE;
  if( err != kOK )
    return err;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Int_t time_cut = 1, pairs_only = 0, nopartner = 0;
  DBRequest request[] = {
    { "timecut",           &time_cut,          kInt,    0, 1 },
    { "pairsonly",         &pairs_only,        kInt,    0, 1 },
    { "nopartner",         &nopartner,         kInt,    0, 1 },
    { 0 }
  };
  err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return err;

  // Set analysis control flags
  SetBit( kDoTimeCut,     time_cut );
  SetBit( kPairsOnly,     pairs_only );
  SetBit( kNoPartner,     nopartner );

  fIsInit = kTRUE;
  return kOK;
}

//_____________________________________________________________________________
Plane* MWDC::MakePlane( const char* name, const char* description,
			THaDetectorBase* parent ) const
{
  // Create an object of the plane class used by this implementation

  return new WirePlane( name, description, parent );
}

//_____________________________________________________________________________
Projection* MWDC::MakeProjection( EProjType type, const char* name,
				  Double_t angle,
				  THaDetectorBase* parent ) const
{
  // Create an object of the projection class used by this implementation.
  // The type of projection determines basically the method with which the
  // Hitpattern will be filled.

  return new ProjectionLR( type, name, angle, parent );
}

//_____________________________________________________________________________
UInt_t MWDC::GetCrateMapDBcols() const
{
  // Return number of columns for the detector crate map in the database.
  // 6 columns means that a 6-th column is present, indicating the resolution
  // of the module, typically a TDC

  return 6;
}

#ifdef MCDATA
//_____________________________________________________________________________
void MWDC::LRPointUpdater::UpdateHit( Podd::MCTrackPoint* pt, Hit* hit,
				      Double_t x ) const
{
  // Updater for MCTrackPoint hit residual aware of LR hits

  assert( pt );

  if( !hit )
    return MCPointUpdater::UpdateHit(pt,hit,x);
  pt->fStatus = 0;
  assert( dynamic_cast<WireHit*>(hit) );
  WireHit* wirehit = static_cast<WireHit*>(hit);
  Double_t dx1 = wirehit->GetPosR() - x;
  Double_t dx2 = wirehit->GetPosL() - x;
  pt->fHitResid = TMath::Min(dx1,dx2);
};

#endif // MCDATA

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

ClassImp(TreeSearch::MWDC)
