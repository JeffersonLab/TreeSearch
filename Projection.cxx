//*-- Author :    Ole Hansen, Jefferson Lab   16-Aug-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Projection.h"
#include "Hitpattern.h"
#include "Plane.h"
#include "THaDetectorBase.h"
#include "PatternTree.h"
#include "PatternGenerator.h"
#include "TreeWalk.h"
#include "Road.h"
#include "Helper.h"
#include "Hit.h"
#include "Tracker.h"   // for Tracker bits

#include "TMath.h"
#include "TString.h"
#include "TBits.h"
#include "TError.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <utility>
#ifdef TESTCODE
#include "TStopwatch.h"
#include <cstring>
#endif

using namespace std;

namespace TreeSearch {

#define ALL(c) (c).begin(), (c).end()

// Need at least 3 planes to do proper fits
static const UInt_t kMinFitPlanes = 3;

// Parameter for angle consistency check in SetAngle (rad)
static const Double_t kAngleTolerance = 1.0 * TMath::DegToRad();

//_____________________________________________________________________________
Projection::Projection( EProjType type, const char* name, Double_t angle,
			THaDetectorBase* parent )
  : THaAnalysisObject( name, name ), fType(type), fNlevels(0),
    fMaxSlope(0.0), fWidth(0.0), fDoElTrackSel(false), fDetector(parent), fPatternTree(0),
    fDummyPlanePattern(0), fFirstPlaneNum(kMaxUInt), fLastPlaneNum(0),
    fMinFitPlanes(kMinFitPlanes), fMaxMiss(0), fRequire1of2(false),
    fPlaneCombos(0), fAltPlaneCombos(0), fMaxPat(kMaxUInt),
    fFrontMaxBinDist(kMaxUInt), fBackMaxBinDist(kMaxUInt), fHitMaxDist(0),
    fConfLevel(1e-3), fHitpattern(0), fRoads(0), fNgoodRoads(0), fNMCRoads(0),
    fRoadCorners(0), fTrkStat(kTrackOK)
{
  // Constructor

  assert( name && parent );

  // angle is the default value. It can be overridden via a database entry.
  SetAngle( angle );

  fTitle.Append(" projection");
  fRoads = new TClonesArray("TreeSearch::Road", 3);
  R__ASSERT(fRoads);

#ifdef TESTCODE
  size_t nbytes = (char*)&t_track - (char*)&n_hits + sizeof(t_track);
  memset( &n_hits, 0, nbytes );
#endif

  ResetBit( kHaveDummies );
}

//_____________________________________________________________________________
Projection::~Projection()
{
  // Destructor

  if( fIsSetup )
    RemoveVariables();
  delete fRoads;
  delete fRoadCorners;
  delete fPatternTree;
  delete fHitpattern;
  if( fAltPlaneCombos != fPlaneCombos )
    delete fAltPlaneCombos;
  delete fPlaneCombos;
}

//_____________________________________________________________________________
void Projection::AddPlane( Plane* pl, Plane* partner )
{
  // Add plane pl (and optional partner plane) to this projection.
  // Sets plane numbers.

  assert( pl );
  assert( !pl->IsDummy() );

  pl->SetPlaneNum( fPlanes.size() );
  fPlanes.push_back( pl );
  // pl->SetProjection( this );
  if( partner ) {
    //    assert( partner->GetZ() > pl->GetZ() ); // planes must be ordered
    assert( partner != pl );
    assert( !partner->IsDummy() );
    partner->SetPlaneNum( fPlanes.size() );
    fPlanes.push_back( partner );
    // partner->SetProjection( this );
  }

  AddDummyPlane( pl, partner );
}

//_____________________________________________________________________________
void Projection::AddDummyPlane( Plane* pl, Plane* partner )
{
  // Add dummy plane pl (and optional partner plane) to this projection.
  // Sets plane numbers.

  assert( pl );

  UInt_t apn = fAllPlanes.size(), appn = 0;
  pl->SetAltPlaneNum( apn );
  fAllPlanes.push_back( pl );
  pl->SetProjection( this );
  if( partner ) {
    assert( partner != pl );
    assert( not (pl->IsDummy() xor partner->IsDummy()) );
    appn = fAllPlanes.size();
    partner->SetAltPlaneNum( appn );
    fAllPlanes.push_back( partner );
    partner->SetProjection( this );
  }
  if( pl->IsDummy() ) {
    SetBit( kHaveDummies );
    SETBIT( fDummyPlanePattern, apn );
    if( partner ) {
      SETBIT( fDummyPlanePattern, appn );
    }
  } else {
    fFirstPlaneNum = TMath::Min( fFirstPlaneNum, apn );
    fLastPlaneNum  = TMath::Max( fLastPlaneNum,  apn );
    if( partner ) {
      fFirstPlaneNum = TMath::Min( fFirstPlaneNum, appn );
      fLastPlaneNum  = TMath::Max( fLastPlaneNum,  appn );
    }
  }
}

//_____________________________________________________________________________
void Projection::Clear( Option_t* )
{
  // Clear event-by-event data

  if( fHitpattern )
    fHitpattern->Clear();

  fRoads->Delete();
  DeleteContainer( fPatternsFound );
  fNgoodRoads = 0;
  fNMCRoads = 0;
  fTrkStat = kTrackOK;

  if( TestBit(kEventDisplay) )
    fRoadCorners->Clear();

#ifdef TESTCODE
  size_t nbytes = (char*)&t_track - (char*)&n_hits + sizeof(t_track);
  memset( &n_hits, 0, nbytes );
#endif
}

//_____________________________________________________________________________
Int_t Projection::Decode( const THaEvData& evdata )
{
  // Decode all planes belonging to this projection

  Int_t sum = 0;
  bool err = false;
  for( vplsiz_t i = 0; i < GetNallPlanes(); ++i ) {
    Plane* pl = fAllPlanes[i];
    Int_t nhits = pl->Decode( evdata );
    //cout<<"   "<<pl->GetName()<<": "<<nhits<<endl;
    if( nhits < 0 ) {
      err = true;
      sum -= nhits;
    } else
      sum += nhits;
  }
  /*/Test code---------------
  bool checkflag = true;
  for( vplsiz_t i = 0; i < GetNallPlanes(); ++i ) {
    Plane* pl = fAllPlanes[i];
    Int_t nhits = pl->GetNhits();
    if(nhits!=1){checkflag = false;break;} 
    //cout<<"   "<<pl->GetName()<<": "<<nhits<<endl;
    TIterator* it = pl->GetHits()->MakeIterator();
    Hit* phit = 0;
    while( (phit = static_cast<Hit*>(it->Next())) ) {
     
    }
    delete it;
  }


  */

  if( err )
    return -sum;

  return sum;
}

//_____________________________________________________________________________
Double_t Projection::GetPlaneZ( UInt_t i ) const
{
  // Return z-position of i-th plane.

  assert( i<fPlanes.size() );
  return fPlanes[i]->GetZ();
}

//_____________________________________________________________________________
void Projection::Reset( Option_t* opt )
{
  // Reset parameters, delete dynamically allocated objects.
  // If opt="FULL", reset everything, including the list of planes.

  fIsInit = kFALSE;
  fMaxSlope = fWidth = 0.0;
  delete fHitpattern; fHitpattern = 0;
  delete fPatternTree; fPatternTree = 0;
  if( fAltPlaneCombos != fPlaneCombos ) {
    delete fAltPlaneCombos; fAltPlaneCombos = 0;
  }
  delete fPlaneCombos; fPlaneCombos = 0;
  delete fRoadCorners; fRoadCorners = 0;
  if( opt and *opt ) {
    TString s(opt);
    if( s.Contains( "FULL", TString::kIgnoreCase )) {
      fPlanes.clear();
      fAllPlanes.clear();
      ResetBit( kHaveDummies );
    }
  }
}

//_____________________________________________________________________________
void Projection::MakePlaneCombos( const vpl_t& planes, TBits*& combos ) const
{
  // Utility function to fill fPlaneCombos and fAltPlaneCombos according
  // to the configuration parameters fMaxMiss, fRequire1of2, and plane
  // required/dummy flags

  assert( combos == 0 );
  UInt_t nc = 1U << planes.size();
  combos = new TBits( nc );
  combos->SetBitNumber( nc-1 );  // Always allow full occupancy
  for( UInt_t i = 1; i <= fMaxMiss; ++i ) {
    UniqueCombo c( planes.size(), i );
    while( c ) {
      // Clear the bit numbers from this combination
      UInt_t bitval = nc-1;
      for( vplsiz_t j = 0; j < c().size(); ++j )
	CLRBIT( bitval, c()[j] );
      // Test if this pattern satisfies other constraints
      vplsiz_t k = 0;
      for( ; k < planes.size(); ++k ) {
	Plane* p = planes[k], *pp;
	// Disallow bit pattern if a required or dummy plane is missing
	if( !TESTBIT(bitval,k) and (p->IsRequired() or p->IsDummy()) )
	  break;
	// If requested, ensure that at least one plane of a plane pair is set
	if( fRequire1of2 and (pp = p->GetPartner()) and !TESTBIT(bitval,k) ) {
	  // Find the index of the partner plane. A bit tedious because calling
	  // GetPlaneNum or GetAltPlaneNum would require us to know which
	  // array of planes we have, and we want to keep things general
	  vplsiz_t j = 0;
	  for( ; j < planes.size(); ++j ) {
	    if( planes[j] == pp )
	      break; // break loop over j
	  }
	  assert( j < planes.size() ); // Partner plane must be in array
	  assert( j != k );            // Plane cannot be partner of itself
	  if( !TESTBIT(bitval,j) )
	    break;  // break loop over k
	}
      }
      assert( bitval < nc );
      if( k == planes.size() ) {
	// No objections were raised (in the loop over k) over this bitval
	combos->SetBitNumber( bitval );
      }
      ++c; // next UniqueCombo
    }
  }
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Projection::Init( const TDatime& date )
{
  // Initialize the Projection object. Called after basic Tracker
  // initialization.
  // Sets up event display support, then continues with standard
  // initialization.

  static const char* const here = "Init";

  Reset();

  // Set up the event display support if requested
  // Set local bit for efficiency
  SetBit( kEventDisplay, fDetector->TestBit(Tracker::kEventDisplay) );
  if( TestBit(kEventDisplay) ) {
    // Create fRoadCorners array needed for event display
    assert( fRoadCorners == 0 );
    try {
      fRoadCorners = new TClonesArray("TreeSearch::Road::Corners", 3);
    }
    catch( bad_alloc ) { fRoadCorners = 0; }
    if( !fRoadCorners or fRoadCorners->IsZombie() ) {
      Error( Here(here), "Out of memory creating event display for "
	     "projection \"%s\". Call expert.", GetName() );
      delete fRoadCorners; fRoadCorners = 0;
      return fStatus = kInitError;
    }
  }

#ifdef MCDATA
  SetBit( kMCdata, fDetector->TestBit(Tracker::kMCdata) );
#endif

  // Standard initialization. This calls ReadDatabase() and DefineVariables()
  THaAnalysisObject::EStatus status = THaAnalysisObject::Init(date);
  if( status )
    return fStatus = status;

  // Determine the width of this projection based on the plane geometry
  assert( fWidth == 0.0 );  // fWidth set to 0.0 in Reset()
  for( UInt_t i = 0; i < GetNallPlanes(); ++i ) {
    Plane* thePlane = fAllPlanes[i];
    assert( thePlane->IsInit() );
    if(thePlane->IsDummy())
      continue;
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
    Double_t s = thePlane->GetStart();
    Double_t d = thePlane->GetPitch();
    //Double_t n = static_cast<Double_t>( thePlane->GetNelem() );
    Double_t lo = s - 0.5*d;
    //Double_t hi = s + (n-0.5)*d;
    Double_t w = TMath::Abs(lo);

    // cout<<thePlane->IsDummy()<<" s: "<<s<<" d: "<<" w: "<<w<<endl;
    // Double_t w = max( TMath::Abs(hi), TMath::Abs(lo) );
    if( w > fWidth )
      fWidth = w;
    /*  cout<<"start: "<<s<<endl;
    cout<<"pitch: "<<d<<endl;
    cout<<"n: "<<n<<endl;
    cout<<"width: "<<w<<endl;getchar();*/
  }
  fWidth *= 2.0;
  if( fWidth < 0.01 ) {
    Error( Here(here), "Error calculating width of projection \"%s\". "
	   "Pitch too small? Fix database.", GetName() );
    return fStatus = kInitError;
  }

  // Determine maxslope, the maximum expected track slope in the projection.
  // width/depth is the maximum geometrically possible slope. It may be
  // further limited by the trigger acceptance, optics, etc.
  Double_t dz = GetZsize();
  assert( dz >= 0.0 );  // else planes not sorted
  if( dz >= 0.01 ) {
    Double_t maxslope = fWidth/dz;
    if( fMaxSlope < 0.01 ) {
      // maxslope not set in the database -> use the calculated value
      fMaxSlope = maxslope;
    } else if( fMaxSlope > maxslope ) {
      // Otherwise limit maxslope to the calculated value
      if( fDebug > 0 ) {
	Warning( Here(here), "For projection \"%s\", maxslope from "
		 "database = %lf exceeds geometric maximum = %lf. "
		 "Using smaller value.", GetName(), fMaxSlope, maxslope );
      }
      fMaxSlope = maxslope;
    }
  } else {
    Error( Here(here), "Error calculating geometric maxslope for plane "
	   "type \"%s\". z-range of planes = %lf too small. Fix database.",
	   GetName(), dz );
    return fStatus = kInitError;
  }

  // Now that the projection's list of planes, width, and maxslope are known,
  // do the level-2 initialization of the projections - load the pattern
  // database and initialize the hitpattern
//   return InitLevel2(date);

// }

// //_____________________________________________________________________________
// THaAnalysisObject::EStatus Projection::InitLevel2( const TDatime& )
// {
  // Level-2 initialization - load pattern database and initialize hitpattern.
  // Requires the planes to be fully initialized.

  // static const char* const here = "InitLevel2";

  // Require at least kMinFitPlanes planes
  // TODO: allow dummy/inactive projections?
  if( GetNplanes() < kMinFitPlanes ) {
    Error( Here(here), "Too few planes of type \"%s\" defined. "
          "Need >= %u, have %u. Fix database.",
	   GetName(), kMinFitPlanes, GetNplanes() );
    return kInitError;
  }

  // No need to set up pattern tree and hitpattern if no tracking requested
  bool doing_tracking = fDetector->TestBit(Tracker::kDoCoarse);
  if( doing_tracking ) {

    vector<Double_t> zpos;
    for( UInt_t i = 0; i < GetNallPlanes(); ++i )
      zpos.push_back( fAllPlanes[i]->GetZ() );
    TreeParam_t tp( fNlevels-1, fWidth, fMaxSlope, zpos );

    if( tp.Normalize() != 0 )
      return fStatus = kInitError;

    // Attempt to read the pattern database from file
    assert( fPatternTree == 0 );
    //TODO: Make the file name
    // const char* filename = "test.tree";
    // fPatternTree = PatternTree::Read( filename, tp );

    // If the tree cannot not be read (or the parameters mismatch), then
    // create it from scratch (takes a few seconds)
    // if( !fPatternTree ) {
      PatternGenerator pg;
      fPatternTree = pg.Generate( tp );
      if( fPatternTree ) {
	// Write the freshly-generated tree to file
	// FIXME: hmmm... we don't necesarily have write permission to DB_DIR
	//       fPatternTree->Write( filename );
      } else
	return fStatus = kInitError;
    // }

    // Set up a hitpattern object with the parameters of this projection
    assert( fHitpattern == 0 );
    try { fHitpattern = MakeHitpattern( *fPatternTree ); }
    catch( bad_alloc ) { fHitpattern = 0; }
    if( !fHitpattern || fHitpattern->IsError() )
      return fStatus = kInitError;
    assert( GetNallPlanes() == fHitpattern->GetNplanes() );

    cout<< "Pattern bin width " <<fHitpattern->GetBinWidth()<<endl;//getchar();

    // Determine maximum search distance (in bins) for combining patterns,
    // separately for front and back planes since they can have different
    // parameters. This is the max distance of bins that can belong to the
    // same hit, provided a hit is present in the plane.
    Plane *front_plane = fPlanes.front(), *back_plane = fPlanes.back();
    Double_t dxf= front_plane->GetMaxLRdist() + 2.0*front_plane->GetResolution();
    Double_t dxb= back_plane->GetMaxLRdist()  + 2.0*back_plane->GetResolution();
    fFrontMaxBinDist = TMath::CeilNint( dxf * fHitpattern->GetBinScale() );
    fBackMaxBinDist  = TMath::CeilNint( dxb * fHitpattern->GetBinScale() );

  } // doing_tracking

  // Special handling of calibration mode: Allow missing hits in calibration
  // planes, and require hits in all other planes
  UInt_t ncalib = 0;
  for( UInt_t k = 0; k < GetNplanes(); ++k ) {
    if( fPlanes[k]->IsCalibrating() )
      ++ncalib;
  }
  if( ncalib > 0 ) {
    Info( Here(here), "Calibrating %d planes in %s-projection.",
	  ncalib, GetName() );
    fMaxMiss = ncalib;
    for( UInt_t k = 0; k < GetNplanes(); ++k ) {
      if( fPlanes[k]->IsCalibrating() ) {
	if( fRequire1of2 and fPlanes[k]->GetPartner() and
	    fPlanes[k]->GetPartner()->IsCalibrating() ) {
	  Error( Here(here), "Cannot require 1 of 2 partner planes and "
		 "calibrate both partners (%s and %s) simultaneously. "
		 "Fix database.", fPlanes[k]->GetName(),
		 fPlanes[k]->GetPartner()->GetName() );
	  return kInitError;
	}
      } else
	fPlanes[k]->SetRequired();
    }
  }

  // Check fMaxMiss (maximum number of missing planes allowed in a hitpattern).
  // This is done here instead of in ReadDatabase because we need GetNplanes().
  assert( kMinFitPlanes <= GetNplanes() ); // ensured above
  UInt_t allowed_maxmiss = GetNplanes()-kMinFitPlanes;
  // There cannot be so many planes missing/calibrating that we don't have at
  // least kMinFitPlanes left
  if( fMaxMiss > allowed_maxmiss ) {
    if( ncalib == 0 ) {
      Error( Here(here), "Illegal number of allowed missing planes = %u. "
	     "Max allowed = %u. Fix database.", fMaxMiss, allowed_maxmiss );
    } else {
      Error( Here(here), "Too many planes in calibration mode, found %u, "
	     "max allowed = %u. Fix database.", ncalib, allowed_maxmiss );
    }
    return kInitError;
  }
  assert( fMaxMiss+kMinFitPlanes <= GetNplanes() );
  fMinFitPlanes = GetNplanes()-fMaxMiss;
  assert( fMinFitPlanes >= kMinFitPlanes );

  // Set up the lookup bitpattern indicating which planes are allowed to have
  // missing hits in fits. The value of the bit pattern of plane hits is used
  // as an index into this table; if the corresponding bit is set, the plane
  // combination is allowed.
  MakePlaneCombos( fPlanes, fPlaneCombos );

  // For TreeSerach, the allowed plane pattern is different in case we have
  // dummy planes. We consider dummy planes to be always required (why else
  // would one configure them?).
  if( TestBit(kHaveDummies) ) {
    assert( GetNplanes() < GetNallPlanes() );
    MakePlaneCombos( fAllPlanes, fAltPlaneCombos );
  } else {
    // If there are no dummy planes, the allowed plane patterns are the same
    // for TreeSearch and for fits.
    fAltPlaneCombos = fPlaneCombos;
  }

  // Determine Chi2 confidence interval limits for the selected CL and the
  // possible degrees of freedom (minfit-2...nplanes-2) of the projection fit
  if( TestBit(kDoChi2) ) {
    fChisqLimits.clear();
    fChisqLimits.resize( GetNplanes()-1, make_pair<Double_t,Double_t>(0,0) );
    for( vec_pdbl_t::size_type dof = kMinFitPlanes-2;
	 dof < fChisqLimits.size(); ++dof ) {
      fChisqLimits[dof].first = TMath::ChisquareQuantile( fConfLevel, dof );
      fChisqLimits[dof].second = TMath::ChisquareQuantile( 1.0-fConfLevel, dof );
    }
  }

  fPatternsFound.reserve( 200 );

  return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t Projection::ReadDatabase( const TDatime& date )
{
  // Read parameters from database

  static const char* const here = "ReadDatabase";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Double_t angle = kBig;
  fHitMaxDist = 0;
  fMaxMiss = 0;
  fMaxPat  = kMaxUInt;
  fConfLevel = 1e-3;
  Int_t req1of2 = 0, disable_chi2 = 0;

  Int_t gbl = Plane::GetDBSearchLevel(fPrefix);
  vector<Int_t> *moduleOrder = new vector<Int_t>;
  const DBRequest request[] = {
    { "angle",           &angle,         kDouble, 0, 1 },
    { "maxslope",        &fMaxSlope,     kDouble, 0, 1, gbl },
    { "search_depth",    &fNlevels,      kUInt},//   0, 0, gbl },
    { "do_eltracksel",   &fDoElTrackSel, kInt,    0, 1, gbl },
    // {"CorrSlope",        &fCorrSlope,    kDouble},
    // {"CorrSlope_min",    &fCorrSlopeMin, kDouble},
    // {"CorrSlope_max",    &fCorrSlopeMax, kDouble},
    // {"CorrIntercept_high", &fCorrInterceptHigh, kDouble},
    // {"CorrIntercept_low",  &fCorrInterceptLow,  kDouble},
    { "module_order",     moduleOrder,   kIntV,   0, 0, gbl},
    { "cluster_maxdist", &fHitMaxDist,   kUInt,   0, 1, gbl },
    { "chi2_conflevel",  &fConfLevel,    kDouble, 0, 1, gbl },
    { "maxmiss",         &fMaxMiss,      kUInt,   0, 1, gbl },
    { "req1of2",         &req1of2,       kInt,    0, 1, gbl },
    { "maxpat",          &fMaxPat,       kUInt,   0, 1, gbl },
    { "disable_chi2",    &disable_chi2,  kInt,    0, 1, gbl },
    { 0 }
  };
  Int_t err = LoadDB( file, date, request, fPrefix );
  if( err ){
    fclose(file);
    return kInitError;
  }
  
  if(fDoElTrackSel){
    cout << "request additional variables for elastic track selection" << endl;
    const DBRequest request_2[] = {
      {"CorrSlope",        &fCorrSlope,    kDouble},
      {"CorrSlope_min",    &fCorrSlopeMin, kDouble},
      {"CorrSlope_max",    &fCorrSlopeMax, kDouble},
      {"CorrIntercept_high", &fCorrInterceptHigh, kDouble},
      {"CorrIntercept_low",  &fCorrInterceptLow,  kDouble},
      { 0 }
    };
    err = LoadDB( file, date, request_2, fPrefix );
    if( err ){
      fclose(file);
      return kInitError;
    }
  }
  
  //cout<<"number of possible module order:"<<moduleOrder->size()<<endl;getchar();
  for(auto &order : *moduleOrder){
    fmoduleOrder[order] = 0;
  }

  fclose(file);
  if( err )
    return kInitError;

  if( fNlevels >= 16 ) {
    Error( Here(here), "Illegal search_depth = %u. Must be < 16. "
	   "Fix database.", fNlevels );
    return kInitError;
  }
  ++fNlevels; // The number of levels is maxdepth+1

  // If angle read, set it, otherwise keep default from call to constructor
  if( angle < kBig )
    SetAngle( angle*TMath::DegToRad() );

  if( fMaxSlope < 0.0 ) {
    Warning( Here(here), "Negative maxslope = %lf makes no sense. "
	     "Using |maxslope|.", fMaxSlope );
    fMaxSlope = -fMaxSlope;
  }

  SetBit( kDoChi2, !disable_chi2 );
  if( TestBit(kDoChi2) ) {
    if( fConfLevel < 0.0 || fConfLevel > 1.0 ) {
      Error( Here(here), "Illegal fit confidence level = %lf. "
	     "Must be 0-1. Fix database.", fConfLevel );
      return kInitError;
    }
  }

  fRequire1of2 = (req1of2 != 0);

  // If any planes defined, update their coordinate offset
  // based on our possibly new angle
  for( vplsiz_t i = 0; i < fPlanes.size(); ++i ) {
    assert( fPlanes[i]->GetProjection() == this ); // else bug in AddPlane
    fPlanes[i]->UpdateOffset();
  }

  fIsInit = kTRUE;
  return kOK;
}

//_____________________________________________________________________________
Int_t Projection::DefineVariables( EMode mode )
{
  // Initialize global variables and lookup table for decoder

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Don't bother with these global variables if no tracking is done
  bool doing_tracking = fDetector->TestBit(Tracker::kDoCoarse);
  if( !doing_tracking ) return kOK;

  // Global variables
  RVarDef vars[] = {
#ifdef TESTCODE
    { "n_hits", "Number of hits used for filling hitpattern", "n_hits" },
    { "n_bins", "Number of bins set in hitpattern", "n_bins" },
    { "n_binhits", "Number of references from bins to hits","n_binhits" },
    { "n_maxhits_bin", "Max number of hits per bin", "maxhits_bin" },
    { "n_test", "Number of pattern comparisons", "n_test"  },
    { "n_pat", "Number of patterns found",   "n_pat"    },
    { "n_roads", "Number of roads before filter",   "n_roads"    },
    { "n_dupl",  "Number of duplicate roads removed",   "n_dupl"    },
    { "n_badfits", "Number of roads found",   "n_badfits"    },
    { "t_treesearch", "Time in TreeSearch (us)", "t_treesearch" },
    { "t_roads", "Time in MakeRoads (us)", "t_roads" },
    { "t_fit", "Time for fitting Roads (us)", "t_fit" },
    { "t_track", "Total time in Track (us)", "t_track" },
    { "rd.nfits", "Number of acceptable fits in road",
                                         "fRoads.TreeSearch::Road.fNfits" },
#endif
    { "nroads","Number of roads (good or bad)",        "GetNroads()"      },
    { "ngood", "Number of good roads",                 "fNgoodRoads"      },
    { "rd.pos",   "Origin of best track (m)",
                                           "fRoads.TreeSearch::Road.fPos" },
    { "rd.slope", "Slope (dx/dz) of best track",
                                         "fRoads.TreeSearch::Road.fSlope" },
    { "rd.chi2",  "Chi2 of best fit",
                                          "fRoads.TreeSearch::Road.fChi2" },
    { "rd.dof",   "Degrees of freedom of best fit",
                                           "fRoads.TreeSearch::Road.fDof" },
    { "rd.pat",   "Bit pattern of plane numbers of hits",
                                  "fRoads.TreeSearch::Road.fPlanePattern" },
    { "rd.good",  "Road has valid data",
                                          "fRoads.TreeSearch::Road.fGood" },
    { "rd.trkstat", "Road reconstruction status",
                                       "fRoads.TreeSearch::Road.fTrkStat" },
    { "rd.fSigRatio", "Road nSignal/nTotal ",
                                       "fRoads.TreeSearch::Road.fSigRatio" },


    { "trkstat", "2D track reconstruction status",  "fTrkStat" },
    { 0 }
  };
  DefineVarsFromList( vars, mode );

  // Additional information about the roads found, for event display
  if( TestBit(kEventDisplay) ) {
    RVarDef vars_evtdisp[] = {
      { "rd.xLL","Lower left corner x coordinate (m)",
	              "fRoadCorners.TreeSearch::Road::Corners.fXLL" },
      { "rd.xLR","Lower right corner x coordinate (m)",
	              "fRoadCorners.TreeSearch::Road::Corners.fXLR" },
      { "rd.zL", "Lower edge z coordinate (m)",
                       "fRoadCorners.TreeSearch::Road::Corners.fZL" },
      { "rd.xUL","Upper left corner x coordinate (m)",
	              "fRoadCorners.TreeSearch::Road::Corners.fXUL" },
      { "rd.xUR","Upper right corner x coordinate (m)",
	              "fRoadCorners.TreeSearch::Road::Corners.fXUR" },
      { "rd.zU", "Upper edge z coordinate (m)",
                       "fRoadCorners.TreeSearch::Road::Corners.fZU" },
      { 0 }
    };
    DefineVarsFromList( vars_evtdisp, mode );
  }

#ifdef MCDATA
  // Additional variables for MC input data
  if( TestBit(kMCdata) ) {
    RVarDef mcvars[] = {
      { "nmcroads","Number of good roads with all MC info", "fNMCRoads"      },
      { "rd.nmcplanes", "Number of planes with hits from MC track",
                                  "fRoads.TreeSearch::Road.fNMCTrackHits" },
      { "rd.mcpat",     "Bit pattern of plane nums w/hits from MC track",
                           "fRoads.TreeSearch::Road.fMCTrackPlanePattern" },
      { "rd.nmcfit",    "Number of planes with hits from fitted MC track",
                               "fRoads.TreeSearch::Road.fNMCTrackHitsFit" },
      { "rd.mcfitpat",  "Bit pattern of planes nums w/fitted MC track hits",
                        "fRoads.TreeSearch::Road.fMCTrackPlanePatternFit" },
      { 0 }
    };
    DefineVarsFromList( mcvars, mode );
  }
#endif

  return 0;
}

//_____________________________________________________________________________
Int_t Projection::FillHitpattern()
{
  // Fill this projection's hitpattern from hits in the planes.
  // Returns the total number of hits processed.

  Int_t ntot = fHitpattern->Fill( fAllPlanes );
  if(fDebug>=3){
    cout << "EF: printing pattern";fHitpattern->Print();
  }

#ifdef TESTCODE
  n_hits = ntot;
  n_bins = fHitpattern->GetBinsSet();
  n_binhits = fHitpattern->GetNhits();
  maxhits_bin = fHitpattern->GetMaxhitBin();
#endif
  
  return ntot;
}

//_____________________________________________________________________________
Int_t Projection::Track()
{
  // Perform tracking in this projection:
  //
  // - match hits to straight-line patterns
  // - combine patterns referring to common sets of hits into one ("Roads")
  // - fit hits in each road
  // - filter roads according to chi^2 and similarity
  //
  // Results in fRoads

  Int_t ret = 0;

  // TreeSearch:
  // Match the hitpattern of the current event against the pattern template
  // database. Results in fPatternsFound.

  assert( fPatternsFound.empty() );
  assert( GetTrackingStatus() == kTrackOK );

#ifdef TESTCODE
  TStopwatch timer, timer_tot;
#endif

  ComparePattern compare( fHitpattern, fAltPlaneCombos, &fPatternsFound,
			  fDummyPlanePattern );
  TreeWalk walk( fNlevels );
  walk( fPatternTree->GetRoot(), compare );

#ifdef VERBOSE
  if( fDebug > 0 ) {
    UInt_t npat = fPatternsFound.size();
    cout << npat << " pattern";
    if( npat!=1 ) cout << "s";
    if( npat > fMaxPat )
      cout << " >>> exceeding limit of " << fMaxPat << ", terminating";
    cout << endl;
  }
#endif
#ifdef TESTCODE
  t_treesearch = 1e6*timer.RealTime();

  n_test = compare.GetNtest();
  n_pat  = fPatternsFound.size();

  timer.Start();
#endif
  
  if( fPatternsFound.empty() ) {
    fTrkStat = kNoPatterns;
    goto quit;
  }
  // Die if too many patterns - noisy event
  if( (UInt_t)fPatternsFound.size() > fMaxPat ) {
    // TODO: keep statistics
    cout<<" Too many patterns: "<<(UInt_t)fPatternsFound.size()<<", giving up..."<<endl;
    fTrkStat = kTooManyPatterns;
    ret = -1;
    goto quit;
  }

  // Combine patterns with common sets of hits into Roads
  MakeRoads();

#ifdef VERBOSE
  if( fDebug > 0 ) {
    if( !fRoads->IsEmpty() ) {
      Int_t nroads = GetNroads();
      cout << nroads << " road";
      if( nroads>1 ) cout << "s";
      cout << endl;
    }
  }
#endif
#ifdef TESTCODE
  n_roads = GetNroads();
#endif

  // This seems to cost much more than it gains
//   // Check for indentical roads or roads that include each other.
//   // Any roads that are eliminated are marked as void.
//   if( GetNroads() > 1 ) {
//     if( RemoveDuplicateRoads() )
//       // Remove empty slots caused by removed duplicate roads
//       fRoads->Compress();
//   }

// #ifdef VERBOSE
//   if( fDebug > 0 ) {
//     if( !fRoads->IsEmpty() ) {
//       Int_t nroads = GetNroads();
//       cout << nroads << " road";
//       if( nroads>1 ) cout << "s";
//       cout << " after filter" << endl;
//     }
//   }
// #endif
#ifdef TESTCODE
  t_roads = 1e6*timer.RealTime();
  timer.Start();
#endif

  // Fit hit positions in the roads to straight lines
  FitRoads();//module order check included
  
  if(fDoElTrackSel){
    RemoveNonElasticTracks();
  }
  

#ifdef TESTCODE
  t_fit   = 1e6*timer.RealTime();
  t_track = 1e6*timer_tot.RealTime();
#endif

#ifdef VERBOSE
  if( fDebug > 0 ) {
    if( !fRoads->IsEmpty() ) {
      Int_t nroads = GetNgoodRoads();
      cout << nroads << " road";
      if( nroads!=1 ) cout << "s";
      cout << " successfully fit" << endl;
    }
  }
#endif

  ret = GetNgoodRoads();

 quit:
#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << "------------ end of projection  " << GetName()
	 << "------------" << endl;
  }
#endif
  return ret;
}


Int_t Projection::RemoveNonElasticTracks(){
    for( UInt_t i = 0; i < GetNroads(); ++i ) {
      Road* rd = static_cast<Road*>(fRoads->UncheckedAt(i));
      assert(rd);
      if(rd -> IsVoid())
	continue;
      if(fCorrSlope * rd->GetSlope() + fCorrInterceptHigh < rd->GetPos()
	 || fCorrSlope * rd->GetSlope() + fCorrInterceptLow > rd->GetPos()
	 || rd->GetSlope() < fCorrSlopeMin
	 || rd->GetSlope() > fCorrSlopeMax){      
	rd -> Void();
	fNgoodRoads--;
      }
    }
  
  return 0;
}


//_____________________________________________________________________________
#ifdef VERBOSE
static void PrintNode( const Node_t& node )
{
  const HitSet& hs = node.second;
  cout << " npl/nhits= " << hs.nplanes << "/" << hs.hits.size()
       << "  hitpos= ";
  UInt_t ipl = 0;
  for( Hset_t::iterator ihit = hs.hits.begin(); ihit != hs.hits.end(); ) {
    assert( (*ihit)->GetPlaneNum() != kMaxUInt );
    while( ipl < (*ihit)->GetPlaneNum() ) { cout << "--/"; ++ipl; }
    bool seq = false;
    do {
      assert( (*ihit)->GetPlaneNum() != kMaxUInt );
      if( seq )  cout << " ";
      cout << (*ihit)->GetPos();
      seq = true;
    } while( ++ihit != hs.hits.end() and (*ihit)->GetPlaneNum() == ipl );
    if( ipl != node.first.link->GetPattern()->GetNbits()-1 ) {
      cout << "/";
      if( ihit == hs.hits.end() )
	cout << "--";
    }
    ++ipl;
  }
  cout << "  pat= ";
  node.first.Print();
}

inline
static void PrintNodeP( const Node_t* node )
{
  PrintNode( *node );
}

#endif

//_____________________________________________________________________________
// Comparison functors used for sorting patterns in MakeRoads()

struct MostPlanes : public binary_function< Node_t*, Node_t*, bool >
{
  bool operator() ( const Node_t* a, const Node_t* b ) const
  {
    // Order patterns by decreasing number of active planes, then
    // decreasing number of hits, then ascending bin numbers
    if( a->second.nplanes != b->second.nplanes )
      return (a->second.nplanes > b->second.nplanes);
    return (a->second.hits.size() != b->second.hits.size()) ?
      (a->second.hits.size() > b->second.hits.size()) :
      (a->first < b->first);
  }
};

struct BinIsLess : public binary_function< Node_t*, Node_t*, bool >
{
  bool operator() ( const Node_t* a, const Node_t* b ) const
  {
    // Order by bin number only
    return ( a->first < b->first );
  }
};

//_____________________________________________________________________________
Int_t Projection::MakeRoads()
{
  // Combine patterns with common sets of hits into Roads.
  //
  // This is the primary de-cloning algorithm. It finds clusters of patterns
  // that share active wires (hits).

  // Sort patterns according to MostPlanes (see above)
  sort( ALL(fPatternsFound), MostPlanes() );

  // Copy patterns to secondary key sorted by bin number only. This key
  // greatly improves lookup speed of potential similar patterns
  typedef set<const Node_t*,BinIsLess> BinOrdNodes_t;
  BinOrdNodes_t nodelookup;

  // Inserting one-by-one with hint is faster than range insert without hint
  // since fPatternsFound is already sorted in a similar order, so the hint
  // is often good
  copy( ALL(fPatternsFound), inserter( nodelookup, nodelookup.end() ));

  assert( fPatternsFound.size() == nodelookup.size() );

#ifdef VERBOSE
  if( fDebug > 2 ) {
    cout << fPatternsFound.size() << " patterns found:" << endl;
    for_each( ALL(fPatternsFound), PrintNodeP );

    cout << "--------------------------------------------" << endl;
    cout << nodelookup.size() << " patterns sorted by bin:" << endl;
    for_each( ALL(nodelookup), PrintNodeP );
  }
#endif

  // Build roads starting with patterns that have the most active planes.
  // These tend to yield the best track candidates.
  for( NodeVec_t::iterator it = fPatternsFound.begin(); it !=
	 fPatternsFound.end(); ++it ) {

    

    const Node_t& nd1 = **it;

    if( nd1.second.used )
      continue;

    // Start a new road with next unused pattern 
    Road* rd = new( (*fRoads)[GetNroads()] ) Road(nd1,this);
    if(fDebug>=3){
      cout << "EF: New road: ";rd->Print();
    }

    // Try to add similar patterns to this road (cf. HitSet::IsSimilarTo)
    // Since only patterns with front bin numbers near the start pattern
    // are candidates, search along the start bin index built above.
    BinOrdNodes_t::iterator jt = nodelookup.find( &nd1 );
    assert( jt != nodelookup.end() );
    assert( (*jt)->first == nd1.first );

    // Test patterns in direction of decreasing front bin number index,
    // beginning with the road start pattern, until they are too far away.

    // need to rescan if fHitMaxDist > 0 because IsInFrontRange may catch
    // more patterns after patterns with new hits have been added
    while( rd->HasGrown() ) {
      rd->ClearGrow();
      // The following runs much slower with a reverse_iterator
      BinOrdNodes_t::iterator jr(jt);
      if( jt != nodelookup.begin() ) {
	--jr;
	while( rd->IsInFrontRange(**jr) ) {
	  if( rd->Add(**jr) ) {
	    // Pattern successfully added
	    // Erase used patterns from the lookup index
	    if( jr == nodelookup.begin() ) {
	      nodelookup.erase( jr );
	      break;
	    } else {
	      nodelookup.erase( jr-- );
	    }
	  } else if( jr != nodelookup.begin() ) {
	    --jr;
	  } else
	    break;
	}
      }
     }
    // Repeat in the forward direction along the index
    rd->SetGrow();
    while( rd->HasGrown() ) {
      rd->ClearGrow();
      BinOrdNodes_t::iterator jf(jt);
      ++jf;
      while( jf != nodelookup.end() and rd->IsInFrontRange(**jf) ) {
	if( rd->Add(**jf) )
	  nodelookup.erase( jf++ );
	else
	  ++jf;
      }
    }
    nodelookup.erase( jt );

    // Update the "used" flags of the road's component patterns
    rd->Finish();
    if(fDebug>=3){
      cout << "EF: road completed: ";rd->Print();
    }
    
    // If event display enabled, export the road's corner coordinates
    if( TestBit(kEventDisplay) ) {
      assert( fRoads->GetLast() == fRoadCorners->GetLast()+1 );
      new( (*fRoadCorners)[fRoads->GetLast()] ) Road::Corners(rd);
    }
  }
  assert( nodelookup.empty() );


#ifdef VERBOSE
  if( fDebug > 2 ) {
    cout << "Generated roads: " << endl;
    for( UInt_t i = 0; i < GetNroads(); ++i ) {
      const Road* rd = GetRoad(i);
      const Road::NodeList_t& ndlst = rd->GetPatterns();
      for_each( ALL(ndlst), PrintNodeP );
      cout << "--------------------------------------------" << endl;
    }
  }
#endif

  return 0;
}

//_____________________________________________________________________________
Bool_t Projection::RemoveDuplicateRoads()
{
  // Check for identical roads or roads that include each other

  // This runs in ~O(N^2) time, but if N>1, it is typically only 2-5.
  bool changed = false, restart = true;
  while( restart ) {
    restart = false;
    for( UInt_t i = 0; i < GetNroads(); ++i ) {
      Road* rd = static_cast<Road*>(fRoads->UncheckedAt(i));
      if( !rd )
	continue;
      for( UInt_t j = i+1; j < GetNroads(); ++j ) {
	Road* rd2 = static_cast<Road*>(fRoads->UncheckedAt(j));
	if( !rd2 )
	  continue;
	if( rd->Include(rd2) ) {
	  fRoads->RemoveAt(j);
	  changed = restart = true;
#ifdef TESTCODE
	  ++n_dupl;
#endif
	} else if( rd2->Include(rd) ) {
	  fRoads->RemoveAt(i);
	  changed = restart = true;
#ifdef TESTCODE
	  ++n_dupl;
#endif
	}
	if( restart )
	  break;
      }
      if( restart )
	break;
    }
  }
  return changed;
}

//_____________________________________________________________________________
Bool_t Projection::FitRoads()
{
  // Fit hits within each road. Store fit parameters with Road.
  // Also, store the hits & positions used by the best fit with Road.
  bool changed = false;

  for( UInt_t i = 0; i < GetNroads(); ++i ) {
    Road* rd = static_cast<Road*>(fRoads->UncheckedAt(i));
    assert(rd);
    if( rd->Fit(fmoduleOrder) ){
      // Count good roads (not void and good fit)
      ++fNgoodRoads;
#ifdef TESTCODE
#ifdef MCDATA
      if(rd->GetPlanePattern()==rd->GetMCTrackPlanePattern())++fNMCRoads;
#endif
#endif
    }
    else {
      changed = true;
#ifdef TESTCODE
      ++n_badfits;
#endif
    }
  }
  if( GetNroads() > 0 && fNgoodRoads == 0 )
    fTrkStat = kFailed2DFits;

  return changed;
}

//_____________________________________________________________________________
void Projection::MakePrefix()
{
  // Set up name prefix for global variables.

  TString basename;
  if( fDetector ) {
    basename = fDetector->GetPrefix();
    basename.Chop();  // delete trailing dot
  } else
    Warning( Here("MakePrefix"), "No parent detector defined. "
	     "Using \"%s\".", GetName() );

  THaAnalysisObject::MakePrefix( basename.Data() );
}

//_____________________________________________________________________________
const char* Projection::GetDBFileName() const
{
  // Return database file name prefix. We want the same database file
  // as our parent detector.

  return fDetector ? fDetector->GetDBFileName() : GetPrefix();
}

//_____________________________________________________________________________
Double_t Projection::GetZsize() const
{
  // Get z_max - z_min of the planes.

  assert( !fPlanes.empty() );

  return fPlanes.back()->GetZ() - fPlanes.front()->GetZ();
}

//_____________________________________________________________________________
void Projection::SetAngle( Double_t angle )
{
  // Set angle of the axis perpendicular to the wires/strips (rad)

  // Ensure that the angle is consistent with the plane type:
  // x and y planes must have angles close to zero or 90 degrees.
  // All other plane types must NOT be close to zero or 90 degrees.
  // Here, "close" means "within kAngleTolerance".
  // Throws bad_angle exception if the restrictions are violated.

  Double_t d = angle -
    TMath::FloorNint( angle/TMath::PiOver2() ) * TMath::PiOver2();
  Bool_t near_zero_or_90 =
    (TMath::Abs(d-TMath::PiOver4()) > TMath::PiOver4()-kAngleTolerance);
  Bool_t x_or_y = (fType == kXPlane or fType == kYPlane);
  if( near_zero_or_90 xor x_or_y ) {
    stringstream msg("Angle = ");
    msg << angle*TMath::RadToDeg()
	<< " too " << (( x_or_y ) ? "far off from" : "close to")
	<< "0 or 90 degrees for projection type \""
	<< kProjParam[fType].name << "\". Fix database.";
    Error( Here("SetAngle"), "%s", msg.str().c_str() );
    throw bad_angle(msg.str());
  }

  // Normalize the angle to -180...180 degrees
  angle = TVector2::Phi_mpi_pi(angle);

  fAxis.Set( TMath::Cos(angle), TMath::Sin(angle) );
}

//_____________________________________________________________________________
void Projection::Print( Option_t* opt ) const
{
  // Print plane type info

  Int_t verbose = 0;
  if( opt ) {
    TString opt_s(opt);
    opt_s.ToLower();
    verbose = opt_s.CountChar('v');
  }

  cout << "Projection:  "
       << GetName()
       << " type="  << (GetType() == kUndefinedType
			? "(undef)" : kProjParam[GetType()].name)
       << " npl="   << GetNplanes()
       << " nlev="  << GetNlevels()
       << " zsize=" << (fPlanes.empty() ? 0.0 : GetZsize())
       << " maxsl=" << GetMaxSlope()
       << " width=" << GetWidth()
       << " angle=" << GetAngle()*TMath::RadToDeg();
  cout << endl;

  if( verbose > 0 ) {
    for( vplsiz_t i = 0; i < fPlanes.size(); ++i ) {
      Plane* pl = fPlanes[i];
      pl->Print(opt);
    }
  }
}

//_____________________________________________________________________________
Hitpattern* Projection::MakeHitpattern( const PatternTree& pt ) const
{
  // Instantiate Hitpattern to be used for this type of projection.
  // The Hitpattern determines how hit information it processed,
  // either as single position measurements or L/R-ambiguous wire hits.

  return new Hitpattern( pt );
}

//_____________________________________________________________________________
EProjType Projection::NameToType( const char* name )
{
  // Return the index corresponding to the given plane name.
  // The comparison is not case-sensitive.

  if( name and *name ) {
    TString s(name);
    s.ToLower();
    for( EProjType type = kTypeBegin; type < kTypeEnd; ++type ) {
      TString ps( kProjParam[type].name );
      ps.ToLower();
      if( s == ps )
	return type;
    }
  }
  return kUndefinedType;
}


//_____________________________________________________________________________
NodeVisitor::ETreeOp
Projection::ComparePattern::operator() ( const NodeDescriptor& nd )
{
  // Test if the pattern from the database that is given by NodeDescriptor
  // is present in the current event's hitpattern

#ifdef TESTCODE
  ++fNtest;
#endif
  // Compute the match pattern and see if it is allowed
  pair<UInt_t,UInt_t> match = fHitpattern->ContainsPattern(nd);
  if( fPlaneCombos->TestBitNumber(match.first)  ) {
    if( nd.depth < fHitpattern->GetNlevels()-1 )
      return kRecurse;

    // Found a match at the bottom of the pattern tree
    Node_t* node = new Node_t;
    node->first = nd;

    // Collect all hits associated with the pattern's bins and save them
    // in the node's HitSet.
    for( UInt_t i = 0; i < fHitpattern->GetNplanes(); ++i ) {
      const vector<Hit*>& hits = fHitpattern->GetHits( i, nd[i] );
      assert( hits.empty() or
	      (hits.front()->GetAltPlaneNum() == i and
	       not hits.front()->GetPlane()->IsDummy()) );
      node->second.hits.insert( ALL(hits) );
    }
    assert( (HitSet::GetAltMatchValue(node->second.hits) xor
	     fDummyPlanePattern) == match.first );
    if( fDummyPlanePattern != 0 ) {
      // If dummy planes are present, then match is given with respect to
      // Plane::GetAltPlaneNum(). We need to calculate the node's
      // plane_pattern and nplanes explicitly wrt to Plane::GetPlaneNum()
      node->second.CalculatePlanePattern();
    } else {
      // No dummy planes, less work :)
      node->second.plane_pattern = match.first;
      node->second.nplanes = match.second;
    }

    // Add the pointer to the new node to the vector of results
    //TODO: stop if max num of patterns reached
    fMatches->push_back( node );
  }
  return kSkipChildNodes;
}

//_____________________________________________________________________________

}  // end namespace TreeSearch

ClassImp(TreeSearch::Projection)

///////////////////////////////////////////////////////////////////////////////
