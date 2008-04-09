//*-- Author :    Ole Hansen, Jefferson Lab   16-Aug-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Projection.h"
#include "Hitpattern.h"
#include "WirePlane.h"
#include "THaDetectorBase.h"
#include "PatternTree.h"
#include "PatternGenerator.h"
#include "TreeWalk.h"
#include "Road.h"
#include "Helper.h"
#include "Hit.h"
#include "MWDC.h"   // for MWDC bits

#include "TMath.h"
#include "TString.h"
#include "TBits.h"
#include "TError.h"

#include <iostream>
#include <sys/time.h>  // for timing
#include <algorithm>
#include <utility>
#ifdef TESTCODE
#include <cstring>
#endif

using namespace std;

namespace TreeSearch {

typedef vector<WirePlane*>::size_type vwsiz_t;
typedef vector<WirePlane*>::iterator  vwiter_t;

//_____________________________________________________________________________
Projection::Projection( Int_t type, const char* name, Double_t angle,
			THaDetectorBase* parent )
  : THaAnalysisObject( name, name ), fType(type), fNlevels(0),
    fMaxSlope(0.0), fWidth(0.0), fDetector(parent), fPatternTree(0),
    fMinFitPlanes(3), fMaxMiss(0), fRequire1of2(true),
    fPlaneCombos(0), fMaxPat(kMaxUInt), fFrontMaxBinDist(kMaxUInt),
    fBackMaxBinDist(kMaxUInt), fHitMaxDist(0), fConfLevel(1.0),
    fHitpattern(0), fRoads(0), fNgoodRoads(0), fRoadCorners(0)
{
  // Constructor

  assert( name && parent );

  // angle is the default value. It can be overridden via a database entry.
  SetAngle( angle );

  fTitle.Append(" projection");
  fRoads = new TClonesArray("TreeSearch::Road", 3);
  R__ASSERT(fRoads);
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
  delete fPlaneCombos;
}

//_____________________________________________________________________________
void Projection::AddPlane( WirePlane* wp, WirePlane* partner )
{
  // Add wire plane wp (and optional partner plane) to this projection. 
  // Sets plane numbers.

  assert(wp);

  wp->SetPlaneNum( fPlanes.size() );
  fPlanes.push_back( wp );
  wp->SetProjection( this );
  if( partner ) {
    assert( partner->GetZ() > wp->GetZ() ); // Planes must be ordered
    partner->SetPlaneNum( fPlanes.size() );
    fPlanes.push_back( partner );
    partner->SetProjection( this );
  }
}

//_____________________________________________________________________________
void Projection::Clear( Option_t* )
{    
  // Clear event-by-event data

  if( fHitpattern )
    fHitpattern->Clear();

  fRoads->Delete();
  fPatternsFound.clear();
  fNgoodRoads = 0;

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
  for( vwsiz_t i = 0; i < GetNplanes(); ++i ) {
    WirePlane* wp = fPlanes[i];
    Int_t nhits = wp->Decode( evdata );
    if( nhits < 0 ) {
      err = true;
      sum -= nhits;
    } else
      sum += nhits;
  }
  if( err )
    return -sum;

  return sum;
}

//_____________________________________________________________________________
Double_t Projection::GetPlaneZ( UInt_t i ) const
{
  // Return the z-position of the i-th wire plane.

  assert( i<fPlanes.size() );
  return fPlanes[i]->GetZ();
}

//_____________________________________________________________________________
void Projection::Reset()
{
  // Reset parameters, clear list of planes, delete Hitpattern
  
  fIsInit = kFALSE;
  fPlanes.clear();
  fMaxSlope = fWidth = 0.0;
  delete fHitpattern; fHitpattern = 0;
  delete fPatternTree; fPatternTree = 0;
  delete fPlaneCombos; fPlaneCombos = 0;
  delete fRoadCorners; fRoadCorners = 0;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Projection::Init( const TDatime& date )
{
  // Initialize the Projection object. Called after MWDC basic initialization.
  // Sets up event display support, then continues with standard 
  // initialization.

  Reset();

  // Set up the event display support if corresponding bit is set in the MWDC
  if( fDetector->TestBit(MWDC::kEventDisplay) ) {
    assert( fRoadCorners == 0 );
    fRoadCorners = new TClonesArray("TreeSearch::Road::Corners", 3);
    R__ASSERT(fRoadCorners);
    // Set local bit to indicate that initialization is done
    SetBit(kEventDisplay);
  }

  // Standard initialization. This calls ReadDatabase() and DefineVariables()
  return THaAnalysisObject::Init(date);
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Projection::InitLevel2( const TDatime& )
{
  // Level-2 initialization - load pattern database and initialize hitpattern.
  // Requires the wire planes to be fully initialized.

  static const char* const here = "InitLevel2";

  vector<Double_t> zpos;
  for( UInt_t i = 0; i < GetNplanes(); ++i )
    zpos.push_back( GetPlaneZ(i) );
  TreeParam_t tp( fNlevels-1, fWidth, fMaxSlope, zpos );
		  
  if( tp.Normalize() != 0 )
    return fStatus = kInitError;

  // Attempt to read the pattern database from file
  assert( fPatternTree == 0 );
  //TODO: Make the file name
  const char* filename = "test.tree";
  fPatternTree = PatternTree::Read( filename, tp );
  
  // If the tree cannot not be read (or the parameters mismatch), then
  // create it from scratch (takes a few seconds)
  if( !fPatternTree ) {
    PatternGenerator pg;
    fPatternTree = pg.Generate( tp );
    if( fPatternTree ) {
      // Write the freshly-generated tree to file
      // FIXME: hmmm... we don't necesarily have write permission to DB_DIR
//       fPatternTree->Write( filename );
    } else 
      return fStatus = kInitError;
  } 

  // Set up a hitpattern object with the parameters of this projection
  assert( fHitpattern == 0 );
  fHitpattern = new Hitpattern( *fPatternTree );
  if( !fHitpattern || fHitpattern->IsError() )
    return fStatus = kInitError;
  assert( GetNplanes() == fHitpattern->GetNplanes() );

  // Determine maximum search distance (in bins) for combining patterns,
  // separately for front and back planes since they can have different
  // parameters.
  // This is the max distance of bins that can belong to the same hit 
  // plus an allowance for extra slope of a pattern if a front/back hit
  // is missing
  WirePlane *front_plane = fPlanes.front(), *back_plane = fPlanes.back();
  Double_t dxf= front_plane->GetMaxLRdist() + 2.0*front_plane->GetResolution();
  Double_t dxb= back_plane->GetMaxLRdist()  + 2.0*back_plane->GetResolution();
  fFrontMaxBinDist = TMath::CeilNint( dxf * fHitpattern->GetBinScale() ) + 2;
  fBackMaxBinDist  = TMath::CeilNint( dxb * fHitpattern->GetBinScale() ) + 2;

  // Check range of fMinFitPlanes (minimum number of planes require for fit)
  // and fMaxMiss (maximum number of missing planes allowed in a hitpattern).
  // This is here instead of in ReadDatabase because we need GetNplanes()
  if( fMinFitPlanes < 3 || fMinFitPlanes > GetNplanes() ) {
    Error( Here(here), "Illegal number of required planes for fitting = %u. "
	   "Must be >= 3 and <= %u. "
	   "Fix database.", fMinFitPlanes, GetNplanes() );
    return kInitError;
  }
  if( fMaxMiss > GetNplanes()-1 ) {
    Error( Here(here), "Illegal number of allowed missing planes = %u. "
	   "Must be <= %u. Fix database.", fMaxMiss, GetNplanes()-1 );
    return kInitError;
  }
  // There cannot be so many planes missing that we don't have at least
  // fMinFitPlanes left
  UInt_t maxmiss = GetNplanes()-fMinFitPlanes;
  if( fMaxMiss > maxmiss ) {
    Warning( Here(here), "Allowed number of missing planes = %u reduced to %u "
	     "to satisfy min_fit_planes = %u.  Fix database.", 
	     fMaxMiss, maxmiss, fMinFitPlanes );
    fMaxMiss = maxmiss;
  }

  // Set up the lookup bitpattern indicating which planes are allowed to have 
  // missing hits. The value of the bit pattern of plane hits is used as an
  // index into this table; if the corresponding bit is set, the plane 
  // combination is allowed.
  assert( fPlaneCombos == 0 );
  UInt_t np = 1U<<GetNplanes();
  fPlaneCombos = new TBits( np );
  fPlaneCombos->SetBitNumber( np-1 );  // Always allow full occupancy
  for( UInt_t i = 1; i <= fMaxMiss; ++i ) {
    UniqueCombo c( GetNplanes(), i );
    while( c ) {
      // Clear the bit numbers from this combination
      UInt_t bitval = np-1;
      for( vector<int>::size_type j = 0; j < c().size(); ++j )
	bitval &= ~( 1U << c()[j] );
      // Test if this pattern satisfies other constraints
      UInt_t k = 0;
      for( ; k < GetNplanes(); ++k ) {
	// Disallow bit pattern if a required plane is missing
	if( 0 == (bitval & (1U<<k)) and fPlanes[k]->IsRequired() )
	  break;
	// If requested, ensure that at least one plane of a plane pair is set
	if( fRequire1of2 and fPlanes[k]->GetPartner() and
	    0 == (bitval & (1U<<k)) and
	    0 == (bitval & (1U<<fPlanes[k]->GetPartner()->GetPlaneNum())) )
	  break;
      }
      assert( bitval < np );
      if( k == GetNplanes() )
	fPlaneCombos->SetBitNumber( bitval );
      ++c;
    }
  }

  // Determine Chi2 confidence interval limits for the selected CL and the
  // possible degrees of freedom (minfit-2...nplanes-2) of the projection fit
  fChisqLimits.clear();
  fChisqLimits.resize( GetNplanes()-1, make_pair<Double_t,Double_t>(0,0) );
  for( vec_pdbl_t::size_type dof = fMinFitPlanes-2;
       dof < fChisqLimits.size(); ++dof ) {
    fChisqLimits[dof].first = TMath::ChisquareQuantile( fConfLevel, dof );
    fChisqLimits[dof].second = TMath::ChisquareQuantile( 1.0-fConfLevel, dof );
  }
  
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
  fMinFitPlanes = 3;
  fMaxMiss = 0;
  fMaxPat  = kMaxUInt;
  fConfLevel = 1e-3;
  Int_t req1of2 = 1;
  const DBRequest request[] = {
    { "angle",           &angle,         kDouble, 0, 1 },
    { "maxslope",        &fMaxSlope,     kDouble, 0, 1, -1 },
    { "search_depth",    &fNlevels,      kUInt,   0, 0, -1 },
    { "cluster_maxdist", &fHitMaxDist,   kUInt,   0, 1, -1 },
    { "min_fit_planes",  &fMinFitPlanes, kUInt,   0, 1, -1 },
    { "chi2_conflevel",  &fConfLevel,    kDouble, 0, 1, -1 },
    { "maxmiss",         &fMaxMiss,      kUInt,   0, 1, -1 },
    { "req1of2",         &req1of2,       kInt,    0, 1, -1 },
    { "maxpat",          &fMaxPat,       kUInt,   0, 1, -1 },
    { 0 }
  };

  Int_t err = LoadDB( file, date, request, fPrefix );
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

  if( fConfLevel < 0.0 || fConfLevel > 1.0 ) {
    Error( Here(here), "Illegal fit confidence level = %lf. "
	   "Must be 0-1. Fix database.", fConfLevel );
    return kInitError;
  }

  fRequire1of2 = (req1of2 != 0);

  fIsInit = kTRUE;
  return kOK;
}

//_____________________________________________________________________________
Int_t Projection::DefineVariables( EMode mode )
{
  // Initialize global variables and lookup table for decoder

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Global variables
  RVarDef vars[] = {
#ifdef TESTCODE
    { "n_hits", "Number of hits used for filling hitpattern", "n_hits" },
    { "n_bins", "Number of bins set in hitpattern", "n_bins" },
    { "n_binhits", "Number of references from bins to hits","n_binhits" },
    { "maxhits_per_bin", "Max number of hits per bin", "maxhits_bin" },
    { "n_test", "Number of pattern comparisons", "n_test"  },
    { "n_pat", "Number of patterns found",   "n_pat"    },
    { "n_roads", "Number of roads before filter",   "n_roads"    },
    { "n_dupl",  "Number of duplicate roads removed",   "n_dupl"    },
    { "n_badfits", "Number of roads found",   "n_badfits"    },
    { "t_treesearch", "Time in TreeSearch (us)", "t_treesearch" },
    { "t_roads", "Time in MakeRoads (us)", "t_roads" },
    { "t_fit", "Time for fitting Roads (us)", "t_fit" },
    { "t_track", "Total time in Track (us)", "t_track" },
#endif
    { "nroads","Number of roads (good or bad)",        "GetNroads()"      },
    { "ngood", "Number of good roads",                 "fNgoodRoads"      },
    { "rd.nfits", "Number of acceptable fits in road",
                                     "fRoads.TreeSearch::Road.GetNfits()" },
    { "rd.pos",   "Origin of best track (m)",
                                           "fRoads.TreeSearch::Road.fPos" },
    { "rd.slope", "Slope (dx/dz) of best track",
                                         "fRoads.TreeSearch::Road.fSlope" },
    { "rd.chi2",  "Chi2 of best fit", 
                                          "fRoads.TreeSearch::Road.fChi2" },
    { "rd.dof",   "Degrees of freedom of best fit",
                                           "fRoads.TreeSearch::Road.fDof" },
    { "rd.good",  "Road has valid data",
                                          "fRoads.TreeSearch::Road.fGood" },
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
 
  return 0;
}

//_____________________________________________________________________________
Int_t Projection::FillHitpattern()
{
  // Fill this projection's hitpattern with hits from the wire planes.
  // Returns the total number of hits processed (where hit pairs on plane
  // and partner plane count as one).

  Int_t ntot = 0;
  for( vwiter_t it = fPlanes.begin(); it != fPlanes.end(); ++it ) {
    WirePlane *plane = *it, *partner = plane->GetPartner();
    // If a plane has a partner (usually with staggered wires), scan them
    // together to resolve some of the L/R ambiguity of the hit positions
    ntot += fHitpattern->ScanHits( plane, partner );
    // If the partner plane was just scanned, don't scan it again
    if( partner ) {
      ++it;
      assert( it != fPlanes.end() );
      assert( *it == partner );
    }
  }
#ifdef TESTCODE
  n_hits = ntot;
  n_bins = fHitpattern->GetBinsSet();
  n_binhits = fHitpattern->GetHitListSize();
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


  // TreeSearch:
  // Match the hitpattern of the current event against the pattern template
  // database. Results in fPatternsFound.

  fPatternsFound.clear();

#ifdef TESTCODE
  struct timeval start, start2, stop, diff;
  gettimeofday( &start, 0 );
  gettimeofday( &start2, 0 );
#endif

  ComparePattern compare( fHitpattern, fPlaneCombos, &fPatternsFound );
  TreeWalk walk( fNlevels );
  walk( fPatternTree->GetRoot(), compare );

#ifdef TESTCODE
  //FIXME: use high-res CPU time timer instead
  gettimeofday(&stop, 0 );
  timersub( &stop, &start2, &diff );
  t_treesearch = 1e6*(Double_t)diff.tv_sec + (Double_t)diff.tv_usec;

  n_test = compare.GetNtest();
  n_pat  = fPatternsFound.size();

  gettimeofday( &start2, 0 );
#endif

  if( fPatternsFound.empty() )
    return 0;
  // Die if too many patterns - noisy event
  if( (UInt_t)fPatternsFound.size() > fMaxPat )
    // TODO: keep statistics
    return -1;

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
  gettimeofday(&stop, 0 );
  timersub( &stop, &start2, &diff );
  t_roads = 1e6*(Double_t)diff.tv_sec + (Double_t)diff.tv_usec;

  gettimeofday( &start2, 0 );
#endif

  // Fit hit positions in the roads to straight lines
  FitRoads();

#ifdef TESTCODE
  gettimeofday(&stop, 0 );
  timersub( &stop, &start2, &diff );
  t_fit = 1e6*(Double_t)diff.tv_sec + (Double_t)diff.tv_usec;
  timersub( &stop, &start, &diff );
  t_track = 1e6*(Double_t)diff.tv_sec + (Double_t)diff.tv_usec;
#endif

#ifdef VERBOSE
  if( fDebug > 0 ) {
    if( !fRoads->IsEmpty() ) {
      Int_t nroads = GetNgoodRoads();
      cout << nroads << " road";
      if( nroads>1 ) cout << "s";
      cout << " successfully fit" << endl;
    }
    cout << "------------ end of projection  " << fName.Data() 
	 << "------------" << endl;
  }
#endif
  return GetNgoodRoads();
}


//_____________________________________________________________________________
#ifdef VERBOSE
static void PrintNode( const Node_t& node )
{
  const HitSet& hs = node.second;
  cout << " npl/nhits= " << hs.nplanes << "/" << hs.hits.size()
       << "  wnums= ";
  UInt_t ipl = 0;
  for( Hset_t::iterator ihit = hs.hits.begin(); ihit != hs.hits.end(); ) {
    while( ipl < (*ihit)->GetPlaneNum() ) { cout << "--/"; ++ipl; }
    bool seq = false;
    do {
      if( seq )  cout << " ";
      cout << (*ihit)->GetWireNum();
      seq = true;
      ++ihit;
    } while( ihit != hs.hits.end() and (*ihit)->GetPlaneNum() == ipl );
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
#endif

//_____________________________________________________________________________
inline
bool Projection::MostPlanes::operator() ( const Node_t& a, 
					  const Node_t& b ) const
{
  // Order patterns by decreasing number of active planes, then
  // decreasing number of hits, then ascending bin numbers
  if( a.second.nplanes > b.second.nplanes ) return true;
  if( a.second.nplanes < b.second.nplanes ) return false;
  if( a.second.hits.size() > b.second.hits.size() ) return true;
  if( a.second.hits.size() < b.second.hits.size() ) return false;
  return ( a.first < b.first );
}

//_____________________________________________________________________________
Int_t Projection::MakeRoads()
{
  // Combine patterns with common sets of hits into Roads.
  //
  // This is the primary de-ghosting algorithm. It finds clusters of patterns
  // that share active wires (hits).

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << fPatternsFound.size() << " patterns found:" << endl;
    for( NodeSet_t::iterator ipat = fPatternsFound.begin(); ipat !=
	   fPatternsFound.end(); ++ipat )
      PrintNode( *ipat );
  }
#endif

  for( NodeSet_t::iterator itn = fPatternsFound.begin(); itn !=
	 fPatternsFound.end(); ) {
    const Node_t& nd1 = *itn;

    // New roads must start with an unused pattern
//     if( nd1.second.used )
//       continue;
    assert( !nd1.second.used );

    // Start a new road with this pattern
    Road* rd = new( (*fRoads)[GetNroads()] ) Road(nd1,this);

    // Try to add the remaining cluster patterns to this road
    NodeSet_t::iterator itn2 = itn;
    // Stop search unless there are any remaining unused patterns
    itn = fPatternsFound.end();
    while( ++itn2 != fPatternsFound.end() ) {
      const Node_t& nd2 = *itn2;
      if( nd2.second.used )
	continue;
      // Test pattern and add it if applicable
      if( rd->Add(nd2) ) {
	// Pattern successfully added
	continue;
      }
      // Either the pattern is out of range or it does not fit in the current
      // cluster, so save current position as next unused pattern
      if( itn == fPatternsFound.end() )
	itn = itn2;
    }
    // Update the "used" flags of the road's component patterns
    rd->Finish();

    // If event display enabled, export the road's corner coordinates
    if( TestBit(kEventDisplay) ) {
      assert( fRoads->GetLast() == fRoadCorners->GetLast()+1 );
      new( (*fRoadCorners)[fRoads->GetLast()] ) Road::Corners(rd);
    }
  }

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << "Generated roads: " << endl;
    for( UInt_t i = 0; i < GetNroads(); ++i ) {
      const Road* rd = GetRoad(i);
      const Road::NodeList_t& ndlst = rd->GetPatterns();
      for( Road::NodeList_t::const_iterator ipat = ndlst.begin(); ipat !=
	     ndlst.end(); ++ipat ) {
	PrintNode( **ipat );
      }
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
    if( rd->IsGood() ) {
      if( rd->Fit() )
	// Count good roads (not void and good fit)
	++fNgoodRoads;
      else {
	// Void roads with bad fits
	rd->Void();
	changed = true;
#ifdef TESTCODE
	++n_badfits;
#endif
      }
    }
  }
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
	     "Using \"%s\".", fName.Data() );

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
  // Set wire angle (rad)
  
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
       << " type="  << GetType()
       << " npl="   << fPlanes.size()
       << " nlev="  << GetNlevels()
       << " zsize=" << GetZsize()
       << " maxsl=" << GetMaxSlope()
       << " width=" << GetWidth()
       << " angle=" << GetAngle()*TMath::RadToDeg();
  cout << endl;

  if( verbose > 0 ) {
    for( vwsiz_t i = 0; i < fPlanes.size(); ++i ) {
      WirePlane* wp = fPlanes[i];
      wp->Print(opt);
    }
  }
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

    // Found a match at the bottom of the pattern tree.
    // Add the pattern descriptor to the set of matches,
    // ordered by Projection::MostPlanes, i.e. most number of active planes
    // first, then most number of hits, then ascending by start bin number.

    // Make new node to insert into the set. (We can't use a map since we
    // want to sort using the HitSet data, too.
    Node_t node = make_pair( nd, HitSet() );
    // Collect all hits associated with the pattern's bins and save them
    // in the node's HitSet.
    for( UInt_t i = 0; i < fHitpattern->GetNplanes(); ++i ) {
      const vector<Hit*>& hits = fHitpattern->GetHits( i, node.first[i] );
      copy( hits.begin(), hits.end(), 
	    inserter( node.second.hits, node.second.hits.end() ));
    }
    assert( HitSet::GetMatchValue(node.second.hits) == match.first );
    node.second.plane_pattern = match.first;
    node.second.nplanes = match.second;

    // Insert the node. This copies the node by value.
    pair<NodeSet_t::iterator,bool> ins = fMatches->insert( node );

    assert(ins.second);  // duplicate matches should never happen
  }
  return kSkipChildNodes;
}

//_____________________________________________________________________________

}  // end namespace TreeSearch

ClassImp(TreeSearch::Projection)

///////////////////////////////////////////////////////////////////////////////
