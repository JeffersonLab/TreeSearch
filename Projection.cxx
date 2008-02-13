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

#include "TMath.h"
#include "TString.h"
#include "TBits.h"

#include <iostream>
#include <sys/time.h>  // for timing
#include <algorithm>
#include <utility>

using namespace std;

namespace TreeSearch {

typedef vector<WirePlane*>::size_type vwsiz_t;
typedef vector<WirePlane*>::iterator  vwiter_t;
typedef vector<Road>::size_type vrsiz_t;

//_____________________________________________________________________________
Projection::Projection( Int_t type, const char* name, Double_t angle,
			THaDetectorBase* parent )
  : THaAnalysisObject( name, name ), fType(type), fNlevels(0),
    fMaxSlope(0.0), fWidth(0.0),
    fHitpattern(0), fPatternTree(0), fDetector(parent),
    fPlaneCombos(0), fLayerCombos(0)
{
  // Constructor

  // angle is the default value. It can be overridden via a database entry.
  SetAngle( angle );

  fTitle.Append(" projection");
}


//_____________________________________________________________________________
Projection::Projection( const Projection& orig )
  : fType(orig.fType), fPlanes(orig.fPlanes), fLayers(orig.fLayers),
    fNlevels(orig.fNlevels), fMaxSlope(orig.fMaxSlope), fWidth(orig.fWidth),
    fSinAngle(orig.fSinAngle), fCosAngle(orig.fCosAngle),
    fHitpattern(0), fPatternTree(0), fDetector(orig.fDetector),
    fPatternsFound(), fRoads()
{
  // Copying

  if( orig.fHitpattern )
    fHitpattern = new Hitpattern(*orig.fHitpattern);
//   if( orig.fPatternTree )
//     fPatternTree = new PatternTree(*orig.fPatternTree);
  if( orig.fPlaneCombos )
    fPlaneCombos = new TBits(*orig.fPlaneCombos);
  if( orig.fLayerCombos )
    fLayerCombos = new TBits(*orig.fLayerCombos);
}

//_____________________________________________________________________________
const Projection& Projection::operator=( const Projection& rhs )
{
  // Assignment

  if( this != &rhs ) {
    fType     = rhs.fType;
    fNlevels    = rhs.fNlevels;
    fPlanes   = rhs.fPlanes;
    fLayers   = rhs.fLayers;
    fMaxSlope = rhs.fMaxSlope;
    fWidth    = rhs.fWidth;
    fSinAngle = rhs.fSinAngle;
    fCosAngle = rhs.fCosAngle;
    fDetector = rhs.fDetector;
    delete fHitpattern;
    if( rhs.fHitpattern )
      fHitpattern = new Hitpattern(*rhs.fHitpattern);
    else
      fHitpattern = 0;
//     if( rhs.fPatternTree )
//       fPatternTree = new PatternTree(*rhs.fPatternTree);
//     else
    fPatternTree = 0;
    fPatternsFound.clear();
    fRoads.clear();
    delete fPlaneCombos;
    if( rhs.fPlaneCombos )
      fPlaneCombos = new TBits(*rhs.fPlaneCombos);
    else
      fPlaneCombos = 0;
    if( rhs.fLayerCombos )
      fLayerCombos = new TBits(*rhs.fLayerCombos);
    else
      fLayerCombos = 0;
  }
  return *this;
}

//_____________________________________________________________________________
Projection::~Projection()
{
  // Destructor

  delete fPatternTree;
  delete fHitpattern;
  delete fPlaneCombos;
  delete fLayerCombos;
}

//_____________________________________________________________________________
void Projection::AddPlane( WirePlane* wp, WirePlane* partner )
{
  // Add wire plane wp (and optional partner plane) to this projection. 
  // Sets plane and layer numbers.

  assert(wp);

  wp->SetPlaneNum( fPlanes.size() );
  fPlanes.push_back( wp );
  if( partner ) {
    assert( partner->GetZ() > wp->GetZ() ); // Planes must be ordered
    partner->SetPlaneNum( fPlanes.size() );
    fPlanes.push_back( partner );
  }

  // Only add the primary plane to the vector of layers; partner planes are
  // linked within the plane objects (get via wp->GetPartner())
  wp->SetLayerNum( fLayers.size() );
  fLayers.push_back( wp );
  if( partner )
    partner->SetLayerNum( wp->GetLayerNum() );
}

//_____________________________________________________________________________
void Projection::Clear( Option_t* opt )
{    
  // Clear event-by-event data

  if( fHitpattern )
    fHitpattern->Clear();

  fPatternsFound.clear();
  fRoads.clear();

#ifdef TESTCODE
  t_treesearch = t_roads = t_track = 0.0;
#endif
}

//_____________________________________________________________________________
Int_t Projection::Decode( const THaEvData& evdata )
{
  // Decode all planes belonging to this projection

  Int_t sum = 0;
  for( vwsiz_t i = 0; i < fLayers.size(); ++i ) {
    WirePlane* wp = fLayers[i];
    Int_t nhits = wp->Decode( evdata );
    sum += nhits;
    if( (wp = wp->GetPartner()) ) {
      Int_t nhits2 = wp->Decode( evdata );
      sum += nhits2;
      if( nhits2 > nhits )
	nhits = nhits2;
    }
    //TODO: nhits now holds the occupancy of this plane or plane pair,
    // use it for occupancy test
  }
  return sum;
}

//_____________________________________________________________________________
void Projection::Reset()
{
  // Reset parameters, clear list of planes, delete Hitpattern
  
  fIsInit = kFALSE;
  fPlanes.clear();
  fLayers.clear();
  fMaxSlope = fWidth = 0.0;
  delete fHitpattern; fHitpattern = 0;
  delete fPatternTree; fPatternTree = 0;
}

//_____________________________________________________________________________
Double_t Projection::GetLayerZ( UInt_t i ) const
{
  // Return the effective z-position of the i-th wire plane layer.  In case
  // of a partnered plane pair, the effective z is halfway between the two
  // planes

  assert( i<fLayers.size() );
  WirePlane* wp = fLayers[i];
  Double_t zpos = wp->GetZ();
  if( wp->GetPartner() )
    zpos = 0.5*( zpos + wp->GetPartner()->GetZ() );

  return zpos;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Projection::InitLevel2( const TDatime& date )
{
  // Level-2 initialization - load pattern database and initialize hitpattern

  TreeParam_t tp;
  tp.maxdepth = fNlevels-1;
  tp.width = fWidth;
  tp.maxslope = fMaxSlope;
  for( UInt_t i = 0; i < GetNlayers(); ++i )
    tp.zpos.push_back( GetLayerZ(i) );

  if( tp.Normalize() != 0 )
    return fStatus = kInitError;

  // Attempt to read the pattern database from file
  delete fPatternTree;
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
  delete fHitpattern;
  fHitpattern = new Hitpattern( *fPatternTree );
  if( !fHitpattern || fHitpattern->IsError() )
    return fStatus = kInitError;

  // Determine maximum search distance (in bins) for combining patterns
  static const Double_t sqrt2 = TMath::Sqrt(2.0);
  fMaxdist.clear();
  assert( GetNlayers() == fHitpattern->GetNplanes() );
  for( UInt_t i = 0; i < GetNlayers(); ++i ) {
    UInt_t dist;
    if( fLayers[i]->GetPartner() ) {
      // Maximum number of bins that a single hit can cover in the midplane
      // of a plane pair
      WirePlane* wpl[2] = { fLayers[i], fLayers[i]->GetPartner() };
      Double_t dx = 0;
      Double_t dz = TMath::Abs( wpl[1]->GetZ() - wpl[0]->GetZ() );
      for( Int_t i = 0; i < 2; ++i ) {
	WirePlane* wp = wpl[i];
	Double_t dx1 = 
	  sqrt2 * wp->GetResolution()  + fMaxSlope * dz  + wp->GetMaxLRdist();
	if( dx1 > dx )
	  dx = dx1;
      }
      dist = TMath::FloorNint( dx * fHitpattern->GetBinScale() );
    } else {
      // Maximum distance of bins that can belong to the same hit in a single
      // plane
      dist = TMath::FloorNint( fHitpattern->GetBinScale() *
        (fLayers[i]->GetMaxLRdist() + 2.0*fLayers[i]->GetResolution()) );
    }    
    fMaxdist.push_back( dist );
  }

  // Set up the lookup table indicating which planes are allowed to have 
  // missing hits. The bit pattern of plane hits is used as an index into
  // this table; if the corresponding bit is set, the plane combo is allowed.

  delete fPlaneCombos;
  UInt_t np = 1U<<GetNplanes();
  fPlaneCombos = new TBits( np );
  // For starters, allow any one plane to be missing
  // TODO: read from database and allow more complex patterns
  UInt_t allbits = np-1;
  fPlaneCombos->SetBitNumber( allbits );
  while(np>>=1)
    fPlaneCombos->SetBitNumber( allbits ^ np );

  delete fLayerCombos;
  UInt_t nl = 1U<<GetNlayers();
  fLayerCombos = new TBits( nl );
  // Allow any one layer to be missing. However, layers made up of 
  // partnered planes must always be present
  //TODO: read from database
  allbits = nl-1;
  fLayerCombos->SetBitNumber( allbits );
//   for( UInt_t i=GetNlayers(); i; ) {
//     nl >>= 1;
//     WirePlane* wp = fLayers[--i];
//     if( !wp->GetPartner() )
//       fLayerCombos->SetBitNumber( allbits ^ nl );
//   }

  return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t Projection::ReadDatabase( const TDatime& date )
{
  // Read parameters from database

  static const char* const here = "ReadDatabase";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  //FIXME: move elsewhere?
  Reset();

  Double_t angle = kBig;
  const DBRequest request[] = {
    { "angle",        &angle,        kDouble, 0, 1 },
    { "maxslope",     &fMaxSlope,    kDouble, 0, 1, -1 },
    { "search_depth", &fNlevels,     kUInt,   0, 0, -1 },
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
    { "t_treesearch", "Time in TreeSearch (us)", "t_treesearch" },
    { "t_roads", "Time in MakeRoads (us)", "t_roads" },
    { "t_track", "Time in Track (us)", "t_track" },
    //    { "", "", "" },
#endif
    { 0 }
  };
  DefineVarsFromList( vars, mode );
  return 0;
}

//_____________________________________________________________________________
Int_t Projection::FillHitpattern()
{
  // Fill this projection's hitpattern with hits from the wire planes.
  // Returns the total number of hits processed (where hit pairs on plane
  // and partner plane count as one).

  Int_t ntot = 0;
  for( vwiter_t it = fLayers.begin(); it != fLayers.end(); ++it ) {
    ntot += fHitpattern->ScanHits( *it, (*it)->GetPartner() );
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

  ComparePattern compare( fHitpattern, fLayerCombos, &fPatternsFound );
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

  // Combine patterns with common sets of hits into Roads
  MakeRoads();

#ifdef TESTCODE
  gettimeofday(&stop, 0 );
  timersub( &stop, &start2, &diff );
  t_roads = 1e6*(Double_t)diff.tv_sec + (Double_t)diff.tv_usec;

  //  n_roads = fRoads.size();

  gettimeofday( &start2, 0 );
#endif

  // RoadFit:
  // Fit hits within each road. Store fit parameters with Road. 
  // Store list of best-fit hits with Road

  //TODO: put in separate routine since the fit will need to be repeated

  // FilterGhosts:
  // Eliminate roads with high chi^2 and roads that appear to be essentially
  // identical to others

  //TODO...

#ifdef TESTCODE
  gettimeofday(&stop, 0 );
  timersub( &stop, &start, &diff );
  t_track = 1e6*(Double_t)diff.tv_sec + (Double_t)diff.tv_usec;
#endif

  return 0;
}


//_____________________________________________________________________________
Int_t Projection::MakeRoads()
{
  // Combine patterns with common sets of hits into Roads.
  //
  // This is the primary de-ghosting algorithm. It groups patterns with
  // common wires (not common hit positions!) in all planes together.

  map<const NodeDescriptor,HitSet>::iterator it1, it2, ref_it1;
  for( it1 = ref_it1 = fPatternsFound.begin(); it1 != fPatternsFound.end(); 
       ++it1 ) {
    pair<const NodeDescriptor,HitSet>& nd1 = *it1;
    assert(nd1.second.used < 3);
    // New roads must contain at least one unused pattern
    if( nd1.second.used )
      continue;
    Road rd(this);
    // Adding the first pattern must make a good match
    if( !rd.Add(nd1) ) {
#ifdef VERBOSE
      cout << ">>>>>>>>> No match, skipped" << endl;
#endif
      continue;
    }
    it2 = ref_it1;
    // Search until end of list or too far right
    while( ++it2 != fPatternsFound.end() and
	   (*it2).first[0] <= nd1.first[0] + fMaxdist[0] ) {
      pair<const NodeDescriptor,HitSet>& nd2 = *it2;
      if( nd1.first[0] <= nd2.first[0] + fMaxdist[0] ) {
	// Try adding unused or partly used pattern to the new road until there
	// are no more patterns that could possibly have common hits.
	if( nd2.second.used < 2 &&
	    it1 != it2 ) // skip the seed in case we run over it
	  rd.Add( nd2 );
      } else {
	// Save last position too far left of it1 (seed of road).
	// This + 1 is where we start the next search.
	ref_it1 = it2;
      }
    }
    // Update the "used" flags of the road's component patterns
    rd.Finish();
    fRoads.push_back(rd);
  }
#ifdef VERBOSE
  if( !fRoads.empty() ) {
    cout << fRoads.size() << " road";
    if( fRoads.size()>1 ) cout << "s";
    cout << endl;
  }
  cout << "------------ end of projection  " << fName.Data() 
       << "------------" << endl;
#endif
  return 0;
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
  
  fSinAngle = TMath::Sin( angle );
  fCosAngle = TMath::Cos( angle );
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
       << " nlay="  << fLayers.size()
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
  UInt_t matchpattern = fHitpattern->ContainsPattern(nd);
  if( fLayerCombos->TestBitNumber(matchpattern)  ) {
    if( nd.depth < fHitpattern->GetNlevels()-1 )
      return kRecurse;

    // Found a match at the bottom of the pattern tree.
    // Add the pattern descriptor to the set of matches,
    // ordered by start bin number (NodeDescriptor::operator<)
    pair<map<const NodeDescriptor,HitSet>::iterator,bool> ins = 
      fMatches->insert( make_pair(nd,HitSet()) );
    assert(ins.second);  // duplicate matches should never happen

    // Retrieve the node that was just inserted.
    pair<const NodeDescriptor,HitSet>& node = *ins.first;

    // Collect all hits associated with the pattern's bins and save them
    // in the node's HitSet
    assert( node.second.hits.empty() );
    for( UInt_t i = fHitpattern->GetNplanes(); i; ) {
      --i;
      const vector<Hit*>& hits = fHitpattern->GetHits( i, node.first[i] );
      copy( hits.begin(), hits.end(), 
	    inserter( node.second.hits, node.second.hits.begin() ));
    }

  }
  return kSkipChildNodes;
}

//_____________________________________________________________________________

}  // end namespace TreeSearch

ClassImp(TreeSearch::Projection)

///////////////////////////////////////////////////////////////////////////////
