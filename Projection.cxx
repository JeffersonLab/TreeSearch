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
#include "TString.h"
#include "PatternTree.h"
#include "PatternGenerator.h"
#include "TreeWalk.h"
#include <iostream>
#include <sys/time.h>  // for timing

using namespace std;

namespace TreeSearch {

typedef vector<WirePlane*>::size_type vwsiz_t;
typedef vector<WirePlane*>::iterator  vwiter_t;

//_____________________________________________________________________________
Projection::Projection( Int_t type, const char* name, Double_t angle,
			THaDetectorBase* parent )
  : THaAnalysisObject( name, name ), fType(type), fNlevels(0),
    fMaxSlope(0.0), fWidth(0.0),
    fHitpattern(0), fPatternTree(0), fDetector(parent)
{
  // Constructor

  // angle is the default value. It can be overridden via a database entry.
  SetAngle( angle );

  fTitle.Append(" projection");
}


//_____________________________________________________________________________
Projection::Projection( const Projection& orig )
  : fType(orig.fType), fPlanes(orig.fPlanes), fNlevels(orig.fNlevels),
    fMaxSlope(orig.fMaxSlope), fWidth(orig.fWidth),
    fSinAngle(orig.fSinAngle), fCosAngle(orig.fCosAngle),
    fHitpattern(0), fPatternTree(0), fDetector(orig.fDetector)
{
  // Copying

  if( orig.fHitpattern )
    fHitpattern = new Hitpattern(*orig.fHitpattern);
//   if( orig.fPatternTree )
//     fPatternTree = new PatternTree(*orig.fPatternTree);
}

//_____________________________________________________________________________
const Projection& Projection::operator=( const Projection& rhs )
{
  // Assignment

  if( this != &rhs ) {
    fType     = rhs.fType;
    fNlevels    = rhs.fNlevels;
    fPlanes   = rhs.fPlanes;
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
  }
  return *this;
}

//_____________________________________________________________________________
Projection::~Projection()
{
  // Destructor

  delete fPatternTree; fPatternTree = 0;
  delete fHitpattern; fHitpattern = 0;
}

//_____________________________________________________________________________
void Projection::AddPlane( WirePlane* wp )
{
  // Add wire plane to this projection

  if( !wp )
    return;

//   const WirePlane* planes[] = { wp, wp->GetPartner(), 0 };
//   const WirePlane** p = planes;
//   while( *p ) {
//     Double_t z = (*p)->GetZ();
//     if( z < fZmin )
//       fZmin = z;
//     if( z > fZmax )
//       fZmax = z;
//     ++p;
//   }
  
  // Only add the primary plane; partner planes are linked within the
  // plane objects
  fPlanes.push_back( wp );
}

//_____________________________________________________________________________
void Projection::Clear( Option_t* opt )
{    
  // Clear event-by-event data

  if( fHitpattern )
    fHitpattern->Clear();

  fPatternsFound.clear();

  //FIXME: TEST
  search_time = 0.0;
}

//_____________________________________________________________________________
Int_t Projection::Decode( const THaEvData& evdata )
{
  // Decode all planes belonging to this projection

  Int_t sum = 0;
  for( vwsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    WirePlane* wp = fPlanes[iplane];
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
  fMaxSlope = fWidth = 0.0;
  delete fHitpattern; fHitpattern = 0;
  delete fPatternTree; fPatternTree = 0;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Projection::InitLevel2( const TDatime& date )
{
  // Level-2 initialization - load pattern database and initialize hitpattern

  TreeParam_t tp;
  tp.maxdepth = fNlevels-1;
  tp.width = fWidth;
  tp.maxslope = fMaxSlope;
  for( vwiter_t it = fPlanes.begin(); it != fPlanes.end(); ++it )
    tp.zpos.push_back( (*it)->GetZ() );

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
  fHitpattern = new Hitpattern( fNlevels, GetNplanes(), fWidth );
  if( !fHitpattern || fHitpattern->IsError() )
    return fStatus = kInitError;

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

  // Register variables in global list
  
  RVarDef vars[] = {
    //{ "n_test",   "Number of pattern comparisons per event", "n_test"     },
    //{ "n_found",  "Number of patterns found per event",      "n_found"    },
    { "t_search", "Search time per event (us)", "search_time"     },
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
  for( vwiter_t it = fPlanes.begin(); it != fPlanes.end(); ++it ) {
    ntot += fHitpattern->ScanHits( *it, (*it)->GetPartner() );
  }
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

  //FIXME: test
  struct timeval start, stop, diff;
  gettimeofday( &start, 0 );

  fPatternsFound.clear();

  ComparePattern compare( fHitpattern, &fPatternsFound );
  TreeWalk walk( fNlevels );
  walk( fPatternTree->GetRoot(), compare );

  //FIXME: test
  gettimeofday(&stop, 0 );
  timersub( &stop, &start, &diff );
  search_time = 1e6*(Double_t)diff.tv_sec + (Double_t)diff.tv_usec;

  cout << fPatternsFound.size() << endl;

  // TreeCombine:
  // Combine patterns with common sets of hits into Roads

  //TODO...

  // RoadFit:
  // Fit hits within each road. Store fit parameters with Road. 
  // Store list of best-fit hits with Road

  //TODO: put in separate routine since the fit will need to be repeated

  // FilterGhosts:
  // Eliminate roads with high chi^2 and roads that appear to be essentially
  // identical to others

  //TODO...

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
      wp = wp->GetPartner();
      if( wp )
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

  // Match?
  //  fProj->n_test++;  //FIXME: test
  if( fHitpattern->ContainsPattern(nd) == fHitpattern->GetNplanes() ) {
    if( nd.depth < fHitpattern->GetNlevels()-1 )
      return kRecurse;

    // Found a match at the bottom of the tree. Add this match to the list
    // of results, ordered by the start bin # (NodeDescriptor::operator<())
    fMatches->insert( nd );
  }
  return kSkipChildNodes;
}

//_____________________________________________________________________________

}  // end namespace TreeSearch

ClassImp(TreeSearch::Projection)

///////////////////////////////////////////////////////////////////////////////
