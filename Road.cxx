//*-- Author :    Ole Hansen, Jefferson Lab   07-Feb-2008

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Road.h"
#include "TreeWalk.h"
#include "Projection.h"
#include "Hitpattern.h"
#include "Hit.h"
#include "WirePlane.h"
#include "TMath.h"
#include "TBits.h"
#include <iostream>
#include <algorithm>

using namespace std;

ClassImp(TreeSearch::Road);

namespace TreeSearch {

// Private class for use while building a Road
struct BuildInfo_t {
  set<Hit*>         fCommonHits;  // Hits common between all patterns
  vector<UShort_t>  fCornerBin;   // Bin numbers defining the corners
                                  // 0=LL, 1=LR, 2=UR, 3=UL 
  //TODO: use fMaxdist[] ??
  BuildInfo_t() : fCornerBin(4,0)
  { fCornerBin[0] = fCornerBin[3] = kMaxUShort; }
};

// Number of points for polygon test
static const vector<double>::size_type kNcorner = 5;

//_____________________________________________________________________________
Road::Road( const Projection* proj ) :
  fProjection(proj), fCornerX(kNcorner), fZL(0), fZU(0),
  fSlope(kBig), fPos(kBig), fChi2(kBig), fGood(false)
#ifdef TESTCODE
  , fSeed(0)
#endif
{
  // Constructor

  fErr[0] = fErr[1] = kBig;

  assert(fProjection);

  fBuild = new BuildInfo_t;
}

//_____________________________________________________________________________
Road::Road( const Road& orig ) :
  fProjection(orig.fProjection),
  fCornerX(orig.fCornerX), fZL(orig.fZL), fZU(orig.fZU),
  fPatterns(orig.fPatterns), fHits(orig.fHits),
  fSlope(orig.fSlope), fPos(orig.fSlope), fChi2(orig.fChi2), fGood(orig.fGood)
#ifdef TESTCODE
  , fSeed(orig.fSeed)
#endif
{
  // Copy constructor

  fErr[0] = orig.fErr[0];
  fErr[1] = orig.fErr[1];

  if( orig.fBuild )
    fBuild = new BuildInfo_t( *orig.fBuild );
  else
    fBuild = 0;
}

//_____________________________________________________________________________
Road& Road::operator=( const Road& rhs )
{
  // Print hit info

  if( this != &rhs ) {
#ifdef TESTCODE
    fSeed = rhs.fSeed;
#endif
    fProjection = rhs.fProjection;
    fCornerX    = rhs.fCornerX;
    fZL         = rhs.fZL;
    fZU         = rhs.fZU;
    fPatterns   = rhs.fPatterns;
    fHits       = rhs.fHits;
    fSlope      = rhs.fSlope;
    fPos        = rhs.fPos;
    fChi2       = rhs.fChi2;
    fErr[0]     = rhs.fErr[0];
    fErr[1]     = rhs.fErr[1];
    fGood       = rhs.fGood;

    delete fBuild;
    if( rhs.fBuild )
      fBuild = new BuildInfo_t( *rhs.fBuild );
    else
      fBuild = 0;
  }
  return *this;
}

//_____________________________________________________________________________
Road::~Road()
{
  // Destructor

  delete fBuild;

}

//_____________________________________________________________________________
#ifdef VERBOSE
static void PrintHits( const set<Hit*>& hits )
{
  //  cout << hits.size() << " hits" << endl;

  for( set<Hit*>::iterator it = hits.begin(); it != hits.end(); ++it ) {
    cout << " ";
    (*it)->Print();
  }
  
}
#endif
//_____________________________________________________________________________
Bool_t Road::CheckMatch( const set<Hit*>& hits ) const
{
  // Match evaluation function. Return true if the numbers of hits in each
  // plane contained in hitcount is sufficient for considering the road
  // cohesive.

  assert(fBuild);

  UInt_t curpat = 0;
  for( set<Hit*>::const_iterator it = hits.begin(); it != hits.end(); ++it )
    curpat |= 1U << ((*it)->GetWirePlane()->GetPlaneNum());

  return fProjection->GetPlaneCombos()->TestBitNumber(curpat);
}

//_____________________________________________________________________________
Bool_t Road::Add( Node_t& nd )
{
  // Check if the hits from the given NodeDescriptor pattern are common
  // with the common hit set already in this road. If so, add the pattern
  // to this road, update the common hits if necessary, and return true.
  // If not, do nothing and return false.
  //
  // Adding only works as long as the road is not yet finished

  if( !fBuild )
    return kFALSE;

#ifdef VERBOSE
  nd.first.Print();
  //nd.first.link->GetPattern()->Print(); nd.first.parent->Print();
  PrintHits(nd.second.hits);
#endif
  bool first = fPatterns.empty();
  if( first ) {
    if( !CheckMatch(nd.second.hits) )
      return kFALSE;
    fBuild->fCommonHits = fHits = nd.second.hits;
#ifdef TESTCODE
    fSeed = &nd;
#endif
  } else {
    set<Hit*> new_commons;
    set_intersection( nd.second.hits.begin(), nd.second.hits.end(),
		      fBuild->fCommonHits.begin(), fBuild->fCommonHits.end(),
		      inserter( new_commons, new_commons.end() ));

#ifdef VERBOSE
    cout << "new/old commons = " << new_commons.size() << " "
	  << fBuild->fCommonHits.size() << endl;
#endif
    assert( new_commons.size() <= fBuild->fCommonHits.size() );
    if( new_commons.size() < fBuild->fCommonHits.size() ) {
      // The set of common hits shrank, so we must check if this would still
      // be a good road
      if( !CheckMatch(new_commons) ) {
	// The new pattern would reduce the set of common hits in the road to
	// too loose a fit, so we reject the new pattern and leave
	// the road as it is
#ifdef VERBOSE
	cout << "failed" << endl;
#endif
	return kFALSE;
      }
      // The new set of common hits is good, so update the build data
      fBuild->fCommonHits.swap( new_commons );
    }
    set<Hit*> new_hits;
    set_union( fHits.begin(), fHits.end(),
	       nd.second.hits.begin(), nd.second.hits.end(),
	       inserter( new_hits, new_hits.begin() ));
#ifdef VERBOSE
    cout << "new/old nhits = " << new_hits.size() << " " 
 	 << fHits.size() << endl;
#endif
    if( new_hits.size() != fHits.size() ) {
      fHits.swap( new_hits );
#ifdef VERBOSE
      PrintHits( fHits );
#endif
    }
  }

  // Save a pointer to this pattern so we can update it later
  fPatterns.push_back(&nd);

  // Expand the road limits if necessary
  for( Int_t i = 0; i < 4; i++ ) {
    UShort_t bin = ( i<2 ) ? nd.first.Start() : nd.first.End(); // lower:upper
    if( i==1 || i==2 )
      // "right" edges
      fBuild->fCornerBin[i] = TMath::Max( bin, fBuild->fCornerBin[i] );
    else
      // "left" edges
      fBuild->fCornerBin[i] = TMath::Min( bin, fBuild->fCornerBin[i] );
  }
#ifdef VERBOSE
  cout << "new npat = " << fPatterns.size() << endl;
  cout << "new left/right = "
       << fBuild->fCornerBin[0]<<" "<< fBuild->fCornerBin[1]<<" "
       << fBuild->fCornerBin[3]<<" "<< fBuild->fCornerBin[2]<< endl;
#endif

  return kTRUE;
}

//_____________________________________________________________________________
void Road::Finish()
{
  // Finish building the road

  assert(fBuild);
  for( list<Node_t*>::iterator it = fPatterns.begin();
       it != fPatterns.end(); ++it ) {
    
    HitSet& hs = (**it).second;
    assert( hs.used < 2 ); // cannot add previously fully used pattern

    // TODO: search only up to first element not in common?
    list<Hit*> not_in_common;
    set_difference( hs.hits.begin(), hs.hits.end(),
		    fBuild->fCommonHits.begin(), fBuild->fCommonHits.end(),
		    back_inserter( not_in_common ) );

    hs.used = not_in_common.empty() ? 2 : 1;
#ifdef VERBOSE
    cout << "used = " << hs.used << " for ";
    (**it).first.Print();
#endif
  }

  // Calculate the corner coordinates
  const UInt_t np1 = fProjection->GetNplanes()-1;
  const UInt_t nl1 = fProjection->GetNlayers()-1;
  const Double_t resL = 2.0*fProjection->GetPlane(0)->GetResolution();
  const Double_t resU = 2.0*fProjection->GetPlane(np1)->GetResolution();
  fCornerX[0] = GetBinX( fBuild->fCornerBin[0] )   - resL;
  fCornerX[1] = GetBinX( fBuild->fCornerBin[1]+1 ) + resL;
  fCornerX[2] = GetBinX( fBuild->fCornerBin[2]+1 ) + resU;
  fCornerX[3] = GetBinX( fBuild->fCornerBin[3] )   - resU;

  const Double_t eps = 1e-3;   // z-shift to include the planes proper
  fZL = fProjection->GetPlaneZ(0) - eps;
  fZU = fProjection->GetPlaneZ(np1) + eps;
  assert( fZL < fZU );

  // Apply correction for difference of layer-z and plane-z
  Double_t dZlayer = fProjection->GetLayerZ(nl1) - fProjection->GetLayerZ(0);
  assert( dZlayer > 1e-3 );
  Double_t slopeL = ( fCornerX[3]-fCornerX[0] ) * (1.0/dZlayer);
  Double_t slopeR = ( fCornerX[2]-fCornerX[1] ) * (1.0/dZlayer);
  Double_t dZL = fProjection->GetLayerZ(0) - fZL;
  Double_t dZU = fZU - fProjection->GetLayerZ(nl1);
  assert( dZL >= 0.0 && dZU >= 0.0 );
  fCornerX[0] -= slopeL * dZL;
  fCornerX[1] -= slopeR * dZL;
  fCornerX[2] += slopeR * dZU;
  fCornerX[3] += slopeL * dZU;

  fCornerX[4] = fCornerX[0];

  assert( fCornerX[0] < fCornerX[1] );
  assert( fCornerX[3] < fCornerX[2] );

  // All done. Put the tools away
  delete fBuild; fBuild = 0;

  return;
}

//_____________________________________________________________________________
Double_t Road::GetBinX( UInt_t bin ) const
{
  // Get X coordinate of left edge of given bin

  const Hitpattern* hpat = fProjection->GetHitpattern();
  return
    static_cast<Double_t>(bin) * hpat->GetBinWidth() - hpat->GetOffset();

}

//_____________________________________________________________________________
void Road::Print( Option_t* opt ) const
{
  // Print road info

}

//_____________________________________________________________________________
Bool_t Road::CollectCoordinates( vector<Point>& points,
				 vector<vector<Point*> >& planepoints )
{
  // Gather hit positions that lie within the Road area

  assert( fCornerX.size() == kNcorner );
  assert( static_cast<UInt_t>(planepoints.size())==fProjection->GetNplanes());

#ifdef VERBOSE
  cout << "Collecting coordinates from: (" << fPatterns.size() << " patterns)"
       << endl;
  PrintHits( fHits );
#ifdef TESTCODE
  cout << "Seed pattern: " << endl;
  fSeed->first.Print();
#endif
#endif
  Double_t zp[kNcorner] = { fZL, fZL, fZU, fZU, fZL };
  UInt_t matchpattern = 0;
  Bool_t good = true;

  // Collect all hit coordinates that lie within the Road
  for( set<Hit*>::iterator it = fHits.begin(); it != fHits.end(); ++it ) {
    const Hit* hit = *it;
    Double_t z = hit->GetZ();
    UInt_t np = hit->GetWirePlane()->GetPlaneNum();
    for( int i = 2; i--; ) {
      Double_t x = i ? hit->GetPosL() : hit->GetPosR();
      if( TMath::IsInside(x, z, kNcorner, &fCornerX[0], zp) ) {
	points.push_back( Point(x,z) );
	planepoints[np].push_back( &points.back() );
	matchpattern |= 1U<<np;
      }
    }
  }
  // Check if this matchpattern is acceptable
  good = fProjection->GetPlaneCombos()->TestBitNumber(matchpattern);

#ifdef VERBOSE
  cout << "Collected:" << endl;
  for( vector<vector<Point*> >::size_type i = planepoints.size(); i; ) {
    --i;
    cout << " pl= " << i;
    if( planepoints[i].empty() )
      cout << " missing";
    else {
      Double_t z = planepoints[i][0]->z;
      cout << " z=" << z << "\t x=";
      for( vector<Point*>::size_type j = 0; j < planepoints[i].size(); ++j ) {
	cout << " " << planepoints[i][j]->x;
	assert( planepoints[i][j]->z == z );
      }
    }
    cout << endl;
  }
  if( !good )
    cout << "REJECTED" << endl;
#endif
  return good;
}

//_____________________________________________________________________________
Bool_t Road::Fit()
{

  fGood = false;
  vector<Point> points;
  vector<vector<Point*> > planepoints( fProjection->GetNplanes() );
  points.reserve( 2*fHits.size() );
  if( !CollectCoordinates(points, planepoints) )
    return false;

  // Determine number of permutations
  // Loop over permutations of hits in the planes
  // Check each permutation for consistency with at least one pattern
  // Fit the permutation
  // Save best Chi2

  return (fGood = true);
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
