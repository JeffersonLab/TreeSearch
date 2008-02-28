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
#include "Helper.h"
#include "TMath.h"
#include "TBits.h"
#include <iostream>
#include <algorithm>
#include <numeric>

using namespace std;

ClassImp(TreeSearch::Road)

namespace TreeSearch {

// Private class for use while building a Road
struct BuildInfo_t {
  Hset_t            fClusterHits;  // Hits common between all patterns
  vector<UShort_t>  fCornerBin;    // Bin numbers defining the corners
                                   // 0=LL, 1=LR, 2=UR, 3=UL 
  BuildInfo_t() : fCornerBin(4,0)
  { fCornerBin[0] = fCornerBin[3] = kMaxUShort; }
};

// Number of points for polygon test
static const vector<double>::size_type kNcorner = 5;
static const UInt_t kMaxNhitCombos = 1000;

typedef vector<Road::Point*> Pvec_t;
typedef Hset_t::iterator siter_t;
#define ALL(c) (c).begin(), (c).end()

//_____________________________________________________________________________
Road::Road( const Projection* proj ) 
  : TObject(), fProjection(proj), fCornerX(kNcorner), fZL(kBig), fZU(kBig), 
    fPos(kBig), fSlope(kBig), fChi2(kBig), fDof(kMaxUInt), fGood(false),
    fBuild(0)
{
  // Constructor

  assert(fProjection);  // Invalid Projection* pointer

  fBuild = new BuildInfo_t;
}

//_____________________________________________________________________________
Road::Road( const Road& orig ) : 
  TObject(orig), fProjection(orig.fProjection),
  fCornerX(orig.fCornerX), fZL(orig.fZL), fZU(orig.fZU),
  fPatterns(orig.fPatterns), fHits(orig.fHits),
  fPos(orig.fPos), fSlope(orig.fSlope), fChi2(orig.fChi2),
  fDof(orig.fDof), fGood(orig.fGood)
{
  // Copy constructor

  if( !orig.fFitData.empty() ) {
    fFitData.reserve( orig.fFitData.size() );
    for( vector<FitResult*>::const_iterator it = orig.fFitData.begin();
	 it != orig.fFitData.end(); ++it )
      fFitData.push_back( new FitResult( **it ));
  }

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
    TObject::operator=(rhs);
    fProjection = rhs.fProjection;
    fCornerX    = rhs.fCornerX;
    fZL         = rhs.fZL;
    fZU         = rhs.fZU;
    fPatterns   = rhs.fPatterns;
    fHits       = rhs.fHits;
    fPos        = rhs.fPos;
    fSlope      = rhs.fSlope;
    fChi2       = rhs.fChi2;
    fDof        = rhs.fDof;
    fGood       = rhs.fGood;
    DeleteContainer( fFitData );
    if( !rhs.fFitData.empty() ) {
      fFitData.reserve( rhs.fFitData.size() );
      for( vector<FitResult*>::const_iterator it = rhs.fFitData.begin();
	   it != rhs.fFitData.end(); ++it )
	fFitData.push_back( new FitResult( **it ));
    }
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

  DeleteContainer( fFitData );
  delete fBuild;

}

//_____________________________________________________________________________
#ifdef VERBOSE
static void PrintHits( const Hset_t& hits )
{
  //  cout << hits.size() << " hits" << endl;

  for( siter_t it = hits.begin(); it != hits.end(); ++it ) {
    cout << " ";
    (*it)->Print();
  }
  
}
#endif
//_____________________________________________________________________________
Bool_t Road::CheckMatch( const Hset_t& hits ) const
{
  // Return true if the hits from the given set either cover all planes
  // or, if planes are missing, the pattern of missing planes is allowed
  // (based on what level of matching the user requests via the database)

  assert(fBuild);  // Not in build mode

  UInt_t curpat = 0;
  for( Hset_t::const_iterator it = hits.begin(); it != hits.end(); ++it )
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

  assert(fBuild);  // Add() again after Finish() not allowed

#ifdef VERBOSE
  cout << "Adding:" << endl;
  nd.first.Print();
  //nd.first.link->GetPattern()->Print(); nd.first.parent->Print();
  PrintHits(nd.second.hits);
#endif
  bool first = fPatterns.empty();
  if( first ) {
    // Check the hits of the first pattern of a new road for missing planes
    if( !CheckMatch(nd.second.hits) )
      return kFALSE;
    // A new road starts out with common hits and all hits that are the same
    fBuild->fClusterHits = fHits = nd.second.hits;
#ifdef VERBOSE
    cout << "New cluster:" << endl;
    PrintHits( fBuild->fClusterHits );
#endif
  } else {
    // Accept this pattern if all its hits are already in the cluster
    if( !includes( ALL(fBuild->fClusterHits), ALL(nd.second.hits),
		   nd.second.hits.key_comp() )) {

      // Check how the new pattern fits into the current cluster

      // Create an instance of the wire number distance comparison function.
      // This allows clustering of wires that are up to maxdist apart.
      Hit::WireDistLess comp( fProjection->GetHitMaxDist() );

      // Get the intersection of the new hits with the hits that are currently
      // clustered together in this road.
      // NB1: Because this comparison function allows for elements to be
      // equivalent even if not equal, the result may contain hits that
      // do not occur in both sets!
      // NB2: We don't use std::set_intersection here because its idea of an
      // intersection of identical elements is not what we need here.
      // NB3: Neither input set is ordered according to comp(), but that is
      // ok here since ordering by Hit::WireNumLess implies ordering by
      // Hit::WireDistLess.
      Hset_t intersection;
      siter_t itnew = nd.second.hits.begin();
      siter_t end1  = nd.second.hits.end();
      siter_t itold = fBuild->fClusterHits.begin();
      siter_t end2  = fBuild->fClusterHits.end();
      while( itnew != end1 and itold != end2 ) {
	if( comp(*itnew, *itold) )
	  ++itnew;
	else if( comp(*itold, *itnew) )
	  ++itold;
	else {
	  intersection.insert( intersection.end(), *itnew );
	  ++itnew;
	}
      }
#ifdef VERBOSE
      cout << "candidate/intersecting hits = " << nd.second.hits.size() << " "
	   << intersection.size() << endl;
      cout << "Intersection:" << endl;
      PrintHits( intersection );
#endif
      assert( intersection.size() <= nd.second.hits.size() );
      if( intersection.size() < nd.second.hits.size() ) {
	if( !CheckMatch(intersection) ) {
	  // The new pattern would reduce the set of common hits in the road to
	  // too loose a match, so we reject the new pattern and leave
	  // the road as it is
#ifdef VERBOSE
	  cout << "no longer a good plane pattern" << endl;
#endif
	  return false;
	}
      }
      // Add the intersection to the prior hits, accumulating the cluster
      fBuild->fClusterHits.insert( ALL(intersection) );
#ifdef VERBOSE
      cout << "New cluster:" << endl;
      PrintHits( fBuild->fClusterHits );
#endif
      // Add the hits of the new pattern to the set of all hits of this road
      Hset_t::size_type oldsize = fHits.size();
      fHits.insert( ALL(nd.second.hits) );
#ifdef VERBOSE
      cout << "new/old nhits = " << fHits.size() << " " 
	   << oldsize << endl;
#endif
      if( fHits.size() != oldsize ) {
#ifdef VERBOSE
	PrintHits( fHits );
#endif
      }
    } 
 #ifdef VERBOSE
    else
      cout << "All hits already in cluster" << endl;
#endif
  }

  // Save a pointer to this pattern so we can update it later
  fPatterns.push_back(&nd);

  // Expand the road limits if necessary
  fBuild->fCornerBin[0] = TMath::Min( nd.first.Start(), fBuild->fCornerBin[0]);
  fBuild->fCornerBin[1] = TMath::Max( nd.first.Start(), fBuild->fCornerBin[1]);
  fBuild->fCornerBin[2] = TMath::Max( nd.first.End(), fBuild->fCornerBin[2] );
  fBuild->fCornerBin[3] = TMath::Min( nd.first.End(), fBuild->fCornerBin[3] );

#ifdef VERBOSE
  cout << "new npat = " << fPatterns.size() << endl;
  cout << "new left/right = "
       << fBuild->fCornerBin[0]<<" "<< fBuild->fCornerBin[1]<<" "
       << fBuild->fCornerBin[3]<<" "<< fBuild->fCornerBin[2]<< endl;
#endif

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t Road::Adopt( const Road* other )
{
  // Widen the road if necessary to include the other road.
  // This only affects the fCornerX data of this road.
  // Returns false and does nothing if the roads are not compatible.

  assert( other && !fBuild );
  if( fBuild )
    return false;

  // Verify that z-limits agree and x-ranges overlap
  const Double_t eps = 1e-6;
  if( TMath::Abs(fZL-other->fZL) > eps or TMath::Abs(fZU-other->fZU) > eps or
      other->fCornerX[1] < fCornerX[0] or fCornerX[1] < other->fCornerX[0] or
      other->fCornerX[2] < fCornerX[3] or fCornerX[2] < other->fCornerX[3] )
    return false;
  
  fCornerX[0] = TMath::Min( fCornerX[0], other->fCornerX[0] );
  fCornerX[1] = TMath::Max( fCornerX[1], other->fCornerX[1] );
  fCornerX[2] = TMath::Max( fCornerX[2], other->fCornerX[2] );
  fCornerX[3] = TMath::Min( fCornerX[3], other->fCornerX[3] );
  fCornerX[4] = fCornerX[0];

  return true;
}

//_____________________________________________________________________________
void Road::Finish()
{
  // Finish building the road

  assert(fBuild);   // Road must be incomplete to be able to Finish()
  for( list<Node_t*>::iterator it = fPatterns.begin();
       it != fPatterns.end(); ++it ) {
    
    HitSet& hs = (**it).second;
    assert( hs.used < 2 ); // cannot add previously fully used pattern

    bool hs_in_cluster = includes( ALL(fBuild->fClusterHits), ALL(hs.hits),
				   hs.hits.key_comp() );
    hs.used = hs_in_cluster ? 2 : 1;
#ifdef VERBOSE
    cout << "used = " << hs.used << " for ";
    (**it).first.Print();
#endif
  }

  // Calculate the corner coordinates. We are building a polygon 
  // with points in the order LL (lower left), LR, UR, UL, LL, as 
  // needed for TMath::IsInside().
  // NB: To include points just at the edge of a road, the width is
  // increased by 2*(average hit position resolution of the plane)
  // on each the left and the right sides.
  // Also, the upper and lower edges are shifted up and down, respectively,
  // so that the points on the planes are guaranteed to be included.
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
  // If the first or last plane are part of a layer, then the x-coordinates
  // obtained from GetBinX above refers to the layer, not the plane. In that
  // case, we project the polygon's sides to the correct z-position (plane z)
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
Bool_t Road::Includes( const Road* other ) const
{
  // Check if all the hits in the given other road are also present in this one

  assert( other && !fBuild );

  return includes( ALL(fHits), ALL(other->fHits), fHits.key_comp() );
}

//_____________________________________________________________________________
Bool_t Road::CollectCoordinates( vector<Pvec_t>& points ) const
{
  // Gather hit positions that lie within the Road area.
  // Return true if the plane occupancy pattern of the selected points 
  // is allowed by Projection::fPlaneCombos, otherwise false.
  // If true is returned, the caller must delete the elements of the 
  // returned vector, otherwise the vector is empty.

  assert( fCornerX.size() == kNcorner );
  DeleteContainerOfContainers( points );

#ifdef VERBOSE
  cout << "Collecting coordinates from: (" << fPatterns.size() << " patterns)"
       << endl;
  PrintHits( fHits );
  cout << "Seed pattern: " << endl;
  fPatterns.front()->first.Print();
#endif
  Double_t zp[kNcorner] = { fZL, fZL, fZU, fZU, fZL };
  Bool_t good = true;

  // Collect the hit coordinates within this Road
  TBits planepattern;
  UInt_t last_np = kMaxUInt;
  for( siter_t it = fHits.begin(); it != fHits.end(); ++it ) {
    Hit* hit = const_cast<Hit*>(*it);
    Double_t z = hit->GetZ();
    UInt_t np = hit->GetWirePlane()->GetPlaneNum();
    for( int i = 2; i--; ) {
      Double_t x = i ? hit->GetPosL() : hit->GetPosR();
      if( TMath::IsInside( x, z, kNcorner, 
			   const_cast<Double_t*>(&fCornerX[0]), zp )) {
	// The hits are sorted by ascending plane number
	if( np != last_np ) {
	  points.push_back( Pvec_t() );
	  planepattern.SetBitNumber(np);
	  last_np = np;
	}
	points.back().push_back( new Point(x, z, hit) );
      }
    }
  }
  // Check if this matchpattern is acceptable
  UInt_t patternvalue = 0;
  assert( planepattern.GetNbytes() <= sizeof(patternvalue) );
  planepattern.Get( &patternvalue );
  good = fProjection->GetPlaneCombos()->TestBitNumber(patternvalue)
    // Need at least MinFitPlanes planes for fitting
    and planepattern.CountBits() >= fProjection->GetMinFitPlanes();
  
#ifdef VERBOSE
  cout << "Collected:" << endl;
  vector<Pvec_t>::iterator ipl = points.begin();
  for( UInt_t i = 0; i<fProjection->GetNplanes(); ++i ) {
    cout << " pl= " << i;
    assert( ipl == points.end() or !(*ipl).empty() );
    if( ipl == points.end() 
	or i != (*ipl).front()->hit->GetWirePlane()->GetPlaneNum() )
      cout << " missing";
    else {
      Double_t z = (*ipl).front()->z;
      cout << " z=" << z << "\t x=";
      for( Pvec_t::iterator it = (*ipl).begin(); it != (*ipl).end(); ++it ) {
	cout << " " << (*it)->x;
	assert( (*it)->z == z );
      }
      ++ipl;
    }
    cout << endl;
  }
#endif
  if( !good ) {
    DeleteContainerOfContainers( points );
#ifdef VERBOSE
    cout << "REJECTED" << endl;
#endif
  }
  return good;
}

//_____________________________________________________________________________
Bool_t Road::Fit()
{
  // Collect hit positions within the Road limits and, if enough points
  // found, fit them to a straight line. If several points found in one
  // or more planes, fit all possible combinations of them.
  // Fit results with acceptable chi2 (Projection::fChisqLimits) are
  // sorted into fFitData, lowest chi2 first.

  fGood = false;


  DeleteContainer( fFitData );
  vector<Pvec_t> points;
  points.reserve( fProjection->GetNplanes() );
  // Collect coordinates of hits that are within the width of the road
  if( !CollectCoordinates(points) )
    return false;

  // Determine number of permutations
  //TODO: protect against overflow
  UInt_t n_combinations = accumulate( ALL(points),
				      (UInt_t)1, SizeMul<Pvec_t>() );
  if( n_combinations > kMaxNhitCombos ) {
    // TODO: keep statistics
    DeleteContainerOfContainers( points );
    return false;
  }

  // Loop over all combinations of hits in the planes
  vector<Pvec_t>::size_type npts = points.size();
  Pvec_t selected;
  fDof = npts-2;
  Projection::pdbl_t chi2_interval = fProjection->GetChisqLimits(fDof);
  Double_t* w = new Double_t[npts];
  for( UInt_t i = 0; i < n_combinations; ++i ) {
    NthCombination( i, points, selected );
    assert( selected.size() == npts );

    // Do linear fit of the points, assuming uncorrelated measurements (x_i)
    // and different resolutions for each point.
    // We fit x = a1 + a2*z (z independent).
    // Notation from: Review of Particle Properties, PRD 50, 1277 (1994)
    Double_t S11 = 0, S12 = 0, S22 = 0, G1 = 0, G2 = 0, chi2 = 0;
    for( Pvec_t::size_type j = 0; j < npts; j++) {
      register Point* p = selected[j];
      register Double_t r = w[j] = 1.0 / ( p->res() * p->res() );
      S11 += r;
      S12 += p->z * r;
      S22 += p->z * p->z * r;
      G1  += p->x * r;
      G2  += p->x * p->z * r;
    }
    Double_t D   = S11*S22 - S12*S12;
    Double_t iD  = 1.0/D;
    Double_t a1  = (G1*S22 - G2*S12)*iD;  // Intercept
    Double_t a2  = (G2*S11 - G1*S12)*iD;  // Slope
    Double_t V11 = S22*iD;                // Intercept error
    Double_t V22 = S11*iD;                // Slope error
    for( Pvec_t::size_type j = 0; j < npts; j++) {
      register Point* p = selected[j];
      Double_t d = a1 + a2*p->z - p->x;
      chi2 += d*d * w[j];
    }

#ifdef VERBOSE
    cout << "Fit:"
	 << " a1 = " << a1 << " (" << V11 << ")"
	 << " a2 = " << a2 << " (" << V22 << ")"
	 << " chi2/dof = " << chi2
	 << endl;
#endif
    // Throw out Chi2's outside of selected confidence interval
    // NB: Obviously, this requires accurate hit resolutions
    //TODO: keep statistics
    //TODO: allow disabling of low cutoff via database switch
    if( chi2 < chi2_interval.first )
      continue;
    if( chi2 > chi2_interval.second )
      continue;
#ifdef VERBOSE
    cout << "ACCEPTED" << endl;
#endif

    // Save fits with acceptable Chi2
    fFitData.push_back( new FitResult(a1,a2,chi2,V11,V22) );
    // Save points used for this fit
    fFitData.back()->fFitCoordinates.swap( selected );

  }
  delete [] w;

  if( fFitData.empty() ) {
    //TODO: keep statistics
    DeleteContainerOfContainers( points );
    return false;
  }

  // Sort fit results by ascending chi2
  sort( ALL(fFitData), FitResult::Chi2IsLess() );

  // Copy best fit results to member variables
  const FitResult* best = fFitData.front();
  fPos   = best->fPos;
  fSlope = best->fSlope;
  fChi2  = best->fChi2;
  
  // Save with the wire planes the hit coordinates used in the good fits 
  for( vector<FitResult*>::size_type ifit = 0;
       ifit < fFitData.size(); ++ifit ) {
    Pvec_t& coord = fFitData[ifit]->fFitCoordinates;
    for( Pvec_t::iterator it = coord.begin(); it != coord.end(); ++it ) {
      Point* p = *it;
      WirePlane* wp = p->hit->GetWirePlane();
      Double_t slope = fFitData[ifit]->fSlope;
      Double_t x_track = fFitData[ifit]->fPos + p->z * slope;
	wp->AddFitCoord( FitCoord( p->hit, this, ifit, p->x, x_track, slope ));
    }
    // At this point, we don't need the saved pointers anymore
    coord.clear();
  }

  DeleteContainerOfContainers( points );
  return (fGood = true);
}

//_____________________________________________________________________________
void Road::Print( Option_t* opt ) const
{
  // Print road info

}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
