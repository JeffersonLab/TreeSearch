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
#include <map>
#include <utility>

using namespace std;

ClassImp(TreeSearch::Road)
ClassImp(TreeSearch::Road::Point)
ClassImp(TreeSearch::Road::Corners)

namespace TreeSearch {

// Private class for use while building a Road
struct BuildInfo_t {
  Hset_t            fClusterHits;  // Hits common between all patterns
  vector<UShort_t>  fCornerBin;    // Bin numbers defining the corners
                                   // 0=LL, 1=LR, 2=UR, 3=UL 
  BuildInfo_t() {
    UShort_t corner_bins[4] = { kMaxUShort, 0, 0, kMaxUShort };
    fCornerBin.assign( corner_bins, corner_bins+4 );
  }
  BuildInfo_t( const Node_t& nd ) : fClusterHits(nd.second.hits) 
  {
    UShort_t corner_bins[4] = { nd.first.Start(), nd.first.Start(),
				nd.first.End(),   nd.first.End() };
    fCornerBin.assign( corner_bins, corner_bins+4 );
  }
};

// Number of points for polygon test
static const size_t kNcorner = 5;
static const UInt_t kMaxNhitCombos = 1000;

typedef Hset_t::iterator siter_t;
#define ALL(c) (c).begin(), (c).end()

//_____________________________________________________________________________
Road::Road( const Projection* proj ) 
  : TObject(), fProjection(proj), fZL(kBig), fZU(kBig), 
    fPos(kBig), fSlope(kBig), fChi2(kBig), fDof(kMaxUInt), fGood(true),
    fBuild(0)
{
  // Construct empty road

  assert(fProjection);  // Invalid Projection* pointer

  memset( fCornerX, 0, kNcorner*sizeof(Double_t) );
  fV[2] = fV[1] = fV[0] = kBig;
  fPoints.reserve( fProjection->GetNplanes() );
  fBuild = new BuildInfo_t;
}

//_____________________________________________________________________________
Road::Road( Node_t& nd, const Projection* proj )
  : TObject(), fProjection(proj), fZL(kBig), fZU(kBig),
    fPos(kBig), fSlope(kBig), fChi2(kBig), fDof(kMaxUInt), fGood(true),
    fPatterns(1,&nd), fHits(nd.second.hits), fBuild(0)
{
  // Construct from pattern

  assert(fProjection);  // Invalid Projection* pointer

  memset( fCornerX, 0, kNcorner*sizeof(Double_t) );
  fV[2] = fV[1] = fV[0] = kBig;
  fPoints.reserve( fProjection->GetNplanes() );
  fBuild = new BuildInfo_t(nd);

#ifdef VERBOSE
  cout << "Adding:" << endl;
  nd.first.Print();
  PrintHits(nd.second.hits);
  cout << "New cluster:" << endl;
  PrintHits( fBuild->fClusterHits );
#endif
}

//_____________________________________________________________________________
Road::Road( const Road& orig ) : 
  TObject(orig), fPatterns(orig.fPatterns), fHits(orig.fHits)
{
  // Copy constructor

  size_t nbytes = (char*)&fGood - (char*)&fProjection + sizeof(fGood);
  memcpy( &fProjection, &orig.fProjection, nbytes );
  
  CopyPointData( orig );

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
    size_t nbytes = (char*)&fGood - (char*)&fProjection + sizeof(fGood);
    memcpy( &fProjection, &rhs.fProjection, nbytes );

    DeleteContainer( fFitData );
    DeleteContainerOfContainers( fPoints );
    CopyPointData( rhs );

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
  DeleteContainerOfContainers( fPoints );
  delete fBuild;

}

//_____________________________________________________________________________
void Road::CopyPointData( const Road& orig )
{
  // Copy fPoints and fFitData. Used by copy c'tor and assignment operator.
  // Creates actual copies of Points because they are managed by the Roads.

  if( orig.fPoints.empty() )
    assert( orig.fFitData.empty() ); // Can't have fit data but no points :-/
  else {
    typedef map<Point*,Point*> Pmap_t;
    Pmap_t xref;
    fPoints.resize( orig.fPoints.size() );
    for( vector<Pvec_t>::size_type i = 0; i < orig.fPoints.size(); ++i ) {
      const Pvec_t& old_planepoints = orig.fPoints[i];
      fPoints[i].reserve( old_planepoints.size() );
      for( Pvec_t::const_iterator it = old_planepoints.begin(); it !=
	     old_planepoints.end(); ++it ) {
	Point* old_point = *it;
	Point* new_point = new Point( *old_point );
	fPoints[i].push_back( new_point );
	// It gets a bit tricky here: To be able to copy the fFitData, which
	// contain pointers to some of the Points in fPoints, we need to keep 
	// track of which old point was copied to which new point.
	pair<Pmap_t::iterator,bool> 
	  ins = xref.insert( make_pair(old_point,new_point) );
	assert( ins.second );  // Duplicate points should never occur
      }
    }

    // Copy fit results
    fFitData.reserve( orig.fFitData.size() );
    for( vector<FitResult*>::const_iterator it = orig.fFitData.begin();
	 it != orig.fFitData.end(); ++it ) {
      // The FitResult copy c'tor does not copy the fFitCoordinates array
      fFitData.push_back( new FitResult( **it ));
      // The copied FitResult must reference the copied Points
      const Pvec_t& old_points = (*it)->GetPoints();
      Pvec_t&       new_points = fFitData.back()->fFitCoordinates;
      new_points.reserve( old_points.size() );
      for( Pvec_t::const_iterator it2 = old_points.begin(); it2 !=
	     old_points.end(); ++it2 ) {
	Point* old_point = *it2;
	Pmap_t::iterator found = xref.find( old_point );
	assert( found != xref.end() );
	new_points.push_back( (*found).second );
      }
    }
  }
}
  
//_____________________________________________________________________________
#ifdef VERBOSE
static void PrintHits( const Hset_t& hits )
{
  //  cout << hits.size() << " hits" << endl;

  for( Hset_t::reverse_iterator it = hits.rbegin(); it != hits.rend(); ++it ) {
    cout << " ";
    (*it)->Print();
  }
  
}
#endif
//_____________________________________________________________________________
inline
Bool_t Road::CheckMatch( const Hset_t& hits ) const
{
  // Return true if the hits from the given set either cover all planes
  // or, if planes are missing, the pattern of missing planes is allowed
  // (based on what level of matching the user requests via the database)

  assert(fBuild);  // This function used only in build mode

  return HitSet::CheckMatch( hits, fProjection->GetPlaneCombos() );
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
  if( fPatterns.empty() ) {
    // Check the hits of the first pattern of a new road for missing planes
    if( !CheckMatch(nd.second.hits) )
      return kFALSE;
    // The sets of all hits and of the common hits in a new road are identical
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
      // Add the intersection to the prior hits, accumulating the cluster.
      // It may seem like nonsense to insert a subset of a set into itself,
      // but note that "intersection" is not a strict subset of fClusterHits.
      // It may contain equivalent hits present in the new pattern but not in
      // the cluster. So the following line actually does something.
      // TODO: could this be done more efficiently?
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
      other->fCornerX[1] + eps < fCornerX[0] or 
      fCornerX[1] + eps < other->fCornerX[0] or
      other->fCornerX[2] + eps < fCornerX[3] or 
      fCornerX[2] + eps < other->fCornerX[3] )
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
  const UInt_t np1 = fProjection->GetNplanes()-1;
  const Double_t resL = 2.0*fProjection->GetPlane(0)->GetResolution();
  const Double_t resU = 2.0*fProjection->GetPlane(np1)->GetResolution();
  fCornerX[0] = GetBinX( fBuild->fCornerBin[0] )   - resL;
  fCornerX[1] = GetBinX( fBuild->fCornerBin[1]+1 ) + resL;
  fCornerX[2] = GetBinX( fBuild->fCornerBin[2]+1 ) + resU;
  fCornerX[3] = GetBinX( fBuild->fCornerBin[3] )   - resU;

  // z-shift to include the planes proper.
  // The upper and lower edges are shifted up and down, respectively,
  // so that the points on the planes are guaranteed to be included.
  const Double_t eps = 1e-3;
  fZL = fProjection->GetPlaneZ(0) - eps;
  fZU = fProjection->GetPlaneZ(np1) + eps;
  assert( fZL < fZU );

  // Adjust x coordinates for the z shift
  Double_t dZplanes = fProjection->GetPlaneZ(np1) - fProjection->GetPlaneZ(0);
  assert( dZplanes > 1e-3 );
  Double_t slopeL = ( fCornerX[3]-fCornerX[0] ) * (1.0/dZplanes);
  Double_t slopeR = ( fCornerX[2]-fCornerX[1] ) * (1.0/dZplanes);
  fCornerX[0] -= slopeL * eps;
  fCornerX[1] -= slopeR * eps;
  fCornerX[2] += slopeR * eps;
  fCornerX[3] += slopeL * eps;

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
Bool_t Road::CollectCoordinates()
{
  // Gather hit positions that lie within the Road area.
  // Return true if the plane occupancy pattern of the selected points 
  // is allowed by Projection::fPlaneCombos, otherwise false.
  // Results are in fPoints.

  DeleteContainerOfContainers( fPoints );

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
    // Prevent duplicate entries for hits with zero drift
    int i = (hit->GetDriftDist() == 0.0) ? 1 : 2;
    do {
      Double_t x = (--i != 0) ? hit->GetPosL() : hit->GetPosR();
      if( TMath::IsInside( x, z, kNcorner, &fCornerX[0], zp )) {
	// The hits are sorted by ascending plane number
	if( np != last_np ) {
	  fPoints.push_back( Pvec_t() );
	  planepattern.SetBitNumber(np);
	  last_np = np;
	}
	fPoints.back().push_back( new Point(x, z, hit) );
      }
    } while( i );
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
  vector<Pvec_t>::reverse_iterator ipl = fPoints.rbegin();
  for( UInt_t i = fProjection->GetNplanes(); i--; ) {
    cout << " pl= " << i;
    assert( ipl == fPoints.rend() or !(*ipl).empty() );
    if( ipl == fPoints.rend() 
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
  if( !good ) {
    cout << "REJECTED" << endl;
  }
#endif

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

  assert( fGood );
  DeleteContainer( fFitData );

  // Collect coordinates of hits that are within the width of the road
  if( !CollectCoordinates() )
    return false;

  // Determine number of permutations
  //TODO: protect against overflow
  UInt_t n_combinations = accumulate( ALL(fPoints),
				      (UInt_t)1, SizeMul<Pvec_t>() );
  if( n_combinations > kMaxNhitCombos ) {
    // TODO: keep statistics
    return false;
  }

  // Loop over all combinations of hits in the planes
  vector<Pvec_t>::size_type npts = fPoints.size();
  Pvec_t selected;
  fDof = npts-2;
  Projection::pdbl_t chi2_interval = fProjection->GetChisqLimits(fDof);
  Double_t* w = new Double_t[npts];
  for( UInt_t i = 0; i < n_combinations; ++i ) {
    NthCombination( i, fPoints, selected );
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
    // Covariance matrix of the fitted parameters
    Double_t V[3] = { S22*iD, -S12*iD, S11*iD };
    for( Pvec_t::size_type j = 0; j < npts; j++) {
      register Point* p = selected[j];
      Double_t d = a1 + a2*p->z - p->x;
      chi2 += d*d * w[j];
    }

#ifdef VERBOSE
    cout << "Fit:"
	 << " a1 = " << a1 << " (" << TMath::Sqrt(V[0]) << ")"
	 << " a2 = " << a2
	 << " chi2 = " << chi2
	 << " ndof = " << fDof
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

#ifdef TESTCODE
    // TESTCODE saves all fits with acceptable Chi2
    fFitData.push_back( new FitResult(a1,a2,chi2,V) );
    fFitData.back()->fFitCoordinates.swap( selected );
#else
    // Production code saves only the best fit
    if( fFitData.empty() || chi2 < fFitData.back()->fChi2 ) {
      if( fFitData.empty() )
	fFitData.push_back( new FitResult(a1,a2,chi2,V) );
      else 
	fFitData.back()->Set(a1,a2,chi2,V);
      // Save points used for this fit
      fFitData.back()->fFitCoordinates.swap( selected );
    }
#endif
  }
  delete [] w;

  if( fFitData.empty() ) {
    //TODO: keep statistics
    return false;
  }

#ifdef TESTCODE
  // Sort fit results by ascending chi2
  sort( ALL(fFitData), FitResult::Chi2IsLess() );

  // Save with the WirePlanes the hit coordinates used in all the good
  // projection (2D) fits 
  for( vector<FitResult*>::size_type ifit = 0; ifit < fFitData.size();
       ++ifit ) {
    Pvec_t& fitpoints = fFitData[ifit]->fFitCoordinates;
    for( Pvec_t::iterator it = fitpoints.begin(); it != fitpoints.end(); 
	 ++it ) {
      Point* p = *it;
      Double_t slope   = fFitData[ifit]->fSlope;
      Double_t x_track = fFitData[ifit]->fPos + p->z * slope;
      FitCoord* fc = p->hit->GetWirePlane()->AddFitCoord
	( FitCoord( p->hit, this, p->x, x_track, slope, kBig, kBig, ifit ));
      p->coord.push_back(fc);
    }
  }
#endif

  // Copy best fit results to member variables
  memcpy( &fPos, &(fFitData.front()->fPos), 6*sizeof(Double_t) );

  return true;
}

//_____________________________________________________________________________
TVector2 Road::Intersect( const Road* other, Double_t z ) const
{
  // Find intersection point in x/y coordinates of the best fit of this
  // road with the best fit of the other in the given z-plane.
  // Both roads must have a good fit and must be from different projections.
  // This function commutes, i.e. a->Intersect(b) == b->Intersect(a)

  assert(other);
  assert(fGood && other->fGood);
  assert(fProjection->GetType() != other->fProjection->GetType());

  Double_t su = fProjection->GetSinAngle();
  Double_t cu = fProjection->GetCosAngle();
  Double_t sv = other->fProjection->GetSinAngle();
  Double_t cv = other->fProjection->GetCosAngle();
  // FIXME: this should be calculated only once
  Double_t inv_denom = 1.0/(sv*cu-su*cv);

  // Standard formulae for the intersection of non-orthogonal coordinates
  // (cf. THaVDCUVTrack.C)
  Double_t x = (GetPos(z) * sv - other->GetPos(z) * su) * inv_denom;
  Double_t y = (other->GetPos(z) * cu - GetPos(z) * cv) * inv_denom;
  return TVector2(x, y);
}

//_____________________________________________________________________________
void Road::Print( Option_t* ) const
{
  // Print road info

  cout << "Road: " << fProjection->GetName() << " type, "
       << fHits.size() << " hits" << endl;
  cout << " TL = (" << fCornerX[3] << "," << fZU << ") "
       << " TR = (" << fCornerX[2] << "," << fZU << ")"
       << endl
       << " BL = (" << fCornerX[0] << "," << fZL << ") "
       << " BR = (" << fCornerX[1] << "," << fZL << ")"
       << endl;

  for( Hset_t::reverse_iterator it = fHits.rbegin(); it != fHits.rend();
       ++it ) {
    cout << " ";
    (*it)->Print();
  }
  if( fGood ) {
    cout << " fitpos = " << fPos << " (" << TMath::Sqrt(fV[0]) << ")"
	 << " slope = " << fSlope
	 << " chi2 = " << fChi2
	 << " ndof = " << fDof
	 << endl;
  }
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
