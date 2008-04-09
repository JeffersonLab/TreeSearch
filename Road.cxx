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
#include "TVector2.h"

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

typedef Hset_t::iterator siter_t;
#define ALL(c) (c).begin(), (c).end()

// Number of points for polygon test
static const size_t kNcorner = 5;
static const UInt_t kMaxNhitCombos = 1000;

//_____________________________________________________________________________
// Private class for building a cluster of patterns
struct BuildInfo_t {
  HitSet            fCluster;      // Copy of start HitSet of the cluster
  vector< pair<UShort_t,UShort_t> > 
                    fLimits; // [nplanes] Min/max bin numbers in each plane

  BuildInfo_t() {}
  BuildInfo_t( const Node_t& node ) : fCluster(node.second) 
  {
    // Construct from given start pattern/hitset
    const NodeDescriptor& nd = node.first;
    assert( fCluster.plane_pattern > 0 and fCluster.nplanes > 0 );
    UInt_t npl = nd.link->GetPattern()->GetNbits();
    fLimits.reserve( npl );
    for( UInt_t i = 0; i < npl; ++i )
      fLimits.push_back( make_pair(nd[i], nd[i]+1) );
  }
  void ExpandWidth( const NodeDescriptor& nd ) {
    // Widen the bin ranges using the bins in the given pattern
    UInt_t npl = nd.link->GetPattern()->GetNbits();
    assert( fLimits.size() == npl or fLimits.empty() );
    if( fLimits.empty() )
      fLimits.assign( npl, make_pair(kMaxUShort, 0) );
    for( UInt_t i = npl; i; ) { --i;
      fLimits[i] = make_pair(TMath::Min(fLimits[i].first,  nd[i]),
			     TMath::Max(fLimits[i].second, UShort_t(nd[i]+1)));
    }
  }
};

//_____________________________________________________________________________
// Helper functor for coordinate conversions
class GetBinX 
  : public unary_function< pair<Double_t,Double_t>, pair<UInt_t,UInt_t> >
{
private:
  const Hitpattern* fHpat;
public:
  explicit GetBinX( const Hitpattern* hpat ) : fHpat(hpat) {}
  pair<Double_t,Double_t> operator() ( const pair<UInt_t,UInt_t>& bin ) const
  {
    // Get X coordinate of left edge of given bin
    Double_t lo = static_cast<Double_t>(bin.first) * fHpat->GetBinWidth() 
      - fHpat->GetOffset();
    Double_t hi = static_cast<Double_t>(bin.second) * fHpat->GetBinWidth() 
      - fHpat->GetOffset();
    return make_pair(lo,hi);
  }
};

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
Road::Road( const Projection* proj ) 
  : TObject(), fProjection(proj), fZL(kBig), fZU(kBig), 
    fPos(kBig), fSlope(kBig), fChi2(kBig), fDof(kMaxUInt), fGood(true),
    fTrack(0), fBuild(0)
{
  // Construct empty road

  assert(fProjection);  // Invalid Projection* pointer

  memset( fCornerX, 0, kNcorner*sizeof(Double_t) );
  fV[2] = fV[1] = fV[0] = kBig;
  fPoints.reserve( fProjection->GetNplanes() );
  fBuild = new BuildInfo_t;
}

//_____________________________________________________________________________
Road::Road( const Node_t& nd, const Projection* proj )
  : TObject(), fProjection(proj), fZL(kBig), fZU(kBig),
    fPos(kBig), fSlope(kBig), fChi2(kBig), fDof(kMaxUInt), fGood(true),
    fTrack(0), fPatterns(1,&nd), fBuild(0)
{
  // Construct from pattern

  assert( fProjection );   // Invalid Projection* pointer
  assert( CheckMatch(nd.second.hits) ); // Start pattern must be a good match

  memset( fCornerX, 0, kNcorner*sizeof(Double_t) );
  fV[2] = fV[1] = fV[0] = kBig;
  fPoints.reserve( fProjection->GetNplanes() );
  fBuild = new BuildInfo_t(nd);

#ifdef VERBOSE
  if( fProjection->GetDebug() > 3 ) {
    cout << "New Road:" << endl;
    nd.first.Print();
    PrintHits(nd.second.hits);
  }
#endif
}

//_____________________________________________________________________________
Road::Road( const Road& orig ) : 
  TObject(orig), fPatterns(orig.fPatterns), fHits(orig.fHits)
{
  // Copy constructor

  size_t nbytes = (char*)&fTrack - (char*)&fProjection + sizeof(fTrack);
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
    size_t nbytes = (char*)&fTrack - (char*)&fProjection + sizeof(fTrack);
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
inline
Bool_t Road::CheckMatch( const Hset_t& hits ) const
{
  // Return true if the hits from the given set either cover all planes
  // or, if planes are missing, the pattern of missing planes is allowed
  // (based on what level of matching the user requests via the database)

  return HitSet::CheckMatch( hits, fProjection->GetPlaneCombos() );
}

//_____________________________________________________________________________
static inline
Bool_t OuterBitsSet( UInt_t p1, UInt_t p2, UInt_t nbits )
{
  // Test if p1 and p2 have the same first bits set and the same last bits set
  if( p1 == 0 or p2 == 0 )
    return false;
  register UInt_t k = 1;
  while( (p1 & k) == 0 and (p2 & k) == 0 )
    k <<= 1;
  if( !(p1 & k) or !(p2 & k) )
    return false;
  k = 1U << (nbits-1);
  while( (p1 & k) == 0 and (p2 & k) == 0 )
    k >>= 1;
  return ( (p1 & k) and (p2 & k) );
}

//_____________________________________________________________________________
inline
Bool_t Road::IsInRange( const NodeDescriptor& nd ) const
{
  // Test if given pattern is within the allowed maximum distance from the
  // current front and back bin ranges of the cluster

  assert( fBuild );
  assert( !fBuild->fLimits.empty() );

  UInt_t fdist = fProjection->GetBinMaxDistF();
  UInt_t bdist = fProjection->GetBinMaxDistB();
  
  return ( nd.Start() + fdist >= fBuild->fLimits.front().first and
	   nd.Start()         <  fBuild->fLimits.front().second + fdist and
	   nd.End()   + bdist >= fBuild->fLimits.back().first and
	   nd.End()           <  fBuild->fLimits.back().second + bdist );
}


//_____________________________________________________________________________
Bool_t Road::Add( const Node_t& nd )
{
  // Check if the hits from the given NodeDescriptor pattern are common
  // with the common hit set already in this road. If so, add the pattern
  // to this road, update the common hits if necessary, and return true.
  // If not, do nothing and return false.
  //
  // Adding can only be done so long as the road is not yet finished

  assert(fBuild);  // Add() after Finish() not allowed

  const HitSet& new_set  = nd.second;
  const Hset_t& new_hits = nd.second.hits;

#ifdef VERBOSE
  if( fProjection->GetDebug() > 3 ) {
    cout << "Adding:" << endl;
    nd.first.Print();
    PrintHits(new_hits);
  }
#endif

  // If hitdist > 0, roads will collect patterns not only identical hits,
  // but also nearby ones (at most hitdist wire numbers apart)
  UInt_t hitdist = fProjection->GetHitMaxDist();

  if( fPatterns.empty() ) {
    // If this is the first pattern, initialize the cluster
    assert( CheckMatch(new_hits) );
    assert( new_set.nplanes > 0 && new_set.plane_pattern > 0 );
    fBuild->fCluster = new_set;
    fBuild->ExpandWidth( nd.first );

#ifdef VERBOSE
    if( fProjection->GetDebug() > 3 ) {
      cout << "New cluster:" << endl;
      PrintHits( new_hits );
    }
#endif
  } 
  else if( IsInRange(nd.first) and 
	   fBuild->fCluster.IsSimilarTo(new_set,hitdist) ) {
    // Accept this pattern if and only if it is a subset of the cluster
    // NB: IsSimilarTo() is a looser match than std::includes(). The new
    // pattern may have extra hits

    // If hitdist > 0, grow the cluster with possible new hits found
    if( hitdist > 0 ) {
      fBuild->fCluster.hits.insert( ALL(new_hits) );
      // Growing clusters are expanded even for lower match levels
      // if they contain hits in the first and last active bins
      if( OuterBitsSet( new_set.plane_pattern, 
			fBuild->fCluster.plane_pattern,
			fProjection->GetNplanes() ))
	fBuild->ExpandWidth( nd.first );
    }
    // Expand the road width, but only for patterns of the same match level.
    // Patterns of lower match level can be artificially wide
    else if( new_set.nplanes == fBuild->fCluster.nplanes )
      fBuild->ExpandWidth( nd.first );

  }
  else {
    // The pattern does not fit into this cluster
#ifdef VERBOSE
    if( fProjection->GetDebug() > 3 )
      cout << "Not a subset of this road" << endl;
#endif
    return false;
  }

  // Save a pointer to the added pattern so we can update it later
  fPatterns.push_back(&nd);

#ifdef VERBOSE
    if( fProjection->GetDebug() > 3 ) 
      cout << "new npat = " << fPatterns.size() << endl;
#endif
  return true;
}

//_____________________________________________________________________________
void Road::Finish()
{
  // Finish building the road

  assert(fBuild);   // Road must be incomplete to be able to Finish()
  for( list<const Node_t*>::iterator it = fPatterns.begin(); it !=
	 fPatterns.end(); ++it ) {
    (**it).second.used = 1;
#ifdef VERBOSE
    if( fProjection->GetDebug() > 3 ) { 
      cout << "used pat = ";
      (**it).first.Print(); 
    }
#endif
  }

  // Save the cluster hits. These are the hits that will be considered for
  // fitting later
  assert( fHits.empty() );
  fHits.swap( fBuild->fCluster.hits );

  // Calculate the corner coordinates fCornerX of a polygon with points in the
  // order LL (lower left), LR, UR, UL, LL, as needed by TMath::IsInside()

  // Convert the bin numbers of the left/right edges to physical coordinates
  vector< pair<Double_t,Double_t> > edgpos;
  UInt_t npl = fProjection->GetNplanes(), lastpl = npl-1;
  edgpos.reserve( npl );
  GetBinX binx( fProjection->GetHitpattern() ); // bin# -> position
  transform( ALL(fBuild->fLimits), back_inserter(edgpos), binx );
  
  // Define the road boundaries as the lines with the slope between the
  // front and back bins, shifted to the left/right to include fully the
  // left/rightmost bin. This may not be optimal, but it is fast, unambiguous,
  // and removes the bias towards the front/back bins.
  TVector2 LL( edgpos.front().first,  fProjection->GetPlaneZ(0) );
  TVector2 LR( edgpos.front().second, fProjection->GetPlaneZ(0) );
  TVector2 UL( edgpos.back().first,   fProjection->GetPlaneZ(lastpl) );
  TVector2 UR( edgpos.back().second,  fProjection->GetPlaneZ(lastpl) );
  Double_t slopeL = (UL.X() - LL.X()) / (UL.Y() - LL.Y());
  Double_t slopeR = (UR.X() - LR.X()) / (UR.Y() - LR.Y());
  Double_t maxdL = 0, maxdR = 0;
  for( UInt_t i = lastpl; i>1; ) { --i;
    Double_t z = fProjection->GetPlaneZ(i);
    // Crossing points in test plane (i) of the L/R front/back connections
    Double_t XL = UL.X() + slopeL * (z - UL.Y());
    Double_t XR = UR.X() + slopeR * (z - UR.Y());
    // Memorize the bins that stick out the most on left/right. Reuse UL/UR
    Double_t lpt = edgpos[i].first, rpt = edgpos[i].second;
    Double_t dL = lpt - XL, dR = rpt - XR;
    if( dL < maxdL ) {
      UL.Set( lpt, z );
      maxdL = dL;
    } 
    if( dR > maxdR ) {
      UR.Set( rpt, z );
      maxdR = dR;
    } 
  }
  // Now UL and UR are points on the left and right road boundary, 
  // respectively. The boundaries share the common slope calculated above.
  
  // To compute the corners of the polygon, shift the z-positions of the first
  // and last planes down and up, respectively, by eps, so that the points 
  // on the planes are guaranteed to be included
  const Double_t eps = 1e-3;
  fZL = fProjection->GetPlaneZ(0) - eps;
  fZU = fProjection->GetPlaneZ(lastpl) + eps;
  assert( fZL < fZU );

  fCornerX[0] = UL.X() + slopeL * (fZL - UL.Y());
  fCornerX[1] = UR.X() + slopeR * (fZL - UR.Y());
  fCornerX[2] = UR.X() + slopeR * (fZU - UR.Y());
  fCornerX[3] = UL.X() + slopeL * (fZU - UL.Y());
  fCornerX[4] = fCornerX[0];

  assert( fCornerX[0] < fCornerX[1] );
  assert( fCornerX[3] < fCornerX[2] );

  // All done. Put the tools away
  delete fBuild; fBuild = 0;

  return;
}

//_____________________________________________________________________________
Bool_t Road::Include( const Road* other )
{
  // Check if this road includes 'other'. Either other's hits are a
  // subset of this road's hit set or other lies entirely within the
  // boundaries of this road.
  // If a match is found, adopt the other road (either widen boundaries
  // or adopt the other's hits)

  static const double eps = 1e-6;

  assert( other and !fBuild and fProjection == other->fProjection );

  if( includes(ALL(fHits), ALL(other->fHits), fHits.key_comp()) ) {
    // Widen the road bundaries
    fCornerX[0] = TMath::Min( fCornerX[0], other->fCornerX[0] );
    fCornerX[1] = TMath::Max( fCornerX[1], other->fCornerX[1] );
    fCornerX[2] = TMath::Max( fCornerX[2], other->fCornerX[2] );
    fCornerX[3] = TMath::Min( fCornerX[3], other->fCornerX[3] );
    fCornerX[4] = fCornerX[0];
    return true;
  }

  if( TMath::Abs(fZL-other->fZL) < eps and
      TMath::Abs(fZU-other->fZU) < eps and
      fCornerX[0] < other->fCornerX[0] + eps and
      other->fCornerX[1] < fCornerX[1] + eps and
      fCornerX[3] < other->fCornerX[3] + eps and
      other->fCornerX[2] < fCornerX[2] + eps ) {
    fHits.insert( ALL(other->fHits) );
    return true;
  }

  return false;
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
  if( fProjection->GetDebug() > 3 ) {
    cout << "Collecting coordinates from: (" << fPatterns.size() 
	 << " patterns)" << endl;
    PrintHits( fHits );
    cout << "Seed pattern: " << endl;
    fPatterns.front()->first.Print();
  }
#endif
  Double_t zp[kNcorner] = { fZL, fZL, fZU, fZU, fZL };
  Bool_t good = true;

  // Collect the hit coordinates within this Road
  TBits planepattern;
  UInt_t last_np = kMaxUInt;
  for( siter_t it = fHits.begin(); it != fHits.end(); ++it ) {
    Hit* hit = const_cast<Hit*>(*it);
    Double_t z = hit->GetZ();
    UInt_t np = hit->GetPlaneNum();
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
  if( fProjection->GetDebug() > 3 ) {
    cout << "Collected:" << endl;
    vector<Pvec_t>::reverse_iterator ipl = fPoints.rbegin();
    for( UInt_t i = fProjection->GetNplanes(); i--; ) {
      cout << " pl= " << i;
      assert( ipl == fPoints.rend() or !(*ipl).empty() );
      if( ipl == fPoints.rend() 
	  or i != (*ipl).front()->hit->GetPlaneNum() )
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
    if( fProjection->GetDebug() > 3 ) 
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
    if( fProjection->GetDebug() > 3 ) cout << "ACCEPTED" << endl;
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
