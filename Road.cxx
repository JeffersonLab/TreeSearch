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
#include <list>
#include <utility>

using namespace std;

ClassImp(TreeSearch::Road);

namespace TreeSearch {

// Private class for use while building a Road
struct BuildInfo_t {
  list<pair<const NodeDescriptor,HitSet>*> fPatterns; // Patterns in this road
  set<Hit*>         fCommonHits;      // Hits common between all patterns
  const Hitpattern* fHitpattern;
  UInt_t            fNlayers;
  const TBits*      fPlaneCombos;     // Bitfield for plane combo lookup
  //TODO: use fMaxdist[] ??
};

//_____________________________________________________________________________
Road::Road( const Projection* proj )
  : fSlope(kBig), fPos(kBig), fChi2(kBig)
{
  // Constructor

  fLeft[0] = fLeft[1] = kMaxUShort;
  fRight[0] = fRight[1] = 0;
  fErr[0] = fErr[1] = kBig;

  assert(proj);

  fBuild = new BuildInfo_t;
  fBuild->fHitpattern  = proj->GetHitpattern();
  fBuild->fNlayers     = proj->GetNlayers();
  fBuild->fPlaneCombos = proj->GetPlaneCombos();

  fNplanes = proj->GetNplanes();

  assert( fBuild->fHitpattern && fBuild->fNlayers && 
	  fNplanes >= fBuild->fNlayers && fBuild->fPlaneCombos );

}

//_____________________________________________________________________________
Road::Road( const Road& orig )
  : fNplanes(orig.fNplanes), fSlope(orig.fSlope), fPos(orig.fSlope),
    fChi2(orig.fChi2)
{
  // Copy constructor

  fLeft[0]  = orig.fLeft[0];  fLeft[1]  = orig.fLeft[1];
  fRight[0] = orig.fRight[0]; fRight[1] = orig.fRight[1];
  fErr[0]   = orig.fErr[0];   fErr[1]   = orig.fErr[1];

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
    fLeft[0]  = rhs.fLeft[0];  fLeft[1]  = rhs.fLeft[1];
    fRight[0] = rhs.fRight[0]; fRight[1] = rhs.fRight[1];
    fSlope  = rhs.fSlope;
    fPos    = rhs.fPos;
    fChi2   = rhs.fChi2;
    fErr[0] = rhs.fErr[0]; fErr[1] = rhs.fErr[1];

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

  return fBuild->fPlaneCombos->TestBitNumber(curpat);
}

//_____________________________________________________________________________
Bool_t Road::Add( pair<const NodeDescriptor,HitSet>& nd )
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
  nd.first.Print(); nd.first.link->GetPattern()->Print();
  nd.first.parent->Print();
  PrintHits(nd.second.hits);
#endif
  bool first = fBuild->fPatterns.empty();
  if( first ) {
    if( !CheckMatch(nd.second.hits) )
      return kFALSE;
    fBuild->fCommonHits = fHits = nd.second.hits;
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
      swap( fBuild->fCommonHits, new_commons );
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
      swap( fHits, new_hits );
#ifdef VERBOSE
      PrintHits( fHits );
#endif
    }
  }

  // Save a pointer to this pattern so we can update it later
  fBuild->fPatterns.push_back(&nd);

  // Expand the road limits if necessary
  fLeft[0]  = TMath::Min( nd.first.Start(), fLeft[0] );
  fLeft[1]  = TMath::Min( nd.first.End(),   fLeft[1] );
  fRight[0] = TMath::Max( nd.first.Start(), fRight[0] );
  fRight[1] = TMath::Max( nd.first.End(),   fRight[1] );

#ifdef VERBOSE
  cout << "new npat = " << fBuild->fPatterns.size() << endl;
  cout << "new left/right = " << fLeft[0]<<" "<<fRight[0]<<" "
       << fLeft[1]<<" "<<fRight[1] << endl;
#endif

  return kTRUE;
}

//_____________________________________________________________________________
void Road::Finish()
{
  // Finish building the road

  assert(fBuild);
  for( list<pair<const NodeDescriptor,HitSet>*>::iterator it = 
	 fBuild->fPatterns.begin(); it != fBuild->fPatterns.end(); ++it ) {
    
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


  // Put the tools away
  // TODO: save npat, nplanes
  delete fBuild; fBuild = 0;

  return;
}

//_____________________________________________________________________________
void Road::Print( Option_t* opt ) const
{
  // Print road info

}

//_____________________________________________________________________________
void Road::CollectCoordinates( UInt_t nplanes, std::vector<Int_t>& hitcount,
			       std::vector<std::vector<Point> >& points )
{
  // 

}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
