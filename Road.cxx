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
#include "TMath.h"
#include "TBits.h"
#include <iostream>
#include <algorithm>

using namespace std;

ClassImp(TreeSearch::Road);

namespace TreeSearch {

// Private class for use while building a Road
struct BuildInfo_t {
  vector<NodeDescriptor*> fPatterns;  // All patterns associated with this road
  set<Hit*>         fCommonHits;      // Hits common between all patterns
  Hitpattern*       fHitpattern;
  UInt_t            fNlayers;
  UInt_t            fNplanes;
  //TODO: add match requirements
  //TODO: use fMaxdist[] ??
};

static const vector<NodeDescriptor*>::size_type kEstNumPat = 20;

typedef vector<Hit*>::const_iterator viter_t;
typedef vector<Hit*>::const_iterator vsiz_t;

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
  fBuild->fHitpattern = proj->GetHitpattern();
  fBuild->fNlayers    = proj->GetNlayers();
  fBuild->fNplanes    = proj->GetNplanes();
  fBuild->fPatterns.reserve( kEstNumPat );

  assert( fBuild->fHitpattern && fBuild->fNlayers && 
	  fBuild->fNplanes >= fBuild->fNlayers );

}

//_____________________________________________________________________________
Road::Road( const Road& orig )
  : fSlope(orig.fSlope), fPos(orig.fSlope), fChi2(orig.fChi2)
{
  // Copy constructor

  fLeft[0]  = orig.fLeft[0];  fLeft[1]  = orig.fLeft[1];
  fRight[0] = orig.fRight[0]; fRight[1] = orig.fRight[1];
  fErr[0]   = orig.fErr[0];   fErr[1]   = orig.fErr[1];

  delete fBuild; fBuild = 0;
  if( orig.fBuild )
    fBuild = new BuildInfo_t( *orig.fBuild );
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

    delete fBuild; fBuild = 0;
    if( rhs.fBuild )
      fBuild = new BuildInfo_t( *rhs.fBuild );
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
static void PrintHits( const vector<Hit*>& hits )
{
  cout << hits.size() << " hits" << endl;

  for( vector<Hit*>::size_type i=0; i<hits.size(); ++i ) {
    cout << " ";
    hits[i]->Print();
  }
  
}

//_____________________________________________________________________________
Bool_t Road::CheckMatch( const set<Hit*>& hits ) const
{
  // Match evaluation function. Return true if the numbers of hits in each
  // plane contained in hitcount is sufficient for considering the road
  // cohesive.

  assert(fBuild);

  TBits curpat( fBuild->fNplanes );
  for( set<Hit*>::const_iterator it = hits.begin(); it != hits.end(); ++it )
    curpat.SetBitNumber( (*it)->GetWirePlane()->GetPlaneNum() );

  // As a simple start, we allow exactly one arbitrary missing plane
  //TODO: allow more general criteria
  const UInt_t kMaxmiss = 1;

  UInt_t nmiss = 0;
  for( UInt_t i = 0; i<fBuild->fNplanes; ++i ) {
    if( !curpat[i] ) {
      ++nmiss;
      if( nmiss > kMaxmiss )
	return false;
    }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t Road::Add( NodeDescriptor& nd )
{
  // Check if the hits from the given NodeDescriptor pattern are common
  // with the common hit set already in this road. If so, add the pattern
  // to this road, update the common hits if necessary, and return true.
  // If not, do nothing and return false.
  //
  // Adding only works as long as the road is not yet finished

  if( !fBuild )
    return kFALSE;

  nd.Print();

  bool first = fBuild->fPatterns.empty();
  // Adding a new pattern: find all the common hits that are also present
  // in the new hit set.

  // First, make a sorted list of all the new hits
  // TODO: Do this only once per pattern, perhaps when we make the patterns
  set<Hit*> new_hits;
  for( UInt_t i = fBuild->fNlayers; i; ) {
    --i;
    const vector<Hit*>& hits = fBuild->fHitpattern->GetHits( i, nd[i] );
    PrintHits(hits);
    copy( hits.begin(), hits.end(), inserter( new_hits, new_hits.begin() ));
  }

  cout << "adding " << new_hits.size() << " hits" <<endl;

  if( first ) {
    if( !CheckMatch(new_hits) )
      return kFALSE;
    swap( fHits, new_hits );
    fBuild->fCommonHits = fHits;
  } else {
    set<Hit*> new_commons;
    set_intersection( new_hits.begin(), new_hits.end(),
		      fBuild->fCommonHits.begin(), fBuild->fCommonHits.end(),
		      inserter( new_commons, new_commons.end() ));

    cout << "new/old commons = " << new_commons.size() << " "
	 << fBuild->fCommonHits.size() << endl;
    
    assert( new_commons.size() <= fBuild->fCommonHits.size() );
    if( new_commons.size() < fBuild->fCommonHits.size() ) {
      // The set of common hits shrank, so we must check if this would still
      // be a good road
      if( !CheckMatch(new_commons) ) {
	// The new pattern would reduce the set of common hits in the road to
	// too loose a fit, so we reject the new pattern and leave
	// the road as it is
	cout << "failed" << endl;
	return kFALSE;
      }
      // The new set of common hits is good, so update the build data
      swap( fBuild->fCommonHits, new_commons );
      set<Hit*> allhits;
      set_union( fHits.begin(), fHits.end(), new_hits.begin(), new_hits.end(),
		 inserter( allhits, allhits.begin() ));
      cout << "new/old nhits = " << allhits.size() << " " 
	   << fHits.size() << endl;
      swap( fHits, allhits );
    }
  }

  // Save a pointer to this pattern so we can update it later
  fBuild->fPatterns.push_back(&nd);
  cout << "new npat = " << fBuild->fPatterns.size() << endl;

  // Expand the road limits if necessary
  assert( nd.link->GetPattern()->GetNbits() == fBuild->fNlayers );
  UInt_t n = fBuild->fNlayers-1;
  fLeft[0]  = TMath::Min( nd[0], fLeft[0] );
  fLeft[1]  = TMath::Min( nd[n], fLeft[1] );
  fRight[0] = TMath::Max( nd[0], fRight[0] );
  fRight[1] = TMath::Max( nd[n], fRight[1] );

  cout << "new left/right = " << fLeft[0]<<" "<<fRight[0]<<" "
       << fLeft[1]<<" "<<fRight[1] << endl;

  return kTRUE;
}

//_____________________________________________________________________________
void Road::Finish()
{
  // Finish building the road

  assert(fBuild);
  for( vector<NodeDescriptor*>::iterator it = fBuild->fPatterns.begin();
       it != fBuild->fPatterns.end(); ++it ) {
    
    NodeDescriptor& nd = **it;
    assert( nd.used < 2 ); // cannot add previously fully used pattern

    // First, make a sorted list of all the new hits
    // TODO: Do this only once per pattern, perhaps when we make the patterns
    set<Hit*> these_hits;
    for( UInt_t i = fBuild->fNlayers; i; ) {
      --i;
      const vector<Hit*>& hits = fBuild->fHitpattern->GetHits( i, nd[i] );
      copy( hits.begin(), hits.end(), 
	    inserter( these_hits, these_hits.begin() ));
    }
    // TODO: search only to first element not in common
    vector<Hit*> not_in_common;
    set_difference( these_hits.begin(), these_hits.end(),
		    fBuild->fCommonHits.begin(), fBuild->fCommonHits.end(),
		    back_inserter( not_in_common ) );

    nd.used = not_in_common.empty() ? 2 : 1;
    cout << "used = " << (int)nd.used << "for ";
    nd.Print();
  }


  // Put the tools away
  delete fBuild; fBuild = 0;

  return;
}

//_____________________________________________________________________________
void Road::Print( Option_t* opt ) const
{
  // Print road info

}

//_____________________________________________________________________________
void CollectCoordinates()
{
  // 

}

#if 0

//_____________________________________________________________________________
MatchCount_t Hitpattern::CommonPlanes( const NodeDescriptor& nd1,
				       const NodeDescriptor& nd2,
				       const vector<UInt_t>& maxdist ) const
{
  // Return the number of physical wire planes (not layers) in which the
  // two patterns nd1 and nd2 share at least one hit

  //TODO: remove this, replaced by Road::Add()

  MatchCount_t matches(0,0);
  bool top = true;
  set<Hit*> sorted_hits;

  for( UInt_t i=fNplanes; i; ) {
    --i;
    UInt_t bin1 = nd1[i];
    UInt_t bin2 = nd2[i];

    // Quit if top bins too far apart
    if( top ) {
      if( (bin1 > bin2 && bin1-bin2 > maxdist[i]) ||
	  (bin1 < bin2 && bin2-bin1 > maxdist[i]) )
	break;
      top = false;
    }

    // Naturally, there are no matches if either hit list is empty
    vsiz_t sz1 = GetHits(i,bin1).size();
    vsiz_t sz2 = GetHits(i,bin2).size();
    if( sz1 == 0 or sz2 == 0 )
      continue;

    // Ensure that hits1 is smaller than hits2 so we have a better chance
    // of saving time later
    if( sz1 > sz2 )
      swap( bin1, bin2 );
    const vector<Hit*>& hits1 = GetHits(i,bin1);
    const vector<Hit*>& hits2 = GetHits(i,bin2);

//     PrintHits( hits1 );
//     PrintHits( hits2 );

    // Determine if we are dealing with a wire plane pair
    vector<Hit*>::const_iterator it = hits1.begin();
    const WirePlane* wp = (*it)->GetWirePlane();
    assert(wp);
    bool doing_pair = (wp->GetPartner() != 0);

    // Naturally, same bins always match
    if( bin1 == bin2 ) {
      ++matches.first;
      ++matches.second;
      if( doing_pair )
	++matches.second;
      continue;
    }
    
    // Fast search for identical elements in the two arrays.
    // This runs in N*log(N) time, like sorting all the hits would.
    // N = hits1.size() + hits2.size() is small anyway.

    // Put hits from the first bin into a sorted container, allowing only
    // unique elements
    sorted_hits.clear();
    copy( hits1.begin(), hits1.end(),
	  inserter( sorted_hits, sorted_hits.begin() ));

    // See if any hits in the second bin match any in the first.
    // This is not totally straightforward with partnered plane pairs.
    // In that case, this really is a parallel search for hits in two planes
    set<Hit*>::iterator found;
    Int_t found_plane = -1;
    for( it = hits2.begin(); it != hits2.end(); ++it ) {
      found = sorted_hits.find( *it );
      if( found != sorted_hits.end() ) {
	if( doing_pair ) {
	  // If this bin is for a plane pair, allow up to two matches	
	  const WirePlane* wp2 = (*found)->GetWirePlane();
	  assert(wp2);
	  if( found_plane < 0 ) {   // no plane found yet
	    found_plane = wp2->GetPlaneNum();
	    ++matches.first;
	    ++matches.second;
	  } else if( (Int_t)wp2->GetPlaneNum() != found_plane ) {
	    assert( wp2->GetPlaneNum() == wp->GetPartner()->GetPlaneNum() );
	    ++matches.second;
	    break;
	  }	    
	} else {
	  // If not a plane pair, we're done once we have the first match
	  ++matches.first;
	  ++matches.second;
	  break;
	}
      }
    }
  }
  return matches;
}
#endif

} // end namespace TreeSearch

///////////////////////////////////////////////////////////////////////////////
