//*-- Author :    Ole Hansen, Jefferson Lab   07-Feb-2008

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Road.h"

#include <iostream>

using namespace std;

ClassImp(TreeSearch::Road);

namespace TreeSearch {

//_____________________________________________________________________________
Road& Road::operator=( const Road& rhs )
{
  // Print hit info

  if( this != &rhs ) {

  }
  return *this;
}

//_____________________________________________________________________________
Bool_t Road::Add( NodeDescriptor* nd )
{
  // Check if the hits from the given NodeDescriptor pattern are common
  // with the common hit set already in this road. If so, add the pattern
  // to this road, update the common hits if necessary, and return true.
  // If not, do nothing and return false.



  return kFALSE;
}

//_____________________________________________________________________________
void Road::Close()
{
  // Finish building the road

  return;
}

//_____________________________________________________________________________
void Road::Print( Option_t* opt ) const
{
  // Print road info

}


#if 0
//_____________________________________________________________________________
// static void PrintHits( const vector<Hit*> hits )
// {
//   cout << hits.size() << " hits" << endl;

//   for( vector<Hit*>::size_type i=0; i<hits.size(); ++i ) {
//     cout << " ";
//     hits[i]->Print();
//   }
  
// }

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
