//*-- Author :    Ole Hansen, Jefferson Lab   11-Feb-2013

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hit.h"
#include "Road.h"

#include <iostream>

using std::cout;
using std::endl;
using std::make_pair;

ClassImp(TreeSearch::Hit)
ClassImp(TreeSearch::MCHit)
ClassImp(TreeSearch::HitSet)

namespace TreeSearch {

//_____________________________________________________________________________
void Hit::Print( Option_t* opt ) const
{
  // Print hit info

  cout << "Hit: plane=" << (fPlane ? fPlane->GetName() : "??")
       << " pos="       << GetPos()
       << " z="         << GetZ()
       << " res="       << GetResolution();
  if( *opt != 'C' )
    cout << endl;
}

//_____________________________________________________________________________
void MCHit::Print( Option_t* ) const
{
  // Print hit info

  Hit::Print("C");
  cout << " MCpos=" << GetMCPos()
       << endl;
}

//_____________________________________________________________________________
Double_t FitCoord::GetChi2() const
{
  // Return chi2 of the fit that used this coordinate

  return fRoad ? fRoad->GetChi2() : kBig;
}

//_____________________________________________________________________________
UInt_t HitSet::GetMatchValue( const Hset_t& hits )
{
  // Return plane occupancy pattern of given hitset

  UInt_t curpat = 0;
  for( Hset_t::const_iterator it = hits.begin(); it != hits.end(); ++it )
    curpat |= 1U << (*it)->GetPlaneNum();

  return curpat;
}

//_____________________________________________________________________________
Bool_t HitSet::IsSimilarTo( const HitSet& tryset, Int_t /* maxdist */ ) const
{
  // If maxdist == 0:
  // Similar to STL includes() algorithm, but allows tryset to have additional
  // hits in a given wire plane if there is at least one included hit in that
  // plane.
  //
  // Example: the following matches, despite the extra hit in "try":
  //   this:  30/   32/40/50/51
  //   try:   --/31 32/40/50/51
  // 
  // Standard includes() implies intersection == set2.
  // This algorithm tests planepattern(intersection) == planepattern(set2)
  //
  // If maxdist > 0:
  // Same as above, but consider hits "equal" not only if they are identical
  // but also if their wire numbers are at most maxdist apart.
  //
  // Example: the following two patters are "similar" with maxdist = 1:
  //   this:  30/32/40/50/51
  //   try:   31/32/40/50/51
  //
  // This mode can be used to build "clusters" of patterns.

  assert( tryset.plane_pattern );

  Hset_t::const_iterator ihits = hits.begin();
  Hset_t::const_iterator ehits = hits.end();
  Hset_t::const_iterator itry  = tryset.hits.begin();
  Hset_t::const_iterator etry  = tryset.hits.end();
  //  Hset_t::key_compare    comp  = hits.key_comp();
  Hit::PosIsLess comp;

  UInt_t intersection_pattern = 0;

  while( ihits != ehits and itry != etry ) {
    if( comp(*itry, *ihits) )
      ++itry;
    else if( comp(*ihits, *itry) )
      ++ihits;
    else {
      intersection_pattern |= 1U << (*itry)->GetPlaneNum();
      ++ihits;
      ++itry;
    }
  }
  return tryset.plane_pattern == intersection_pattern;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

