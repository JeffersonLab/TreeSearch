//*-- Author :    Ole Hansen, Jefferson Lab   11-Feb-2013

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::WireHit                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "WireHit.h"
#include "TimeToDistConv.h"

#include <iostream>

using std::cout;
using std::endl;

ClassImp(TreeSearch::WireHit)
ClassImp(TreeSearch::MCWireHit)
ClassImp(TreeSearch::HitSet)

namespace TreeSearch {

//_____________________________________________________________________________
void WireHit::Print( Option_t* opt ) const
{
  // Print hit info

  cout << "Hit: wire=" << GetWireNum()
       << "/" << (fPlane ? fPlane->GetName() : "??")
       << " wpos="     << GetPos()
       << " z="        << GetZ()
       << " res="      << GetResolution()
       << " time="     << GetDriftTime()
       << " drift="    << GetDriftDist();
  if( *opt != 'C' )
    cout << endl;
}

//_____________________________________________________________________________
Double_t WireHit::ConvertTimeToDist( Double_t slope )
{
  // Convert drift time to drift distance. 'slope' is the approximate
  // slope of the track.
  // Updates the internal variables fPosL and fPosR.
  // Must be called before doing analysis of drift chamber hits.

#ifdef NDEBUG
  WirePlane* wp = static_cast<WirePlane*>(fPlane);
#else
  WirePlane* wp = dynamic_cast<WirePlane*>(fPlane);
  assert(wp);
#endif
  Double_t dist = wp->GetTTDConv()->ConvertTimeToDist(fTime, slope);
  fPosL = fPos-dist;
  fPosR = fPos+dist;
  return dist;
}

//_____________________________________________________________________________
void MCWireHit::Print( Option_t* ) const
{
  // Print hit info

  WireHit::Print("C");
  cout << " MCpos=" << GetMCPos()
       << endl;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

