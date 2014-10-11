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
ClassImp(TreeSearch::HitSet)
#ifdef MCDATA
ClassImp(TreeSearch::MCWireHit)
#endif

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
UInt_t WireHit::GetNumPos() const
{
  // Number of positions available from this hit. Normally 2 for wire hits,
  // unless the drift time is zero.

  return (GetDriftDist() == 0.0) ? 1 : 2;
}

//_____________________________________________________________________________
Double_t WireHit::GetPosI( UInt_t i ) const
{
  // Get the i-th position. i can be 0 or 1, where 0 ==> fPosR, 1 == fPosL

  assert( i<2 );
  return (i != 0) ? GetPosL() : GetPosR();
}

//_____________________________________________________________________________
Double_t WireHit::ConvertTimeToDist( Double_t slope )
{
  // Convert drift time to drift distance. 'slope' is the approximate
  // slope of the track.
  // Updates the internal variables fPosL and fPosR.
  // Must be called before doing analysis of drift chamber hits.

  assert( dynamic_cast<WirePlane*>(fPlane) );
  WirePlane* wp = static_cast<WirePlane*>(fPlane);
  Double_t dist = wp->GetTTDConv()->ConvertTimeToDist(fTime, slope);
  fPosL = fPos-dist;
  fPosR = fPos+dist;
  return dist;
}

#ifdef MCDATA
//_____________________________________________________________________________
void MCWireHit::Print( Option_t* ) const
{
  // Print hit info

  WireHit::Print("C");
  MCPrint();
}
#endif

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

