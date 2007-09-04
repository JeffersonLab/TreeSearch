///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TimeToDistConv                                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TimeToDistConv.h"
#include "TMath.h"

ClassImp(TreeSearch::TimeToDistConv)
ClassImp(TreeSearch::LinearTTD)

namespace TreeSearch {

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// LinearTTD                                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
LinearTTD::LinearTTD() : fDriftVel(kBig)
{
  // Constructor.
  fNparam = 1;
}

//_____________________________________________________________________________
Double_t LinearTTD::GetParameter( UInt_t i ) const
{
  // Get i-th parameter

  return i==0 ? GetDriftVel() : kBig;
}

//_____________________________________________________________________________
Int_t LinearTTD::SetParameters( const vector<double>& parameters )
{
  // Set parameters. The first element of parameter array is interpreted as
  // the drift velocity in m/s. Further elements are ignored.

  if( parameters.empty() )
    return -1;

  fDriftVel = parameters[0];
  return 0;
}

//_____________________________________________________________________________
Double_t LinearTTD::ConvertTimeToDist( Double_t time, Double_t slope ) const
{
  // Time in s. Return distance in m. slope is used to project the distance
  // of closest approach onto the wire plane (1/cos correction).
  
  // Don't bother with very small slopes
  if( TMath::Abs(slope) < 1e-2 )
    return fDriftVel * time;
  
  // NB: 1/cos = sqrt(1+tan^2)
  return fDriftVel * time *  TMath::Sqrt( 1.0 + slope*slope );
}


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch
