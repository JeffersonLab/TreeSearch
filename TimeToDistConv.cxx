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

//FIXME: flesh out these skeleton functions

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// LinearTTD                                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
LinearTTD::LinearTTD( Double_t drift_velocity ) : fDriftVel(drift_velocity)
{
  // Constructor. Set drift velocity (m/s).
}

//_____________________________________________________________________________
LinearTTD::~LinearTTD()
{
  // Destructor
}

//_____________________________________________________________________________
Double_t LinearTTD::ConvertTimeToDist( Double_t time, Double_t slope )
{
  // Time in s. Return distance in m. 
  
  // NB: 1/cos = sqrt(1+tan^2)
  return fDriftVel * time * TMath::Sqrt( 1.0 + slope*slope );
}


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch
