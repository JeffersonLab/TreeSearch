///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TimeToDistConv                                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TimeToDistConv.h"


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
  
  // OK, this _is_ simple... 
  // A slightly better version could apply a slope correction.
 
  return fDriftVel * time;
}


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch
