#ifndef ROOT_TreeSearch_TimeToDistConv
#define ROOT_TreeSearch_TimeToDistConv

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::TimeToDistConv                                                //
//                                                                           //
// Base class for algorithms for converting drift time into  drift distance  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

namespace TreeSearch {

  class TimeToDistConv {

  public:
    virtual ~TimeToDistConv() {}

    virtual Double_t ConvertTimeToDist( Double_t time, Double_t slope ) = 0;

  protected:

    TimeToDistConv() {}
    TimeToDistConv( const TimeToDistConv& );
    TimeToDistConv& operator=( const TimeToDistConv& );

    ClassDef(TimeToDistConv,0)     // Drift time to diatance converter ABC
  };


  //___________________________________________________________________________
  // LinearTTD
  //
  // Simple linear conversion of drit time (s) and drift distance (m).
  // dist = drift_velocity * time.
  
  class LinearTTD : public TimeToDistConv {

  public:
    LinearTTD( Double_t drift_velocity );
    virtual ~LinearTTD();

    virtual Double_t ConvertTimeToDist( Double_t time, Double_t slope );

    Double_t GetDriftVel() { return fDriftVel; }

protected:

    Double_t fDriftVel;   // Drift velocity (m/s)

    ClassDef(LinearTTD,0)  // Linear drift time-to-distance converter
  };

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
