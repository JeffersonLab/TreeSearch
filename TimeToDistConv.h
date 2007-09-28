#ifndef ROOT_TreeSearch_TimeToDistConv
#define ROOT_TreeSearch_TimeToDistConv

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::TimeToDistConv                                                //
//                                                                           //
// Base class for algorithms for converting drift time to drift distance     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>

using std::vector;

namespace TreeSearch {

  extern const Double_t kBig;

  class TimeToDistConv {

  public:
    virtual ~TimeToDistConv() {}

    virtual Double_t ConvertTimeToDist( Double_t time, 
					Double_t slope ) const = 0;
            Int_t    GetNparam() const { return fNparam; }
    virtual Double_t GetParameter( UInt_t i ) const { return kBig; }
    virtual Int_t    SetParameters( const vector<double>& param ) { return 0; }

  protected:

    TimeToDistConv() : fNparam(0) {}
    TimeToDistConv( const TimeToDistConv& );
    TimeToDistConv& operator=( const TimeToDistConv& );

    Int_t  fNparam;    // Number of parameters

    ClassDef(TimeToDistConv,0)     // Drift time to diatance converter ABC
  };


  //___________________________________________________________________________
  // LinearTTD
  //
  // Simple linear conversion of drift time (s) and drift distance (m).
  
  class LinearTTD : public TimeToDistConv {

  public:
    LinearTTD();
    virtual ~LinearTTD() {}

    virtual Double_t ConvertTimeToDist( Double_t time, Double_t slope ) const;
    virtual Int_t    SetParameters( const vector<double>& param );

            Double_t GetDriftVel() const { return fDriftVel; }
    virtual Double_t GetParameter( UInt_t i ) const;

protected:

    Double_t fDriftVel;   // Drift velocity (m/s)

    ClassDef(LinearTTD,0)  // Linear drift time-to-distance converter
  };

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
