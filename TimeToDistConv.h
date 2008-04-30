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
            UInt_t   GetNparam() const { return fNparam; }
    virtual Double_t GetParameter( UInt_t ) const { return kBig; }
    virtual Int_t    SetParameters( const vector<double>& ) { return 0; }

  protected:

    TimeToDistConv( UInt_t npar = 0 ) : fNparam(npar) {}
    TimeToDistConv( const TimeToDistConv& rhs ) : fNparam(rhs.fNparam) {}
    TimeToDistConv& operator=( const TimeToDistConv& rhs )
    {
      if( this != &rhs ) {
	fNparam = rhs.fNparam;
      }
      return *this;
    }

    UInt_t  fNparam;    // Number of parameters

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

  //___________________________________________________________________________
  // TanhFitTTD
  //
  // Conversion of drift time (s) to distance using fit to tanh function
  
  class TanhFitTTD : public TimeToDistConv {

  public:
    TanhFitTTD();
    virtual ~TanhFitTTD() {}

    virtual Double_t ConvertTimeToDist( Double_t time, Double_t slope ) const;
    virtual Double_t GetParameter( UInt_t i ) const;
    virtual Int_t    SetParameters( const vector<double>& param );

protected:

    Double_t fDriftVel;   // Drift velocity (m/s)
    Double_t fC0;         // c0 (m)
    Double_t fC2;         // c2 (m/s^2)
    Double_t fT0;         // t0 (s)
    Double_t fInvC0;      // 1/c0 (1/m) for efficiency

    ClassDef(TanhFitTTD,0)  // TanhFit drift time-to-distance converter
  };

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
