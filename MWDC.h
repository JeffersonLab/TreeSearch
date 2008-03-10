#ifndef ROOT_TreeSearch_MWDC
#define ROOT_TreeSearch_MWDC

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaTrackingDetector.h"
#include "THaDetMap.h"
#include "Types.h"
#include "TMatrixDSym.h"
#include <vector>

class THaTrack;
class THaBenchmark;
class TClonesArray;
class THashTable;

using std::vector;

namespace TreeSearch {
  class WirePlane;
  class Projection;
  class Road;

  extern const Double_t kBig;

  class MWDC : public THaTrackingDetector {
    friend class WirePlane;

  public:
    MWDC( const char* name, const char* description = "", 
	  THaApparatus* app = 0 );
    virtual ~MWDC();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    //    virtual Int_t   End(THaRunBase* run);
    virtual EStatus Init( const TDatime& date );
    virtual Int_t   CoarseTrack( TClonesArray& tracks );
    virtual Int_t   FineTrack( TClonesArray& tracks );
    virtual Int_t   DefineVariables( EMode mode = kDefine );
    virtual void    Print(const Option_t* opt) const;
    virtual void    SetDebug( Int_t level );
    void            EnableBenchmarks( Bool_t b = kTRUE );

    void            EnableEventDisplay( Bool_t enable = true );
    Double_t        GetRefTime( UInt_t i ) const 
    { return (i<(UInt_t)fRefMap->GetSize()) ? fRefTime[i] : kBig; }

    EProjType       NameToType( const char* name );

#ifdef TESTCODE
    Int_t           GetEvNum() const { return fEvNum; }
    vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }
    vector<TreeSearch::Projection*>& GetListOfProjections() { return fProj; }
#endif
    // Analysis control flags. Set via database.
    enum {
      kDoTimeCut  = BIT(14), // Use TDC time cut in decoder (defined in planes)
      kPairsOnly  = BIT(15), // Accept only pairs of hits from plane pairs
      kMCdata     = BIT(16), // Assume input is Monte Carlo data
      kNoPartner  = BIT(17), // Never partner wire planes
      k3dFastMatch= BIT(18), // Use fast 3D matching algorithm (auto detected)
      kEventDisplay = BIT(19)  // Support event display
    };

  protected:
    typedef vector<Road*> Rvec_t;

    vector<WirePlane*>   fPlanes;  // Wire planes
    vector<Projection*>  fProj;    // Plane projections

    THaDetMap*      fRefMap;    // Map of reference channels for VME readout
    vector<float>   fRefTime;   // [fRefMap->GetSize()] ref channel data

    THashTable*     fCrateMap;  // Map of MWDC DAQ modules

    // Paremeters for 3D projection matching
    Double_t        f3dMatchvalScalefact; // Correction for fast 3D matchval
    Double_t        f3dMatchCut;          // Maximum allowed 3D match error

    // Event data
#ifdef TESTCODE
    Int_t           fEvNum;     // Current event number (for diagnostics)
#endif
    static Double_t CalcChisquare( const Rvec_t& roads,
				   const vector<Double_t>& coef );
    static Int_t    FitTrack( const Rvec_t& roads, vector<Double_t>& coef,
			      Double_t& chi2, TMatrixDSym* coef_covar = 0 );

    // Helper functions for getting DAQ module parameters - used by Init
    UInt_t   LoadDAQmodel( THaDetMap::Module* m ) const;
    Double_t LoadDAQresolution( THaDetMap::Module* m ) const;
    UInt_t   GetDAQnchan( THaDetMap::Module* m ) const;

    // Podd interface
    virtual Int_t   ReadDatabase( const TDatime& date );

    ClassDef(MWDC,0)   // Tree search reconstruction of BigBite MWDCs
  };

}

///////////////////////////////////////////////////////////////////////////////

#endif
