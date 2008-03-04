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
    virtual Int_t   FindVertices(TClonesArray&) {return 0;}
    virtual void    Print(const Option_t* opt) const;
    virtual void    SetDebug( Int_t level );
    void            EnableBenchmarks( Bool_t b = kTRUE );

    Double_t        GetRefTime( UInt_t i ) const 
    { return (i<(UInt_t)fRefMap->GetSize()) ? fRefTime[i] : kBig; }

    EProjType       NameToType( const char* name );

    //FIXME: for testing
    vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }
    vector<TreeSearch::Projection*>& GetListOfProjections() { return fProj; }

    // Analysis control flags. Set via database. Defaults are usually fine.
    enum {
      kDoTimeCut  = BIT(14), // Use TDC time cut in decoder (defined in planes)
      kPairsOnly  = BIT(15), // Accept only pairs of hits from plane pairs
      kMCdata     = BIT(16), // Assume input is Monte Carlo data
      kNoPartner  = BIT(17)  // Never partner wire planes
    };

  protected:

    typedef vector<TreeSearch::Road*> Rvec_t;

    vector<WirePlane*>   fPlanes;  // Wire planes
    vector<Projection*>  fProj;    // Plane projections

    THaDetMap*      fRefMap;    // Map of reference channels for VME readout
    vector<float>   fRefTime;   // [fRefMap->GetSize()] ref channel data

    THashTable*     fCrateMap;  // Map of MWDC DAQ modules

    // Paremeters for 3D projection matching
    Bool_t          f3dSimpleMatch;       // Use simplified algorithm
    Double_t        f3dMatchvalScalefact; // Correction for simple 3D matchval
    Double_t        f3dMatchCut;          // Maximum allowed 3D match error

    static  Double_t CalcChisquare( const Rvec_t& rds,
				    const vector<Double_t>& coef );
    static  Int_t   FitTrack( const Rvec_t& roads, vector<Double_t>& coef,
			      Double_t& chi2 );
    virtual Int_t   ReadDatabase( const TDatime& date );

    // Helper functions for getting DAQ module parameters - used by Init
    UInt_t   LoadDAQmodel( THaDetMap::Module* m ) const;
    Double_t LoadDAQresolution( THaDetMap::Module* m ) const;
    UInt_t   GetDAQnchan( THaDetMap::Module* m ) const;

    ClassDef(MWDC,0)   // Tree search reconstruction of BigBite MWDCs
  };

}

///////////////////////////////////////////////////////////////////////////////

#endif
