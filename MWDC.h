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
#include <utility>
#include <set>
#include <list>
#include <cassert>

class THaTrack;
class THaBenchmark;
class TClonesArray;
class THashTable;

using std::vector;

namespace TreeSearch {
  class WirePlane;
  class Projection;
  class Road;
  class ThreadCtrl;  // Defined in implementation

  typedef std::vector<Road*> Rvec_t;
  typedef std::set<Road*>    Rset_t;

  extern const Double_t kBig;

  class MWDC : public THaTrackingDetector {
  public:
    MWDC( const char* name, const char* description = "", 
	  THaApparatus* app = 0 );
    virtual ~MWDC();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );
    virtual Int_t   CoarseTrack( TClonesArray& tracks );
    virtual Int_t   FineTrack( TClonesArray& tracks );
    virtual Int_t   DefineVariables( EMode mode = kDefine );
    virtual void    Print(const Option_t* opt) const;
    virtual void    SetDebug( Int_t level );

    void            EnableEventDisplay( Bool_t enable = true );
    const pdbl_t&   GetChisqLimits( UInt_t i ) const;
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
      kEventDisplay = BIT(19), // Support event display
      kDoCoarse   = BIT(20), // Do coarse tracking (if disabled, decode only)
      kDoFine     = BIT(21)  // Do fine tracking (implies kDoCoarse)
    };

  protected:
    friend class WirePlane;
    class TrackFitWeight;
    struct FitRes_t {
      vector<Double_t> coef;
      Double_t matchval;
      Double_t chi2;
      Int_t    ndof;
      Rvec_t*  roads;
      FitRes_t() : roads(0) {}
    };

    typedef std::vector<WirePlane*>  Wpvec_t;
    typedef std::vector<Projection*> Prvec_t;

    Wpvec_t        fPlanes;      // Wire planes
    Prvec_t        fProj;        // Plane projections
    Wpvec_t        fCalibPlanes; // Planes in calibration mode

    THaDetMap*     fRefMap;      // Map of reference channels for VME readout
    vector<float>  fRefTime;     // [fRefMap->GetSize()] ref channel data

    THashTable*    fCrateMap;    // Map of MWDC DAQ modules

    // Multithread support
    Int_t          fMaxThreads;  // Maximum simultaneously active threads
    ThreadCtrl*    fThreads;     //! Thread controller

    // Paremeters for 3D projection matching
    Double_t       f3dMatchvalScalefact; // Correction for fast 3D matchval
    Double_t       f3dMatchCut;          // Maximum allowed 3D match error

    // Track fit cut parameters
    Int_t          fMinNdof;     // Minimum number of points in fit - 4
    vec_pdbl_t     fChisqLimits; // lo/hi onfidence interval limits on Chi2

    // Event data
    Int_t          fFailNhits; // Too many hits in wire plane(s)
    Int_t          fFailNpat;  // Too many treesearch patterns found
#ifdef TESTCODE
    UInt_t         fNcombos;   // Number of road combinations tried
    UInt_t         fN3dFits;   // Number of track fits done (=good road combos)
    Int_t          fEvNum;     // Current event number
    Double_t       t_track, t_3dmatch, t_3dfit, t_coarse; // times in us
#endif

    void      Add3dMatch( const Rvec_t& selected, Double_t matchval,
			  std::list<std::pair<Double_t,Rvec_t> >& combos_found,
			  Rset_t& unique_found ) const;
    void      FindNearestHits( WirePlane* wp, const THaTrack* track,
			       const Rvec_t& roads ) const;
    void      FitErrPrint( Int_t err ) const;
    Int_t     FitTrack( const Rvec_t& roads, vector<Double_t>& coef,
			Double_t& chi2, TMatrixDSym* coef_covar = 0 ) const;
    template< typename Action > static
    Action    ForAllTrackPoints( const Rvec_t& roads, 
				 const vector<Double_t>& coef, Action action );
    UInt_t    MatchRoads( vector<Rvec_t>& roads,
			  std::list<std::pair<Double_t,Rvec_t> >& combos_found,
			  Rset_t& unique_found );
			  
    THaTrack* NewTrack( TClonesArray& tracks, const FitRes_t& fit_par );
    Bool_t    PassTrackCuts( const FitRes_t& fit_par ) const;

    // Helper functions for getting DAQ module parameters - used by Init
    UInt_t    LoadDAQmodel( THaDetMap::Module* m ) const;
    Double_t  LoadDAQresolution( THaDetMap::Module* m ) const;
    UInt_t    GetDAQnchan( THaDetMap::Module* m ) const;

    // Podd interface
    virtual Int_t   ReadDatabase( const TDatime& date );

    ClassDef(MWDC,0)   // Tree search reconstruction of BigBite MWDCs
  };

  //___________________________________________________________________________
  inline const pdbl_t& MWDC::GetChisqLimits( UInt_t i ) const
  { 
    assert( i < fChisqLimits.size() );
    return fChisqLimits[i];
  }

}

///////////////////////////////////////////////////////////////////////////////

#endif
