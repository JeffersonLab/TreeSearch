//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 31-Jan-2013
//
#ifndef ROOT_TreeSearch_Tracker
#define ROOT_TreeSearch_Tracker

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Tracker                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaTrackingDetector.h"
#include "THaDetMap.h"
#include "Types.h"
#include <vector>
#include <utility>
#include <set>
#include <list>
#include <cassert>
#include <exception>
#include "TMatrixDSym.h"
#include "TRotation.h"
#ifndef MCDATA
#include "THaEvData.h"
#else
#include "SimDecoder.h"
#endif

class THaTrack;
class THaBenchmark;
class TClonesArray;
class THashTable;
class TClass;
class TSeqCollection;

using std::vector;

namespace TreeSearch {
  class Plane;
  class Projection;
  class Road;
  class Hit;
  class ThreadCtrl;  // Defined in implementation

  typedef std::vector<Road*> Rvec_t;
  typedef std::set<Road*>    Rset_t;
  typedef std::vector<UInt_t> vec_uint_t;

  extern const Double_t kBig;

  class Tracker : public THaTrackingDetector {
  public:
    Tracker( const char* name, const char* description = "",
	     THaApparatus* app = 0 );
    virtual ~Tracker();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );
    virtual Int_t   CoarseTrack( TClonesArray& tracks );
    virtual Int_t   FineTrack( TClonesArray& tracks );
    virtual Int_t   DefineVariables( EMode mode = kDefine );
    virtual void    Print(const Option_t* opt) const;
    virtual void    SetDebug( Int_t level );

    virtual Int_t   Begin( THaRunBase* r=0 );
    virtual Int_t   End( THaRunBase* r=0 );

    void            EnableEventDisplay( Bool_t enable = true );
    const pdbl_t&   GetChisqLimits( UInt_t i ) const;
    const TRotation& GetRotation()     const { return fRotation; }
    const TRotation& GetInvRotation()  const { return fInvRot; }
    Bool_t          IsRotated()        const { return fIsRotated; }

    // Analysis control flags. Set via database.
    enum {
#ifdef MCDATA
      kMCdata        = BIT(16), // Assume input is Monte Carlo data
#endif
      k3dFastMatch   = BIT(18), // Use fast 3D matching algo (auto detected)
      kEventDisplay  = BIT(19), // Support event display
      kDoCoarse      = BIT(20), // Do coarse tracking (if unset, decode only)
      kDoFine        = BIT(21), // Do fine tracking (implies kDoCoarse)
      kDoChi2        = BIT(22), // Apply chi2 cut to 3D tracks
      kProjTrackToZ0 = BIT(23)  // Project tracks to global z = 0
    };

#ifdef TESTCODE
    Int_t           GetEvNum() const { return fEvNum; }
    vector<TreeSearch::Plane*>& GetListOfPlanes() { return fPlanes; }
    vector<TreeSearch::Projection*>& GetListOfProjections() { return fProj; }
#endif

    enum ETrackingStatus {
      kTrackOK             = 0, // At least one seemingly good track found
      // Decode
      kTooManyRawHits      = 1, // Raw occupancy too high
      // Coarse Track
      kNoTrackingRequested = 2, // CoarseTracking turned off
      kProjTrackError      = 3, // Failure in Projection::Track()
      kTooFewProj          = 4, // Too few projections with tracks
      kFailed3DMatch       = 5, // No roads matched in 3D
      kFailedTrackCuts     = 6, // All road combos failed ndof & chi2 cuts
      kFailedOptimalN      = 7, // Failed 3D de-ghosting algorithm
      // MatchRoads
      kTooManyRoadCombos   = 8, // Overflow in MatchRoads
      kNoRoadCombos        = 9  // Product of road vector sizes = 0 (bug?)
    };
    ETrackingStatus GetTrackingStatus() const { return fTrkStat; }

#ifdef MCDATA
    // Bad configuration exception, may be thrown by Decode
    class bad_config : public std::runtime_error {
    public:
      bad_config( const std::string& what_arg )
	: std::runtime_error(what_arg) {}
    };
#endif

  protected:
    friend class Plane;
    class TrackFitWeight;
    struct FitRes_t {
      vector<Double_t> coef;
      Double_t matchval;
      Double_t chi2;
      Int_t    ndof;
      Rvec_t*  roads;
      FitRes_t() : roads(0) {}
    };

    typedef std::vector<Plane*>  Rpvec_t;
    typedef std::vector<Projection*> Prvec_t;

    // Planes and projections
    Rpvec_t        fPlanes;           // Readout planes
    Prvec_t        fProj;             // Plane projections
    Rpvec_t        fCalibPlanes;      // Planes in calibration mode

    // Configuration and parameters
    THashTable*    fCrateMap;         // Map of DAQ modules
    Double_t       fMinProjAngleDiff; // Min coord axis angle diff required
    TRotation      fRotation;         // Rotation Tracker -> lab frame
    TRotation      fInvRot;           // Rotation lab frame -> Tracker
    Bool_t         fIsRotated;        // Tracker frame is rotated
    Bool_t         fAllPartnered;     // All planes have partners

    // Multithread support
    UInt_t         fMaxThreads;       // Maximum simultaneously active threads
    ThreadCtrl*    fThreads;          //! Thread controller

    // Parameters for 3D projection matching
    UInt_t         fMinReqProj;  // Minimum # proj required for 3D match
    Double_t       f3dMatchvalScalefact; // Correction for fast 3D matchval
    Double_t       f3dMatchCut;  // Maximum allowed 3D match error
    vec_uint_t     f3dIdx;       // Lookup table proj index -> fast 3d index

    // Track fit cut parameters
    Int_t          fMinNdof;     // Minimum number of points in fit-4
    vec_pdbl_t     fChisqLimits; // lo/hi confidence interval limits on Chi2

    // Event-by-event data
    ETrackingStatus fTrkStat;    // Reconstruction status

    // Only needed for TESTCODE, but kept for binary compatibility
    UInt_t         fNcombos;     // # of road combinations tried
    UInt_t         fN3dFits;     // # of track fits done (=good road combos)
    Int_t          fEvNum;       // Current event number
    Double_t       t_track, t_3dmatch, t_3dfit, t_coarse; // times in us
    UInt_t         fN3dMatchs;   // # of hits (amp corr) matched for 3D track reconstruction

    void      Add3dMatch( const Rvec_t& selected, Double_t matchval,
			  std::list<std::pair<Double_t,Rvec_t> >& combos_found,
			  Rset_t& unique_found ) const;
    void      FitErrPrint( Int_t err ) const;
    Int_t     FitTrack( const Rvec_t& roads, vector<Double_t>& coef,
			Double_t& chi2, TMatrixDSym* coef_covar = 0 ) const;
    template< typename Action > static
    Action    ForAllTrackPoints( const Rvec_t& roads,
				 const vector<Double_t>& coef, Action action );
    THaTrack* NewTrack( TClonesArray& tracks, const FitRes_t& fit_par );
    Bool_t    PassTrackCuts( const FitRes_t& fit_par ) const;

    UInt_t    MatchRoadsGeneric( vector<Rvec_t>& roads, UInt_t ncombos,
		   std::list<std::pair<Double_t,Rvec_t> >& combos_found,
		   Rset_t& unique_found );

    UInt_t    MatchRoadsFast3D( vector<Rvec_t>& roads, UInt_t ncombos,
	           std::list<std::pair<Double_t,Rvec_t> >& combos_found,
	           Rset_t& unique_found );

    // Virtualization of the tracker class, specialized Trackers may/must
    // override
    virtual Plane* MakePlane( const char* name, const char* description = "",
			      THaDetectorBase* parent = 0 ) const = 0;
    virtual Projection* MakeProjection( EProjType type, const char* name,
					Double_t angle,
					THaDetectorBase* parent ) const;
    virtual UInt_t GetCrateMapDBcols() const = 0;

    virtual UInt_t MatchRoads( vector<Rvec_t>& roads,
	         std::list<std::pair<Double_t,Rvec_t> >& combos_found,
		 Rset_t& unique_found );
    virtual UInt_t MatchRoadsImpl( vector<Rvec_t>& roads, UInt_t ncombos,
                 std::list<std::pair<Double_t,Rvec_t> >& combos_found,
                 Rset_t& unique_found ) = 0;

    virtual THaAnalysisObject::EStatus PartnerPlanes() = 0;

    virtual Int_t NewTrackCalc( Int_t idx, THaTrack* newTrack,
				const TVector3& pos, const TVector3& dir,
				const FitRes_t& fit_par );


    // Helper functions for getting DAQ module parameters - used by Init
    UInt_t    LoadDAQmodel( THaDetMap::Module* m ) const;
    Double_t  LoadDAQresolution( THaDetMap::Module* m ) const;
    UInt_t    GetDAQnchan( THaDetMap::Module* m ) const;

    // Podd interface
    virtual Int_t  ReadDatabase( const TDatime& date );

#ifdef MCDATA
    // Monte Carlo data support
    class MCPointUpdater;
    const Podd::SimDecoder* fMCDecoder; //! MC data decoder (if kMCdata)
    MCPointUpdater* fMCPointUpdater;    //! MC track point update (if kMCdata)
    Bool_t          fChecked;           // Configuration checked to be good
    vec_uint_t      fMCHitBits;         // Per-track MC hit planepattern
    vec_uint_t      f3dMatchBits;       // Per-track MC hit planepattern

    virtual Hit*  FindHitForMCPoint( Podd::MCTrackPoint* pt,
				     MCPointUpdater* updater ) const;
    virtual THaTrack* FindTrackForMCPoint( Podd::MCTrackPoint* pt,
					   TClonesArray& tracks,
					   MCPointUpdater* updater ) const;
    virtual Int_t FitMCPoints( Podd::MCTrack* mctrk ) const;

    // Helper class to use with FindHitForMCPoint & FindTrackForMCPoint
    class MCPointUpdater {
    public:
      virtual ~MCPointUpdater() {}
      virtual void UpdateHit( Podd::MCTrackPoint* pt, Hit* hit,
			      Double_t x ) const;
      virtual void UpdateTrack( Podd::MCTrackPoint* pt, Plane* pl,
				THaTrack* track, Double_t x ) const;
    };
#endif

  private:
    // Shared between Init and ReadDatabase
    Int_t    fDBmaxmiss;
    Double_t fDBconf_level;

    ClassDef(Tracker,0)   // Tracking system analyzed using TreeSearch reconstruction
  };

  //___________________________________________________________________________
  inline const pdbl_t& Tracker::GetChisqLimits( UInt_t i ) const
  {
    assert( i < fChisqLimits.size() );
    return fChisqLimits[i];
  }

}

///////////////////////////////////////////////////////////////////////////////

#endif
