#ifndef ROOT_TreeSearch_Projection
#define ROOT_TreeSearch_Projection

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaAnalysisObject.h"
#include "TreeWalk.h"   // for NodeVisitor
#include "Hit.h"        // for Node_t
#include "Types.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TVector2.h"
#include <vector>
#include <set>
#include <cassert>
#include <functional>
#include <list>
#include <exception>
#include <string>
#include <unordered_map>

class THaDetectorBase;
class TBits;

namespace TreeSearch {

  class Hitpattern;
  class PatternTree;
  class Road;
  class Plane;

  typedef std::vector<Plane*>            vpl_t;
  typedef std::vector<Plane*>::size_type vplsiz_t;
  typedef std::vector<Plane*>::iterator  vpliter_t;

  class Projection : public THaAnalysisObject {
  public:

    Projection( EProjType type, const char* name, Double_t angle,
		THaDetectorBase* parent );
    Projection() : fType(kUndefinedType), fDetector(0), fPatternTree(0),
		   fPlaneCombos(0), fAltPlaneCombos(0), fHitpattern(0),
		   fRoads(0), fRoadCorners(0) {} // ROOT RTTI
    virtual ~Projection();

    void            AddPlane( Plane* pl, Plane* partner = 0 );
    void            AddDummyPlane( Plane* pl, Plane* partner = 0 );
    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );
    // EStatus         InitLevel2( const TDatime& date );
    virtual void    Print( Option_t* opt="" ) const;
    void            Reset( Option_t* opt="" );

    Int_t           FillHitpattern();
    Int_t           Track();
    Int_t           MakeRoads();
    Int_t           RemoveNonElasticTracks(); //remove false tracks and non elasticTracks base on x(y)/x'(y') correlation

    static EProjType NameToType( const char* name );

    Double_t        GetAngle()        const;
    const TVector2& GetAxis()         const { return fAxis; }
    UInt_t          GetBinMaxDistB()  const { return fBackMaxBinDist; }
    UInt_t          GetBinMaxDistF()  const { return fFrontMaxBinDist; }
    const pdbl_t&   GetChisqLimits( UInt_t i ) const;
    Double_t        GetConfLevel()    const { return fConfLevel; }
    Double_t        GetCosAngle()     const { return fAxis.X(); }
    UInt_t          GetHitMaxDist()   const { return fHitMaxDist; }
    Hitpattern*     GetHitpattern()   const { return fHitpattern; }
    Double_t        GetMaxSlope()     const { return fMaxSlope; }
    UInt_t          GetMinFitPlanes() const { return fMinFitPlanes; }
    UInt_t          GetNgoodRoads()   const { return fNgoodRoads; }
    UInt_t          GetNlevels()      const { return fNlevels; }
    UInt_t          GetNpatterns()    const;
    UInt_t          GetNplanes()      const { return (UInt_t)fPlanes.size(); }
    UInt_t          GetNroads()       const;
    TBits*          GetPlaneCombos()  const { return fPlaneCombos; }
    Plane*          GetPlane ( UInt_t plane ) const { return fPlanes[plane]; }
    Double_t        GetPlaneZ( UInt_t plane ) const;
    Road*           GetRoad  ( UInt_t i )     const;
    Double_t        GetSinAngle()     const { return fAxis.Y(); }
    EProjType       GetType()         const { return fType; }
    Double_t        GetWidth()        const { return fWidth; }
    Double_t        GetZsize()        const;

    // Dummy plane support
    UInt_t          GetDummyPlanePattern() const { return fDummyPlanePattern; }
    UInt_t          GetFirstPlaneNum()     const { return fFirstPlaneNum; }
    UInt_t          GetLastPlaneNum()      const { return fLastPlaneNum; }

    void            SetPatternTree( PatternTree* pt ) { fPatternTree = pt; }

    const vpl_t&    GetListOfPlanes() const { return fPlanes; }

    Bool_t          DoingChisqTest() const  { return TestBit(kDoChi2); }

    // Analysis control flags
    enum {
      kEventDisplay = BIT(14), // Support event display
      kHaveDummies  = BIT(15), // Dummy planes present
      kDoChi2       = BIT(22)  // Apply chi2 cut to 2D fits
#ifdef MCDATA
    , kMCdata       = BIT(23)  // Assume input is Monte Carlo data
#endif
    };

    // Helper functors for STL algos
    struct ByType
      : public std::binary_function< Projection*, Projection*, bool >
    {
      bool operator() ( const Projection* a, const Projection* b ) const
      {
	assert( a && b );
	return( a->GetType() < b->GetType() );
      }
    };
    class TypeEquals : public std::unary_function< Projection*, bool > {
    public:
      TypeEquals( EProjType t ) : type(t) {}
      bool operator() ( const Projection* p ) const
      { assert(p); return ( p->GetType() == type ); }
    private:
      EProjType type;
    };

    // Bad angle exception, may be thrown by SetAngle and constructors
    class bad_angle : public std::invalid_argument {
    public:
      bad_angle( const std::string& what_arg )
	: std::invalid_argument(what_arg) {}
    };

    enum ETrackingStatus {
      kTrackOK             = 0,
      // Track
      kNoPatterns          = 1, // TreeSearch found no patterns
      kTooManyPatterns     = 2, // TreeSearch found too many patterns
      // FitRoads
      kFailed2DFits        = 3, // No roads with good fits
    };
    ETrackingStatus GetTrackingStatus() const { return fTrkStat; }

  protected:
    typedef std::vector<Node_t*> NodeVec_t;

    // Configuration
    EProjType        fType;          // Projection type (u,v,x,y...)
    vpl_t            fPlanes;     // Active (non-dummy) planes in this projection
    vpl_t            fAllPlanes;  // All planes in this projection (incl dummies)
    UInt_t           fNlevels;       // Number of levels of search tree
    Double_t         fMaxSlope;      // Maximum physical track slope (0=perp)
    Double_t         fWidth;         // Width of tracking region (m)
    
    //TrackSlopeSelection
    Bool_t           fDoElTrackSel;// turn on or off track preselection based on their slope
    std::unordered_map<Int_t, Int_t>         fmoduleOrder;   // 
    Double_t         fCorrSlope;
    Double_t         fCorrSlopeMin;
    Double_t         fCorrSlopeMax;
    Double_t         fCorrInterceptHigh;
    Double_t         fCorrInterceptLow;

    TVector2         fAxis;          // Projection axis, normal to strips
    THaDetectorBase* fDetector;      //! Parent detector
    PatternTree*     fPatternTree;   // Precomputed template database

    UInt_t           fDummyPlanePattern; // Bitpattern of dummy plane numbers
    UInt_t           fFirstPlaneNum; // Idx of first active plane in fAllPlanes
    UInt_t           fLastPlaneNum;  // Idx of last active plane in fAllPlanes

    // Plane occupancy parameters
    UInt_t           fMinFitPlanes;  // Min num of planes required for fitting
    UInt_t           fMaxMiss;       // Allowed number of missing planes
    Bool_t           fRequire1of2;   // Require hit in at least one plane of a pair
    TBits*           fPlaneCombos;   // Allowed plane occupancy patterns
    TBits*           fAltPlaneCombos;// Allowed plane patterns including dummies

    // Road construction control
    UInt_t           fMaxPat;        // Sanity cut on number of patterns
    UInt_t           fFrontMaxBinDist; // Max pattern dist in front plane
    UInt_t           fBackMaxBinDist;  // Max pattern dist in back plane
    UInt_t           fHitMaxDist;    // Max allowed distance between hits for
                                     // clustering patterns into roads

    // Fit statistics cut parameters
    Double_t         fConfLevel;     // Requested confidence level for chi2 cut
    vec_pdbl_t       fChisqLimits;   // lo/hi onfidence interval limits on Chi2

    // Event-by-event results
    Hitpattern*      fHitpattern;    // Hitpattern of current event
    NodeVec_t        fPatternsFound; // Patterns found by TreeSearch
    TClonesArray*    fRoads;         // Roads found by MakeRoads
    UInt_t           fNgoodRoads;    // Good roads in fRoads
    UInt_t           fNMCRoads;      // Good roads with all their hits from MC in fRoads
    TClonesArray*    fRoadCorners;   // Road corners, for event display
    ETrackingStatus  fTrkStat;       // Reconstruction status

    // Statistics (only needed for TESTCODE, but kept for binary compatibility)
    UInt_t n_hits, n_bins, n_binhits, maxhits_bin;
    UInt_t n_test, n_pat, n_roads, n_dupl, n_badfits;
    Double_t t_treesearch, t_roads, t_fit, t_track;

    Bool_t  FitRoads();
    Bool_t  RemoveDuplicateRoads();
    void    SetAngle( Double_t a );
    UInt_t  GetNallPlanes() const { return (UInt_t)fAllPlanes.size(); }

    virtual Hitpattern* MakeHitpattern( const PatternTree& ) const;
    virtual void MakePlaneCombos( const vpl_t& planes, TBits*& combos ) const;

    // Podd interface
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );
    virtual const char* GetDBFileName() const;
    virtual void MakePrefix();

    // NodeVisitor class for comparing patterns in the tree with the
    // hitpattern. Matches represent candidates for track roads and are
    // added to the list of roads for further analysis
    class ComparePattern : public NodeVisitor {
    public:
      ComparePattern( const Hitpattern* hitpat, const TBits* combos,
		      NodeVec_t* matches, UInt_t dummypattern = 0 )
	: fHitpattern(hitpat), fPlaneCombos(combos), fMatches(matches),
	  fDummyPlanePattern(dummypattern)
#ifdef TESTCODE
	, fNtest(0)
#endif
      { assert(fHitpattern && fPlaneCombos && fMatches); }
      virtual ETreeOp operator() ( const NodeDescriptor& nd );
#ifdef TESTCODE
      UInt_t GetNtest() const { return fNtest; }
#endif
    private:
      const Hitpattern* fHitpattern;   // Hitpattern to compare to
      const TBits*      fPlaneCombos;  // Allowed plane occupancy patterns
      NodeVec_t*        fMatches;      // Set of matching patterns
      UInt_t            fDummyPlanePattern;  // Dummy plane # bitpattern
#ifdef TESTCODE
      UInt_t fNtest;  // Number of pattern comparisons
#endif
    };

  private:
    // Prevent default copying, assignment
    Projection( const Projection& orig );
    const Projection& operator=( const Projection& rhs );

    ClassDef(Projection,0)  // A track projection plane
  };

  //___________________________________________________________________________
  inline
  Double_t Projection::GetAngle() const
  {
    // Return axis angle in rad, normalized to [-pi,pi]

    return TMath::ATan2( fAxis.Y(), fAxis.X() );
  }

  //___________________________________________________________________________
  inline const pdbl_t& Projection::GetChisqLimits( UInt_t i ) const
  {
    assert( i < fChisqLimits.size() );
    return fChisqLimits[i];
  }

  //___________________________________________________________________________
  inline
  UInt_t Projection::GetNpatterns() const
  {
    return static_cast<UInt_t>( fPatternsFound.size() );
  }

  //___________________________________________________________________________
  inline
  UInt_t Projection::GetNroads() const
  {
    // Return total number of roads found (including voided ones)
    assert( fRoads && fRoads->GetLast()+1 >= 0  );
    return static_cast<UInt_t>( fRoads->GetLast()+1 );
  }

  //___________________________________________________________________________
  inline
  Road* Projection::GetRoad( UInt_t i ) const
  {
    // Get i-th road (i=0..GetNroads()-1) found in this projection.
    // The pointer should always be non-zero, but check with Road::IsGood()
    // whether the road is good

    assert( fRoads && i < GetNroads() );
    return (Road*)fRoads->UncheckedAt(i); //static_cast would require Road.h
  }

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
