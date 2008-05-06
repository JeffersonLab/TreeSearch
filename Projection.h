#ifndef ROOT_TreeSearch_Projection
#define ROOT_TreeSearch_Projection

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaAnalysisObject.h"
#include "TreeWalk.h"  // for NodeVisitor, Node_t
#include "Types.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TVector2.h"
#include <vector>
#include <set>
#include <cassert>
#include <functional>
#include <list>

class THaDetectorBase;
class TBits;

namespace TreeSearch {

  class Hitpattern;
  class PatternTree;
  class WirePlane;
  class Road;

  class Projection : public THaAnalysisObject {
  public:

    Projection( EProjType type, const char* name, Double_t angle,
		THaDetectorBase* parent );
    Projection() : fDetector(0), fPatternTree(0), fPlaneCombos(0),
		   fHitpattern(0), fRoads(0), fRoadCorners(0) {} // ROOT RTTI
    virtual ~Projection();

    void            AddPlane( WirePlane* wp, WirePlane* partner = 0 );
    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );
    EStatus         InitLevel2( const TDatime& date );
    virtual void    Print( Option_t* opt="" ) const;
    void            Reset();

    Int_t           FillHitpattern();
    Int_t           Track();
    Int_t           MakeRoads();


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
    WirePlane*      GetPlane ( UInt_t plane ) const { return fPlanes[plane]; }
    Double_t        GetPlaneZ( UInt_t plane ) const;
    Road*           GetRoad  ( UInt_t i )     const;
    Double_t        GetSinAngle()     const { return fAxis.Y(); }
    EProjType       GetType()         const { return fType; }
    Double_t        GetWidth()        const { return fWidth; }
    Double_t        GetZsize()        const;

    void            SetMaxSlope( Double_t m ) { fMaxSlope = m; }
    void            SetPatternTree( PatternTree* pt ) { fPatternTree = pt; }
    void            SetWidth( Double_t width ) { fWidth = width; }
    
#ifdef TESTCODE
    std::vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }
#endif

    // Analysis control flags
    enum {
      kEventDisplay = BIT(14) // Support event display
    };

  protected:
    typedef std::vector<Node_t*> NodeVec_t;

    // Configuration
    EProjType        fType;          // Type of plane (u,v,x,y...)
    std::vector<WirePlane*> fPlanes; // Wire planes in this projection
    UInt_t           fNlevels;       // Number of levels of search tree
    Double_t         fMaxSlope;      // Maximum physical track slope (0=perp)
    Double_t         fWidth;         // Width of tracking region (m)
    TVector2         fAxis;          // Nominal projection axis normal to wires
    THaDetectorBase* fDetector;      //! Parent detector
    PatternTree*     fPatternTree;   // Precomputed template database

    // Plane occupancy parameters
    UInt_t           fMinFitPlanes;  // Min num of planes required for fitting
    UInt_t           fMaxMiss;       // Allowed number of missing planes
    Bool_t           fRequire1of2;   // Require one plane of a plane pair
    TBits*           fPlaneCombos;   // Allowed plane occupancy patterns

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
    TClonesArray*    fRoadCorners;   // Road corners, for event display

#ifdef TESTCODE
    // Statistics
    UInt_t n_hits, n_bins, n_binhits, maxhits_bin;
    UInt_t n_test, n_pat, n_roads, n_dupl, n_badfits;
    Double_t t_treesearch, t_roads, t_fit, t_track;
#endif

    Bool_t  FitRoads();
    Bool_t  RemoveDuplicateRoads();
    void    SetAngle( Double_t a );

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
		      NodeVec_t* matches )
	: fHitpattern(hitpat), fPlaneCombos(combos), fMatches(matches)
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
    // Return wire angle in rad, normalized to [-pi,pi]
    
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
