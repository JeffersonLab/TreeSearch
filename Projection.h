#ifndef ROOT_TreeSearch_Projection
#define ROOT_TreeSearch_Projection

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaAnalysisObject.h"
#include "TreeWalk.h"
#include "Hit.h"
#include "TMath.h"
#include "TClonesArray.h"
#include <vector>
#include <map>
#include <cassert>
#include <utility>

class THaDetectorBase;
class TBits;

namespace TreeSearch {

  class Hitpattern;
  class PatternTree;
  class WirePlane;
  class Road;

  class Projection : public THaAnalysisObject {
  public:

    typedef std::pair<Double_t,Double_t> pdbl_t;

    Projection( Int_t type, const char* name, Double_t angle,
		THaDetectorBase* parent = 0 );
    virtual ~Projection();

    void            AddPlane( WirePlane* wp, WirePlane* partner = 0 );
    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    EStatus         InitLevel2( const TDatime& date );
    virtual void    Print( Option_t* opt="" ) const;
    void            Reset();

    Int_t           FillHitpattern();
    Int_t           Track();
    Int_t           MakeRoads();


    Double_t        GetAngle() const;
    const pdbl_t&   GetChisqLimits( UInt_t i ) const;
    Double_t        GetConfLevel() const { return fConfLevel; }
    Double_t        GetCosAngle() const { return fCosAngle; }
    UInt_t          GetHitMaxDist() const { return fHitMaxDist; }
    Hitpattern*     GetHitpattern() const { return fHitpattern; }
    TBits*          GetLayerCombos() const { return fLayerCombos; }
    WirePlane*      GetLayer( UInt_t layer )  const { return fLayers[layer]; }
    Double_t        GetLayerZ( UInt_t layer ) const;
    Double_t        GetMaxSlope() const { return fMaxSlope; }
    UInt_t          GetMinFitPlanes() const { return fMinFitPlanes; }
    UInt_t          GetNlevels()  const { return fNlevels; }
    UInt_t          GetNlayers()  const { return (UInt_t)fLayers.size(); }
    UInt_t          GetNpatterns() const;
    UInt_t          GetNplanes()  const { return (UInt_t)fPlanes.size(); }
    Int_t           GetNroads()   const { return fRoads->GetLast()+1; }
    UInt_t          GetBinMaxDist() const { return fBinMaxDist; }
    TBits*          GetPlaneCombos() const { return fPlaneCombos; }
    WirePlane*      GetPlane( UInt_t plane )  const { return fPlanes[plane]; }
    Double_t        GetPlaneZ( UInt_t plane ) const;
    Double_t        GetSinAngle() const { return fSinAngle; }
    Int_t           GetType()  const { return fType; }
    Double_t        GetWidth() const { return fWidth; }
    Double_t        GetZsize() const;

    void            SetMaxSlope( Double_t m ) { fMaxSlope = m; }
    void            SetPatternTree( PatternTree* pt ) { fPatternTree = pt; }
    void            SetWidth( Double_t width ) { fWidth = width; }
    
    //FIXME: for testing
//  std::vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }
//  std::vector<TreeSearch::WirePlane*>& GetListOfLayers() { return fLayers; }

  protected:
    typedef std::vector<pdbl_t> vec_pdbl_t;

    Int_t           fType;        // Type of plane (u,v,x,y...)
    std::vector<WirePlane*> fPlanes; // Wire planes in this projection
    std::vector<WirePlane*> fLayers; // Effective detector planes (wp pairs)
    UInt_t          fNlevels;     // Number of levels of search tree
    Double_t        fMaxSlope;    // Maximum physical track slope (0=perp)
    Double_t        fWidth;       // Width of tracking region (m)
    Double_t        fSinAngle;    // Sine of wire angle
    Double_t        fCosAngle;    // Cosine of wire angle
    UInt_t          fHitMaxDist;  // Max allowed distance between hits for
                                  // clustering patterns into roads
    UInt_t          fBinMaxDist;  // Search distance for MakeRoads

    Hitpattern*     fHitpattern;  // Hitpattern of current event
    PatternTree*    fPatternTree; // Precomputed template database

    THaDetectorBase* fDetector;    //! Parent detector

    typedef std::map<const NodeDescriptor,HitSet> NodeMap_t;
    NodeMap_t        fPatternsFound; // Patterns found by TreeSearch
    TClonesArray*    fRoads;       // Roads found by MakeRoads

    TBits*           fPlaneCombos; // Allowed plane combos with missing hits
    TBits*           fLayerCombos; // Allowed layer combos with missing hits

    UInt_t           fMinFitPlanes;// Min number of planes required for fitting
    Double_t         fConfLevel;   // Requested confidence level for chi2 cut
    vec_pdbl_t       fChisqLimits; // lo/hi onfidence interval limits on Chi2

#ifdef TESTCODE
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
		      std::map<const NodeDescriptor,HitSet>* matches )
	: fHitpattern(hitpat), fLayerCombos(combos), fMatches(matches)
#ifdef TESTCODE
	, fNtest(0)
#endif
      { assert(fHitpattern && fLayerCombos && fMatches); }
      virtual ETreeOp operator() ( const NodeDescriptor& nd );
#ifdef TESTCODE
      UInt_t GetNtest() const { return fNtest; }
#endif
    private:
      const Hitpattern* fHitpattern;   // Hitpattern to compare to
      const TBits*      fLayerCombos;  // Allowed combos of missing layers
      // Set of matching patterns
      std::map<const NodeDescriptor,HitSet>* fMatches;
#ifdef TESTCODE
      UInt_t fNtest;  // Number of pattern comparisons
#endif
    };

  private:
    // Prevent default construction, copying, assignment
    Projection() {};
    Projection( const Projection& orig );
    const Projection& operator=( const Projection& rhs );

    ClassDef(Projection,0)  // A track projection plane
  };

  //___________________________________________________________________________
  inline
  Double_t Projection::GetAngle() const
  {
    // Return wire angle in rad, normalized to [-pi,pi]
    Double_t a = TMath::ASin(fSinAngle);
    if( fCosAngle < 0.0 )
      return (fSinAngle > 0.0) ? TMath::TwoPi() - a : -TMath::TwoPi() - a;
  
    return a;
  }

  //___________________________________________________________________________
  inline const Projection::pdbl_t& Projection::GetChisqLimits( UInt_t i ) const
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

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
