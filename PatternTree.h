#ifndef ROOT_TreeSearch_PatternTree
#define ROOT_TreeSearch_PatternTree

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Pattern.h"
#include "TreeWalk.h"
#include <vector>
#include <map>
#include <iostream>
#include <cassert>

using std::vector;

namespace TreeSearch {

  class TreeParam_t {
  public:
    TreeParam_t( UInt_t maxdepth, Double_t width, Double_t maxslope,
		 const vector<Double_t>& zpos ) 
      : fMaxdepth(maxdepth), fNormalized(false), fWidth(width),
	fMaxslope(maxslope), fZpos(zpos) {}
    Int_t Normalize();
    UInt_t   maxdepth() const { return fMaxdepth; }
    Double_t width()    const { return fWidth; }
    Double_t maxslope() const { return fMaxslope; }
    const vector<Double_t>& zpos() const { return fZpos; }
  private: 
    UInt_t    fMaxdepth;    // Depth of tree
    Bool_t    fNormalized;  // maxslope and zpos are normalized
    Double_t  fWidth;       // Physical detector width (needed in Hitpattern)
    Double_t  fMaxslope;    // Max slope (dx/dz) of tracks
    vector<Double_t> fZpos; // z-positions of pattern planes
  };

  class PatternTree {
  public:
    PatternTree( const TreeParam_t& param,
		 UInt_t nPatterns = 0, UInt_t nLinks = 0 );
    virtual ~PatternTree();
    // TODO: copy c'tor, assignment (see below)
    // TODO: implement Read()

    static PatternTree* Read( const char* filename, const TreeParam_t& param );

    void   Print( Option_t* opt="", std::ostream& os = std::cout );
    Int_t  Write( const char* filename );

    Bool_t IsOK()       const { return fParamOK; }
    UInt_t GetNlevels() const { return fParameters.maxdepth()+1; }
    UInt_t GetNplanes() const { return fParameters.zpos().size(); }
    const TreeParam_t& GetParameters() const { return fParameters; }
    Link*  GetRoot() { return fLinks.empty() ? 0 : &fLinks.front(); }
    Double_t GetWidth() const { return fParameters.width(); }

    // Copy an arbitrary tree into the PatternTree array structures
    class CopyPattern : public NodeVisitor {
    public:
      CopyPattern( PatternTree* tree ) : fTree(tree) { assert(fTree); }
      virtual ETreeOp operator() ( const NodeDescriptor& nd );

    private:
      PatternTree* fTree;    // Tree object to fill
      std::map<Pattern*,Int_t> fMap;   // Index map for serializing 
      void AddChild( Pattern* parent, Pattern* child, Int_t type );
    };
    friend class CopyPattern;

  private:
    
    typedef vector<Pattern>::size_type vpsz_t;
    typedef vector<Link>::size_type vlsz_t;
    typedef vector<UShort_t>::size_type vsiz_t;

    TreeParam_t      fParameters; // Tree parameters (levels, width, depth)
    Bool_t           fParamOK;    // Flag: Parameters are tested valid

    vector<Pattern>  fPatterns;   // Array of all patterns
    vector<Link>     fLinks;      // Array of all links
    vector<UShort_t> fBits;       // Array of all pattern bits
 
    // Variables for unserializing the tree

    vpsz_t           fNpat;       // Current pattern count 
    vlsz_t           fNlnk;       // Current link count
    vsiz_t           fNbit;       // Current bit count

    // Disallow copying and assignment for now. The vectors can NOT be copied 
    // directly since they contain pointers to the other vectors' elements!
    PatternTree( const PatternTree& orig );
    const PatternTree& operator=( const PatternTree& rhs );

    ClassDef(PatternTree,0)   // Precomputed template database
  };

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
