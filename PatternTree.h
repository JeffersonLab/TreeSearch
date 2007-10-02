#ifndef ROOT_TreeSearch_PatternTree
#define ROOT_TreeSearch_PatternTree

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>

namespace TreeSearch {

  //  class Pattern;
  class Link;

  struct TreeParam_t {
    UInt_t    maxdepth;
    Double_t  width;
    Double_t  maxslope;
    std::vector<Double_t > zpos;
  };

  class PatternTree {
  public:
    PatternTree( const TreeParam_t& param,
		 UInt_t nPatterns = 0, UInt_t nLinks = 0 );
    virtual ~PatternTree();
 
    Int_t  AddPattern( const Link* patt );
    Int_t  Read( const char* filename );

    const  TreeParam_t& GetParameters() const { return fParameters; }
    const  Link*        GetRoot()       const { return fRoot; }
    static Int_t        Normalize( TreeParam_t& param );

  private:
    
    TreeParam_t  fParameters; // Tree parameters (levels, width, depth)
    Link*        fRoot;       // Root node 

    ClassDef(PatternTree,0)   // Precomputed template database
  };

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
