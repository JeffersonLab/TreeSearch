#ifndef ROOT_TreeSearch_PatternTree
#define ROOT_TreeSearch_PatternTree

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>
#include <string>

namespace TreeSearch {

  class Pattern;

  struct TreeParam_t {
    UInt_t    maxdepth;
    Double_t  width;
    Double_t  maxslope;
    std::vector<Double_t > zpos;
  };

  class PatternTree {
  public:
    PatternTree( const TreeParam_t& param );
    PatternTree( const TreeParam_t& param, UInt_t nPatterns, UInt_t nLinks );
    virtual ~PatternTree();
 
    Int_t  AddPattern(int);
    Int_t  Read( const std::string& filename );

    const  TreeParam_t& GetParameters() const { return fParameters; }
    const  Pattern*     GetRoot()       const { return fRoot; }
    static Int_t        Normalize( TreeParam_t& param );

  private:
    
    TreeParam_t  fParameters; // Tree parameters (levels, width, depth)
    Pattern*     fRoot;       // Root node 

    ClassDef(PatternTree,0)   // Precomputed template database
  };


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
