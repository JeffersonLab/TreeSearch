#ifndef ROOT_TreeSearch_PatternGenerator
#define ROOT_TreeSearch_PatternGenerator

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternGenerator                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Pattern.h"
#include "PatternTree.h"
#include "TreeWalk.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>

using std::vector;

namespace TreeSearch {

  class PatternGenerator {

  public:
    PatternGenerator();
    virtual ~PatternGenerator();

    PatternTree* Generate( TreeParam_t parameters );
    PatternTree* Generate( UInt_t maxdepth, Double_t detector_width, 
			   const char* zpos, Double_t maxslope );

    void  Print( Option_t* opt="", std::ostream& os = std::cout ) const;
    Int_t Write( const char* filename );

  private:

    struct Statistics_t {
      UInt_t nPatterns, nLinks, nBytes, MaxChildListLength, MaxHashDepth,
	nHashBytes;
      ULong64_t nAllPatterns;
      Double_t  BuildTime;
    };

    UInt_t         fNlevels;     // Number of levels of the tree (0-nlevels-1)
    UInt_t         fNplanes;     // Number of hitpattern planes
    Double_t       fMaxSlope;    // Max allowed slope, normalized units (0-1)
    vector<double> fZ;           // z positions of planes, normalized (0-1)

    vector<Link*>  fHashTable;   // Hashtab for indexing patterns during build
    Statistics_t   fStats;       // Tree statistics
    TreeWalk       fTreeWalk;    // TreeWalk functor

    enum EOperation { kDelete, kResetRefIndex };
    void     AddHash( Pattern* pat );
    void     CalcStatistics();
    void     DoTree( EOperation op );
    Pattern* Find( const Pattern& pat );
    bool     LineCheck( const Pattern& pat );
    void     MakeChildNodes( Pattern* parent, UInt_t depth );
    bool     TestSlope( const Pattern& pat, UInt_t depth );

    ClassDef(PatternGenerator,0)   // Generator for pattern template database

  }; // end class PatternGenerator

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
