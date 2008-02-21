#ifndef ROOT_TreeSearch_PatternGenerator
#define ROOT_TreeSearch_PatternGenerator

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternGenerator                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Pattern.h"
#include "PatternTree.h"
#include <vector>

using std::vector;

namespace TreeSearch {

  class PatternGenerator {
    //    friend class Test_PatternGenerator;
  public:
    PatternGenerator();
    virtual ~PatternGenerator();

    PatternTree* Generate( TreeParam_t parameters );
    PatternTree* Generate( UInt_t maxdepth, Double_t detector_width, 
			   const char* zpos, Double_t maxslope );

    struct Statistics_t {
      UInt_t nPatterns, nLinks, nBytes, MaxChildListLength, nHashBytes;
      ULong64_t nAllPatterns;
      Double_t  BuildTime;
    };

    Pattern* GetRoot() const { return fHashTable[0].fPattern; }
    const Statistics_t& GetStatistics() const { return fStats; }

    void  Print( Option_t* opt="", std::ostream& os = std::cout ) const;

  private:

    class HashNode {
      friend class PatternGenerator;
    private:
      Pattern* fPattern;    // Bit pattern treenode
      UInt_t   fMinDepth;   // Minimum valid depth for this pattern (<=16)
      void     UsedAtDepth( UInt_t depth ) {
	if( depth < fMinDepth ) fMinDepth = depth;
      }
    public:
      HashNode( Pattern* pat = 0 ) : fPattern(pat), fMinDepth(kMaxUInt) {}
      Pattern* GetPattern() const { return fPattern; }
    };

    UInt_t         fNlevels;     // Number of levels of the tree (0-nlevels-1)
    UInt_t         fNplanes;     // Number of hitpattern planes
    Double_t       fMaxSlope;    // Max allowed slope, normalized units (0-1)
    vector<double> fZ;           // z positions of planes, normalized (0-1)

    vector<HashNode> fHashTable; // Hashtab for indexing patterns during build
    Statistics_t   fStats;       // Tree statistics

    HashNode* AddHash( Pattern* pat );
    void      CalcStatistics();
    void      ClearStatistics();
    void      DeleteTree();
    HashNode* Find( const Pattern& pat );
    UInt_t    Hash( const Pattern& pat ) const;
    bool      LineTest( const Pattern& pat ) const;
    void      MakeChildNodes( HashNode* parent, UInt_t depth );
    bool      SlopeTest( const Pattern& pat, UInt_t depth ) const;

    ClassDef(PatternGenerator,0)   // Generator for pattern template database

  }; // end class PatternGenerator

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
