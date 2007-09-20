#ifndef ROOT_TreeSearch_PatternGenerator
#define ROOT_TreeSearch_PatternGenerator

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternGenerator                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Pattern.h"
#include <vector>

using std::vector;

namespace TreeSearch {

  class PatternTree;
  class ListNode;

  class PatternGenerator {
  public:
    PatternGenerator();
    virtual ~PatternGenerator();

    PatternTree* Generate( UInt_t nlevels, Double_t detector_width, 
			   const vector<double>& zpos, Double_t maxslope );

    void Print( Option_t* opt="" );

  private:
    UInt_t         fNlevels;     // Number of levels of the tree (0-nlevels-1)
    UInt_t         fNplanes;     // Number of hitpattern planes
    Double_t       fMaxSlope;    // Max allowed slope, normalized units (0-1)
    vector<double> fZ;           // z positions of planes, normalized (0-1)

    vector<ListNode*> fHashTable;// Hashtab for indexing patterns during build

    void     AddHash( Pattern* pat );
    Pattern* Find( const Pattern& pat );
    bool     TestSlope( const Pattern& pat, UInt_t depth );
    bool     LineCheck( const Pattern& pat );
    void     MakeChildNodes( Pattern* parent, UInt_t depth );

    enum EOperation { kDelete, kResetRefIndex, kCountPatterns,
		      kCountChildNodes, kMaxChildlistLength, 
		      kMaxHashDepth, kBytesRequired };
    Int_t    DoTree( EOperation op );

    // Utility class for iterating over child patterns
    class ChildIter {
    private:
      const Pattern fParent;
      Pattern   fChild;
      Int_t     fCount;
      Int_t     fType;
    public:
      ChildIter( const Pattern& parent ) 
	: fParent(parent), fChild(parent), fType(0) { reset(); }
      ChildIter&      operator++();
      const ChildIter operator++(int) { 
	ChildIter clone(*this);
	++(*this); 
	return clone;
      }
      Pattern& operator*()            { return fChild; }
      operator bool()  const { return (fCount >= 0); }
      bool     operator!()      const { return !((bool)*this); }
      Int_t    type()           const { return fType; }
      void     reset() { 
	fCount = 1<<fParent.GetNbits();
	++(*this);
      }
    };

    ClassDef(PatternGenerator,0)   // Generator for pattern template database

  }; // end class PatternGenerator

///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
