#ifndef ROOT_TreeSearch_PatternTree
#define ROOT_TreeSearch_PatternTree

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>
#include <cassert>
#include <cstring>

using std::vector;

namespace TreeSearch {

  extern const Double_t kBig;

  class ListNode;

  class Pattern {
    friend class PatternTree;
  private:
    UShort_t*  fBits;        // [fNbits] Bit pattern array
    ListNode*  fChild;       // Linked list of child patterns
    UChar_t    fNbits;       // Bit count
    UChar_t    fMinDepth;    // Minimum valid depth for this pattern
    UChar_t    fMaxDepth;    // Maximum valid depth (used during building)
    Int_t      fRefIndex;    // Reference index for serializing the tree
  public:
    Pattern( UInt_t nplanes );
    Pattern( const Pattern& orig );
    Pattern& operator=( const Pattern& rhs );
    //    Pattern( UShort_t* bitloc );
    ~Pattern();
    bool operator==( const Pattern& rhs ) const {
      assert( fNbits == rhs.fNbits );
      return (0 == memcmp( fBits, rhs.fBits, fNbits*sizeof(UShort_t)));
    }
    bool operator!=( const Pattern& rhs ) const { return !(*this == rhs ); }
    ListNode*  AddChild( Pattern& child, Int_t type );
    Pattern*   FindChild( const Pattern& pat, Int_t type );
    UShort_t*  GetBits()            { return fBits; }
    UInt_t     GetMinDepth() const  { return fMinDepth; }
    UInt_t     GetWidth()    const  { return fBits[fNbits-1]-fBits[0]; }
    UInt_t     GetNbits()    const  { return fNbits; }
    Int_t      Hash() const { return fNbits>0 ? fBits[fNbits-1]-fBits[1] : 0; }
    UShort_t&  operator[](UInt_t i) { assert(i<fNbits); return fBits[i]; }
    UShort_t   operator[](UInt_t i) const { assert(i<fNbits); return fBits[i];}
    void       UsedAtDepth( UInt_t depth );
  };

  class ListNode {
    friend class PatternTree;
    friend class Pattern;
  private:
    Pattern*   fPattern;    // Bit pattern treenode
    ListNode*  fNext;       // Next list element
    Int_t      fOp;         // Operation to be applied to pattern
  public:
    ListNode( Pattern* pat, ListNode* next, Int_t op )
      : fPattern(pat), fNext(next), fOp(op) {}
    Pattern*   GetPattern() const { return fPattern; }
    ListNode*  Next() const { return fNext; }
  };

  class PatternTree {
  public:
    PatternTree( UInt_t depth, const vector<double>& zpos );
    virtual ~PatternTree();

    void AddHash( Pattern* pat );
    void DeleteBuildTree();
    Pattern* Find( const Pattern& pat );
    bool TestSlope( const Pattern& pat, UInt_t depth );
    bool LineCheck( const Pattern& pat );
    Int_t Init();

    void MakeChildNodes( Pattern& parent, UInt_t depth );

    Double_t  GetMaxSlope() const { return fMaxSlope; }

    void      SetMaxSlope( Double_t slope ) { fMaxSlope = slope; }


  private:
    UInt_t         fDepth;       // Depth of tree (levels = 0 - depth-1)
    UInt_t         fNplanes;     // Number of hitpattern planes
    Double_t       fMaxSlope;    // Max allowed slope, normalized units (0-1)
    vector<double> fZ;           // z positions of planes, normalized (0-1)

    ListNode*      fRoot;        // Root node 

    vector<ListNode*> fHashTable;// Hashtable of patterns (for Build)
    UInt_t         fNpatterns;   // Pattern count

    ClassDef(PatternTree,0)   // Precomputed template database
  };

  //___________________________________________________________________________
  inline
  void Pattern::UsedAtDepth( UInt_t depth )
  {
    // Mark pattern as having been used at given depth

    if( depth < fMinDepth ) fMinDepth = depth;
    if( depth > fMaxDepth ) fMaxDepth = depth;
  }


///////////////////////////////////////////////////////////////////////////////

}  // end namespace TreeSearch

#endif
