#ifndef TreeSearch__TreeWalk
#define TreeSearch__TreeWalk

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::TreeWalk                                                      //
//                                                                           //
// Generic function object to traverse a PatternTree.                        //
//                                                                           //
// TreeWalk implements an internal iterator pattern [E. Gamma et al.,        //
// "Design Patterns", Addison-Wesley, 1995] that applies generic operation   //
// to each tree element.  The operation itself represents a simplified form  //
// of a visitor pattern [ibid.], where the simplification lies in the fact   //
// that the tree contains only a single class of nodes instead of many.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Node.h"
#include <fstream>
#include <iostream>
#include <map>

namespace TreeSearch {

  //___________________________________________________________________________
  // Base class for "Visitors" to the pattern tree nodes
  class NodeVisitor {
  public:
    enum ETreeOp { kRecurse, kRecurseUncond, kSkipChildNodes, kError };

    virtual ETreeOp operator() ( const NodeDescriptor& nd ) = 0;
    virtual ~NodeVisitor() {}
  protected:
    // Functions to allow derived classes to modify Links and Patterns.
    // (needed since C++ friend access doesn't extend to derived classes)
    static void SetLinkPattern( Link* link, Pattern* pattern );
    static void SetLinkType( Link* link, Int_t type );
    static void SetLinkNext( Link* link, Link* next );
    static void SetPatternChild( Pattern* pat, Link* link );
  };


  //___________________________________________________________________________
  // The actual tree iterator class
  class TreeWalk {
  private:
    UInt_t   fNlevels;  // Number of levels in tree
  public:
    TreeWalk( UInt_t nlevels = 0 ) : fNlevels(nlevels) {}
    virtual ~TreeWalk() {}
    void SetNlevels( UInt_t n ) { fNlevels = n; }
    
    NodeVisitor::ETreeOp 
    operator() ( Link* link, NodeVisitor& op, Pattern* parent = 0,
		 UInt_t depth = 0, UInt_t shift = 0,
		 Bool_t mirrored = false ) const;
    ClassDef(TreeWalk, 0)  // Generic traversal function for a PatternTree
  };


  //___________________________________________________________________________
  // TreeSearch::WritePattern
  // TreeSearch::CountPattern
  // TreeSearch::PrintPattern
  //
  // Common "visitor" objects for tree traversal operations. The object's
  // operator() is applied to each tree node. These can be used with the
  // TreeWalk iterator

  //___________________________________________________________________________
  // Write patterns to binary file
  class WritePattern : public NodeVisitor {
  public:
    WritePattern( const char* filename, size_t index_size = sizeof(Int_t) );
    virtual ~WritePattern() { delete os; }
    virtual ETreeOp operator() ( const NodeDescriptor& nd );

  private:
    std::ofstream* os;        // Output file stream
    size_t    fIdxSiz;        // Byte size of the pattern count (1, 2, or 4)
    std::map<Pattern*,Int_t> fMap; // Index map for serializing 

    // Because of the pointer member, disallow copying and assignment
    WritePattern( const WritePattern& orig );
    WritePattern& operator=( const WritePattern& rhs );
  };

  //___________________________________________________________________________
  // Count unique patterns (including shifts)
  class CountPattern : public NodeVisitor {
  public:
    CountPattern() : fCount(0) {}
    virtual ETreeOp operator() ( const NodeDescriptor& )
    { fCount++; return kRecurse; }
    ULong64_t GetCount() const { return fCount; }

  private:
    ULong64_t    fCount;    // Pattern count
  };

  //___________________________________________________________________________
  // Pretty print to output stream (and count) all actual patterns
  class PrintPattern : public NodeVisitor {
  public:
    PrintPattern( std::ostream& ostr = std::cout, bool dump = false ) 
      : os(ostr), fCount(0), fDump(dump) {}
    virtual ETreeOp operator() ( const NodeDescriptor& nd );
    ULong64_t GetCount() const { return fCount; }

  private:
    std::ostream&  os;        // Destinaton stream
    ULong64_t      fCount;    // Pattern count
    Bool_t         fDump;     // Dump mode (one line per pattern)
  };

  /////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
