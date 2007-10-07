#ifndef TreeSearch__NodeVisitor
#define TreeSearch__NodeVisitor

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::NodeVisitor                                                   //
//                                                                           //
// Common "visitor" objects for tree traversal operations. The object's      //
// operator() is applied to each tree node. These should be used with the    //
// TreeWalk iterator (see TreeWalk.h)                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TreeWalk.h"
#include <fstream>
#include <iomanip>
#include <map>

namespace TreeSearch {

//___________________________________________________________________________
// Base class for "Visitor" classes to the tree nodes.
//   class NodeVisitor {
//     virtual ETreeOp operator() ( const NodeDescriptor& nd ) = 0;
//   };

// Write patterns to binary file
class WritePattern {
public:
  WritePattern( const char* filename, size_t index_size = sizeof(Int_t) );
  ~WritePattern() { delete os; }
  TreeWalk::ETreeOp operator() ( const NodeDescriptor& nd );

private:
  std::ofstream* os;        // Output file stream
  size_t    fIdxSiz;        // Byte size of the pattern count (1, 2, or 4)
  std::map<Pattern*,Int_t> fMap; // Index map for serializing 

  // Because of the pointer member, disallow copying and assignment
  WritePattern( const WritePattern& orig );
  WritePattern& operator=( const WritePattern& rhs );
};

// Count unique patterns (including shifts)
class CountPattern {
public:
  CountPattern() : fCount(0) {}
  TreeWalk::ETreeOp operator() ( const NodeDescriptor& nd )
  { fCount++; return TreeWalk::kRecurse; }
  ULong64_t GetCount() const { return fCount; }

private:
  ULong64_t    fCount;    // Pattern count
};

// Pretty print to output stream (and count) all actual patterns
class PrintPattern {
public:
  PrintPattern( std::ostream& ostr = std::cout, bool dump = false ) 
    : os(ostr), fCount(0), fDump(dump) {}
  TreeWalk::ETreeOp operator() ( const NodeDescriptor& nd );
  ULong64_t GetCount() const { return fCount; }

private:
  std::ostream&  os;        // Destinaton stream
  ULong64_t      fCount;    // Pattern count
  Bool_t         fDump;     // Dump mode (one line per pattern)
};

/////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
