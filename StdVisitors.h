#ifndef TreeSearch__StdVisitors
#define TreeSearch__StdVisitors

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::WritePattern                                                  //
// TreeSearch::CountPattern                                                  //
// TreeSearch::PrintPattern                                                  //
//                                                                           //
// Common "visitor" objects for tree traversal operations. The object's      //
// operator() is applied to each tree node. These can be used with the       //
// TreeWalk iterator (see TreeWalk.h)                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "NodeVisitor.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>

namespace TreeSearch {

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
  virtual ETreeOp operator() ( const NodeDescriptor& nd )
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
