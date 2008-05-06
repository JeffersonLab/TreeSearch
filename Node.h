#ifndef TreeSearch__Node
#define TreeSearch__Node

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Node                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Pattern.h"
#include <utility>

namespace TreeSearch {

  //___________________________________________________________________________
  // Structure to describe a tree node (an actual pattern, including shifts
  // and mirroring) 
  struct NodeDescriptor {
    Link*    link;     // Linked-list node pointing to a base pattern
    Pattern* parent;   // Parent node
    UShort_t shift;    // Shift of the base pattern to its actual position
    Bool_t   mirrored; // Pattern is mirrored
    UChar_t  depth;    // Current recursion depth

    NodeDescriptor( Link* ln, Pattern* p, UShort_t shft, Bool_t mir, 
		    UChar_t dep )
      : link(ln), parent(p), shift(shft), mirrored(mir), depth(dep)
    { assert(ln && ln->GetPattern()); }
    NodeDescriptor() {}
    ~NodeDescriptor() {}

    UShort_t Start() const { return shift; }
    UShort_t End()   const { return (*this)[link->GetPattern()->GetNbits()-1];}
    void     Print() const;

    // operator[] returns actual bit value in the i-th plane
    UShort_t  operator[](UInt_t i) const {
      if( i == 0 ) return shift;
      Pattern& pat = *(link->GetPattern());
      if( mirrored ) return shift - pat[i];
      else           return shift + pat[i];
    }
    // Comparison operators
    bool operator<( const NodeDescriptor& rhs ) const {
      if( shift < rhs.shift ) return true;
      if( shift > rhs.shift ) return false;
      UInt_t nb = link->GetPattern()->GetNbits();
      for( UInt_t i = 1; i < nb; ++i ) {
	if( (*this)[i] < rhs[i] )  return true;
	if( (*this)[i] > rhs[i] )  return false;
      }
      return false; // Patterns are equal
    }
    bool operator<=( const NodeDescriptor& rhs ) const {
      if( shift < rhs.shift ) return true;
      if( shift > rhs.shift ) return false;
      UInt_t nb = link->GetPattern()->GetNbits();
      for( UInt_t i = 1; i < nb; ++i ) {
	if( (*this)[i] < rhs[i] )  return true;
	if( (*this)[i] > rhs[i] )  return false;
      }
      return true; // Patterns are equal
    }
    bool operator>( const NodeDescriptor& rhs ) const {
      return !operator<=(rhs);
    }
    bool operator>=( const NodeDescriptor& rhs ) const {
      return !operator<(rhs);
    }
    bool operator==( const NodeDescriptor& rhs ) const {
      return ( shift == rhs.shift && mirrored == rhs.mirrored &&
	       *(link->GetPattern()) == *(rhs.link->GetPattern()) );
    }
    bool operator!=( const NodeDescriptor& rhs ) const {
      return !operator==(rhs);
    }
  };    

  class HitSet;
  typedef std::pair<NodeDescriptor,HitSet> Node_t;

  /////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
