#ifndef ROOT_TreeSearch_EProjType
#define ROOT_TreeSearch_EProjType

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::EProjType                                                     //
//                                                                           //
// Enumerator for wire projection types (U, V, X) upported by the BigBite    //
// MWDC tracking detector class.                                             //
//                                                                           //
// Type Y is currently left out to avoid unnecessary empty vector elements   //
// and loop iterations.                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

namespace TreeSearch {

  // Types of projections
  enum EProjType { kUndefinedType = -1, kUPlane, kTypeBegin = kUPlane,
		   kVPlane, kXPlane, kTypeEnd };

  // Overloaded operator++ allows us to iterate over the EProjType enum.
  inline
  EProjType& operator++( EProjType& e )
  { switch(e) { 
    case kUndefinedType: return e=kUPlane;
    case kUPlane: return e=kVPlane;
    case kVPlane: return e=kXPlane;
    case kXPlane: return e=kTypeEnd;
    case kTypeEnd: default: return e=kUndefinedType;
    }
  }
  inline
  const EProjType operator++( EProjType& e, int )
  { EProjType r(e); ++e; return r; }

}

// Other common types used by TreeSearch

#include <utility>
#include <vector>

namespace TreeSearch {

  typedef std::pair<double,double>    pdbl_t;
  typedef std::vector<pdbl_t>         vec_pdbl_t;

}

///////////////////////////////////////////////////////////////////////////////

#endif
