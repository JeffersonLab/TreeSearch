#ifndef ROOT_TreeSearch_EProjType
#define ROOT_TreeSearch_EProjType

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::EProjType                                                     //
//                                                                           //
// Enumerator for wire projection types (U, V, X, Y) supported by the        //
// GEM tracking detector class.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdexcept>

namespace TreeSearch {

  // Types of projections
  enum EProjType { kUndefinedType = -1, kUPlane, kTypeBegin = kUPlane,
		   kVPlane, kXPlane, kYPlane, kTypeEnd };

  // Overloaded prefix operator++ allows us to iterate over the EProjType enum.
  inline
  EProjType& operator++( EProjType& e )
  {
    switch(e) {
    case kUndefinedType: return e=kUPlane;
    case kUPlane: return e=kVPlane;
    case kVPlane: return e=kXPlane;
    case kXPlane: return e=kYPlane;
    case kYPlane: return e=kTypeEnd;
    case kTypeEnd: default:
      throw std::range_error("EProjType out of range trying to increment "
			     "kTypeEnd");
    }
  }
  // Postfix operator++
  inline
  EProjType operator++( EProjType& e, int )
  { EProjType r(e); ++e; return r; }

  // Projection default parameters
  struct ProjectionParam {
    const char*   name;
    double        angle;   // default angle, can be overridden from database
  };
#ifndef __CINT__
  static const ProjectionParam kProjParam[] = {
    { "u", -60.0 },
    { "v",  60.0 },
    { "x",   0.0 },
    { "y",  90.0 }
  };
#endif
}

///////////////////////////////////////////////////////////////////////////////

#endif
