#ifndef ROOT_TreeSearch_Types
#define ROOT_TreeSearch_Types

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Types.h                                                                   //
//                                                                           //
// Common data types used by TreeSearch                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "EProjType.h"

#include <utility>
#include <vector>

namespace TreeSearch {

  typedef std::pair<double,double>    pdbl_t;
  typedef std::vector<pdbl_t>         vec_pdbl_t;

}

///////////////////////////////////////////////////////////////////////////////

#endif
