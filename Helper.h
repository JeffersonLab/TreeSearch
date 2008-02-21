#ifndef ROOT_TreeSearch_Helper
#define ROOT_TreeSearch_Helper

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Helper                                                        //
//                                                                           //
// Helper classes and functions                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>
#include <cassert>
#include <functional>

namespace TreeSearch {

  //___________________________________________________________________________
  template< typename VectorElem > void 
  NthPermutation( UInt_t n, 
		  std::vector<std::vector<VectorElem> >& vec,
		  std::vector<VectorElem*>& selected )
  {
    // Get the n-th permutation of the elements in vec and
    // put result in selected. selected[k] is one of the
    // vec[k].size() elements in the k-th plane.
    // The output vector the_points must be resized to vec.size()
    // before calling this function.

    assert( !vec.empty() and selected.size() == vec.size() );

    UInt_t vsiz = vec.size(), k;
    for( UInt_t j = vsiz; j--; ) {
      typename std::vector<VectorElem>::size_type npt = vec[j].size();
      assert(npt);
      if( npt == 1 )
	k = 0;
      else {
	k = n % npt;
	n /= npt;
      }
      // Put pointer to selected element in result
      selected[j] = &vec[j][k];
    }
  }

  //___________________________________________________________________________
  template< typename VectorElem > struct SizeMul :
    public std::binary_function< typename std::vector<VectorElem>::size_type, 
				 std::vector<VectorElem>, 
				 typename std::vector<VectorElem>::size_type >
  {
    // Get the product of the sizes of the vectors in a vector, ignoring 
    // empty ones
    typedef typename std::vector<VectorElem>::size_type vsiz_t;
    vsiz_t operator() ( vsiz_t val, const std::vector<VectorElem>& vec ) const 
    { return ( vec.empty() ? val : val * vec.size() ); }
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
