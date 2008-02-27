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
  NthCombination( UInt_t n, 
		  std::vector<std::vector<VectorElem> >& vec,
		  std::vector<VectorElem>& selected )
  {
    // Get the n-th permutation of the elements in vec and
    // put result in selected. selected[k] is one of the
    // vec[k].size() elements in the k-th plane.

    assert( !vec.empty() );

    UInt_t vsiz = vec.size(), k;
    selected.resize( vsiz );
    for( UInt_t j = vsiz; j--; ) {
      typename std::vector<VectorElem>::size_type npt = vec[j].size();
      assert(npt);
      if( npt == 1 )
	k = 0;
      else {
	k = n % npt;
	n /= npt;
      }
      // Copy the selected element
      selected[j] = vec[j][k];
    }
  }

  //___________________________________________________________________________
  template< typename VectorElem > struct SizeMul :
    public std::binary_function< typename std::vector<VectorElem>::size_type, 
				 std::vector<VectorElem>, 
				 typename std::vector<VectorElem>::size_type >
  {
    // Function object to get the product of the sizes of the vectors in a 
    // vector, ignoring empty ones
    typedef typename std::vector<VectorElem>::size_type vsiz_t;
    vsiz_t operator() ( vsiz_t val, const std::vector<VectorElem>& vec ) const
    { return ( vec.empty() ? val : val * vec.size() ); }
  };

  //___________________________________________________________________________
  // Iterator over unique combinations of k out of N elements. The element 
  // numbers are accessible in a vector<int> of size k.
  class UniqueCombo {
    typedef std::vector<int> vint_t;
  public:
    UniqueCombo( int N, int k ) : fN(N), fGood(true)
    { assert( 0<k && k<=N ); for( int i=0; i<k; ++i ) fCurrent.push_back(i); }
    // Default copy and assignment are fine

    UniqueCombo& operator++()
    { fGood = recursive_plus( fCurrent.size()-1 ); return *this; }
    const UniqueCombo operator++(int)
    { UniqueCombo clone(*this); ++*this; return clone; }
    const vint_t& operator()() const { return fCurrent; }
    const vint_t& operator*()  const { return fCurrent; }
    bool operator==( const UniqueCombo& rhs ) const
    { return ( fN == rhs.fN and fCurrent == rhs.fCurrent ); }
    bool operator!=( const UniqueCombo& rhs ) const { return !(*this==rhs); }
    operator bool() const  { return fGood; }
    bool operator!() const { return !((bool)*this); }
    
  private:
    bool recursive_plus( vint_t::size_type pos ) {
      if( fCurrent[pos] < fN+int(pos-fCurrent.size()) ) {
	++fCurrent[pos];
	return true;
      }
      if( pos==0 )
	return false;
      if( recursive_plus(pos-1) ) {
	fCurrent[pos] = fCurrent[pos-1]+1;
	return true;
      }
      return false;
    }
    int    fN;
    vint_t fCurrent;
    bool   fGood;
  };

  //___________________________________________________________________________
  template< typename Container >
  void DeleteContainer( Container& c )
  {
    // Delete all elements of given container of pointers
    for( typename Container::iterator it = c.begin(); it != c.end(); ++it ) {
      delete *it;
    }
    c.clear();
  }

  //___________________________________________________________________________
  template< typename ContainerOfContainers >
  void DeleteContainerOfContainers( ContainerOfContainers& cc )
  {
    // Delete all elements of given container of containers of pointers
    for( typename ContainerOfContainers::iterator it = cc.begin();
	 it != cc.end(); ++it ) {
      DeleteContainer( *it );
    }
    cc.clear();
  }
  
///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
