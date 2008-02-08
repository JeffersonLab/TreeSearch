//*-- Author :    Ole Hansen, Jefferson Lab   07-Feb-2008

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Road.h"

#include <iostream>

using namespace std;

ClassImp(TreeSearch::Road);

namespace TreeSearch {

//_____________________________________________________________________________
Road& Road::operator=( const Road& rhs )
{
  // Print hit info

  if( this != &rhs ) {

  }
  return *this;
}

//_____________________________________________________________________________
void Road::Print( Option_t* opt ) const
{
  // Print road info

}


//_____________________________________________________________________________

} // end namespace TreeSearch

///////////////////////////////////////////////////////////////////////////////
