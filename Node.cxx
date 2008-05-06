///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Node                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Node.h"
#include <iostream>

using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
void NodeDescriptor::Print() const
{
  // Print bit contents of the pattern of this NodeDescriptor

  cout << "(" << (UInt_t)depth << "/" << (mirrored ? "-" : "+") << ") ";
  for( UInt_t i=0; i<link->GetPattern()->GetNbits(); i++ ) {
    cout << (*this)[i] << " ";
  }
  cout << endl;
}


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
