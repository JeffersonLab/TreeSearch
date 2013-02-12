//*-- Author :    Ole Hansen, Jefferson Lab   05-Mar-2010

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hit                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GEMHit.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(TreeSearch::GEMHit)
ClassImp(TreeSearch::MCGEMHit)

namespace TreeSearch {

//_____________________________________________________________________________
void GEMHit::Print( Option_t* opt ) const
{
  // Print hit info

  cout << "Hit: plane=" << (fPlane ? fPlane->GetName() : "??")
       << " pos="       << GetPos()
       << " z="         << GetZ()
       << " res="       << GetResolution()
       << " size="      << GetSize()
       << " type="      << GetType();
  if( *opt != 'C' )
    cout << endl;
}

//_____________________________________________________________________________
void MCGEMHit::Print( Option_t* ) const
{
  // Print hit info

  GEMHit::Print("C");
  cout << " MCpos=" << GetMCPos()
       << endl;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

