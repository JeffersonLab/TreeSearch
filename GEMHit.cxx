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
#ifdef MCDATA
ClassImp(TreeSearch::MCGEMHit)
#endif

namespace TreeSearch {

//_____________________________________________________________________________
void GEMHit::Print( Option_t* opt ) const
{
  // Print hit info

  Hit::Print("C");
  cout << " size="  << GetSize()
       << " type="  << GetType()
       << " time="  << GetPeaktime();

  if( *opt != 'C' )
    cout << endl;
}

#ifdef MCDATA
//_____________________________________________________________________________
void MCGEMHit::Print( Option_t* ) const
{
  // Print hit info

  GEMHit::Print("C");
  MCPrint();
}
  // Int_t MCGEMHit::Compare(const TObject *obj) const
  // {
  //    if(((MCGEMHit*)obj)->GetPos()==GetPos()) return 0;
  //  }
#endif



///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

