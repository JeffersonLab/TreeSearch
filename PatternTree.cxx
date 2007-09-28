///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternTree.h"
#include <iostream>
using namespace std;

ClassImp(TreeSearch::PatternTree)

namespace TreeSearch {

//_____________________________________________________________________________
PatternTree::PatternTree( const TreeParam_t& parameters )
  : fParameters(parameters), fRoot(0)
{
  // Constructor. 

  Normalize(fParameters);
}

//_____________________________________________________________________________
  PatternTree::PatternTree( const TreeParam_t& parameters, UInt_t nPatterns,
			    UInt_t nLinks )
  : fParameters(parameters), fRoot(0)
{
  // Constructor. 

  Normalize(fParameters);
}

//_____________________________________________________________________________
PatternTree::~PatternTree()
{
  // Destructor

  // TODO: delete tree
}

//_____________________________________________________________________________
Int_t PatternTree::AddPattern(int i)
{

  cout << "AddPattern:" <<  i << endl;
  return 0;
}

//_____________________________________________________________________________
Int_t PatternTree::Normalize( TreeParam_t& param )
{
  // Normalize parameter values to width = 1 and z_max = 1.


  return 0;
}

//_____________________________________________________________________________

//_____________________________________________________________________________

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
