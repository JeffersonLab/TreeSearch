///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::PatternTree                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PatternTree.h"
#include <iostream>
#include <cstdlib>
using namespace std;

ClassImp(TreeSearch::PatternTree)

namespace TreeSearch {

// Tree node with extra fields, used for building the pattern tree
class FullNode : public TreeNode {
public:
protected:
  UShort_t    fMaxDepth;
  UShort_t    fNbits;
  Int_t       fRefCount;
};

//_____________________________________________________________________________
PatternTree::PatternTree()
{

}

//_____________________________________________________________________________
PatternTree::~PatternTree()
{

}

//_____________________________________________________________________________
void PatternTree::Make()
{
  
  int size = sizeof(TreeNode)+2*sizeof(UShort_t);
  div_t isiz = div(size,sizeof(int*));
  if( isiz.rem )
    isiz.quot++;
  TreeNode* tn = (TreeNode*)new int*[isiz.quot];
  new(tn) TreeNode;
  for(int i=0; i<4; i++)
    tn->fBit[i] = 2*i+1;
  for(int i=0; i<4; i++)
    cout << tn->fBit[i] << endl;
    
  delete [] tn;

}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch
