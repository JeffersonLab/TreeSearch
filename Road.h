#ifndef ROOT_TreeSearch_Road
#define ROOT_TreeSearch_Road

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Road                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

namespace TreeSearch {

  class NodeDescriptor;

  class Road {

  public:
    Road( NodeDescriptor* nd ) {}
    Road( const Road& ) {}
    Road& operator=( const Road& );
    virtual ~Road() {}

    Int_t Add( NodeDescriptor* nd ) { return 0;}
    Int_t Close() { return 0; }
    UInt_t GetNpat() const { return 0; };

    void Print( Option_t* opt="" ) const;

  protected:

    ClassDef(Road,1)  // Region containing track candidate hits 
  };


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
