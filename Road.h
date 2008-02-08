#ifndef ROOT_TreeSearch_Road
#define ROOT_TreeSearch_Road

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Road                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>

using std::vector;

namespace TreeSearch {

  class NodeDescriptor; // Defined in TreeWalk.h
  class BuildInfo_t;    // Defined in implementation

  class Road {

  public:
    explicit Road( NodeDescriptor* nd ) {}
    Road( const Road& ) {}
    Road& operator=( const Road& );
    virtual ~Road() {}

    Bool_t Add( NodeDescriptor* nd );
    void   Close();

    void Print( Option_t* opt="" ) const;

  protected:

    // Bin numbers defining the corners
    UShort_t  fStart[2];
    UShort_t  fEnd[2];

    // Fit results
    Double_t  fSlope;
    Double_t  fPos;
    Double_t  fChi2;
    Double_t  fErr[2];

    // Data used while building
    BuildInfo_t* fBifo;
    //TODO:
    // array of hit positions within the road (x/z coordinates)
    // function to collect these
    // array of NodeDescriptors* collected
    // array of common hits
    // array of counters of common hits per plane
    // match evaluation function (deals with missing planes, 
    //     1 out of 2 handling etc)
    

    ClassDef(Road,1)  // Region containing track candidate hits 
  };


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
