#ifndef ROOT_TreeSearch_Road
#define ROOT_TreeSearch_Road

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Road                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <vector>
#include <set>

namespace TreeSearch {

  class Hit;
  class Projection;
  class NodeDescriptor; // Defined in TreeWalk.h
  class BuildInfo_t;    // Defined in implementation

  class Road {

  public:
    explicit Road( const Projection* proj );
    Road( const Road& );
    Road& operator=( const Road& );
    virtual ~Road();

    Bool_t Add( const NodeDescriptor& nd );
    void   Finish();

    void Print( Option_t* opt="" ) const;

  protected:

    // Corrdinates of hit positions for fitting
    struct Point {
      Double_t x, z;
    };

    // Bin numbers defining the corners
    UShort_t  fLeft[2];   // Left corner bin, 0=lower, 1=upper
    UShort_t  fRight[2];

    std::set<Hit*>     fHits;   // All hits collected from this road's patterns
    std::vector<Point> fPoints; // Hit positions within the road

    // Fit results
    Double_t  fSlope;
    Double_t  fPos;
    Double_t  fChi2;
    Double_t  fErr[2];

    // Data used while building
    BuildInfo_t* fBuild;

    Bool_t CheckMatch( const std::set<Hit*>& hits ) const;
    void   CollectCoordinates();

    ClassDef(Road,1)  // Region containing track candidate hits 
  };


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
