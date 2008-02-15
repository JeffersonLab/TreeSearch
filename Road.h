#ifndef ROOT_TreeSearch_Road
#define ROOT_TreeSearch_Road

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Road                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <set>
#include <utility>
#include <vector>
#include <list>

namespace TreeSearch {

  class Projection;
  class NodeDescriptor;
  class Hit;
  class HitSet;
  class BuildInfo_t;    // Defined in implementation

  class Road {

  public:
    typedef std::pair<const NodeDescriptor,HitSet> Node_t;

    explicit Road( const Projection* proj );
    Road( const Road& );
    Road& operator=( const Road& );
    virtual ~Road();

    Bool_t Add( Node_t& nd );
    void   Finish();
    Bool_t Fit();
    Bool_t IsGood() const { return fGood; }
    void   Print( Option_t* opt="" ) const;

  protected:

    // Coordinates of hit positions, for track fitting
    struct Point {
      Double_t x, z;
      Point() : x(0), z(0) {}
      Point( Double_t _x, Double_t _z ) : x(_x), z(_z) {}
    };

    const Projection* fProjection;     // Projection that this Road belongs to

    std::vector<Double_t> fCornerX;    // x positions of corners
    Double_t              fZL, fZU;    // z +/- eps of first/last plane 


    std::list<Node_t*>     fPatterns;   // Patterns in this road
    std::set<Hit*>         fHits;       // All hits linked to the patterns

    // Fit results
    Double_t  fSlope;
    Double_t  fPos;
    Double_t  fChi2;
    Double_t  fErr[2];

    Bool_t    fGood;     // Road successfully built and fit

    // Data used while building
    BuildInfo_t* fBuild;

#ifdef TESTCODE
    Node_t*   fSeed;
#endif

    Bool_t   CheckMatch( const std::set<Hit*>& hits ) const;
    Bool_t   CollectCoordinates( std::vector<Point>& points,
		       std::vector<std::vector<Point*> >& planepoints);
    Double_t GetBinX( UInt_t bin ) const;

    ClassDef(Road,1)  // Region containing track candidate hits 
  };


///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
