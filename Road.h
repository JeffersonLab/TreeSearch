#ifndef ROOT_TreeSearch_Road
#define ROOT_TreeSearch_Road

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Road                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hit.h"
#include <set>
#include <utility>
#include <vector>
#include <list>
#include <functional>
#include <cassert>

namespace TreeSearch {

  class Projection;
  class NodeDescriptor;
  class Hit;
  class HitSet;
  class BuildInfo_t;    // Defined in implementation

  class Road : public TObject {

  public:
    typedef std::pair<const NodeDescriptor,HitSet> Node_t;

    explicit Road( const Projection* proj );
    Road() {} // For internal ROOT use
    Road( const Road& );
    Road& operator=( const Road& );
    virtual ~Road();

    Bool_t         Add( Node_t& nd );
    Bool_t         Adopt( const Road* other );
    virtual Int_t  Compare( const TObject* obj ) const;
    void           Finish();
    Bool_t         Fit();
    Double_t       GetChi2( UInt_t ifit ) const;
    UInt_t         GetNfits() const { return (UInt_t)fFitData.size(); }
    Bool_t         Includes( const Road* other ) const;
    Bool_t         IsGood() const { return fGood; }
    virtual Bool_t IsSortable () const { return kTRUE; }
    virtual void   Print( Option_t* opt="" ) const;

    // Coordinates of hit positions, for track fitting
    struct Point {
      Point() : x(0), hit(0) {}
      Point( Double_t _x, Double_t _z, Hit* _hit ) 
	: x(_x), z(_z), hit(_hit) { assert(hit); }
      Double_t res() const { return hit->GetResolution(); }
      Double_t x;    // Selected x coordinates
      Double_t z;    // z coordinate
      Hit*     hit;  // Associated hit (stored in WirePlane)
    };

  protected:

    const Projection*     fProjection; //! Projection that this Road belongs to

    std::vector<Double_t> fCornerX;    // x positions of corners
    Double_t              fZL, fZU;    // z +/- eps of first/last plane 


    std::list<Node_t*>    fPatterns;   // Patterns in this road
    Hset_t                fHits;       // All hits linked to the patterns
    
    // Fit results
    struct FitResult {
      Double_t  fPos, fSlope, fChi2, fPosErr, fSlopeErr;
      std::vector<Point*>   fFitCoordinates;
      FitResult( Double_t pos, Double_t slope, Double_t chi2,
		 Double_t pos_err, Double_t slope_err )
	: fPos(pos), fSlope(slope), fChi2(chi2),
	  fPosErr(pos_err), fSlopeErr(slope_err) {}
      FitResult() {}
      // Sort fit results by ascending chi2
      bool operator<( const FitResult& rhs ) const 
      { return ( fChi2 < rhs.fChi2 ); }

      struct Chi2IsLess
	: public std::binary_function< FitResult*, FitResult*, bool >
      {
	bool operator() ( const FitResult* a, const FitResult* b ) const
	{ assert(a&&b); return ( a->fChi2 < b->fChi2 ); }
      };
    };

    std::vector<FitResult*> fFitData; // Good fit results, sorted by chi2

    // Best fit results (copy of fFitData.begin() for global variable access)
    Double_t  fPos;      // Track origin
    Double_t  fSlope;    // Track slope
    Double_t  fChi2;     // Chi2 of fit
    UInt_t    fDof;      // Degrees of freedom of fit (nhits-2)

    Bool_t    fGood;     // Road successfully built and fit


    BuildInfo_t* fBuild; //! Working data for building

    Bool_t   CheckMatch( const Hset_t& hits ) const;
    Bool_t   CollectCoordinates( std::vector< std::vector<Point> >& ) const;
    void     DeleteFitResults();
    Double_t GetBinX( UInt_t bin ) const;

    ClassDef(Road,1)  // Region containing track candidate hits and fit results
  };

  //___________________________________________________________________________
  inline
  Int_t Road::Compare( const TObject* obj ) const 
  {
    // Used for sorting Roads in a TClonesArray or similar.
    // A Road is "less than" another if the chi2 of its best fit is smaller.
    // Returns -1 if this is smaller than rhs, 0 if equal, +1 if greater.

    // Require identical classes of objects
    assert( obj && IsA() == obj->IsA() );

    //TODO: take fDof into account & compare statistical significance
    if( fChi2 < static_cast<const Road*>(obj)->fChi2 ) return -1;
    if( fChi2 > static_cast<const Road*>(obj)->fChi2 ) return  1;
    return 0;
  }

  //---------------------------------------------------------------------------
  inline
  void Road::DeleteFitResults()
  {
    // Delete all elements of fFitData

    for( std::vector<FitResult*>::iterator it = fFitData.begin();
	 it != fFitData.end(); ++it )
      delete *it;
    fFitData.clear();
  }
  
  //---------------------------------------------------------------------------
  inline
  Double_t Road::GetChi2( UInt_t ifit ) const
  {
    // Return chi2 of i-th fit
    assert( ifit < GetNfits() );
    return fFitData[ifit]->fChi2;
  }

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
