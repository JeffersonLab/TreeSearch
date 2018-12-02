#ifndef ROOT_TreeSearch_Road
#define ROOT_TreeSearch_Road

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Road                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Hit.h"
#include "TVector2.h"
#include <set>
#include <utility>
#include <vector>
#include <list>
#include <functional>
#include <cassert>
#include <cstring>
#include <unordered_map>
#include <utility>

class THaTrack;

using std::vector;

namespace TreeSearch {

  class Projection;
  class BuildInfo_t;    // Defined in implementation

  class Road : public TObject {

  public:
    //_________________________________________________________________________
    // Coordinates of hit positions, for track fitting
    struct Point {
      Point() : x(0), hit(0) {}
      Point( Double_t _x, Double_t _z, Hit* _hit )
	: x(_x), z(_z), hit(_hit) { assert(hit); }
      virtual ~Point() {}
      Double_t res() const { return hit->GetResolution(); }
      Double_t x;    // Selected x coordinates
      Double_t z;    // z coordinate
      Hit*     hit;  // Associated hit (stored in Plane)
      ClassDef(Point,0)
    };

    typedef std::vector<Road::Point*>  Pvec_t;
    typedef std::list<const Node_t*>   NodeList_t;

    //_________________________________________________________________________
    // For global variable access/event display
    friend class Corners;
    class Corners : public TObject {
    public:
      explicit Corners( Road* rd )
	: fXLL(rd->fCornerX[0]), fXLR(rd->fCornerX[1]), fZL(rd->fZL),
	  fXUL(rd->fCornerX[3]), fXUR(rd->fCornerX[2]), fZU(rd->fZU) {}
      Corners() {}  // For ROOT RTTI
      virtual ~Corners() {}
    private:
      Double_t fXLL;  // Lower left corner x coord
      Double_t fXLR;  // Lower right corner x coord
      Double_t fZL;   // Lower edge z coord
      Double_t fXUL;  // Upper left corner x coord
      Double_t fXUR;  // Upper right corner x coord
      Double_t fZU;   // Upper edge z coord
      ClassDef(Corners,0)
    };

    //_________________________________________________________________________
    explicit Road( const Projection* proj );
    Road( const Node_t& nd, const Projection* proj );
    Road() : fProjection(0), fTrack(0), fBuild(0) {} // For internal ROOT use
    Road( const Road& );
    Road& operator=( const Road& );
    virtual ~Road();

    Bool_t         Add( const Node_t& nd );
    void           ClearGrow() { fGrown = false; }
    virtual Int_t  Compare( const TObject* obj ) const;
    void           Finish();
    Bool_t         Fit(std::unordered_map<Int_t, Int_t> &moduleOrder);
    Double_t       GetChi2()    const { return fChi2; }
    UInt_t         GetNdof()    const { return fDof; }
    const Hset_t&  GetHits()    const { return fHits; }
    const Pvec_t&  GetPoints()  const { return fFitCoord; }
    const vector<std::pair<Pvec_t, Double_t>>& GetAllCandidatePoints() const { return fAllFitCoord; }
    Double_t       GetPos()     const { return fPos; }
    Double_t       GetPos( Double_t z ) const { return fPos + z*fSlope; }
    Double_t       GetPosErrsq( Double_t z ) const;
    UInt_t         GetPlanePattern()  const { return fPlanePattern; }
    const Projection* GetProjection() const { return fProjection; }
    Double_t       GetSlope()   const { return fSlope; }
    THaTrack*      GetTrack()   const { return fTrack; }
    Bool_t         HasGrown()   const { return fGrown; }
    Bool_t         Include( const Road* other );
    TVector2       Intersect( const Road* other, Double_t z ) const;
    Bool_t         IsGood()     const { return fGood; }
    Bool_t         IsInBackRange( const Node_t& nd ) const;
    Bool_t         IsInFrontRange( const Node_t& nd ) const;
    Bool_t         IsInRange( const Node_t& nd ) const;
    virtual Bool_t IsSortable() const { return kTRUE; }
    Bool_t         IsVoid()     const { return !fGood; }
    virtual void   Print( Option_t* opt="" ) const;
    void           SetGrow() { fGrown = true; }
    void           SetTrack( THaTrack* track ) { fTrack = track; }
    void           SetSigRatio(double ratio) { fSigRatio = ratio; }
    void           Void() { fGood = false; }
    Bool_t   UpdateFitCoord(const Pvec_t& fitCoord);

#ifdef MCDATA
    // MC truth data, for diagnosing tracking failures
    UInt_t GetNMCTrackHits()        const {return fNMCTrackHits;}        // # planes with hits from MC track
    UInt_t GetMCTrackPlanePattern()  const {return fMCTrackPlanePattern;} // Planepattern of MC track hits
    UInt_t GetNMCTrackHitsFit()       const {return fNMCTrackHitsFit;}     // # MC track hits used in best fit
    UInt_t GetMCTrackPlanePatternFit() const {return fMCTrackPlanePatternFit;} // Planepattern of fitted MC track hits
#endif

#ifdef VERBOSE
    const NodeList_t& GetPatterns() const { return fPatterns; }
#endif

    struct PosIsLess
      : public std::binary_function< Road*, Road*, bool >
    {
      bool operator() ( const Road* a, const Road* b ) const
      { return ( a->GetPos() < b->GetPos() ); }
    };

    class PosIsNear {
    public:
      PosIsNear( Double_t tolerance ) : fTol(tolerance) {}
      bool operator() ( const Road* rd, Double_t pos ) const
      { return ( rd->GetPos() + fTol < pos ); }
      bool operator() ( Double_t pos, const Road* rd ) const
      { return ( pos + fTol < rd->GetPos() ); }
    private:
      Double_t fTol;
    };

    enum ETrackingStatus {
      kTrackOK              = 0,
      kTooFewPlanesWithHits = 1, // Too few planes with hits
      kTooManyHitCombos     = 2, // Too many hit combinations to fit
      kNoGoodFit            = 3, // No fit with good chi2
    };
    ETrackingStatus GetTrackingStatus() const { return fTrkStat; }

  protected:

    NodeList_t     fPatterns;   // Patterns in this road
    Hset_t         fHits;       // All hits linked to the patterns
    vector<Pvec_t> fPoints;     // All hit coordinates within road [nplanes][]
    Pvec_t         fFitCoord;   // fPoints used in best fit [nplanes]
    vector<std::pair<Pvec_t, Double_t>> fAllFitCoord; // all valid fPoints combinations [nthCombination][nplanes]
    UInt_t         fPlanePattern; // Bitpattern of planes in best fit

    const Projection* fProjection; //! Projection that this Road belongs to

    Double_t       fCornerX[5]; // x positions of corners
    Double_t       fZL;         // z-eps of first plane
    Double_t       fZU;         // z+eps of last plane

#ifdef MCDATA
    // MC truth data, for diagnosing tracking failures
    UInt_t         fNMCTrackHits;        // # planes with hits from MC track
    UInt_t         fMCTrackPlanePattern; // Planepattern of MC track hits
    UInt_t         fNMCTrackHitsFit;     // # MC track hits used in best fit
    UInt_t         fMCTrackPlanePatternFit; // Planepattern of fitted MC track hits
#endif

    // Best fit results
    Double_t       fPos;        // Track origin
    Double_t       fSlope;      // Track slope
    Double_t       fChi2;       // Chi2 of fit
    Double_t       fV[3];       // Covar matrix of param (V11, V12=V21, V22)
    UInt_t         fDof;        // Degrees of freedom of fit (nhits-2)

    Bool_t         fGood;       // Road successfully fit
    THaTrack*      fTrack;      // The lowest-chi2 3D track using this road

    BuildInfo_t*   fBuild;      //! Working data for building
    Bool_t         fGrown;      //! Add() added hits in front or back plane
    ETrackingStatus fTrkStat;   // Reconstruction status
    Double_t       fSigRatio;   // number of "signal" hits over total hits in this road, only valid in MC runs 


    // Only needed for TESTCODE
    UInt_t         fNfits;      // Statistics: num fits with acceptable chi2

    Bool_t   CheckMatch( const Hset_t& hits ) const;
    Bool_t   CollectCoordinates();
    

  private:
    void     CopyPointData( const Road& orig );

    ClassDef(Road,0)  // Region containing track candidate hits and fit results
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
  Double_t Road::GetPosErrsq( Double_t z ) const
  {
    // Return square of uncertainty in x = a1+z2*z for best fit (in m^2)

    return fV[0] + 2.0*fV[1]*z + fV[2]*z*z;
  }

  //---------------------------------------------------------------------------
  inline
  Bool_t Road::IsInRange( const Node_t& nd ) const
  {
    // Test if given pattern is within the allowed maximum distance from the
    // current front and back bin ranges of the cluster

    return ( IsInFrontRange(nd) and IsInBackRange(nd) );
  }

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch


#endif
