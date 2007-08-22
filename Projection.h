#ifndef ROOT_TreeSearch_Projection
#define ROOT_TreeSearch_Projection

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include <vector>
#include <string>

using std::vector;
using std::string;

namespace TreeSearch {

  class Hitpattern;
  class PatternTree;
  class WirePlane;

  class Projection {
  public:

    Projection( Int_t type, const char* name, Double_t angle );
    Projection( const Projection& orig );
    Projection& operator=( const Projection& rhs );
    virtual ~Projection();

    void  Reset();
    Int_t Init( UInt_t depth = 0 );

    void           AddPlane( WirePlane* wp );
    Int_t          GetType() const { return fType; }
    const char*    GetName() const { return fName.c_str(); }
    Double_t       GetMaxSlope() const { return fMaxSlope; }
    Double_t       GetDepth() const { return fZmax - fZmin; }
    Double_t       GetZmin()  const { return fZmin; }
    Double_t       GetZmax()  const { return fZmax; }
    Double_t       GetWidth() const { return fWidth; }
    Double_t       GetAngle() const;
    Double_t       GetSinAngle() const { return fSinAngle; }
    Double_t       GetCosAngle() const { return fCosAngle; }
    Hitpattern*    GetHitpattern() const { return fHitpattern; }

    void           SetMaxSlope( Double_t m ) { fMaxSlope = m; }
    void           SetPatternTree( PatternTree* pt ) { fPatternTree = pt; }
    void           SetWidth( Double_t width ) { fWidth = width; }

    //FIXME: for testing
    vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }

  private:
    Int_t               fType;        // Type of plane (u,v,x,y...)
    string              fName;        // Type name ('u', 'v', etc.)
    vector<WirePlane*>  fPlanes;      // Wire planes of this type
    Double_t            fMaxSlope;    // Maximum physical track slope (0=perp)
    Double_t            fWidth;       // Width of tracking region (m)
    Double_t            fZmin;        // z-position of first wire plane
    Double_t            fZmax;        // z-position of last wire plane
    Double_t            fSinAngle;    // Sine of wire angle
    Double_t            fCosAngle;    // Cosine of wire angle

    Hitpattern*         fHitpattern;  // Hitpattern of current event
    PatternTree*        fPatternTree; // Precomputed template database

    ClassDef(Projection,1)  // A track projection plane
  };

  //___________________________________________________________________________
  inline
  Double_t Projection::GetAngle() const {
    // Return wire angle in rad, normalized to [-pi,pi]
    Double_t a = TMath::ASin(fSinAngle);
    if( fCosAngle < 0.0 )
      return (fSinAngle > 0.0) ? TMath::TwoPi() - a : -TMath::TwoPi() - a;
  
    return a;
  }

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
