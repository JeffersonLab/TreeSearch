//*-- Author :    Ole Hansen, Jefferson Lab   16-Aug-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Projection.h"
#include "Hitpattern.h"
#include "WirePlane.h"

namespace TreeSearch {

//_____________________________________________________________________________
Projection::Projection( Int_t type, const char* name, Double_t angle ) 
  //, const PatternTree* pt )
  : fType(type), fName(name?name:""), fMaxSlope(0.0), fWidth(0.0), 
    fHitpattern(NULL), fPatternTree(NULL)
    //    fPatternTree(pt)
{
  // Constructor

  fSinAngle = TMath::Sin(angle);
  fCosAngle = TMath::Cos(angle);

}


//_____________________________________________________________________________
Projection::Projection( const Projection& orig )
  : fType(orig.fType), fName(orig.fName), fPlanes(orig.fPlanes), 
    fMaxSlope(orig.fMaxSlope), fWidth(orig.fWidth),
    fSinAngle(orig.fSinAngle), fCosAngle(orig.fCosAngle),
    fHitpattern(NULL), fPatternTree(NULL)
{
  // Copying

  if( orig.fHitpattern )
    fHitpattern = new Hitpattern(*orig.fHitpattern);
//   if( orig.fPatternTree )
//     fHitpattern = new PatternTree(*orig.fPatternTree);
}

//_____________________________________________________________________________
Projection& Projection::operator=( const Projection& rhs )
{
  // Assignment

  if( this != &rhs ) {
    fType     = rhs.fType;
    fName     = rhs.fName;
    fPlanes   = rhs.fPlanes;
    fMaxSlope = rhs.fMaxSlope;
    fWidth    = rhs.fWidth;
    fSinAngle = rhs.fSinAngle;
    fCosAngle = rhs.fCosAngle;
    delete fHitpattern;
    if( rhs.fHitpattern )
      fHitpattern = new Hitpattern(*rhs.fHitpattern);
    else
      fHitpattern = NULL;
//     if( rhs.fPatternTree )
//       fPatternTree = new PatternTree(*rhs.fPatternTree);
//     else
//       fPatternTree = NULL;
  }
  return *this;
}

//_____________________________________________________________________________
Projection::~Projection()
{
  // Destructor

  //  delete fPatternTree;
  delete fHitpattern;
}

//_____________________________________________________________________________
void Projection::Reset()
{
  // Reset parameters, clear list of planes, delete Hitpattern
  
  fPlanes.clear();
  fMaxSlope = fWidth = 0.0;
  delete fHitpattern; fHitpattern = NULL;
  //  delete fPatternTree; fPatternTree = NULL;
}

//_____________________________________________________________________________
Int_t Projection::Init( UInt_t depth )
{
  // Initialize hitpattern

  if( depth == 0 && !fPatternTree )
    return -1;
  Int_t nplanes = fPlanes.size();
  if( nplanes == 0 )
    return -2;

  delete fHitpattern;
  
  if( depth > 0 )
    fHitpattern = new Hitpattern( depth, nplanes, fWidth );
  //FIXME:
//   else
//     fHitpattern = new Hitpattern( *fPatternTree );

  return 0;
}

//_____________________________________________________________________________

}  

ClassImp(TreeSearch::Projection)

///////////////////////////////////////////////////////////////////////////////
