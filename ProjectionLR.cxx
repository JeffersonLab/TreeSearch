//*-- Author :    Ole Hansen, Jefferson Lab   11-Feb-2013

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::ProjectionLR                                                  //
//                                                                           //
// A collection of horizontal drift chamber wire planes with a given wire    //
// direction (u,v,x etc.). Uses Hitpattern with extra support for hits       //
// with left/right ambiguity and staggered wire planes to resolve some       //
// of the ambiguity.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "ProjectionLR.h"
#include "HitpatternLR.h"

using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
ProjectionLR::ProjectionLR( EProjType type, const char* name, Double_t angle,
			    THaDetectorBase* parent )
  : Projection( type, name, angle, parent )
{
  // Constructor
}

//_____________________________________________________________________________
ProjectionLR::~ProjectionLR()
{
  // Destructor

}

//_____________________________________________________________________________
Hitpattern* ProjectionLR::MakeHitpattern( const PatternTree& pt )
{
  // Instantiate a HitpatternLR that interprets hits as L/R-ambiguous wire hits.

  return new HitpatternLR( pt );
}

//_____________________________________________________________________________

}  // end namespace TreeSearch

ClassImp(TreeSearch::ProjectionLR)

///////////////////////////////////////////////////////////////////////////////
