//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 11-Feb-2013
//
#ifndef ROOT_TreeSearch_ProjectionLR
#define ROOT_TreeSearch_ProjectionLR

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::ProjectionLR                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Projection.h"

namespace TreeSearch {

  class ProjectionLR : public Projection {
  public:

    ProjectionLR( EProjType type, const char* name, Double_t angle,
		  THaDetectorBase* parent );
    virtual ~ProjectionLR();

    virtual Hitpattern* MakeHitpattern( const PatternTree& ) const;

  private:
    // Prevent default copying, assignment
    ProjectionLR( const Projection& orig );
    const ProjectionLR& operator=( const Projection& rhs );

    ClassDef(ProjectionLR,0)  // A track projection plane with LR hitpattern
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
