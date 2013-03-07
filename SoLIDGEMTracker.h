//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 06-Feb-2012
//
#ifndef ROOT_SoLID_GEMTracker
#define ROOT_SoLID_GEMTracker

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SoLID::GEMTracker                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GEMTracker.h"
#include "TString.h"

namespace SoLID {

  using TreeSearch::Plane;

  class GEMTracker : public TreeSearch::GEMTracker {
  public:
    GEMTracker( const char* name, const char* description = "", 
		THaApparatus* app = 0 );
    virtual ~GEMTracker();

  protected:

    virtual Plane* MakePlane( const char* name, const char* description = "",
			      THaDetectorBase* parent = 0 ) const;

    virtual THaAnalysisObject::EStatus PartnerPlanes();

     // Podd interface
    virtual const char* GetDBFileName() const;
    virtual void  MakePrefix();
    virtual Int_t ReadGeometry( FILE* file, const TDatime& date,
				Bool_t required = kTRUE );

    // Geometry
    Double_t  fPhi;          // Phi rotation of this sector (rad)

    // Configuration
    TString   fDBPrefix;     // Safe storage for database file name prefix

    ClassDef(GEMTracker,0)   // Collection of GEM trackers in one SoLID sector
  };

}

///////////////////////////////////////////////////////////////////////////////

#endif
