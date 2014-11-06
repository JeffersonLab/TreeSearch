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

    void SetSectorNumber( Int_t i ) { fSector = i; }

  protected:

    virtual Plane* MakePlane( const char* name, const char* description = "",
			      THaDetectorBase* parent = 0 ) const;

    virtual THaAnalysisObject::EStatus PartnerPlanes();
    virtual Int_t NewTrackCalc( Int_t idx, THaTrack* newTrack,
				const TVector3& pos, const TVector3& dir );

     // Podd interface
    virtual const char* GetDBFileName() const;
    virtual void  MakePrefix();
    virtual Int_t ReadGeometry( FILE* file, const TDatime& date,
				Bool_t required = kTRUE );

    // Geometry
    Int_t     fSector;       // Sector number (for analysis convenience)
    Double_t  fPhi;          // Phi rotation of this sector (rad)

    // Configuration
    TString   fDBPrefix;     // Safe storage for database file name prefix

#ifdef MCDATA
    virtual Int_t FitMCPoints( Podd::MCTrack* mctrk ) const;
#endif

    ClassDef(GEMTracker,0)   // Collection of GEM trackers in one SoLID sector
  };

}

///////////////////////////////////////////////////////////////////////////////

#endif
