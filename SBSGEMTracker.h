//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 06-Feb-2012
//
#ifndef ROOT_SBS_GEMTracker
#define ROOT_SBS_GEMTracker

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBS::GEMTracker                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GEMTracker.h"
#include "TString.h"

namespace SBS {

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
				const TVector3& pos, const TVector3& dir,
				const FitRes_t& fit_par );

     // Podd interface
    virtual const char* GetDBFileName() const;
    virtual void  MakePrefix();
    virtual Int_t ReadGeometry( FILE* file, const TDatime& date,
				Bool_t required = kTRUE );

    // Geometry
    Int_t     fSector;     // Sector number (for analysis convenience)
    //Parameters for the chamber localisation
    Double_t  fDMag;       // distance of Spectrometer "Zero" from the Hall Center 
    Double_t  fXOffset;    // X offset of the plane/sector wrt the spectrometer central ray
    Double_t  fThetaH;     // Horizontal Spectrometer rotation (wrt hall pivot)
    Double_t  fThetaV;     // Vertical Spectrometer rotation (wrt Spec "Zero")
    
    // Configuration
    TString   fDBPrefix;     // Safe storage for database file name prefix

#ifdef MCDATA
    virtual Int_t FitMCPoints( Podd::MCTrack* mctrk ) const;
#endif

    ClassDef(GEMTracker,0)   // Collection of GEM trackers in one SBS sector
  };

}

///////////////////////////////////////////////////////////////////////////////

#endif
