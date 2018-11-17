//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 11-Jan-2010
//
#ifndef ROOT_TreeSearch_GEMTracker
#define ROOT_TreeSearch_GEMTracker

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::GEMTracker                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Tracker.h"

namespace TreeSearch {

  class GEMTracker : public Tracker {
  public:
    GEMTracker( const char* name, const char* description = "",
		THaApparatus* app = 0 );
    virtual ~GEMTracker();

    virtual EStatus Init( const TDatime& date );

    // Analysis control flags. Set via database.
    enum {
      kDoADCCut   = BIT(14), // Use ADC cut in decoder (defined in planes)
      k3dCorrAmpl = BIT(17), // Use amplitude correation 3D matching algorithm
    };

  protected:

    // Parameters for 3D projection matching
    UInt_t         fNPlanePartnerPairs;  // Number of plane pairs available for matching
    UInt_t         fMaxCorrMismatches;   // Max # planes w/o amplitude match
    Double_t       fMaxCorrNsigma;       // Amplitude correlation cutoff (#sig)
    Double_t       fminAdcToCheck;       // skip adc corr check if one projection has adc less than this value and the other doesn't have hit

    UInt_t MatchRoadsCorrAmpl( vector<Rvec_t>& roads, UInt_t ncombos,
			       std::list<std::pair<Double_t,Rvec_t> >& combos_found,
			       Rset_t& unique_found );

    virtual UInt_t GetCrateMapDBcols() const;

    virtual Plane* MakePlane( const char* name, const char* description = "",
			      THaDetectorBase* parent = 0 ) const;
    virtual UInt_t MatchRoadsImpl( vector<Rvec_t>& roads, UInt_t ncombos,
				   std::list<std::pair<Double_t,Rvec_t> >& combos_found,
				   Rset_t& unique_found );

    virtual THaAnalysisObject::EStatus PartnerPlanes();

     // Podd interface
    virtual Int_t  ReadDatabase( const TDatime& date );

#ifdef MCDATA
    class GEMPointUpdater : public MCPointUpdater {
    public:
      virtual void UpdateHit( Podd::MCTrackPoint* pt, Hit* hit,
			      Double_t x ) const;
    };
#endif

    ClassDef(GEMTracker,0)   // TreeSearch track reconstruction for GEM trackers
  };

}

///////////////////////////////////////////////////////////////////////////////

#endif
