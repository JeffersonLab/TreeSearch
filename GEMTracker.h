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
    UInt_t         fMaxCorrMismatches;   // Max # planes w/o amplitude match
    Double_t       fMaxCorrNsigma;       // Amplitude correlation cutoff (#sig)

    // ===> TODO: candidates for virtual functions
    // void      FindNearestHits( ReadoutPlane* wp, const THaTrack* track,
    // 			       const Rvec_t& roads ) const;

    // UInt_t    MatchRoads( vector<Rvec_t>& roads,
    // 			  std::list<std::pair<Double_t,Rvec_t> >& combos_found,
    // 			  Rset_t& unique_found );
			  
    // Podd interface
    virtual Int_t   ReadDatabase( const TDatime& date );

    virtual TClass* GetPlaneClass() const;

    ClassDef(GEMTracker,0)   // Tree search track reconstruction for GEM trackers
  };

}

///////////////////////////////////////////////////////////////////////////////

#endif
