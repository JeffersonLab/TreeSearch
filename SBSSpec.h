#ifndef ROOT_SBS_SBSSpec
#define ROOT_SBS_SBSSpec

//////////////////////////////////////////////////////////////////////////
//
// SBS::SBSSpec
//
// SBS spectrometer class containing GEM trackers
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h"
#include "TVector2.h"

#ifndef KBIG
#define KBIG 1e38
#endif

namespace SBS {

  // Helper class for holding additional track data not in THaTrack
  class SBSTrackInfo : public TObject {
  public:
    SBSTrackInfo()
      : fSector(-1), fX_T(KBIG), fY_T(KBIG), fXdir_T(KBIG), fYdir_T(KBIG)
      , f3dMatchBits(0), fN3dMatch(0)
#ifdef MCDATA
      , fMCHitBits(0), fNMCHits(0)
#endif
    {}
    SBSTrackInfo( Int_t sector, const TVector3& point,
		  const TVector3& direction)
      : TObject(), fSector(sector), fX_T(point.X()), fY_T(point.Y()),
      fXdir_T(direction.X()/direction.Z()), fYdir_T(direction.Y()/direction.Z())
      , f3dMatchBits(0), fN3dMatch(0)
#ifdef MCDATA
      , fMCHitBits(0), fNMCHits(0)
#endif
    {}

    Int_t    Get3DMatchBits() const     { return f3dMatchBits; }
    Int_t    GetN3DMatch()    const     { return fN3dMatch; }
    void     Set3DMatchBits( Int_t bits ) { f3dMatchBits = bits; }
    void     SetN3DMatch( Int_t nmatch )  { fN3dMatch = nmatch; }
   
#ifdef MCDATA
    Int_t    GetMCHitBits() const       { return fMCHitBits; }
    Int_t    GetNMCHits()   const       { return fNMCHits; }
    void     SetMCHitBits( Int_t bits ) { fMCHitBits = bits; }
    void     SetNMCHits( Int_t nhits )  { fNMCHits = nhits; }
#endif

  protected:

    // Track coordinates in the SBS transport coordinates.
    Int_t     fSector;     // Sector where this track was reconstructed
    Double_t  fX_T;        // X in transport coordinates 
    Double_t  fY_T;        // Y in transport coordinates 
    Double_t  fXdir_T;     // X direction in transport coordinates 
    Double_t  fYdir_T;     // Y direction transport coordinates 
    Int_t     f3dMatchBits; // Bitpattern of 3D hit ampcorr matches
    Int_t     fN3dMatch;   // number of 3D hit ampcorr matches
    
#ifdef MCDATA
    // Diagnostic info derived from Monte Carlo truth data
    Int_t     fMCHitBits;  // Bitpattern of plane #s w/MC hit used by this track
    Int_t     fNMCHits;    // Number of MC hits in this track (# bits set in fMCHitBits)
#endif

    ClassDef(SBSTrackInfo,0) // SBS track coordinates
  };

  class SBSSpec : public THaSpectrometer {

  public:
    SBSSpec( const char* name, const char* description, UInt_t nsectors );
    virtual ~SBSSpec();

    virtual void      Clear( Option_t* opt="" );
    virtual EStatus   Init( const TDatime& run_time );
    virtual Int_t     FindVertices( TClonesArray& tracks );
    virtual Int_t     TrackCalc();

    TClonesArray*     GetTrackInfo() { return fSBSTrackInfo; }

  protected:

    TClonesArray*     fSBSTrackInfo;   // SBS-specific per-track info

    virtual Int_t     DefineVariables( EMode mode = kDefine );
    virtual Int_t     ReadRunDatabase( const TDatime& date );

    ClassDef(SBSSpec,0) // SBS spectrometer
  };

} // end namespace SBS

#endif

