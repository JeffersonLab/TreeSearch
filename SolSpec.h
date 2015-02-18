#ifndef ROOT_SoLID_SolSpec
#define ROOT_SoLID_SolSpec

//////////////////////////////////////////////////////////////////////////
//
// SoLID::SolSpec
//
// SoLID spectrometer class containing GEM trackers
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h"
#include "TVector2.h"

#ifndef KBIG
#define KBIG 1e38
#endif

namespace SoLID {

  // Helper class for holding additional track data not in THaTrack
  class SolTrackInfo : public TObject {
  public:
    SolTrackInfo()
      : fSector(-1), fRtransv(KBIG), fTheta(KBIG), fPhi(KBIG),
	fPhiRot(KBIG), fThetaDir(KBIG), fPhiDir(KBIG), fPhiDirRot(KBIG)
#ifdef MCDATA
      , fMCHitBits(0), fNMCHits(0)
#endif
    {}
    SolTrackInfo( Int_t sector, const TVector3& point,
		  const TVector3& direction, Double_t phi0 )
      : TObject(), fSector(sector), fRtransv(point.Perp()),
	fTheta(point.Theta()), fPhi(point.Phi()),
	fPhiRot(TVector2::Phi_mpi_pi(fPhi-phi0)),
	fThetaDir(direction.Theta()), fPhiDir(direction.Phi()),
	fPhiDirRot(TVector2::Phi_mpi_pi(fPhiDir-phi0))
#ifdef MCDATA
      , fMCHitBits(0), fNMCHits(0)
#endif
    {}

#ifdef MCDATA
    Int_t    GetMCHitBits() const       { return fMCHitBits; }
    Int_t    GetNMCHits()   const       { return fNMCHits; }
    void     SetMCHitBits( Int_t bits ) { fMCHitBits = bits; }
    void     SetNMCHits( Int_t nhits )  { fNMCHits = nhits; }
#endif

  protected:

    // Track coordinates in the cylindrical system SoLID uses.
    // Aziumthal angles are measured wrt to the x-axis, which is the centerline
    // of sector 0. Positive angles mean clockwise rotation when looking
    // downstream along the beam.
    Int_t     fSector;     // Sector where this track was reconstructed
    Double_t  fRtransv;    // Transverse distance of track point from beam [m]
    Double_t  fTheta;      // Polar angle of track point wrt target center [rad]
    Double_t  fPhi;        // Azimuth of track point [rad]
    Double_t  fPhiRot;     // Azimuth of track point wrt to sector center [rad]
    Double_t  fThetaDir;   // Polar angle of track direction [rad]
    Double_t  fPhiDir;     // Azimuth of track direction [rad]
    Double_t  fPhiDirRot;  // Azimuth of track dir wrt to sector center [rad]

#ifdef MCDATA
    // Diagnostic info derived from Monte Carlo truth data
    Int_t     fMCHitBits;  // Bitpattern of plane #s w/MC hit used by this track
    Int_t     fNMCHits;    // Number of MC hits in this track (# bits set in fMCHitBits)
#endif

    ClassDef(SolTrackInfo,0) // SoLID track coordinates
  };

  class SolSpec : public THaSpectrometer {

  public:
    SolSpec( const char* name, const char* description, UInt_t nsectors );
    virtual ~SolSpec();

    virtual void      Clear( Option_t* opt="" );
    virtual EStatus   Init( const TDatime& run_time );
    virtual Int_t     FindVertices( TClonesArray& tracks );
    virtual Int_t     TrackCalc();

    TClonesArray*     GetTrackInfo() { return fSolTrackInfo; }

  protected:

    TClonesArray*     fSolTrackInfo;   // SoLID-specific per-track info

    virtual Int_t     DefineVariables( EMode mode = kDefine );
    virtual Int_t     ReadRunDatabase( const TDatime& date );

    ClassDef(SolSpec,0) // SoLID spectrometer
  };

} // end namespace SoLID

#endif

