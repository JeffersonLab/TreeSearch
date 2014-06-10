//*-- Author :    Ole Hansen, Jefferson Lab   06-Feb-2012

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// SoLID::GEMPlane                                                          //
//                                                                          //
// A 1-dimensional readout plane of a GEM chamber.                          //
// Use two of these objects for a 2-d readout plane.                        //
//                                                                          //
// This version specifies the geometry in cylindrical coordinates,          //
// as needed for SoLID                                                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "SoLIDGEMPlane.h"
#include "SoLIDGEMTracker.h"

using namespace std;

namespace SoLID {

//_____________________________________________________________________________
GEMPlane::GEMPlane( const char* name, const char* description,
		    THaDetectorBase* parent )
  : TreeSearch::GEMPlane(name,description,parent)
{
  // Constructor

  assert( dynamic_cast<SoLID::GEMTracker*>(fTracker) );

}

//_____________________________________________________________________________
GEMPlane::~GEMPlane()
{
  // Destructor.

  // if( fIsSetup )
  //   RemoveVariables();

}

//_____________________________________________________________________________
// void GEMPlane::Clear( Option_t* opt )
// {
//   // Clear event-by-event data (hits)

//   TreeSearch::GEMPlane::Clear(opt);

// }

// //_____________________________________________________________________________
// Int_t GEMPlane::Decode( const THaEvData& evData )
// {
//   // Decode and pre-process raw data.
//   //
//   // This routine determines hit positions and integrated ADC amplitudes from
//   // the raw GEM strip ADC data

//   return TreeSearch::GEMPlane::Decode( evData );
// }

//_____________________________________________________________________________
// Int_t GEMPlane::DefineVariables( EMode mode )
// {
//   // initialize global variables

//   if( mode == kDefine && fIsSetup ) return kOK;
//   fIsSetup = ( mode == kDefine );

//   return kOK;
// }

//_____________________________________________________________________________
// Int_t GEMPlane::ReadDatabase( const TDatime& date )
// {
//   // Read database

//   //  static const char* const here = "ReadDatabase";

//   // Read the database for the base class, quit if error
//   Int_t status = TreeSearch::GEMPlane::ReadDatabase(date);
//   if( status != kOK )
//     return status;

//   fIsInit = true;
//   return kOK;
// }

//___________________________________________________________________________
Bool_t GEMPlane::Contains( Double_t x, Double_t y ) const
{
  // Check if the given point is within the active area of this plane.
  // x/y are plane coordinates, i.e. relative to fOrigin of this GEMPlane.
  //
  // The active area of this type of GEM plane is a section of a ring.

  // Convert the point to lab coordinates (without the sector's phi rotation)
  // NB: x/y are in the Tracker frame, so they already include fOrigin
  const TVector3& tracker_origin = fTracker->GetOrigin();
  Double_t xs = x + tracker_origin.X();
  Double_t ys = y + tracker_origin.Y();
  Double_t r2 = xs*xs + ys*ys;
  if( r2 <= fRmin2 or r2 >= fRmax2 ){
    return false;
  }
  // Get the point's phi angle between -pi and pi and check if it falls
  // in the angular range of this ring segment. Recall:
  // - By definition, the x-axis of the plane frame has phi = 0
  // - Phi angles have already been normalized to the range -pi to pi
  // - fPhiMin < fPhiMax is assured by construction
  Double_t phi = TMath::ATan2( ys, xs );
  return ( fPhiMin < phi and phi < fPhiMax );
}

//_____________________________________________________________________________
Int_t GEMPlane::ReadGeometry( FILE* file, const TDatime& date,
			      Bool_t /* required */ )
{
  // Read basic geometry for a SoLID GEM readout plane defined in cylindrical
  // coordinates.
  //
  // Read: rmin, rmax, dphi, phioff, z, dz
  // Units: [m] or [deg]
  //
  // rmin and rmax are the inner and outer radii of the ring delimiting the
  // plane.
  //
  // dphi is the full angular width of the active area.
  // The central phi of the plane is defined by the enclosing Tracker detector.
  // The plane does not need to know it and can assume phi = 0 without loss
  // of generality.
  //
  // phioff is an optional offset to be applied to the central phi defining the
  // active area. The angular limits are -dphi/2 + phioff and +dphi/2 + phioff.
  // The phi offset is meant to represent a shift of the aperture, like a
  // collimator, not a real rotation of the readout plane/strip angles.
  //
  // NOTE: phioff DOES NOT rotate the strips. The strip angles are defined
  // by the projection that this plane belongs to. They are relative to
  // the x-axis (phi = 0).
  //
  // z is the z-position of this plane in absolute coordinates IN THE LAB FRAME
  //
  // dz is the plane thickness. It is currently unused.

  static const char* const here = "ReadGeometry";

  Double_t dphi, phioff = 0.0, z, dz = 1e-3;
  Int_t gbl = GetDBSearchLevel(fPrefix);
  DBRequest request[] = {
    // "position"-like parameters
    { "rmin",   &fRmin,      kDouble, 0, 0, 0, "inner radius [m]" },
    { "rmax",   &fRmax,      kDouble, 0, 0, 0, "outer radius [m]" },
    { "z",      &z,          kDouble, 0, 0, 0,
                                           "z-position of readout plane [m]" },
    // "size"-like parameters
    { "dphi",   &dphi,       kDouble, 0, 0, gbl,
                                   "full angular width of active area [deg]" },
    { "phioff", &phioff,     kDouble, 0, 1, 0,
      "offset rotation of sector wrt central phi of sector (optional) [deg]" },
    { "dz",     &dz,         kDouble, 0, 1, gbl,  "extent in z (unused) [m]" },
    { 0 }
  };
  Int_t err = LoadDB( file, date, request, fPrefix );
  if( err )
    return err;

  // Sanity checks
  if( fRmin < 0 or fRmax < 0 ) {
    Error( Here(here), "rmin and rmax must be positive. Read %.1lf and %.1lf. "
	   "Fix database.", fRmin, fRmax );
    return kInitError;
  }
  if( fRmin >= fRmax ) {
    Error( Here(here), "rmin = %.1lf must be less than rmax = %.1lf. "
	   "Fix database.", fRmin, fRmax );
    return kInitError;
  }
  // Limit the opening angle to simplify subsequent calculations
  if( dphi <= 0 or dphi >= 90. ) {
    Error( Here(here), "Illegal value for dphi = %.1lf. Must be > 0 and <= "
	   "90 degrees. Fix database.", dphi );
    return kInitError;
  }
  // Ensure that the offset is sane, i.e. don't let any angles get larger than
  // +/- 90 degrees
  Double_t max_off = 90.0 - 0.5*dphi;
  if( TMath::Abs(phioff) >= max_off ) {
    Error( Here(here), "Illegal value for phioff = %.1lf. Must be between "
	   "%.1lf and %.1lf degrees. Fix database.",
	   phioff, -max_off, max_off );
    return kInitError;
  }
  if( dz <= 0 ) {
    Error( Here(here), "Illegal value for dz = %lf. Must be positive. "
	   "Fix database.", dz );
    return kInitError;
  }

  // Keep angles in rad
  dphi   *= TMath::DegToRad();
  phioff *= TMath::DegToRad();

  // Compute derived quantities
  Double_t phi2 = 0.5*dphi;
  fRmin2  = fRmin*fRmin;
  fRmax2  = fRmax*fRmax;
  fPhiMin = -phi2 + phioff;
  fPhiMax =  phi2 + phioff;
  assert( fPhiMin < fPhiMax );
  assert( fPhiMin > -TMath::PiOver2() );
  assert( fPhiMax < TMath::PiOver2() );

  // Define the origin of the plane in the same way as in libsolgem: It is
  // the center of the bounding box of the ring segment, calculated WITHOUT
  // any phi offset rotation.
  Double_t xmin = fRmin * TMath::Cos(phi2), xmax = fRmax;
  fOrigin.SetXYZ( 0.5*(xmin+xmax), 0.0, z );

  // If there is a phi offset, rotate the origin by the offset angle, so it
  // stays in the center of the actual (rotated) active area.
  // This origin is still in the lab frame. SoLID::GEMTracker::Init later
  // subtracts the origin of the first plane (smallest z) from all the planes.
  Double_t xs = -kBig, ys = -kBig;
  if( phioff != 0.0 ) {
    fOrigin.RotateZ(phioff);

    // fSize is not used by SoLID::GEMTrackers; Contains(x,y) uses the actual
    // shape of the active area since it is not rectangular.
    // Calculate fSize anyway for completeness, defined as the smallest symmetric
    // bounding box of the ring section around fOrigin AFTER rotation by phioff.
    TVector2 L( TMath::Cos(fPhiMin), TMath::Sin(fPhiMin) );
    TVector2 R( TMath::Cos(fPhiMax), TMath::Sin(fPhiMax) );
    TVector2 C[4] = { L, R, L, R };  // = TL, TR, BL, BR
    TVector2 org = fOrigin.XYvector();
    if ( fPhiMin < 0 && fPhiMax > 0 ) C[0].Set( 1.0, C[0].Y() );
    for( Int_t i = 0; i < 4; ++i ) {
      C[i] *= ( i<2 ) ? fRmax : fRmin;
      C[i] -= org;
      Double_t xa = TMath::Abs(C[i].X());
      Double_t ya = TMath::Abs(C[i].Y());
      if( xa > xs ) xs = xa;
      if( ya > ys ) ys = ya;
    }
  } else {
    // Simplified case of no phi offset where everything is symmetrical
    xs = 0.5 * (xmax-xmin);
    ys = fRmax * TMath::Sin(phi2);
  }
  fSize[0] = xs;  // Note: defined as half-width
  fSize[1] = ys;
  fSize[2] = dz;

  return kOK;
}

//_____________________________________________________________________________

}

ClassImp(SoLID::GEMPlane)

///////////////////////////////////////////////////////////////////////////////

