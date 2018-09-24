//*-- Author :    Ole Hansen, Jefferson Lab   06-Feb-2012

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// SBS::GEMPlane                                                          //
//                                                                          //
// A 1-dimensional readout plane of a GEM chamber.                          //
// Use two of these objects for a 2-d readout plane.                        //
//                                                                          //
// This version specifies the geometry in cylindrical coordinates,          //
// as needed for SBS                                                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "SBSGEMPlane.h"
#include "SBSGEMTracker.h"




using namespace std;

namespace SBS {

//_____________________________________________________________________________
GEMPlane::GEMPlane( const char* name, const char* description,
		    THaDetectorBase* parent )
  : TreeSearch::GEMPlane(name,description,parent)
{
  // Constructor

  assert( dynamic_cast<SBS::GEMTracker*>(fTracker) );

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
  Int_t GEMPlane::ReadDatabase( const TDatime& date )
 {
   // Read database

   //  static const char* const here = "ReadDatabase";

   // Read the database for the base class, quit if error
   Int_t status = TreeSearch::GEMPlane::ReadDatabase(date);
   if( status != kOK )
     return status;

   fIsInit = true;
   return kOK;
   //add subGEMPlane here
   //newModule = new SBS::GEMModule(name, description, parent);
   //fModules.push_back(newModule); //needs only Geo information, can just be an struct/map
 }

//___________________________________________________________________________
Bool_t GEMPlane::Contains( Double_t x, Double_t y ) const
{
  return ( abs(x) < fDX/2.0 and abs(y) < fDY/2.0 );
}

//_____________________________________________________________________________
Int_t GEMPlane::ReadGeometry( FILE* file, const TDatime& date,
			      Bool_t /* required */ )
{
  // A GEM plane is basically just a box, characterized by an extension in x (dx) and y (dy)
  // and by five parameters for its location: 
  // dmag, the distance of the magnet to the target, 
  // which is taken as a refence for the spectrometer coordinates 
  // d0, the distance of the box center to the spectrometer reference, 
  // X_offset, the shift in X_transport of the plane wrt the central ray of the spectrometer, 
  // and two angles: one of horizontal rotation (thetaH), which translates the SBS angle,
  // one of vertical rotation (thetaV), translating the "bending" of the SBS wrt the xOz plane.
  // The "Zero" of the spectrometer is defined as: ( -dmag*sin(thetaH), 0, dmag*cos(thetaH) ); 
  //
  // Warning: the x direction of the box is actually in the -y direction of the lab, 
  //          the y direction of the box being in the + x direction of the lab 
  // NB: the plane also uses xOffset, the translation in X (in transport coordinates) of the chamber.
  // This variable is read by SBSGEMTracker.
  
  // In this routine, only D0 is read. It is the parameter which will allow to order the planes in the tracker.
  // The other parameters are read in class SBSGEMTracker 
  
  // the box is centered on the "central ray" which is then rotated by thetaH and thetaV.
  // To illustrate this, let's take a dumb example: 
  // if thetaH = thetaV = 0; then the center of the box would be located on the Z axis, 
  // and the box coverage on the x (resp y) direction would be from -dx/2 to +dx/2 
  // (resp -dy/2 to +dy/2)
  
  // The rotation to place the box in the lab is made the following way.
  // first the box is rotated by thetaH wrt y direction (i.e. x, z, modified, y conserved)
  // then, the box is rotated by thetaV wrt x' direction obtained at the previous step 
  // (i.e. y, z' modified, x' conserved)
  // lastly, the frame is rotated by 90 degrees wrt z" direction (i.e. x', y', modified, z conserved)
  
  // The box z location is to be understood as the box location on the z" axis 
  // obtained with the rotation defined before, and is made wrt x, y = 0, 0.
  
  
  static const char* const here = "ReadGeometry";

  Double_t depth;
  //Int_t gbl = GetDBSearchLevel(fPrefix);
  DBRequest request[] = {
    {"d0",          &fD0,           kDouble, 0, 1},
    {"dx",          &fDX,           kDouble, 0, 1},
    {"dy",          &fDY,           kDouble, 0, 1},
    {"dz",          &depth,         kDouble, 0, 1},
    {"nmodule",     &nModule,       kInt,    0, 1},
    { 0 }
  };

  Int_t err = LoadDB( file, date, request, fPrefix );
  if( err )
    return err;
  

  for (int j=1; j<nModule+1; j++){
    mOffsets[j]=0;
    Double_t offset_temp;
    DBRequest module_request[] = {
      {"offset",          &offset_temp,           kDouble, 0, 1},
      { 0 }
    };
    ostringstream module_prefix(fPrefix, ios_base::ate);
    module_prefix<<j<<".";
    cout<<module_prefix.str().c_str()<<endl;
    int err = LoadDB( file, date, module_request, module_prefix.str().c_str());
    if( err ) exit(2);
    mOffsets[j-1] = offset_temp;
    //cout<<offset_temp<<" "<<mOffsets[j]<<endl;getchar();     
  }

  // Sanity checks
  if( fDX < 0 or fDY < 0 ) {
    Error( Here(here), "dx and dy must be positive. Read %.1lf and %.1lf. "
  	   "Fix database.", fDX, fDY );
    return kInitError;
  }
  if( fD0 < 0 ) {
    Error( Here(here), "d0  must be positive. Read %.1lf. "
  	   "Fix database.", fD0 );
    return kInitError;
  }
  if( depth <= 0 ) {
    Error( Here(here), "Illegal value for dz = %lf. Must be positive. "
  	   "Fix database.", depth );
    return kInitError;
  }
  
  // double x0, y0, z0;
  // x0 = -dmag*sin(fThetaH)*1.0e3;;
  
  fOrigin.SetXYZ(0.0, 0.0, fD0);
  
  fSize[0] = fDX;
  fSize[1] = fDY;
  fSize[2] = depth;
  cout<<"SBS:GEMPlane readgeo done"<<endl;
  return kOK;
}

//_____________________________________________________________________________

}

ClassImp(SBS::GEMPlane)

///////////////////////////////////////////////////////////////////////////////

