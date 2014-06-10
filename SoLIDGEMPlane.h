//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 06-Feb-2012
//
#ifndef ROOT_SoLID_GEMPlane
#define ROOT_SoLID_GEMPlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SoLID::GEMPlane                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GEMPlane.h"

namespace SoLID {

  class GEMPlane : public TreeSearch::GEMPlane {

  public:
    GEMPlane( const char* name, const char* description = "",
	      THaDetectorBase* parent = 0 );
    GEMPlane() {} // For ROOT RTTI
    virtual ~GEMPlane();

    //    virtual void    Clear( Option_t* opt="" );
    //    virtual Int_t   Decode( const THaEvData& );
    //    virtual void    Print( Option_t* opt="" ) const;

    virtual Bool_t  Contains( Double_t x, Double_t y ) const;

  protected:

    // Podd interface
    //    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t ReadGeometry( FILE* file, const TDatime& date,
				Bool_t required = kTRUE );
    //    virtual Int_t DefineVariables( EMode mode = kDefine );

    // Geometry, configuration
    Double_t  fRmin;         // Inner radius (m)
    Double_t  fRmax;         // Outer radius (m)
    Double_t  fPhiMin;       // Lower phi angle limit (rad)
    Double_t  fPhiMax;       // Upper phi angle limit (rad)

    // Calculated parameters, for efficiency
    Double_t  fRmin2;        // Inner radius squared (m^2)
    Double_t  fRmax2;        // Outer radius squared (m^2)

    ClassDef(GEMPlane,0)     // GEM readout plane in cylindrical coordinates
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace SoLID


#endif
