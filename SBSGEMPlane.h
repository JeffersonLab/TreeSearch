//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 06-Feb-2012
//
#ifndef ROOT_SBS_GEMPlane
#define ROOT_SBS_GEMPlane

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBS::GEMPlane                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GEMPlane.h"
#include <map>
#include <sstream>

namespace SBS {

  class GEMPlane : public TreeSearch::GEMPlane {

  public:
    GEMPlane( const char* name, const char* description = "",
	      THaDetectorBase* parent = 0 );
    GEMPlane() {} // For ROOT RTTI
    virtual ~GEMPlane();

    // virtual void    Clear( Option_t* opt="" );
    // virtual Int_t   Decode( const THaEvData& );
    // virtual void    Print( Option_t* opt="" ) const;
    
    virtual Bool_t  Contains( Double_t x, Double_t y ) const;
    virtual Double_t GetModuleOffsets(Int_t module) {
      return mOffsets.find(module)->second;//mOffsets[module];
      //might not matter, but who knows
    };
  protected:

    // Podd interface
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t ReadGeometry( FILE* file, const TDatime& date,
				Bool_t required = kTRUE );
    //    virtual Int_t DefineVariables( EMode mode = kDefine );

    // Geometry, configuration
    Double_t  fD0;         // distance of Plane from spectrometer "Zero"
    Double_t  fDX;      // Coverage in X_T
    Double_t  fDY;      // Coverage in Y_T
    Int_t  nModule;  //number of modules in this Plane
    std::map<Int_t, Double_t> mOffsets; //map of module offsets

    ClassDef(GEMPlane,0)     // GEM readout plane in cylindrical coordinates
  };

///////////////////////////////////////////////////////////////////////////////

} // end namespace SBS


#endif
