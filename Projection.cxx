//*-- Author :    Ole Hansen, Jefferson Lab   16-Aug-2007

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Projection.h"
#include "Hitpattern.h"
#include "WirePlane.h"
#include "THaDetectorBase.h"
#include <iostream>

using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
  Projection::Projection( Int_t type, const char* name, Double_t angle,
			  THaDetectorBase* parent )
  //, const PatternTree* pt )
  : THaAnalysisObject( name, name ), fType(type), fDepth(0),
    fMaxSlope(0.0), fWidth(0.0), fZmin(kBig), fZmax(-kBig), 
    fHitpattern(NULL), fPatternTree(NULL), fDetector(parent)
    //    fPatternTree(pt)
{
  // Constructor

  // angle is the default value. It can be overridden via a database entry.
  SetAngle( angle );

  fTitle.Append(" projection");
}


//_____________________________________________________________________________
Projection::Projection( const Projection& orig )
  : fType(orig.fType), fPlanes(orig.fPlanes), fDepth(orig.fDepth),
    fMaxSlope(orig.fMaxSlope), fWidth(orig.fWidth), fZmin(orig.fZmin),
    fZmax(orig.fZmax), fSinAngle(orig.fSinAngle), fCosAngle(orig.fCosAngle),
    fHitpattern(NULL), fPatternTree(NULL), fDetector(orig.fDetector)
{
  // Copying

  if( orig.fHitpattern )
    fHitpattern = new Hitpattern(*orig.fHitpattern);
//   if( orig.fPatternTree )
//     fHitpattern = new PatternTree(*orig.fPatternTree);
}

//_____________________________________________________________________________
Projection& Projection::operator=( const Projection& rhs )
{
  // Assignment

  if( this != &rhs ) {
    fType     = rhs.fType;
    fDepth    = rhs.fDepth;
    fPlanes   = rhs.fPlanes;
    fMaxSlope = rhs.fMaxSlope;
    fWidth    = rhs.fWidth;
    fZmin     = rhs.fZmin;
    fZmax     = rhs.fZmax;
    fSinAngle = rhs.fSinAngle;
    fCosAngle = rhs.fCosAngle;
    fDetector = rhs.fDetector;
    delete fHitpattern;
    if( rhs.fHitpattern )
      fHitpattern = new Hitpattern(*rhs.fHitpattern);
    else
      fHitpattern = NULL;
//     if( rhs.fPatternTree )
//       fPatternTree = new PatternTree(*rhs.fPatternTree);
//     else
//       fPatternTree = NULL;
  }
  return *this;
}

//_____________________________________________________________________________
Projection::~Projection()
{
  // Destructor

  //  delete fPatternTree;
  delete fHitpattern;
}

//_____________________________________________________________________________
void Projection::AddPlane( WirePlane* wp )
{
  // Add wire plane to this projection

  if( !wp )
    return;

  const WirePlane* planes[] = { wp, wp->GetPartner(), 0 };
  const WirePlane** p = planes;
  while( *p ) {
    Double_t z = (*p)->GetZ();
    if( z < fZmin )
      fZmin = z;
    if( z > fZmax )
      fZmax = z;
    ++p;
  }
  
  fPlanes.push_back( wp );
}

//_____________________________________________________________________________
void Projection::SetAngle( Double_t angle )
{
  // Set wire angle (rad)
  
  fSinAngle = TMath::Sin( angle );
  fCosAngle = TMath::Cos( angle );
}

//_____________________________________________________________________________
void Projection::Reset()
{
  // Reset parameters, clear list of planes, delete Hitpattern
  
  fPlanes.clear();
  fMaxSlope = fWidth = 0.0;
  fZmin = kBig;
  fZmax = -kBig;
  delete fHitpattern; fHitpattern = NULL;
  //  delete fPatternTree; fPatternTree = NULL;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Projection::Init( const TDatime& date )
{
  // Initialize hitpattern

  //  static const char* const here = "Init";

  // Initialize ourselves. This calls our ReadDatabase().
  EStatus status = THaAnalysisObject::Init(date);
  if( status )
    return fStatus = status;

//   if( depth == 0 && !fPatternTree )
//     return -1;
//   Int_t nplanes = fPlanes.size();
//   if( nplanes == 0 )
//     return -2;

//   delete fHitpattern;
//   fHitpattern = NULL;
  
//   if( depth > 0 )
//     fHitpattern = new Hitpattern( depth, nplanes, fWidth );
  //Fixme:
//   else
//     fHitpattern = new Hitpattern( *fPatternTree );

  return fStatus = kOK;
}

//_____________________________________________________________________________
void Projection::MakePrefix()
{
  // Set up name prefix for global variables. 

  TString basename;
  if( fDetector ) {
    basename = fDetector->GetPrefix();
    basename.Chop();  // delete trailing dot
  } else
    Warning( Here("MakePrefix"), "No parent detector defined. "
	     "Using \"%s\".", fName.Data() );

  THaAnalysisObject::MakePrefix( basename.Data() );
}

//_____________________________________________________________________________
const char* Projection::GetDBFileName() const
{
  // Return database file name prefix. We want the same database file
  // as our parent detector.

  return fDetector ? fDetector->GetDBFileName() : GetPrefix();
}

//_____________________________________________________________________________
Int_t Projection::ReadDatabase( const TDatime& date )
{
  static const char* const here = "ReadDatabase";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Double_t angle = kBig;

  DBRequest request[] = {
    { "angle",        &angle,        kDouble, 0, 1 },
    { "maxslope",     &fMaxSlope,    kDouble, 0, 1, -1 },
    { "search_depth", &fDepth,       kUInt,   0, 0, -1 },
    { 0 }
  };

  Int_t err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return kInitError;

  if( fDepth > 20 ) {
    Error( Here(here), "Illegal search_depth = %u. Must be <= 20. "
	     "Fix database.", fDepth );
    return kInitError;
  }

  if( angle < kBig )
    SetAngle( angle*TMath::DegToRad() );

  if( fMaxSlope < 0.0 ) {
    Warning( Here(here), "Negative maxslope = %lf makes no sense. "
	     "Using |maxslope|.", fMaxSlope );
    fMaxSlope = -fMaxSlope;
  }

  return kOK;
}

//_____________________________________________________________________________
void Projection::Print( Option_t* opt ) const
{    
  // Print plane type info

  cout << "Projection:  " 
       << GetName()
       << " type=" << GetType()
       << " npl="  << fPlanes.size()
       << " depth=" << GetDepth()
       << " zsize=" << GetZsize()
       << " maxsl=" << GetMaxSlope()
       << " width=" << GetWidth()
       << " angle=" << GetAngle()*TMath::RadToDeg();
  cout << endl;

  //TODO: print plane details
}

//_____________________________________________________________________________

}  // end namespace TreeSearch

ClassImp(TreeSearch::Projection)

///////////////////////////////////////////////////////////////////////////////
