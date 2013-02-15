//*-- Author :    Ole Hansen, Jefferson Lab   08-Feb-2013

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TreeSearch::Plane                                                        //
//                                                                          //
// Base class for a 1-dimensional tracker sensor plane.                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "Plane.h"
#include "Hit.h"
#include "Tracker.h"
#include "Projection.h"

#include "THaDetMap.h"
#include "THaEvData.h"
#include "TClonesArray.h"
#include "TH1.h"

#include "TClass.h"

#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;

namespace TreeSearch {

//_____________________________________________________________________________
Plane::Plane( const char* name, const char* description,
	      THaDetectorBase* parent )
  : THaSubDetector(name,description,parent), fPlaneNum(kMaxUInt),
    fDefinedNum(kMaxUInt), fType(kUndefinedType), fStart(0), 
    fPitch(0), fCoordOffset(0), fPartner(0), fProjection(0),
    fResolution(0), fMaxHits(kMaxUInt), fHits(0), fFitCoords(0),
    fHitMap(0)
{
  // Constructor

  static const char* const here = "Plane";

  assert( name && parent );

  fTracker = dynamic_cast<Tracker*>( GetMainDetector() );
  assert( fTracker );

  try {
    fFitCoords = new TClonesArray("TreeSearch::FitCoord", 20 );
  }
  catch( std::bad_alloc ) {
    Error( Here(here), "Out of memory constructing plane %s. "
	   "Call expert.", name );
    MakeZombie();
    return;
  }

#ifdef TESTCODE
  // No diagnostic histograms by default
  ResetBit(kDoHistos);
#endif
}

//_____________________________________________________________________________
Plane::~Plane()
{
  // Destructor.

  // Histograms in fHist should be automatically deleted by ROOT when 
  // the output file is closed

  delete fFitCoords;
  delete fHits;
}

//_____________________________________________________________________________
FitCoord* Plane::AddFitCoord( const FitCoord& coord )
{ 
  // Add given fit coordinate data to this plane's array of fit coordinates

  return
    new( (*fFitCoords)[GetNcoords()] ) FitCoord(coord);
}

//_____________________________________________________________________________
Int_t Plane::Begin( THaRunBase* /* run */ )
{
  // Book diagnostic histos if not already done

  assert( fIsInit && fNelem > 0 );

#ifdef TESTCODE
  if( TestBit(kDoHistos) ) {
    //TODO: check if exist
    string hname, htitle;
    hname   = fName + "_hitmap";
    htitle  = fName + " hitmap";
    fHitMap = new TH1F( hname.c_str(), htitle.c_str(), fNelem, 0, fNelem );
  }
#endif
  return 0;
}

//_____________________________________________________________________________
void Plane::Clear( Option_t* opt )
{    
  // Clear event-by-event data (hits)

  assert( fIsInit );

  fHits->Clear(opt);
  fFitCoords->Clear(opt);
}

//_____________________________________________________________________________
Int_t Plane::End( THaRunBase* /* run */ )
{
  // Write diagnostic histograms to file

  assert( fIsInit );

#ifdef TESTCODE
  if( TestBit(kDoHistos) ) {
    assert( fHitMap );
    fHitMap->Write();
  }
#endif
  return 0;
}

//_____________________________________________________________________________
Int_t Plane::ReadDatabaseCommon( const TDatime& date )
{
  // Read database

  static const char* const here = "ReadDatabase";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Set defaults
  TString plane_type;
  Int_t required = 0, do_histos = 0;
  fMaxHits   = kMaxUInt;

  // Read fOrigin (plane position) and fSize. fOrigin is the chamber
  // position relative to the Tracker reference frame - which in turn can be
  // offset as a whole. Thus, with respect to some absolute frame
  // (whatever it may be), each Plane is positioned at 
  // fOrigin(Tracker) + fOrigin(plane)
  Int_t status = ReadGeometry( file, date, kTRUE );

  if( status == kOK ) {
    vector<Int_t>* detmap = 0;
    try {
      // Oddly, putting this container on the stack may cause a stack overflow
      detmap = new vector<Int_t>;

      const DBRequest request[] = {
	{ "detmap",         detmap,           kIntV },
	{ "xp.res",         &fResolution,     kDouble,  0, 0, -1 },
	{ "maxhits",        &fMaxHits,        kUInt,    0, 1, -1 },
	{ "required",       &required,        kInt,     0, 1 },
	{ "type",           &plane_type,      kTString, 0, 1 },
	{ "description",    &fTitle,          kTString, 0, 1 },
	{ "partner",        &fPartnerName,    kTString, 0, 1 },
	{ "do_histos",      &do_histos,       kInt,     0, 1, -1 },
	{ 0 }
      };

      status = LoadDB( file, date, request, fPrefix );

      if( status == kOK ) {
	UInt_t flags = TestBit(kHaveRefChans) ? THaDetMap::kFillRefChan : 0;
	// Parse the detector map of the data channels
	if( FillDetMap( *detmap, flags, here ) == 0 )
	  status = kInitError;
      }
      delete detmap;
    }
    // Catch exceptions here so that we can close the file and clean up
    catch(...) {
      delete detmap;
      fclose(file);
      throw;
    }
  }

  // Finished reading the database
  fclose(file);
  if( status != kOK )
    return status;

  // Retrieve DAQ module parameters for our crateslots
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fDetMap->GetModule(imod);
    fTracker->LoadDAQmodel(d);
    fTracker->LoadDAQresolution(d);
    if( TestBit(kTDCReadout) )
      d->MakeTDC();
    else
      d->MakeADC();
    UInt_t nchan = fTracker->GetDAQnchan(d);
    if( d->hi >= nchan ) {
      Error( Here(here), "Detector map channel out of range for module "
          "cr/sl/lo/hi = %u/%u/%u/%u. Must be < %u. Fix database.",
          d->crate, d->slot, d->lo, d->hi, nchan );
      return kInitError;
    }
  }

  // Determine the type of this plane. If the optional plane type variable is
  // not given, use the first character of the plane name.
  TString name = plane_type.IsNull() ? fName[0] : plane_type[0];
  fType = Projection::NameToType( name );
  if( fType == kUndefinedType ) {
    TString names;
    for( EProjType i = kTypeBegin; i < kTypeEnd; ++i ) {
      names += kProjParam[i].name;
      if( i+1 != kTypeEnd ) 
        names += " ";
    }
    Error( Here(here), "Unsupported plane type \"%s\". Must be one of "
        "%s. Fix database.", name.Data(), names.Data() );
    return kInitError;
  }

  if( required )
    SetBit( kIsRequired );

#ifdef TESTCODE
  if( do_histos )
    SetBit( kDoHistos );
#endif

  // Update fCoordOffset based on the orientation of fProjection (if defined)
  UpdateOffset();

  return kOK;
}
//_____________________________________________________________________________
void Plane::SetPartner( Plane* p )
{    
  // Partner this plane with plane 'p'. The meaning of "partner" depends on the
  // type of detector:
  // For wire chambers:
  //  -> a plane immediately next to this one with staggered wires
  // For 2-D strip readouts (GEMs):
  //  -> the plane that represents the other readout coordinate

  if( p )
    p->fPartner = this;
  else if( fPartner ) {
    assert( this == fPartner->fPartner );
    fPartner->fPartner = 0;
  }
  fPartner = p;

  return;
}

//_____________________________________________________________________________
void Plane::UpdateOffset()
{
  // Update fCoordOffset based on the orientation of the associated projection
  // and this plane's fOrigin.

  if( fProjection ) {
    // Calculate the offset to be applied to the sensor positions if the
    // chamber has a position offset. This will translate the hit positions
    // into the Tracker frame, which is the reference frame for tracking
    TVector2 off( fOrigin.X(), fOrigin.Y() );
    fCoordOffset = off * fProjection->GetAxis();
  } else
    fCoordOffset = 0;
}

//_____________________________________________________________________________
void Plane::Print( Option_t* ) const
{    
  // Print plane info

  cout << IsA()->GetName() << ":  #" << GetPlaneNum() << " "
       << GetName() << " \"" << GetTitle() << "\"\t"
       << fNelem << " channels\t"
       << "z = " << GetZ();
  if( fPartner ) {
    cout << "\t partner = " 
	 << fPartner->GetName();
    //	 << fPartner->GetTitle();
  }
  cout << endl;
}

//_____________________________________________________________________________
// Int_t Plane::Compare( const TObject* obj ) const 
// {
//   // Used to sort planes in a TCollection/TList by z-position

//   // Fail if comparing to some other type of object
//   assert( obj && IsA() == obj->IsA() );

//   if( obj == this )
//     return 0;

//   const Plane* other = static_cast<const Plane*>( obj );
//   if( GetZ() < other->GetZ() ) return -1;
//   if( GetZ() > other->GetZ() ) return  1;
//   return 0;
// }

//_____________________________________________________________________________
Double_t Plane::GetMaxSlope() const
{ 
  return fProjection ? fProjection->GetMaxSlope() : kBig;
}

//_____________________________________________________________________________

}

ClassImp(TreeSearch::Plane)

///////////////////////////////////////////////////////////////////////////////

