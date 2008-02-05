#ifndef ROOT_TreeSearch_MWDC
#define ROOT_TreeSearch_MWDC

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaTrackingDetector.h"
#include "Projection.h"
#include "THaDetMap.h"
#include <vector>

class THaTrack;
class THaBenchmark;
class TClonesArray;
class THashTable;

using std::vector;

namespace TreeSearch {
  class WirePlane;
  class Projection;

  extern const Double_t kBig;

  // Types of projections supported by the MWDC. See also operator++ below
  enum EProjType { kUndefinedType = -1, kUPlane, kTypeBegin = kUPlane,
		   kVPlane, kXPlane, kTypeEnd };

  class MWDC : public THaTrackingDetector {
    friend class WirePlane;

  public:
    MWDC( const char* name, const char* description = "", 
	  THaApparatus* app = 0 );
    virtual ~MWDC();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    //    virtual Int_t   End(THaRunBase* run);
    virtual EStatus Init( const TDatime& date );
    virtual Int_t   CoarseTrack( TClonesArray& tracks );
    virtual Int_t   FineTrack( TClonesArray& tracks );
    virtual Int_t   DefineVariables( EMode mode = kDefine );
    virtual Int_t   FindVertices(TClonesArray&) {return 0;}
    virtual void    Print(const Option_t* opt) const;
    virtual void    SetDebug( Int_t level );
    void            EnableBenchmarks( Bool_t b = kTRUE );

    Double_t        GetRefTime( UInt_t i ) const 
    { return (i<(UInt_t)fRefMap->GetSize()) ? fRefTime[i] : kBig; }

    EProjType       NameToType( const char* name );

    //FIXME: for testing
    vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }
    vector<TreeSearch::Projection*>& GetListOfProjections() { return fProj; }

    // Analysis control flags. Set via database. Defaults are usually fine.
    enum {
      kDoTimeCut  = BIT(14), // Use TDC time cut in decoder (defined in planes)
      kPairsOnly  = BIT(15), // Accept only pairs of hits from plane pairs
      kMCdata     = BIT(16), // Assume input is Monte Carlo data
      kNoPartner  = BIT(17)  // Never partner wire planes
    };

  protected:

    vector<WirePlane*>   fPlanes;  // Wire planes
    vector<Projection*>  fProj;    // Plane projections

    THaDetMap*      fRefMap;    // Map of reference channels for VME readout
    vector<float>   fRefTime;   // [fRefMap->GetSize()] ref channel data

    THashTable*     fCrateMap;  // Map of MWDC DAQ modules

    virtual Int_t   ReadDatabase( const TDatime& date );

    // Helper functions for getting DAQ module parameters - used by Init
    UInt_t   LoadDAQmodel( THaDetMap::Module* m ) const;
    Double_t LoadDAQresolution( THaDetMap::Module* m ) const;
    UInt_t   GetDAQnchan( THaDetMap::Module* m ) const;

    ClassDef(MWDC,0)   // Tree search reconstruction of BigBite MWDCs
  };

  // Overloaded operator++ allows us to iterate over the EProjType enum.
  inline
  EProjType& operator++( EProjType& e )
  { switch(e) { 
    case kUndefinedType: return e=kUPlane;
    case kUPlane: return e=kVPlane;
    case kVPlane: return e=kXPlane;
    case kXPlane: return e=kTypeEnd;
    case kTypeEnd: default: return e=kUndefinedType;
    }
  }
  inline
  const EProjType operator++( EProjType e, int )
  { EProjType r(e); ++e; return r; }

}

///////////////////////////////////////////////////////////////////////////////

#endif
