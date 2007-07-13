#ifndef ROOT_TreeSearch_MWDC
#define ROOT_TreeSearch_MWDC

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaTrackingDetector.h"

class THaTrack;
class THaBenchmark;
class TClonesArray;
class TList;

namespace TreeSearch {
  class WirePlane;
  class PlaneGroup;
  class Hit;

  class MWDC : public THaTrackingDetector {

  public:
    MWDC( const char* name, const char* description = "", 
	  THaApparatus* app = NULL );
    virtual ~MWDC();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual Int_t   End(THaRunBase* run);
    virtual EStatus Init( const TDatime& date );
    virtual Int_t   CoarseTrack( TClonesArray& tracks );
    virtual Int_t   FineTrack( TClonesArray& tracks );
    virtual Int_t   DefineVariables( EMode mode = kDefine );
    virtual Int_t   FindVertices(TClonesArray&) {return 0;}
    virtual void    Print(const Option_t* opt) const;
    virtual void    SetDebug( Int_t level );
    void            EnableBenchmarks( Bool_t b = kTRUE );

    TList*          GetListOfPlanes() { return fPlanes; }

  protected:

  

    virtual Int_t  ReadDatabase( const TDatime& date );


    TList*  fPlanes;   // Wire planes

    ClassDef(MWDC,0)   // Tree search reconstruction of BigBite MWDCs
  };
}

///////////////////////////////////////////////////////////////////////////////

#endif
