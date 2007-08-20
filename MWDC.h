#ifndef ROOT_TreeSearch_MWDC
#define ROOT_TreeSearch_MWDC

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::MWDC                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaTrackingDetector.h"
#include "Projection.h"
#include <vector>

class THaTrack;
class THaBenchmark;
class TClonesArray;

using std::vector;

namespace TreeSearch {
  class WirePlane;
  class Projection;

  extern const Double_t kBig;

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

    Double_t        GetUangle() const { return fUangle; }
    Double_t        GetVangle() const { return fVangle; }

    //FIXME: for testing
    vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }
    vector<TreeSearch::Projection*>& GetListOfProjections() { return fProj; }

    Projection::EProjType MakePlaneType( const char* name, Double_t& type );

  protected:

    vector<WirePlane*>   fPlanes;  // Wire planes
    vector<Projection*>  fProj;    // Plane projections

    Double_t  fUangle;   // angle of U wires wrt X-axis (-180 .. +180 deg)
    Double_t  fVangle;   // angle of V wires wrt X-axis (-180 .. +180 deg)

    virtual Int_t  ReadDatabase( const TDatime& date );

    ClassDef(MWDC,0)   // Tree search reconstruction of BigBite MWDCs
  };
}

///////////////////////////////////////////////////////////////////////////////

#endif
