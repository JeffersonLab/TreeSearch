#ifndef ROOT_THaTsMWDC
#define ROOT_THaTsMWDC

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaTsMWDC                                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaTrackingDetector.h"

class THaTsMWDCPlane;
class THaTsMWDCGroup;
class THaTsMWDCHit;
class THaTrack;
class THaBenchmark;
class TClonesArray;

class THaTsMWDC : public THaTrackingDetector {

public:
  THaTsMWDC( const char* name, const char* description = "", 
	     THaApparatus* a = NULL );
  virtual ~THaTsMWDC();

  virtual Int_t   Decode( const THaEvData& );
  virtual EStatus Init( const TDatime& date );
  virtual Int_t   CoarseTrack( TClonesArray& tracks );
  virtual Int_t   FineTrack( TClonesArray& tracks );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   FindVertices(TClonesArray&) {return 0;}
  virtual void    Print(const Option_t* opt) const;
  virtual void    SetDebug( Int_t level );
          void    EnableBenchmarks( Bool_t b = kTRUE );


protected:


  Int_t  ReadDatabase( const TDatime& date );
  void   Clear( Option_t* opt="" );
  Int_t  End(THaRunBase* run);


  ClassDef(THaTsMWDC,0)   // Tree search reconstruction of BigBite MWDCs
};

////////////////////////////////////////////////////////////////////////////////

#endif
