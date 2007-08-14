#ifndef ROOT_TreeSearch_Hitpattern
#define ROOT_TreeSearch_Hitpattern

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TBits.h"
#include <cstring>

namespace TreeSearch {

  // Modified version of TBits that clears the bits without deleting the
  // reserved memory
  class Bits : public TBits {
  public:
    Bits( UInt_t nbits = 8 ) : TBits(nbits) {}
    Bits( const TBits& orig ) : TBits(orig) {}
    Bits& operator=( const TBits& rhs ) {
      TBits::operator=(rhs); return *this;
    }
    virtual ~Bits() {}
    void FastClear() { memset(fAllBits,0,fNbytes); }
    virtual void Clear( Option_t* opt="" ) { if (fAllBits) FastClear(); }
    ClassDef(Bits,1)  // Bit container with memory-neutral Clear method
  };

  class PatternTree;
  class WirePlane;

  class Hitpattern {

  public:
    //    Hitpattern( const PatternTree& pt );
    Hitpattern( UInt_t depth, UInt_t nplanes, Double_t width );
    virtual ~Hitpattern();

    Bits*  GetRow( UInt_t i ) const { return i<fNplanes ? fPattern[i] : 0; }
    Bits** GetPattern() const { return fPattern; }
    Double_t GetWidth() const { return fWidth; }

    Int_t  SetPoints( WirePlane* hitsA, WirePlane* hitsB,
		      UInt_t plane );

    void   Clear( Option_t* opt="" );
    void   Print( Option_t* opt="" ) const;

//FIXME: add Draw() (=event display)

  protected:

    UInt_t   fDepth;    // Depth of the pattern tree
    UInt_t   fNplanes;  // Number of wire planes contained in the pattern
    Double_t fWidth;    // Physical width of region encompassing hitpattern (m)
    Bits**   fPattern;  // [fNPlanes] pattern at all fDepth resolutions

    ClassDef(Hitpattern,0)  // Wire chamber hitpattern at multiple resolutions
  };
}

///////////////////////////////////////////////////////////////////////////////

#endif
