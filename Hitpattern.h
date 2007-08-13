#ifndef ROOT_TreeSearch_Hitpattern
#define ROOT_TreeSearch_Hitpattern

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Hitpattern                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//FIXME: add Draw() (=event display)

#include "TBits.h"
#include <cstring>

class TSeqCollection;

namespace TreeSearch {

  // Modified version of TBits that clears the bits without deleting the
  // reserved memory
  class Bits : public TBits {
  public:
    Bits( UInt_t nbits ) : TBits(nbits) {}
    virtual ~Bits() {}
    virtual void Clear( Option_t* opt="" ) { 
      if (fAllBits) memset(fAllBits,0,fNbytes);
    }
    ClassDef(Bits,1)  // Bit container with memory-neutral Clear method
  };

  class PatternTree;

  class Hitpattern {

  public:
    //    Hitpattern( const PatternTree& pt );
    Hitpattern( UInt_t depth, UInt_t nplanes, Double_t width );
    virtual ~Hitpattern();

    Int_t SetPoints( TSeqCollection* hitsA, TSeqCollection* hitsB,
		     UInt_t plane );

    void Clear( Option_t* opt="" );
    void Print( Option_t* opt="" ) const;

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
