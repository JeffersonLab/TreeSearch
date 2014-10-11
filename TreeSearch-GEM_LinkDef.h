#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace TreeSearch;

#pragma link C++ class TreeSearch::GEMTracker+;
#pragma link C++ class TreeSearch::GEMPlane+;
#pragma link C++ class TreeSearch::GEMHit+;
#ifdef MCDATA
#pragma link C++ class TreeSearch::MCGEMHit+;
#endif

#endif
