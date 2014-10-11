#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace TreeSearch;

#pragma link C++ class TreeSearch::MWDC+;
#pragma link C++ class TreeSearch::WirePlane+;
#pragma link C++ class TreeSearch::WireHit+;
#ifdef MCDATA
#pragma link C++ class TreeSearch::MCWireHit+;
#endif
#pragma link C++ class TreeSearch::HitpatternLR+;
#pragma link C++ class TreeSearch::ProjectionLR+;
#pragma link C++ class TreeSearch::TimeToDistConv+;
#pragma link C++ class TreeSearch::LinearTTD+;
#pragma link C++ class TreeSearch::TanhFitTTD+;

#endif
