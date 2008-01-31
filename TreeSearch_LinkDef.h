#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace TreeSearch;

#pragma link C++ class TreeSearch::MWDC+;
#pragma link C++ class TreeSearch::WirePlane+;
#pragma link C++ class TreeSearch::Hit+;
#pragma link C++ class TreeSearch::MCHit+;
#pragma link C++ class TreeSearch::HitPairIter+;
#pragma link C++ class TreeSearch::TimeToDistConv+;
#pragma link C++ class TreeSearch::LinearTTD+;
#pragma link C++ class TreeSearch::Bits+;
#pragma link C++ class TreeSearch::Hitpattern+;
#pragma link C++ class TreeSearch::Projection+;
#pragma link C++ class TreeSearch::PatternTree+;
#pragma link C++ class TreeSearch::PatternGenerator+;
#pragma link C++ class TreeSearch::PatternGenerator::Statistics_t+;
#pragma link C++ class TreeSearch::TreeWalk+;
#pragma link C++ class TreeSearch::NodeDescriptor+;
#pragma link C++ class TreeSearch::TreeParam_t+;

#pragma link C++ class std::pair<TObject*,TObject*>;
#pragma link C++ class std::vector<TreeSearch::WirePlane*>;
#pragma link C++ class std::vector<TreeSearch::Projection*>;

#endif
