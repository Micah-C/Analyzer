#ifndef AnalyzeTopTagger_h
#define AnalyzeTopTagger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeTopTagger
{
public:
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

   AnalyzeTopTagger();
   ~AnalyzeTopTagger(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, std::string filetag = "", bool isQuiet = false);
   void     InitHistos();
   void     WriteHistos();

};

#endif
