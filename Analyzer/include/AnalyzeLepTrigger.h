#ifndef AnalyzeLepTrigger_h
#define AnalyzeLepTrigger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeLepTrigger 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    
    AnalyzeLepTrigger();
    ~AnalyzeLepTrigger(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
    bool passTriggerGeneral( std::vector<std::string>& myTriggerVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass );
    void printTriggerList( const std::vector<std::string>& TriggerNames );
    void doesTriggerExist( std::vector<std::string>& myTrigTestVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass );
    
};

#endif