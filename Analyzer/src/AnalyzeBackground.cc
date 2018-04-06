#define AnalyzeBackground_cxx
#include "Analyzer/Analyzer/include/AnalyzeBackground.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <iostream>

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"

void AnalyzeBackground::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    std::vector<std::string> jettypes {"pt30", "pt45"};
    for(std::string jettype : jettypes)
    {
        std::string base = "h_njets_" + jettype;
        my_histos.emplace(base,new TH1D(base.c_str(),base.c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l",new TH1D( (base+"_0l").c_str(),(base+"_0l").c_str(),17,0,17));
        my_histos.emplace(base + "_1l",new TH1D((base+"_1l").c_str(),(base+"_1l").c_str(),15,0,15));
        my_histos.emplace(base + "_2l",new TH1D((base+"_2l").c_str(),(base+"_2l").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b",new TH1D((base+"_0l_g1b").c_str(),(base+"_0l_g1b").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b",new TH1D((base+"_1l_g1b").c_str(),(base+"_1l_g1b").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b",new TH1D((base+"_2l_g1b").c_str(),(base+"_2l_g1b").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b_ht500",new TH1D((base+"_0l_g1b_ht500").c_str(),(base+"_0l_g1b_ht500").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b_ht500",new TH1D((base+"_1l_g1b_ht500").c_str(),(base+"_1l_g1b_ht500").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b_ht500",new TH1D((base+"_2l_g1b_ht500").c_str(),(base+"_2l_g1b_ht500").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b_g1t",new TH1D((base+"_0l_g1b_g1t").c_str(),(base+"_0l_g1b_g1t").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b_g1t",new TH1D((base+"_1l_g1b_g1t").c_str(),(base+"_1l_g1b_g1t").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b_g1t",new TH1D((base+"_2l_g1b_g1t").c_str(),(base+"_2l_g1b_g1t").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b_g1t_ht500",new TH1D((base+"_0l_g1b_g1t_ht500").c_str(),(base+"_0l_g1b_g1t_ht500").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b_g1t_ht500",new TH1D((base+"_1l_g1b_g1t_ht500").c_str(),(base+"_1l_g1b_g1t_ht500").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b_g1t_ht500",new TH1D((base+"_2l_g1b_g1t_ht500").c_str(),(base+"_2l_g1b_g1t_ht500").c_str(),15,0,15));

        my_histos.emplace(base + "_0l_g1b_2t_ht500",new TH1D((base+"_0l_g1b_2t_ht500").c_str(),(base+"_0l_g1b_2t_ht500").c_str(),17,0,17));

        my_histos.emplace(base + "_1l_g1b_mbl",new TH1D((base+"_1l_g1b_mbl").c_str(),(base+"_1l_g1b_mbl").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt1",new TH1D((base+"_1l_g1b_mbl_bdt1").c_str(),(base+"_1l_g1b_mbl_bdt1").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt2",new TH1D((base+"_1l_g1b_mbl_bdt2").c_str(),(base+"_1l_g1b_mbl_bdt2").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt3",new TH1D((base+"_1l_g1b_mbl_bdt3").c_str(),(base+"_1l_g1b_mbl_bdt3").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt4",new TH1D((base+"_1l_g1b_mbl_bdt4").c_str(),(base+"_1l_g1b_mbl_bdt4").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher1",new TH1D((base+"_1l_g1b_mbl_fisher1").c_str(),(base+"_1l_g1b_mbl_fisher1").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher2",new TH1D((base+"_1l_g1b_mbl_fisher2").c_str(),(base+"_1l_g1b_mbl_fisher2").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher3",new TH1D((base+"_1l_g1b_mbl_fisher3").c_str(),(base+"_1l_g1b_mbl_fisher3").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher4",new TH1D((base+"_1l_g1b_mbl_fisher4").c_str(),(base+"_1l_g1b_mbl_fisher4").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_g1t",new TH1D((base+"_1l_g1b_mbl_g1t").c_str(),(base+"_1l_g1b_mbl_g1t").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_0t",new TH1D((base+"_1l_g1b_mbl_0t").c_str(),(base+"_1l_g1b_mbl_0t").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_1t1",new TH1D((base+"_1l_g1b_mbl_1t1").c_str(),(base+"_1l_g1b_mbl_1t1").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_1t2",new TH1D((base+"_1l_g1b_mbl_1t2").c_str(),(base+"_1l_g1b_mbl_1t2").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_1t3",new TH1D((base+"_1l_g1b_mbl_1t3").c_str(),(base+"_1l_g1b_mbl_1t3").c_str(),15,0,15));
        
        // For Z->ll control region
        my_histos.emplace(base + "_2l_onZ",new TH1D((base+"_2l_onZ").c_str(),(base+"_2l_onZ").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b",new TH1D((base+"_2l_onZ_g1b").c_str(),(base+"_2l_onZ_g1b").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl",new TH1D((base+"_2l_onZ_g1b_nombl").c_str(),(base+"_2l_onZ_g1b_nombl").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_g1t",new TH1D((base+"_2l_onZ_g1b_g1t").c_str(),(base+"_2l_onZ_g1b_g1t").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt1",new TH1D((base+"_2l_onZ_g1b_nombl_bdt1").c_str(),(base+"_2l_onZ_g1b_nombl_bdt1").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt2",new TH1D((base+"_2l_onZ_g1b_nombl_bdt2").c_str(),(base+"_2l_onZ_g1b_nombl_bdt2").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt3",new TH1D((base+"_2l_onZ_g1b_nombl_bdt3").c_str(),(base+"_2l_onZ_g1b_nombl_bdt3").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt4",new TH1D((base+"_2l_onZ_g1b_nombl_bdt4").c_str(),(base+"_2l_onZ_g1b_nombl_bdt4").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher1",new TH1D((base+"_2l_onZ_g1b_fisher1").c_str(),(base+"_2l_onZ_g1b_fisher1").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher2",new TH1D((base+"_2l_onZ_g1b_fisher2").c_str(),(base+"_2l_onZ_g1b_fisher2").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher3",new TH1D((base+"_2l_onZ_g1b_fisher3").c_str(),(base+"_2l_onZ_g1b_fisher3").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher4",new TH1D((base+"_2l_onZ_g1b_fisher4").c_str(),(base+"_2l_onZ_g1b_fisher4").c_str(),15,0,15));
    }
}

void AnalyzeBackground::Loop(NTupleReader& tr, double weight, int maxevents, std::string filetag, bool isQuiet)
{
    while(tr.getNextEvent())
    {
        const int& ntops_3jet = tr.getVar<int>("ntops_3jet");
        const int& ntops_2jet = tr.getVar<int>("ntops_2jet");
        const int& ntops_1jet = tr.getVar<int>("ntops_1jet");
        const std::string& runtype = tr.getVar<std::string>("runtype");
        const std::vector<TLorentzVector>& Jets         = tr.getVec<TLorentzVector>("Jets");
        const std::vector<double>&      Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const std::vector<std::string>& TriggerNames           = tr.getVec<std::string>("TriggerNames");
        const std::vector<int>&         TriggerPass            = tr.getVec<int>("TriggerPass");
        const TopTaggerResults* ttr         = tr.getVar<TopTaggerResults*>("ttr");
        const std::vector<TopObject*>& tops = ttr->getTops(); 

        const std::vector<TLorentzVector>& BJets = tr.getVec<TLorentzVector>("BJets");
        const int& NBJets = tr.getVar<int>("NBJets");
        const std::vector<TLorentzVector>& BJets_pt45 = tr.getVec<TLorentzVector>("BJets_pt45");
        const int& NBJets_pt45 = tr.getVar<int>("NBJets_pt45");
        const std::vector<TLorentzVector>& GoodMuons = tr.getVec<TLorentzVector>("GoodMuons");
        const std::vector<int>& GoodMuonsCharge      = tr.getVec<int>("GoodMuonsCharge");
        const std::vector<TLorentzVector>& GoodElectrons = tr.getVec<TLorentzVector>("GoodElectrons");
        const std::vector<int>& GoodElectronsCharge      = tr.getVec<int>("GoodElectronsCharge");
        const std::vector<TLorentzVector>& GoodLeptons = tr.getVec<TLorentzVector>("GoodLeptons");
        const int& NGoodLeptons = tr.getVar<int>("NGoodLeptons");
        const double& HT_trigger = tr.getVar<double>("HT_trigger");
        const double& Mbl        = tr.getVar<double>("Mbl");

        const bool& fisher_bin1 = tr.getVar<bool>("fisher_bin1");
        const bool& fisher_bin2 = tr.getVar<bool>("fisher_bin2");
        const bool& fisher_bin3 = tr.getVar<bool>("fisher_bin3");
        const bool& fisher_bin4 = tr.getVar<bool>("fisher_bin4");
        const bool& bdt_bin1    = tr.getVar<bool>("bdt_bin1");
        const bool& bdt_bin2    = tr.getVar<bool>("bdt_bin2");
        const bool& bdt_bin3    = tr.getVar<bool>("bdt_bin3");
        const bool& bdt_bin4    = tr.getVar<bool>("bdt_bin4");

        if (maxevents > 0 && tr.getEvtNum() >= maxevents) break;        

        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        
        // Exclude events with MadGraph HT > 100 from the DY inclusive sample
        if(filetag == "DYJetsToLL_M-50_Incl")
        {
            const double& madHT = tr.getVar<double>("madHT");
            if (madHT > 100) continue;
        }

        // Make sure event weight is not 0 for data
        double eventweight = 1.;
        if(runtype != "Data")
        {
            eventweight = tr.getVar<double>("Weight");
        }        
        // ------------------------------
        // -- Trigger for data
        // ------------------------------
      
        bool passTriggerAllHad = PassTriggerAllHad(TriggerNames, TriggerPass);
        bool passTriggerMuon = PassTriggerMuon(TriggerNames, TriggerPass);
        bool passTriggerElectron = PassTriggerElectron(TriggerNames, TriggerPass);
        if (runtype == "Data")
        {
            if (filetag == "Data_JetHT")
            {
                if(!passTriggerAllHad) continue;
                else
                {
                    passTriggerMuon = false;
                    passTriggerElectron = false;
                }
            }
            if (filetag == "Data_SingleMuon")
            { 
                if(!passTriggerMuon) continue;
                else
                {
                    passTriggerAllHad = false;
                    passTriggerElectron = false;
                }
            }
            if (filetag == "Data_SingleElectron")
            {
                if(!passTriggerElectron) continue;
                else
                {
                    passTriggerAllHad = false;
                    passTriggerMuon = false;
                }
            }
        }

        // -------------------------------
        // -- Basic event selection stuff
        // -------------------------------

        // Check whether event would pass the trigger requirement
        bool passTrigger = true;
        int rec_njet_pt45(0) ;
        int rec_njet_pt30(0) ;
        for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) 
        {
            TLorentzVector jlv( Jets.at(rji) ) ;
            if (abs(jlv.Eta()) > 2.4) continue;
            if ( jlv.Pt() > 30 ) 
            {
                rec_njet_pt30++;
            }
            if ( jlv.Pt() > 45 ) 
            {
                rec_njet_pt45++ ;
            }
        } 
        if ( !( HT_trigger>500 && rec_njet_pt45>=6 ) ) 
        {
            passTrigger = false;
        }

        bool passTrigger0l = false;
        bool passTrigger1l = false;
        bool passTrigger2l = false;
        if(runtype == "Data")
        {
            if (GoodMuons.size() > 0)
            {
                passTrigger1l = passTriggerMuon && (filetag == "Data_SingleMuon");
                passTrigger2l = passTriggerMuon && (filetag == "Data_SingleMuon");
            } 
            else if (GoodElectrons.size() > 0)
            {
                passTrigger1l = passTriggerElectron && (filetag == "Data_SingleElectron");
                passTrigger2l = passTriggerElectron && (filetag == "Data_SingleElectron");
            }
            else
            {
                passTrigger0l = passTriggerAllHad && (filetag == "Data_JetHT");
            }
        } 
        else 
        {
            passTrigger0l = true;
            passTrigger1l = true;
            passTrigger2l = true;
        }

        bool passBaseline0l = NGoodLeptons==0 && rec_njet_pt45>=6 && HT_trigger > 500 && NBJets_pt45 >= 1;
        bool passBaseline1l = NGoodLeptons==1 && rec_njet_pt30>=6 ;
        bool passBaseline2l = NGoodLeptons==2;

        bool passNtop = tops.size() >= 1;
        bool passNb = NBJets >= 1;
        bool onZ = false;
        bool passMbl_2l = false;
        if ( (GoodMuons.size() == 2) && (GoodMuonsCharge[0] != GoodMuonsCharge[1]) )
        {
            double mll = (GoodMuons[0] + GoodMuons[1]).M();
            if( mll > 81.2 && mll < 101.2)
                onZ = true;          

            // check whether a bl pair passes the M(b,l) cut
            for (TLorentzVector myb : BJets)
            {
                double mass_bl_1 = (GoodMuons[0] + myb).M();
                if(mass_bl_1 < 180 && mass_bl_1 > 30)
                    passMbl_2l = true;
                double mass_bl_2 = (GoodMuons[1] + myb).M();
                if(mass_bl_2 < 180 && mass_bl_2 > 30)
                    passMbl_2l = true;
            }
        } 
        else if ( (GoodElectrons.size() == 2) && (GoodElectronsCharge[0] != GoodElectronsCharge[1]) )
        {
            double mll = (GoodElectrons[0] + GoodElectrons[1]).M();
            if( mll > 81.2 && mll < 101.2)
                onZ = true;  
            // check whether a bl pair passes the M(b,l) cut
            for (TLorentzVector myb : BJets)
            {
                double mass_bl_1 = (GoodElectrons[0] + myb).M();
                if(mass_bl_1 < 180 && mass_bl_1 > 30)
                    passMbl_2l = true;
                double mass_bl_2 = (GoodElectrons[1] + myb).M();
                if(mass_bl_2 < 180 && mass_bl_2 > 30)
                    passMbl_2l = true;
            }
        }
        
        // -------------------------
        // -- Check r(j) behavior --
        // -------------------------

        int njets_rj = 0;
        std::vector<std::string> jettypes {"pt30", "pt45"};
        for(std::string jettype : jettypes)
        {
            if(jettype == "pt30")
                njets_rj = rec_njet_pt30;
            else if(jettype == "pt45")
                njets_rj = rec_njet_pt45;
            std::string base = "h_njets_" + jettype;

            my_histos[base]->Fill(njets_rj, eventweight);
            if(NGoodLeptons == 0 && passTrigger0l)
                my_histos[base+"_0l"]->Fill(njets_rj, eventweight);
            else if(NGoodLeptons == 1 && passTrigger1l)
                my_histos[base+"_1l"]->Fill(njets_rj, eventweight);
            else if(NGoodLeptons == 2 && passTrigger2l)
                my_histos[base+"_2l"]->Fill(njets_rj, eventweight);
            if (passNb)
            {
                if(NGoodLeptons == 0 && passTrigger0l)
                    my_histos[base+"_0l_g1b"]->Fill(njets_rj, eventweight);
                else if(NGoodLeptons == 1 && passTrigger1l)
                    my_histos[base+"_1l_g1b"]->Fill(njets_rj, eventweight);
                else if(NGoodLeptons == 2 && passTrigger2l)
                    my_histos[base+"_2l_g1b"]->Fill(njets_rj, eventweight);
              
                if (HT_trigger > 500)
                {
                    if(NGoodLeptons == 0 && passTrigger0l)
                        my_histos[base+"_0l_g1b_ht500"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 1 && passTrigger1l)
                        my_histos[base+"_1l_g1b_ht500"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 2 && passTrigger2l)
                        my_histos[base+"_2l_g1b_ht500"]->Fill(njets_rj, eventweight);
                }
              
                if (passNtop)
                {
                    if(NGoodLeptons == 0 && passTrigger0l)
                        my_histos[base+"_0l_g1b_g1t"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 1 && passTrigger1l)
                        my_histos[base+"_1l_g1b_g1t"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 2 && passTrigger2l)
                        my_histos[base+"_2l_g1b_g1t"]->Fill(njets_rj, eventweight);
                  
                    if (HT_trigger > 500)
                    {
                        if(NGoodLeptons == 0 && passTrigger0l)
                        {
                            my_histos[base+"_0l_g1b_g1t_ht500"]->Fill(njets_rj, eventweight);
                            if(tops.size() == 2)
                                my_histos[base+"_0l_g1b_2t_ht500"]->Fill(njets_rj, eventweight);
                        }
                        else if(NGoodLeptons == 1 && passTrigger1l)
                            my_histos[base+"_1l_g1b_g1t_ht500"]->Fill(njets_rj, eventweight);
                        else if(NGoodLeptons == 2 && passTrigger2l)
                            my_histos[base+"_2l_g1b_g1t_ht500"]->Fill(njets_rj, eventweight);
                    }
                }
            }

            // Dedicated 1l region
            if(NGoodLeptons == 1 && passNb && passTrigger1l)
            {
                bool passMtop = Mbl < 180 && Mbl > 30;
                
                if(passMtop)
                {
                    my_histos[base+"_1l_g1b_mbl"]->Fill(njets_rj, eventweight);

                    if(bdt_bin1)
                        my_histos[base+"_1l_g1b_mbl_bdt1"]->Fill(njets_rj, eventweight);
                    if(bdt_bin2)
                        my_histos[base+"_1l_g1b_mbl_bdt2"]->Fill(njets_rj, eventweight);
                    if(bdt_bin3)
                        my_histos[base+"_1l_g1b_mbl_bdt3"]->Fill(njets_rj, eventweight);
                    if(bdt_bin4)
                        my_histos[base+"_1l_g1b_mbl_bdt4"]->Fill(njets_rj, eventweight);

                    if(fisher_bin1)
                        my_histos[base+"_1l_g1b_mbl_fisher1"]->Fill(njets_rj, eventweight);
                    if(fisher_bin2)
                        my_histos[base+"_1l_g1b_mbl_fisher2"]->Fill(njets_rj, eventweight);
                    if(fisher_bin3)
                        my_histos[base+"_1l_g1b_mbl_fisher3"]->Fill(njets_rj, eventweight);
                    if(fisher_bin4)
                        my_histos[base+"_1l_g1b_mbl_fisher4"]->Fill(njets_rj, eventweight);

                    if(passNtop)
                    {
                        my_histos[base+"_1l_g1b_mbl_g1t"]->Fill(njets_rj, eventweight);

                        if(tops.size() == 1 && ntops_1jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t1"]->Fill(njets_rj, eventweight);
                        if(tops.size() == 1 && ntops_2jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t2"]->Fill(njets_rj, eventweight);
                        if(tops.size() == 1 && ntops_3jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t3"]->Fill(njets_rj, eventweight);
                    }
                    else 
                    {
                        my_histos[base+"_1l_g1b_mbl_0t"]->Fill(njets_rj, eventweight);
                    }
                }
            }
          
            // Now for the Z->ll region
            // h_njets_2l_onZ_g1b_g1t
            if (onZ && passTrigger2l)
            {
                my_histos[base+"_2l_onZ"]->Fill(njets_rj, eventweight);
                if(passNb)
                {
                    my_histos[base+"_2l_onZ_g1b"]->Fill(njets_rj, eventweight);

                    if(fisher_bin1)
                        my_histos[base+"_2l_onZ_g1b_fisher1"]->Fill(njets_rj, eventweight);
                    if(fisher_bin2)
                        my_histos[base+"_2l_onZ_g1b_fisher2"]->Fill(njets_rj, eventweight);
                    if(fisher_bin3)
                        my_histos[base+"_2l_onZ_g1b_fisher3"]->Fill(njets_rj, eventweight);
                    if(fisher_bin4)
                        my_histos[base+"_2l_onZ_g1b_fisher4"]->Fill(njets_rj, eventweight);
                  
                    if(!passMbl_2l)
                    {
                        my_histos[base+"_2l_onZ_g1b_nombl"]->Fill(njets_rj, eventweight);
                        if(bdt_bin1)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt1"]->Fill(njets_rj, eventweight);
                        if(bdt_bin2)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt2"]->Fill(njets_rj, eventweight);
                        if(bdt_bin3)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt3"]->Fill(njets_rj, eventweight);
                        if(bdt_bin4)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt4"]->Fill(njets_rj, eventweight);
                    }
                  
                    if(passNtop)
                    {
                        my_histos[base+"_2l_onZ_g1b_g1t"]->Fill(njets_rj, eventweight);
                    }
                }
            }
        }
    }
}


bool AnalyzeBackground::PassTriggerGeneral(std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    bool passTrigger = false;
    for(unsigned int i=0; i<TriggerNames.size(); ++i)
    {
        if(TriggerPass.at(i) != 1)
            continue;
        std::string trigname = TriggerNames.at(i);
        if( std::any_of(mytriggers.begin(), mytriggers.end(), [&] (std::string s) { return trigname.find(s)!=std::string::npos; }) )
        {
            passTrigger = true;
            break;
        }
    }
    return passTrigger;

}

bool AnalyzeBackground::PassTriggerAllHad(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {
        //"HLT_PFHT1050", // 2017 trigger
        //"HLT_PFHT900"
        //"HLT_PFHT380_SixPFJet32_DoublePFBTagCSV", // 2017 trigger
        //"HLT_PFHT430_SixPFJet40_PFBTagCSV", // 2017 trigger
        "HLT_PFHT450_SixJet40_BTagCSV",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
            };
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

bool AnalyzeBackground::PassTriggerMuon(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

bool AnalyzeBackground::PassTriggerElectron(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_Ele27_WPTight_Gsf"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

void AnalyzeBackground::WriteHistos()
{
    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
