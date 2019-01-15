#define Analyze1Lep_cxx
#include "Analyzer/Analyzer/include/Analyze1Lep.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/Framework/include/Utility.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

Analyze1Lep::Analyze1Lep() : initHistos(false)
{
}

void Analyze1Lep::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos, 
                             const std::vector<TH2DInfo>& hist2DInfos,  const std::vector<TH2DProfileInfo>& hist2DProfileInfos)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("fwm2_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm2_top6_1l_ge7j_ge1b","fwm2_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm3_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm3_top6_1l_ge7j_ge1b","fwm3_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm4_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm4_top6_1l_ge7j_ge1b","fwm4_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm5_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm5_top6_1l_ge7j_ge1b","fwm5_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm6_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm6_top6_1l_ge7j_ge1b","fwm6_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm7_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm7_top6_1l_ge7j_ge1b","fwm7_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm8_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm8_top6_1l_ge7j_ge1b","fwm8_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm9_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm9_top6_1l_ge7j_ge1b","fwm9_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm10_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm10_top6_1l_ge7j_ge1b","fwm10_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev0_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("jmt_ev0_top6_1l_ge7j_ge1b","jmt_ev0_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev1_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("jmt_ev1_top6_1l_ge7j_ge1b","jmt_ev1_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev2_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("jmt_ev2_top6_1l_ge7j_ge1b","jmt_ev2_top6_1l_ge7j_ge1b", 50, 0, 1 ) );

    for(unsigned int i = 1; i <= 7 ; i++) //Bad hard code
    {
        my_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b",  std::make_shared<TH1D>(("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b",   std::make_shared<TH1D>(("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 20, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 40, 0, 200));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 40, 0, 200));
    }

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            my_histos.emplace(hInfo.name+mycut.first, 
                              std::make_shared<TH1D>((hInfo.name+mycut.first).c_str(),(hInfo.name+mycut.first).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
        }

        for(const auto& h2dInfo : hist2DInfos)
        {
            my_2d_histos.emplace(h2dInfo.name+mycut.first, 
                                 std::make_shared<TH2D>((h2dInfo.name+mycut.first).c_str(),(h2dInfo.name+mycut.first).c_str(), 
                                                        h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
        }

        for(const auto& h2dProfile : hist2DProfileInfos)
        {
            my_2d_tp_histos.emplace(h2dProfile.name+mycut.first,
                                    std::make_shared<TProfile2D>((h2dProfile.name+mycut.first).c_str(),(h2dProfile.name+mycut.first).c_str(), 
                                                                 h2dProfile.nBinsX, h2dProfile.lowX, h2dProfile.highX, h2dProfile.nBinsY, h2dProfile.lowY, h2dProfile.highY, h2dProfile.lowZ, h2dProfile.highZ));
        }
    }
}

void Analyze1Lep::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& MET                  = tr.getVar<double>("MET");
        const auto& METPhi               = tr.getVar<double>("METPhi");
        const auto& ntops                = tr.getVar<int>("ntops");
        const auto& ntops_3jet           = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet           = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet           = tr.getVar<int>("ntops_1jet");
        const auto& runtype              = tr.getVar<std::string>("runtype");     
        const auto& filetag              = tr.getVar<std::string>("filetag");
        const auto& NJet                 = tr.getVar<int>("NJets");
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30");
        const auto& Jets                 = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt30        = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodJets_pt45       = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodBJets_pt45      = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodMuons           = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons       = tr.getVar<int>("NGoodElectrons");
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons");
        const auto& GoodLeptons          = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& fisher_val           = tr.getVar<double>("fisher_val");
        const auto& HT_trigger_pt30      = tr.getVar<double>("HT_trigger_pt30");
        const auto& JetID                = tr.getVar<bool>("JetID");
        const auto& passTrigger          = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC        = tr.getVar<bool>("passTriggerMC");
        const auto& passMadHT            = tr.getVar<bool>("passMadHT");
        const auto& passBaseline1l_Good  = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline1e1m     = tr.getVar<bool>("passBaseline1e1m_Good");
        const auto& Mbl                  = tr.getVar<double>("Mbl");
        const auto& MblVec               = tr.getVec<double>("MblVec");
        const auto& passBlind            = tr.getVar<bool>("passBlindLep_Good");            
        const auto& deepESM_val          = tr.getVar<double>("deepESM_val");
        const auto& deepESM_bin1         = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESM_bin2         = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESM_bin3         = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESM_bin4         = tr.getVar<bool>("deepESM_bin4");
        const auto& deepESM_binNum       = tr.getVar<int>("deepESM_binNum");
        const auto& fwm2_top6            = tr.getVar<double>("fwm2_top6");
        const auto& fwm3_top6            = tr.getVar<double>("fwm3_top6");
        const auto& fwm4_top6            = tr.getVar<double>("fwm4_top6");
        const auto& fwm5_top6            = tr.getVar<double>("fwm5_top6");
        const auto& fwm6_top6            = tr.getVar<double>("fwm6_top6");
        const auto& fwm7_top6            = tr.getVar<double>("fwm7_top6");
        const auto& fwm8_top6            = tr.getVar<double>("fwm8_top6");
        const auto& fwm9_top6            = tr.getVar<double>("fwm9_top6");
        const auto& fwm10_top6           = tr.getVar<double>("fwm10_top6");
        const auto& jmt_ev0_top6         = tr.getVar<double>("jmt_ev0_top6");
        const auto& jmt_ev1_top6         = tr.getVar<double>("jmt_ev1_top6");
        const auto& jmt_ev2_top6         = tr.getVar<double>("jmt_ev2_top6");
        const auto& Jets_cm_top6         = tr.getVec<TLorentzVector>("Jets_cm_top6");

        // ------------------------
        // -- Print event number
        // ------------------------       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if(tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() );

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight = 1.0;
        double eventweight = 1.0;
        double leptonweight = 1.0;
        double PileupWeight = 1.0;
        double bTagWeight = 1.0;
        double htDerivedweight = 1.0;
        double topPtScaleFactor = 1.0;
        double prefiringScaleFactor = 1.0;
        double FSRUp = 1.0;
        double FSRDown = 1.0;
        double FSRUp_2 = 1.0;
        double FSRDown_2 = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonweight = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedweight = tr.getVar<double>("htDerivedweight");
            topPtScaleFactor = tr.getVar<double>("topPtScaleFactor");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            FSRUp = tr.getVar<double>("PSweight_FSRUp");
            FSRDown = tr.getVar<double>("PSweight_FSRDown");
            FSRUp_2 = tr.getVar<double>("PSweight_FSRUp_2");
            FSRDown_2 = tr.getVar<double>("PSweight_FSRDown_2");
            
            weight *= eventweight*leptonweight*bTagWeight*htDerivedweight*prefiringScaleFactor;
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------
        bool pass_0l              = NGoodLeptons == 0;
        bool pass_1l              = NGoodLeptons == 1;
        bool pass_ht              = HT_trigger_pt30 > 300;
        bool pass_MBL             = (50 < Mbl && Mbl < 250);
        bool pass_1e_1m = false;
        if(NGoodLeptons == 2)
            pass_1e_1m = GoodLeptons[0].first != GoodLeptons[1].first;
        bool pass_njet_pt30       = NGoodJets_pt30 >= 7;
        bool pass_njet_pt30_1btag = NGoodBJets_pt30 >= 1;
        bool pass_njet_pt30_2btag = NGoodBJets_pt30 >= 2;
        bool pass_1e              = false;
        bool pass_1m              = false;
        if(pass_1l)
        { 
            if(GoodLeptons[0].first == "e") pass_1e = true;
            if(GoodLeptons[0].first == "m") pass_1m = true;
        }
        bool pass_5to6njet_pt30 = (NGoodJets_pt30 == 5 || NGoodJets_pt30 == 6);

        bool passBaseline1l_AllJets = JetID              &&
                                 HT_trigger_pt30 > 300   &&
                                 passMadHT               &&
                                 passTrigger             &&
                                 passTriggerMC           &&
                                 passBlind               &&
                                 NGoodBJets_pt30 >= 1    &&
                                 (50 < Mbl && Mbl < 250) &&
            (
                ((runtype != "Data" || filetag.find("Data_SingleMuon")     != std::string::npos) && NGoodMuons == 1     && NGoodElectrons == 0)
                                                                                     ||
                ((runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) && NGoodElectrons == 1 && NGoodMuons == 0)
            );
        
        // ---------------------------
        // --    1 Top selection
        // ---------------------------
        //exactly selections
        bool pass_1t    = ntops==1;
        bool pass_1t_d1 = pass_1t && deepESM_bin1, pass_1t_d2 = pass_1t && deepESM_bin2, pass_1t_d3 = pass_1t && deepESM_bin3, pass_1t_d4 = pass_1t && deepESM_bin4;

        bool pass_1t1 = pass_1t && ntops_1jet==1;
        bool pass_1t2 = pass_1t && ntops_2jet==1;
        bool pass_1t3 = pass_1t && ntops_3jet==1;
        bool pass_1t2or3 = pass_1t2 or pass_1t3;

        bool pass_1t1_d1 = pass_1t1 && deepESM_bin1, pass_1t1_d2 = pass_1t1 && deepESM_bin2, pass_1t1_d3 = pass_1t1 && deepESM_bin3, pass_1t1_d4 = pass_1t1 && deepESM_bin4;
        bool pass_1t2_d1 = pass_1t2 && deepESM_bin1, pass_1t2_d2 = pass_1t2 && deepESM_bin2, pass_1t2_d3 = pass_1t2 && deepESM_bin3, pass_1t2_d4 = pass_1t2 && deepESM_bin4;
        bool pass_1t3_d1 = pass_1t3 && deepESM_bin1, pass_1t3_d2 = pass_1t3 && deepESM_bin2, pass_1t3_d3 = pass_1t3 && deepESM_bin3, pass_1t3_d4 = pass_1t3 && deepESM_bin4;
        bool pass_1t2or3_d1 = pass_1t2or3 && deepESM_bin1, pass_1t2or3_d2 = pass_1t2or3 && deepESM_bin2, pass_1t2or3_d3 = pass_1t2or3 && deepESM_bin3, pass_1t2or3_d4 = pass_1t2or3 && deepESM_bin4;

        //at least selections
        bool pass_ge1t    = ntops>=1;
        bool pass_ge1t_d1 = pass_ge1t && deepESM_bin1, pass_ge1t_d2 = pass_ge1t && deepESM_bin2, pass_ge1t_d3 = pass_ge1t && deepESM_bin3, pass_ge1t_d4 = pass_ge1t && deepESM_bin4;

        bool pass_ge1t1 = pass_ge1t && ntops_1jet>=1;
        bool pass_ge1t2 = pass_ge1t && ntops_2jet>=1;
        bool pass_ge1t3 = pass_ge1t && ntops_3jet>=1;

        bool pass_ge1t1_d1 = pass_ge1t1 && deepESM_bin1, pass_ge1t1_d2 = pass_ge1t1 && deepESM_bin2, pass_ge1t1_d3 = pass_ge1t1 && deepESM_bin3, pass_ge1t1_d4 = pass_ge1t1 && deepESM_bin4;
        bool pass_ge1t2_d1 = pass_ge1t2 && deepESM_bin1, pass_ge1t2_d2 = pass_ge1t2 && deepESM_bin2, pass_ge1t2_d3 = pass_ge1t2 && deepESM_bin3, pass_ge1t2_d4 = pass_ge1t2 && deepESM_bin4;
        bool pass_ge1t3_d1 = pass_ge1t3 && deepESM_bin1, pass_ge1t3_d2 = pass_ge1t3 && deepESM_bin2, pass_ge1t3_d3 = pass_ge1t3 && deepESM_bin3, pass_ge1t3_d4 = pass_ge1t3 && deepESM_bin4;

        // ------------------------------------------------
        // --    Tempory Home of the W+Jets Control Region
        // ------------------------------------------------
        double mT = -1;
        bool pass_mT = false;
        bool pass_Mbl_all = false;
        if(pass_1l)
        {
            TLorentzVector metLV;
            metLV.SetPtEtaPhiE(MET,0,METPhi,MET);
            mT = utility::calcMT(GoodLeptons[0].second, metLV);
            pass_mT = mT > 50 && mT < 110;
        
            for(int i = 0; i < Jets.size(); ++i)
            {
                if(!GoodJets_pt30[i]) continue;
                double mbl = (GoodLeptons[0].second+Jets[i]).M();
                if(mbl > 50 && mbl < 250)
                {
                    pass_Mbl_all = true;
                    break;
                }
            }            
        }

        const auto& NGoodBJets_pt30_loose = tr.getVar<int>("NGoodBJets_pt30_loose");        
        bool pass_0b_loose = NGoodBJets_pt30_loose == 0;

        bool passBaseline1l_WCR = JetID                   &&
                                  passMadHT               &&
                                  passTrigger             &&
                                  passTriggerMC           &&
                                  HT_trigger_pt30 > 300   &&
                                  pass_mT                 && 
                                  pass_0b_loose           && 
                                  !pass_Mbl_all           && 
                                  MET>30                  &&
            (
                ((runtype != "Data" || filetag.find("Data_SingleMuon")     != std::string::npos) && NGoodMuons == 1     && NGoodElectrons == 0)
                                                                                     ||
                ((runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) && NGoodElectrons == 1 && NGoodMuons == 0)
            );

        bool evenEvent = tr.getEvtNum() % 2 == 0;

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map_1l 
        {
            {""                                , true                                                                     },
            {"_1l"                             , pass_1l && pass_ht                                                       },
            {"_1l_ge7j"                        , pass_1l && pass_ht && pass_njet_pt30 && JetID                            },
            {"_1l_ge1b"                        , pass_1l && pass_ht && pass_njet_pt30_1btag && JetID                      },
            {"_1l_ge2b"                        , pass_1l && pass_ht && pass_njet_pt30_2btag && JetID                      },
            {"_1l_0b_ge300ht_50to110mt_ge30MET", passBaseline1l_WCR                                                       },
            {"_1l_0b_ge300ht_50to110mt_ge30MET_even", passBaseline1l_WCR && evenEvent                                          },
            {"_1l_0b_ge300ht_50to110mt_ge30MET_odd" , passBaseline1l_WCR && !evenEvent                                         },
            {"_1e_1m_ge2b_le5j"                , passBaseline1e1m                                                         },
            {"_1l_1t"                          , pass_1l && pass_ht && pass_1t && JetID                                   },
            {"_1l_ge1t"                        , pass_1l && pass_ht && pass_ge1t && JetID                                 },
            {"_1l_ge7j_ge2b"                   , passBaseline1l_AllJets && pass_njet_pt30_2btag                           },
            {"_1l_ge7j_ge1b_noMbl"             , pass_1l && pass_ht && pass_njet_pt30 && pass_njet_pt30_1btag && JetID    },
            {"_1l_ge7j_ge1b"                   , passBaseline1l_Good                                                      },                         
            {"_1e_ge7j_ge1b"                   , passBaseline1l_Good && pass_1e                                           },                         
            {"_1m_ge7j_ge1b"                   , passBaseline1l_Good && pass_1m                                           },                         
            {"_1l_ge7j_ge1b_d1"                , passBaseline1l_Good && deepESM_bin1                                      },                         
            {"_1l_ge7j_ge1b_d2"                , passBaseline1l_Good && deepESM_bin2                                      },                         
            {"_1l_ge7j_ge1b_d3"                , passBaseline1l_Good && deepESM_bin3                                      },                         
            {"_1l_ge7j_ge1b_d4"                , passBaseline1l_Good && deepESM_bin4                                      },                         
            {"_1l_1j_ge1b"                     , passBaseline1l_AllJets && NGoodJets_pt30 == 1                            },
            {"_1l_2j_ge1b"                     , passBaseline1l_AllJets && NGoodJets_pt30 == 2                            },
            {"_1l_3j_ge1b"                     , passBaseline1l_AllJets && NGoodJets_pt30 == 3                            },
            {"_1l_4j_ge1b"                     , passBaseline1l_AllJets && NGoodJets_pt30 == 4                            },
            {"_1l_5j_ge1b"                     , passBaseline1l_AllJets && NGoodJets_pt30 == 5                            },
            {"_1l_6j_ge1b"                     , passBaseline1l_AllJets && NGoodJets_pt30 == 6                            },
            {"_1l_7j_ge1b"                     , passBaseline1l_Good && NGoodJets_pt30 == 7                               },
            {"_1l_8j_ge1b"                     , passBaseline1l_Good && NGoodJets_pt30 == 8                               },
            {"_1l_9j_ge1b"                     , passBaseline1l_Good && NGoodJets_pt30 == 9                               },
            {"_1l_10j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 10                              },
            {"_1l_11j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 11                              },
            {"_1l_12j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 12                              },
            {"_1l_13j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 13                              },
            {"_1l_14j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 14                              },
            {"_1l_15j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 15                              },
            {"_1l_ge7j_ge1b_1t"                , passBaseline1l_Good && pass_1t                                           },
            {"_1l_ge7j_ge1b_1t_d1"             , passBaseline1l_Good && pass_1t_d1                                        },
            {"_1l_ge7j_ge1b_1t_d2"             , passBaseline1l_Good && pass_1t_d2                                        }, 
            {"_1l_ge7j_ge1b_1t_d3"             , passBaseline1l_Good && pass_1t_d3                                        }, 
            {"_1l_ge7j_ge1b_1t_d4"             , passBaseline1l_Good && pass_1t_d4                                        },
            {"_1l_ge7j_ge1b_ge1t"              , passBaseline1l_Good && pass_ge1t                                         },
            {"_1l_ge7j_ge1b_ge1t_d1"           , passBaseline1l_Good && pass_ge1t_d1                                      },
            {"_1l_ge7j_ge1b_ge1t_d2"           , passBaseline1l_Good && pass_ge1t_d2                                      }, 
            {"_1l_ge7j_ge1b_ge1t_d3"           , passBaseline1l_Good && pass_ge1t_d3                                      }, 
            {"_1l_ge7j_ge1b_ge1t_d4"           , passBaseline1l_Good && pass_ge1t_d4                                      },
            {"_1l_5j_ge1b_htCorr"              , passBaseline1l_AllJets && NGoodJets_pt30 == 5                            },
            {"_1l_6j_ge1b_htCorr"              , passBaseline1l_AllJets && NGoodJets_pt30 == 6                            },
            {"_1l_7j_ge1b_htCorr"              , passBaseline1l_Good && NGoodJets_pt30 == 7                               },
            
            {"_1l_ge7j_ge1b_1t1"               , passBaseline1l_Good && pass_1t1                                          },
            {"_1l_ge7j_ge1b_1t2"               , passBaseline1l_Good && pass_1t2                                          },
            {"_1l_ge7j_ge1b_1t3"               , passBaseline1l_Good && pass_1t3                                          },
            {"_1l_ge7j_ge1b_1t2or3"            , passBaseline1l_Good && pass_1t2or3                                       },
            {"_1l_ge7j_ge1b_1t1_d1"            , passBaseline1l_Good && pass_1t1_d1                                       },
            {"_1l_ge7j_ge1b_1t1_d2"            , passBaseline1l_Good && pass_1t1_d2                                       },
            {"_1l_ge7j_ge1b_1t1_d3"            , passBaseline1l_Good && pass_1t1_d3                                       },
            {"_1l_ge7j_ge1b_1t1_d4"            , passBaseline1l_Good && pass_1t1_d4                                       },
            {"_1l_ge7j_ge1b_1t2_d1"            , passBaseline1l_Good && pass_1t2_d1                                       },
            {"_1l_ge7j_ge1b_1t2_d2"            , passBaseline1l_Good && pass_1t2_d2                                       },
            {"_1l_ge7j_ge1b_1t2_d3"            , passBaseline1l_Good && pass_1t2_d3                                       },
            {"_1l_ge7j_ge1b_1t2_d4"            , passBaseline1l_Good && pass_1t2_d4                                       },
            {"_1l_ge7j_ge1b_1t3_d1"            , passBaseline1l_Good && pass_1t3_d1                                       },
            {"_1l_ge7j_ge1b_1t3_d2"            , passBaseline1l_Good && pass_1t3_d2                                       },
            {"_1l_ge7j_ge1b_1t3_d3"            , passBaseline1l_Good && pass_1t3_d3                                       },
            {"_1l_ge7j_ge1b_1t3_d4"            , passBaseline1l_Good && pass_1t3_d4                                       },
            {"_1l_ge7j_ge1b_1t2or3_d1"         , passBaseline1l_Good && pass_1t2or3_d1                                    },
            {"_1l_ge7j_ge1b_1t2or3_d2"         , passBaseline1l_Good && pass_1t2or3_d2                                    },
            {"_1l_ge7j_ge1b_1t2or3_d3"         , passBaseline1l_Good && pass_1t2or3_d3                                    },
            {"_1l_ge7j_ge1b_1t2or3_d4"         , passBaseline1l_Good && pass_1t2or3_d4                                    },
            {"_1l_ge7j_ge1b_ge1t1"             , passBaseline1l_Good && pass_ge1t1                                        },
            {"_1l_ge7j_ge1b_ge1t2"             , passBaseline1l_Good && pass_ge1t2                                        },
            {"_1l_ge7j_ge1b_ge1t3"             , passBaseline1l_Good && pass_ge1t3                                        },
            {"_1l_ge7j_ge1b_ge1t1_d1"          , passBaseline1l_Good && pass_ge1t1_d1                                     },
            {"_1l_ge7j_ge1b_ge1t1_d2"          , passBaseline1l_Good && pass_ge1t1_d2                                     },
            {"_1l_ge7j_ge1b_ge1t1_d3"          , passBaseline1l_Good && pass_ge1t1_d3                                     },
            {"_1l_ge7j_ge1b_ge1t1_d4"          , passBaseline1l_Good && pass_ge1t1_d4                                     },
            {"_1l_ge7j_ge1b_ge1t2_d1"          , passBaseline1l_Good && pass_ge1t2_d1                                     },
            {"_1l_ge7j_ge1b_ge1t2_d2"          , passBaseline1l_Good && pass_ge1t2_d2                                     },
            {"_1l_ge7j_ge1b_ge1t2_d3"          , passBaseline1l_Good && pass_ge1t2_d3                                     },
            {"_1l_ge7j_ge1b_ge1t2_d4"          , passBaseline1l_Good && pass_ge1t2_d4                                     },
            {"_1l_ge7j_ge1b_ge1t3_d1"          , passBaseline1l_Good && pass_ge1t3_d1                                     },
            {"_1l_ge7j_ge1b_ge1t3_d2"          , passBaseline1l_Good && pass_ge1t3_d2                                     },
            {"_1l_ge7j_ge1b_ge1t3_d3"          , passBaseline1l_Good && pass_ge1t3_d3                                     },
            {"_1l_ge7j_ge1b_ge1t3_d4"          , passBaseline1l_Good && pass_ge1t3_d4                                     },                                                 
        };

        std::vector<TH1DInfo> histInfos = {
            {    "h_njets",         20,   0.0,   20.0},
            {"blind_njets",         20,   0.0,   20.0},
            {    "h_ntops",         10,   0.0,   10.0},
            {"blind_ntops",         10,   0.0,   10.0},
            {    "h_nb",            10,   0.0,   10.0},
            {"blind_nb",            10,   0.0,   10.0},
            {    "h_fisher",       200,  -0.5,    0.5},
            {"blind_fisher",       200,  -0.5,    0.5},
            {    "h_deepESM",      200,   0.0,    1.0},
            {"blind_deepESM",      200,   0.0,    1.0},
            {    "h_deepESMMerged",  4,   0.5,    4.5},
            {"blind_deepESMMerged",  4,   0.5,    4.5},
            {    "h_ht",           300,   0.0, 3000.0},
            {"blind_ht",           300,   0.0, 3000.0},
            {    "h_mbl",          300,   0.0,  300.0},
            {"blind_mbl",          300,   0.0,  300.0},
            {    "h_lPt",          200,   0.0, 2000.0},
            {"blind_lPt",          200,   0.0, 2000.0},
            {    "h_lEta",         100,  -6.0,    6.0},
            {"blind_lEta",         100,  -6.0,    6.0},
            {    "h_jPt",          200,   0.0, 2000.0},
            {"blind_jPt",          200,   0.0, 2000.0},
            {    "h_jEta",         100,  -6.0,    6.0},
            {"blind_jEta",         100,  -6.0,    6.0},
            {    "h_allMbl",       300,   0.0,  300.0},            
            {"blind_allMbl",       300,   0.0,  300.0},
            {"h_weight",               200,  -5.0, 5.0},
            {"h_leptonweight",         200,  -5.0, 5.0},
            {"h_PileupWeight",         200,  -5.0, 5.0},
            {"h_bTagWeight",           200,  -5.0, 5.0},
            {"h_htDerivedweight",      200,  -5.0, 5.0},
            {"h_prefiringScaleFactor", 200,  -5.0, 5.0},
            {"h_FSRUp",                200,  -5.0, 5.0},
            {"h_FSRUp_2",              200,  -5.0, 5.0},
            {"h_FSRDown",              200,  -5.0, 5.0},
            {"h_FSRDown_2",            200,  -5.0, 5.0},
        };

        std::vector<TH2DInfo> hist2DInfos = {
            {    "h_njets_fisher",  15, 0,   15, 200, -0.5, 0.5},
            {"blind_njets_fisher",  15, 0,   15, 200, -0.5, 0.5},
            {    "h_njets_deepESM", 15, 0,   15, 200,  0.0, 1.0},
            {"blind_njets_deepESM", 15, 0,   15, 200,  0.0, 1.0},
            {    "h_njets_mbl",     15, 0,   15, 300,  0.0, 300},
            {"blind_njets_mbl",     15, 0,   15, 300,  0.0, 300},
            {    "h_ht_deepESM",   300, 0, 3000, 200,  0.0, 1.0},
            {"blind_ht_deepESM",   300, 0, 3000, 200,  0.0, 1.0},
        };

        std::vector<TH2DProfileInfo> hist2DProfileInfos = {
            {"h_njets_deepESMMerged_preFireSF", 15,0,15, 4,0.5,4.5, 0,1.0},
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map_1l, histInfos, hist2DInfos, hist2DProfileInfos);
            initHistos = true;
        }

        // Global cuts
        if ( !(passTriggerMC && passTrigger && passMadHT && passBlind) ) continue;

        for(auto& kv : cut_map_1l)
        {
            if(kv.second)
            {
                double w = weight;
                if(kv.first.find("50to110mt") != std::string::npos) w = eventweight*leptonweight*bTagWeight*prefiringScaleFactor;
                my_histos["h_njets"               +kv.first]->Fill(NGoodJets_pt30, w);
                my_histos["h_ntops"               +kv.first]->Fill(ntops, w);
                my_histos["h_nb"                  +kv.first]->Fill(NGoodBJets_pt30, w);
                my_histos["h_fisher"              +kv.first]->Fill(fisher_val, w);
                my_histos["h_deepESM"             +kv.first]->Fill(deepESM_val, w);
                my_histos["h_deepESMMerged"       +kv.first]->Fill(deepESM_binNum, w);
                my_histos["h_ht"                  +kv.first]->Fill(HT_trigger_pt30, w);
                my_histos["h_mbl"                 +kv.first]->Fill(Mbl, w);
                my_histos["h_weight"              +kv.first]->Fill(weight, w);
                my_histos["h_leptonweight"        +kv.first]->Fill(leptonweight, w);
                my_histos["h_PileupWeight"        +kv.first]->Fill(PileupWeight, w);
                my_histos["h_bTagWeight"          +kv.first]->Fill(bTagWeight, w);
                my_histos["h_htDerivedweight"     +kv.first]->Fill(htDerivedweight, w);
                my_histos["h_prefiringScaleFactor"+kv.first]->Fill(prefiringScaleFactor, w);
                my_histos["h_FSRUp"               +kv.first]->Fill(FSRUp);
                my_histos["h_FSRDown"             +kv.first]->Fill(FSRDown);
                my_histos["h_FSRUp_2"             +kv.first]->Fill(FSRUp_2);
                my_histos["h_FSRDown_2"           +kv.first]->Fill(FSRDown_2);
                for(const auto& l : GoodLeptons)
                {
                    my_histos["h_lPt"+kv.first]->Fill(l.second.Pt(), w);
                    my_histos["h_lEta"+kv.first]->Fill(l.second.Eta(), w);
                }
                for(const auto& mbl : MblVec)
                {
                    my_histos["h_allMbl"+kv.first]->Fill(mbl, w);
                }
                for(int j = 0; j < Jets.size(); j++)
                {
                    if(!GoodJets_pt30[j]) continue;
                    my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                    my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                }
                my_2d_histos["h_njets_fisher"+kv.first]->Fill(NGoodJets_pt30, fisher_val, w);
                my_2d_histos["h_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, w);
                my_2d_histos["h_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, w);
                my_2d_histos["h_ht_deepESM"+kv.first]->Fill(HT_trigger_pt30, deepESM_val, w);
                my_2d_tp_histos["h_njets_deepESMMerged_preFireSF"+kv.first]->Fill(NJet, deepESM_binNum, prefiringScaleFactor, w);

                if ( NGoodJets_pt30 <= 7 )
                {
                    my_histos["blind_njets"         +kv.first]->Fill(NGoodJets_pt30, w);
                    my_histos["blind_ntops"         +kv.first]->Fill(ntops, w);
                    my_histos["blind_nb"            +kv.first]->Fill(NGoodBJets_pt30, w);
                    my_histos["blind_fisher"        +kv.first]->Fill(fisher_val, w);
                    my_histos["blind_deepESM"       +kv.first]->Fill(deepESM_val, w);
                    my_histos["blind_deepESMMerged" +kv.first]->Fill(deepESM_binNum, w);
                    my_histos["blind_ht"            +kv.first]->Fill(HT_trigger_pt30, w);
                    my_histos["blind_mbl"           +kv.first]->Fill(Mbl, w);
                    for(const auto l : GoodLeptons)
                    {
                        my_histos["blind_lPt"+kv.first]->Fill(l.second.Pt(), w);
                        my_histos["blind_lEta"+kv.first]->Fill(l.second.Eta(), w);
                    }
                    for(const auto& mbl : MblVec)
                    {
                        my_histos["blind_allMbl"+kv.first]->Fill(mbl, w);
                    }
                    for(int j = 0; j < Jets.size(); j++)
                    {
                        if(!GoodJets_pt30[j]) continue;
                        my_histos["blind_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                        my_histos["blind_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                    }
                    my_2d_histos["blind_njets_fisher"+kv.first]->Fill(NGoodJets_pt30, fisher_val, w);
                    my_2d_histos["blind_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, w);
                    my_2d_histos["blind_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, w);
                    my_2d_histos["blind_ht_deepESM"+kv.first]->Fill(HT_trigger_pt30, deepESM_val, w);
                }
            }
        }

        // Fill histos for deepESM and Fisher training
        if(passBaseline1l_Good)
        {
            my_histos["fwm2_top6_1l_ge7j_ge1b"]->Fill(fwm2_top6, weight);
            my_histos["fwm3_top6_1l_ge7j_ge1b"]->Fill(fwm3_top6, weight);
            my_histos["fwm4_top6_1l_ge7j_ge1b"]->Fill(fwm4_top6, weight);
            my_histos["fwm5_top6_1l_ge7j_ge1b"]->Fill(fwm5_top6, weight);
            my_histos["fwm6_top6_1l_ge7j_ge1b"]->Fill(fwm6_top6, weight);
            my_histos["fwm7_top6_1l_ge7j_ge1b"]->Fill(fwm7_top6, weight);
            my_histos["fwm8_top6_1l_ge7j_ge1b"]->Fill(fwm8_top6, weight);
            my_histos["fwm9_top6_1l_ge7j_ge1b"]->Fill(fwm9_top6, weight);
            my_histos["fwm10_top6_1l_ge7j_ge1b"]->Fill(fwm10_top6, weight);
            my_histos["jmt_ev0_top6_1l_ge7j_ge1b"]->Fill(jmt_ev0_top6, weight);
            my_histos["jmt_ev1_top6_1l_ge7j_ge1b"]->Fill(jmt_ev1_top6, weight);
            my_histos["jmt_ev2_top6_1l_ge7j_ge1b"]->Fill(jmt_ev2_top6, weight);

            for(unsigned int i = 0; i < Jets_cm_top6.size(); i++)
            {
                my_histos["Jet_cm_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Pt()), weight);
                my_histos["Jet_cm_eta_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Eta()), weight);
                my_histos["Jet_cm_phi_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Phi()), weight);
                my_histos["Jet_cm_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).M()), weight);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).Pt(), weight);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).Eta(), weight);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).Phi(), weight);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).M(), weight);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Pt(), weight);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Eta(), weight);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Phi(), weight);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).M(), weight);
            }
        }
    } // end of event loop
}

void Analyze1Lep::WriteHistos(TFile* outfile)
{
    outfile->cd();
    
    for(const auto& p : my_histos) 
    {
        p.second->Write();
    }
    
    for(const auto& p : my_2d_histos) 
    {
        p.second->Write();
    }

    for(const auto& p : my_tp_histos) 
    {
        p.second->Write();
    }

    for(const auto& p : my_2d_tp_histos)
    {
        p.second->Write();
    }
    
    for(const auto& p : my_efficiencies) 
    {
        p.second->Write();
    }    
}
