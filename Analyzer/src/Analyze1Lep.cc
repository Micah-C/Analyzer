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

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

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
    my_histos.emplace("GoodLeptons_pt_1_1l_ge7j_ge1b", std::make_shared<TH1D>("GoodLeptons_pt_1_1l_ge7j_ge1b","GoodLeptons_pt_1_1l_ge7j_ge1b", 150, 0, 1500 ) );
    my_histos.emplace("GoodLeptons_eta_1_1l_ge7j_ge1b", std::make_shared<TH1D>("GoodLeptons_eta_1_1l_ge7j_ge1b","GoodLeptons_eta_1_1l_ge7j_ge1b", 100, -6, 6 ) );
    my_histos.emplace("GoodLeptons_phi_1_1l_ge7j_ge1b", std::make_shared<TH1D>("GoodLeptons_phi_1_1l_ge7j_ge1b","GoodLeptons_phi_1_1l_ge7j_ge1b", 80, -4, 4 ) );
    my_histos.emplace("GoodLeptons_m_1_1l_ge7j_ge1b", std::make_shared<TH1D>("GoodLeptons_m_1_1l_ge7j_ge1b","GoodLeptons_m_1_1l_ge7j_ge1b", 20, 0, 200 ) );

    for(unsigned int i = 1; i <= 7 ; i++) //Bad hard code
    {
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta00_05_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta00_05_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta00_05_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta05_10_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta05_10_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta05_10_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta10_15_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta10_15_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta10_15_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta15_20_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta15_20_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta15_20_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta20_Inf_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta20_Inf_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta20_Inf_p").c_str(), 150, 0, 1500, 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta00_05_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta00_05_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta00_05_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta05_10_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta05_10_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta05_10_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta10_15_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta10_15_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta10_15_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta15_20_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta15_20_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta15_20_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta20_Inf_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta20_Inf_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_eta20_Inf_n").c_str(), 150, 0, 1500, 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_eta_pt_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 100, -6, 6, 150, 0, 1500 ));
        my_2d_histos.emplace("Jet_cm_eta_m_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 100, -6, 6, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_phi_pt_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 80, -4, 4, 150, 0, 1500 ));
        my_2d_histos.emplace("Jet_cm_phi_m_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 80, -4, 4, 200, 0, 200 ));

        my_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b",  std::make_shared<TH1D>(("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b",   std::make_shared<TH1D>(("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 40, 0, 200));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 40, 0, 200));
    }

    my_histos.emplace("fwm2_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm2_top6_1l_ge7j_ge1b_NoHTweight","fwm2_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm3_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm3_top6_1l_ge7j_ge1b_NoHTweight","fwm3_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm4_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm4_top6_1l_ge7j_ge1b_NoHTweight","fwm4_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm5_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm5_top6_1l_ge7j_ge1b_NoHTweight","fwm5_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm6_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm6_top6_1l_ge7j_ge1b_NoHTweight","fwm6_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm7_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm7_top6_1l_ge7j_ge1b_NoHTweight","fwm7_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm8_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm8_top6_1l_ge7j_ge1b_NoHTweight","fwm8_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm9_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm9_top6_1l_ge7j_ge1b_NoHTweight","fwm9_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("fwm10_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("fwm10_top6_1l_ge7j_ge1b_NoHTweight","fwm10_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev0_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("jmt_ev0_top6_1l_ge7j_ge1b_NoHTweight","jmt_ev0_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev1_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("jmt_ev1_top6_1l_ge7j_ge1b_NoHTweight","jmt_ev1_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev2_top6_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("jmt_ev2_top6_1l_ge7j_ge1b_NoHTweight","jmt_ev2_top6_1l_ge7j_ge1b_NoHTweight", 50, 0, 1 ) );
    my_histos.emplace("GoodLeptons_pt_1_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("GoodLeptons_pt_1_1l_ge7j_ge1b_NoHTweight","GoodLeptons_pt_1_1l_ge7j_ge1b_NoHTweight", 150, 0, 1500 ) );
    my_histos.emplace("GoodLeptons_eta_1_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("GoodLeptons_eta_1_1l_ge7j_ge1b_NoHTweight","GoodLeptons_eta_1_1l_ge7j_ge1b_NoHTweight", 100, -6, 6 ) );
    my_histos.emplace("GoodLeptons_phi_1_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("GoodLeptons_phi_1_1l_ge7j_ge1b_NoHTweight","GoodLeptons_phi_1_1l_ge7j_ge1b_NoHTweight", 80, -4, 4 ) );
    my_histos.emplace("GoodLeptons_m_1_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>("GoodLeptons_m_1_1l_ge7j_ge1b_NoHTweight","GoodLeptons_m_1_1l_ge7j_ge1b_NoHTweight", 20, 0, 200 ) );

    for(unsigned int i = 1; i <= 7 ; i++) //Bad hard code
    {

        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_p").c_str(), 150, 0, 1500, 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_n").c_str(), 150, 0, 1500, 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_eta_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_eta_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_eta_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 100, -6, 6, 150, 0, 1500 ));
        my_2d_histos.emplace("Jet_cm_eta_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_eta_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_eta_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 100, -6, 6, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_phi_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_phi_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_phi_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 80, -4, 4, 150, 0, 1500 ));
        my_2d_histos.emplace("Jet_cm_phi_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_phi_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_phi_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 80, -4, 4, 200, 0, 200 ));

        my_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight",  std::make_shared<TH1D>(("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>(("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH1D>(("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight",   std::make_shared<TH1D>(("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b_NoHTweight").c_str(), 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(), 200, 0.0, 1.0, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(), 200, 0.0, 1.0, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(), 200, 0.0, 1.0, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b_NoHTweight").c_str(), 200, 0.0, 1.0, 40, 0, 200));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(), 20, 0, 20, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(), 20, 0, 20, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(), 20, 0, 20, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(),("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b_NoHTweight").c_str(), 20, 0, 20, 40, 0, 200));
    }


    my_histos.emplace("fwm2_top6_QCDCR", std::make_shared<TH1D>("fwm2_top6_QCDCR","fwm2_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm3_top6_QCDCR", std::make_shared<TH1D>("fwm3_top6_QCDCR","fwm3_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm4_top6_QCDCR", std::make_shared<TH1D>("fwm4_top6_QCDCR","fwm4_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm5_top6_QCDCR", std::make_shared<TH1D>("fwm5_top6_QCDCR","fwm5_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm6_top6_QCDCR", std::make_shared<TH1D>("fwm6_top6_QCDCR","fwm6_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm7_top6_QCDCR", std::make_shared<TH1D>("fwm7_top6_QCDCR","fwm7_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm8_top6_QCDCR", std::make_shared<TH1D>("fwm8_top6_QCDCR","fwm8_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm9_top6_QCDCR", std::make_shared<TH1D>("fwm9_top6_QCDCR","fwm9_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("fwm10_top6_QCDCR", std::make_shared<TH1D>("fwm10_top6_QCDCR","fwm10_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev0_top6_QCDCR", std::make_shared<TH1D>("jmt_ev0_top6_QCDCR","jmt_ev0_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev1_top6_QCDCR", std::make_shared<TH1D>("jmt_ev1_top6_QCDCR","jmt_ev1_top6_QCDCR", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev2_top6_QCDCR", std::make_shared<TH1D>("jmt_ev2_top6_QCDCR","jmt_ev2_top6_QCDCR", 50, 0, 1 ) );

    for(unsigned int i = 1; i <= 7 ; i++) //Bad hard code
    {
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta00_05_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta00_05_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta00_05_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta05_10_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta05_10_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta05_10_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta10_15_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta10_15_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta10_15_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta15_20_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta15_20_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta15_20_p").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta20_Inf_p", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta20_Inf_p").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta20_Inf_p").c_str(), 150, 0, 1500, 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta00_05_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta00_05_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta00_05_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta05_10_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta05_10_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta05_10_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta10_15_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta10_15_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta10_15_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta15_20_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta15_20_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta15_20_n").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta20_Inf_n", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta20_Inf_n").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR_eta20_Inf_n").c_str(), 150, 0, 1500, 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR", std::make_shared<TH2D>(("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_pt_m_"+std::to_string(i)+"_QCDCR").c_str(), 150, 0, 1500, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_eta_pt_"+std::to_string(i)+"_QCDCR", std::make_shared<TH2D>(("Jet_cm_eta_pt_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_eta_pt_"+std::to_string(i)+"_QCDCR").c_str(), 100, -6, 6, 150, 0, 1500 ));
        my_2d_histos.emplace("Jet_cm_eta_m_"+std::to_string(i)+"_QCDCR", std::make_shared<TH2D>(("Jet_cm_eta_m_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_eta_m_"+std::to_string(i)+"_QCDCR").c_str(), 100, -6, 6, 200, 0, 200 ));
        my_2d_histos.emplace("Jet_cm_phi_pt_"+std::to_string(i)+"_QCDCR", std::make_shared<TH2D>(("Jet_cm_phi_pt_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_phi_pt_"+std::to_string(i)+"_QCDCR").c_str(), 80, -4, 4, 150, 0, 1500 ));
        my_2d_histos.emplace("Jet_cm_phi_m_"+std::to_string(i)+"_QCDCR", std::make_shared<TH2D>(("Jet_cm_phi_m_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_phi_m_"+std::to_string(i)+"_QCDCR").c_str(), 80, -4, 4, 200, 0, 200 ));

        my_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_QCDCR",  std::make_shared<TH1D>(("Jet_cm_pt_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_QCDCR").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_QCDCR", std::make_shared<TH1D>(("Jet_cm_eta_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_QCDCR").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_QCDCR", std::make_shared<TH1D>(("Jet_cm_phi_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_QCDCR").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_QCDCR",   std::make_shared<TH1D>(("Jet_cm_m_"+std::to_string(i)+"_QCDCR").c_str(),("Jet_cm_m_"+std::to_string(i)+"_QCDCR").c_str(), 200, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_deepESM_QCDCR", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_deepESM_QCDCR").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_deepESM_QCDCR").c_str(), 200, 0.0, 1.0, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_deepESM_QCDCR", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_deepESM_QCDCR").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_deepESM_QCDCR").c_str(), 200, 0.0, 1.0, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_deepESM_QCDCR", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_deepESM_QCDCR").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_deepESM_QCDCR").c_str(), 200, 0.0, 1.0, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_deepESM_QCDCR", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_deepESM_QCDCR").c_str(),("Jet_cm_m_"+std::to_string(i)+"_deepESM_QCDCR").c_str(), 200, 0.0, 1.0, 40, 0, 200));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_njets_QCDCR", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_njets_QCDCR").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_njets_QCDCR").c_str(), 20, 0, 20, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_njets_QCDCR", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_njets_QCDCR").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_njets_QCDCR").c_str(), 20, 0, 20, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_njets_QCDCR", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_njets_QCDCR").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_njets_QCDCR").c_str(), 20, 0, 20, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_njets_QCDCR", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_njets_QCDCR").c_str(),("Jet_cm_m_"+std::to_string(i)+"_njets_QCDCR").c_str(), 20, 0, 20, 40, 0, 200));
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

    my_histos.emplace( "h_cutFlow", std::make_shared<TH1D>("h_cutFlow", "h_cutFlow", 9,0,9));    
}

void Analyze1Lep::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    
    std::string suffix = "";
    while( tr.getNextEvent() )
    {
        const auto& MET                       = tr.getVar<double>("MET");
        const auto& METPhi                    = tr.getVar<double>("METPhi");
        const auto& ntops                     = tr.getVar<int>("ntops"+suffix);
        const auto& runtype                   = tr.getVar<std::string>("runtype");     
        const auto& filetag                   = tr.getVar<std::string>("filetag");
        const auto& Jets                      = tr.getVec<TLorentzVector>("Jets"+suffix);
        const auto& GoodJets_pt30             = tr.getVec<bool>("GoodJets_pt30"+suffix);
        const auto& NonIsoMuonJets_pt30       = tr.getVec<bool>("NonIsoMuonJets_pt30"+suffix);
        const auto& NJet                      = tr.getVar<int>("NJets");
        const auto& NGoodJets_pt30            = tr.getVar<int>("NGoodJets_pt30"+suffix);
        const auto& NNonIsoMuonJets_pt30      = tr.getVar<int>("NNonIsoMuonJets_pt30"+suffix);
        const auto& NGoodBJets_pt30           = tr.getVar<int>("NGoodBJets_pt30"+suffix);
        const auto& NGoodMuons                = tr.getVar<int>("NGoodMuons"+suffix);
        const auto& NGoodElectrons            = tr.getVar<int>("NGoodElectrons"+suffix);
        const auto& NGoodLeptons              = tr.getVar<int>("NGoodLeptons"+suffix);
        const auto& NVtx                      = tr.getVar<int>("NVtx");
        const auto& GoodElectrons             = tr.getVec<bool>("GoodElectrons"+suffix);
        const auto& Electrons                 = tr.getVec<TLorentzVector>("Electrons");
        const auto& Muons                     = tr.getVec<TLorentzVector>("Muons");
        const auto& GoodMuons                 = tr.getVec<bool>("GoodMuons"+suffix);
        const auto& Electrons_MiniIso         = tr.getVec<double>("Electrons_MiniIso");
        const auto& Electrons_passIso         = tr.getVec<bool>("Electrons_passIso");
        const auto& Muons_MiniIso             = tr.getVec<double>("Muons_MiniIso");
        const auto& fixedGridRhoFastjetAll    = tr.getVar<double>("fixedGridRhoFastjetAll");
        const auto& GoodLeptons               = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons"+suffix);
        const auto& GoodNonIsoMuons           = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodNonIsoMuons"+suffix);
        const auto& HT_trigger_pt30           = tr.getVar<double>("HT_trigger_pt30"+suffix);
        const auto& HT_NonIsoMuon_pt30        = tr.getVar<double>("HT_NonIsoMuon_pt30"+suffix);
        const auto& JetID                     = tr.getVar<bool>("JetID"+suffix);
        const auto& correct2018Split          = tr.getVar<bool>("correct2018Split"+suffix);
        const auto& passTrigger               = tr.getVar<bool>("passTrigger"+suffix);
        const auto& passTriggerMC             = tr.getVar<bool>("passTriggerMC"+suffix);
        const auto& passMETFilters            = tr.getVar<bool>("passMETFilters"+suffix);
        const auto& passMadHT                 = tr.getVar<bool>("passMadHT"+suffix);
        const auto& passBaseline1l_Good       = tr.getVar<bool>("passBaseline1l_Good"+suffix);
        const auto& passBaseline1e1m          = tr.getVar<bool>("passBaseline1e1m_Good"+suffix);
        const auto& passBaselineGoodOffline1l = tr.getVar<bool>("passBaselineGoodOffline1l"+suffix);
        const auto& passBaseline1l_NonIsoMuon = tr.getVar<bool>("passBaseline1l_NonIsoMuon"+suffix);
        const auto& passHEMVeto               = tr.getVar<bool>("passHEMVeto"+suffix);
        const auto& Mbl                       = tr.getVar<double>("Mbl"+suffix);
        const auto& MblVec                    = tr.getVec<double>("MblVec"+suffix);
        const auto& passBlind                 = tr.getVar<bool>("passBlindLep_Good"+suffix);            
        const auto& deepESM_val               = tr.getVar<double>("deepESM_val"+suffix);
        const auto& deepESM_valNonIsoMuon     = tr.getVar<double>("deepESM_valNonIsoMuon"+suffix);
        const auto& deepESM_bin1              = tr.getVar<bool>("deepESM_bin1"+suffix);
        const auto& deepESM_bin2              = tr.getVar<bool>("deepESM_bin2"+suffix);
        const auto& deepESM_bin3              = tr.getVar<bool>("deepESM_bin3"+suffix);
        const auto& deepESM_bin4              = tr.getVar<bool>("deepESM_bin4"+suffix);
        const auto& deepESM_binNum            = tr.getVar<int>("deepESM_binNum"+suffix);
        const auto& fwm2_top6                 = tr.getVar<double>("fwm2_top6"+suffix);
        const auto& fwm3_top6                 = tr.getVar<double>("fwm3_top6"+suffix);
        const auto& fwm4_top6                 = tr.getVar<double>("fwm4_top6"+suffix);
        const auto& fwm5_top6                 = tr.getVar<double>("fwm5_top6"+suffix);
        const auto& fwm6_top6                 = tr.getVar<double>("fwm6_top6"+suffix);
        const auto& fwm7_top6                 = tr.getVar<double>("fwm7_top6"+suffix);
        const auto& fwm8_top6                 = tr.getVar<double>("fwm8_top6"+suffix);
        const auto& fwm9_top6                 = tr.getVar<double>("fwm9_top6"+suffix);
        const auto& fwm10_top6                = tr.getVar<double>("fwm10_top6"+suffix);
        const auto& jmt_ev0_top6              = tr.getVar<double>("jmt_ev0_top6"+suffix);
        const auto& jmt_ev1_top6              = tr.getVar<double>("jmt_ev1_top6"+suffix);
        const auto& jmt_ev2_top6              = tr.getVar<double>("jmt_ev2_top6"+suffix);
        const auto& Jets_cm_top6              = tr.getVec<TLorentzVector>("Jets_cm_top6"+suffix);
        const auto& eventCounter              = tr.getVar<int>("eventCounter");

        // ------------------------
        // -- Print event number
        // ------------------------       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if(tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() );

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight=1.0, weightNoHT=1.0, weightQCDCR=1.0, weightNoBTag=1.0;
        double pdfUp = 1.0, pdfDown = 1.0, pdfNom = 1.0, weightPDFDown = 1.0, weightPDFUp = 1.0, weightPDFNom = 1.0;
        double htUp = 1.0, htDown = 1.0, weightHTDown = 1.0, weightHTUp = 1.0;
        double scaleUp = 1.0, scaleDown = 1.0, scaleNom = 1.0, weightScaleDown = 1.0, weightScaleUp = 1.0, weightScaleNom = 1.0;
        double eventweight=1.0, leptonweight=1.0, bTagWeight=1.0, prefiringScaleFactor=1.0, pileupWeight=1.0, htDerivedweight=1.0;
        double bTagWeightUp = 1.0, bTagWeightDown = 1.0, weightBtagDown = 1.0, weightBtagUp = 1.0;
        double prefiringScaleFactorUp=1.0, prefiringScaleFactorDown=1.0;
        double weightNoLepton=1.0, weightNoPU=1.0, weightPUdown=1.0, weightPUup=1.0, pileupWeightUp=1.0, pileupWeightDown=1.0;
        double weightNoPrefire=1.0, weightPrefireDown=1.0, weightPrefireUp=1.0;
        double isrUp=1.0, isrDown=1.0, fsrUp=1.0, fsrDown=1.0;
        double weightISRUp=1.0, weightISRDown=1.0, weightFSRUp=1.0, weightFSRDown=1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF"+suffix);
            const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF"+suffix);
            const auto& muNonIso     = tr.getVar<double>("totNonIsoMuonSF"+suffix);
            leptonweight = eleLepWeight*muLepWeight;
            
            pileupWeight = tr.getVar<double>("puWeightCorr"+suffix);
            pileupWeightUp = tr.getVar<double>("puSysUpCorr"+suffix);
            pileupWeightDown = tr.getVar<double>("puSysDownCorr"+suffix);

            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central"+suffix);
            bTagWeightUp = tr.getVar<double>("bTagSF_EventWeightSimple_Up"+suffix);
            bTagWeightDown = tr.getVar<double>("bTagSF_EventWeightSimple_Down"+suffix);

            htDerivedweight = tr.getVar<double>("htDerivedweight"+suffix);
            htUp = tr.getVar<double>("htScaleUp"+suffix);
            htDown = tr.getVar<double>("htScaleDown"+suffix);

            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor"+suffix);
            prefiringScaleFactorUp = tr.getVar<double>("prefiringScaleFactorUp"+suffix);
            prefiringScaleFactorDown = tr.getVar<double>("prefiringScaleFactorDown"+suffix);

            fsrUp = tr.getVar<double>("PSweight_FSRUp_2"+suffix);
            fsrDown = tr.getVar<double>("PSweight_FSRDown_2"+suffix);
            isrUp = tr.getVar<double>("PSweight_ISRUp_2"+suffix);
            isrDown = tr.getVar<double>("PSweight_ISRDown_2"+suffix);
            
            pdfUp = tr.getVar<double>("PDFweightUp"+suffix);
            pdfDown = tr.getVar<double>("PDFweightDown"+suffix);
            pdfNom = tr.getVar<double>("PDFweightNom"+suffix);
 
            scaleUp = tr.getVar<double>("scaleWeightUp"+suffix);
            scaleDown = tr.getVar<double>("scaleWeightDown"+suffix);
            scaleNom = tr.getVar<double>("scaleWeightNom"+suffix);

            weightQCDCR *= eventweight*muNonIso*prefiringScaleFactor*pileupWeight;
            weightNoHT *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight;
            weightNoLepton *= eventweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
            weightNoBTag *= eventweight*leptonweight*prefiringScaleFactor*pileupWeight*htDerivedweight;

            weightNoPU *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*htDerivedweight;
            weightPUdown *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*htDerivedweight*pileupWeightDown;
            weightPUup *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*htDerivedweight*pileupWeightUp;

            weightBtagDown *= eventweight*leptonweight*bTagWeightDown*prefiringScaleFactor*htDerivedweight*pileupWeight;
            weightBtagUp *= eventweight*leptonweight*bTagWeightUp*prefiringScaleFactor*htDerivedweight*pileupWeight;

            weightHTDown *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*htDown*pileupWeight;
            weightHTUp *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*htUp*pileupWeight;

            weightNoPrefire *= eventweight*leptonweight*bTagWeight*pileupWeight*htDerivedweight;
            weightPrefireDown *= eventweight*leptonweight*bTagWeight*prefiringScaleFactorDown*htDerivedweight*pileupWeight;
            weightPrefireUp *= eventweight*leptonweight*bTagWeight*prefiringScaleFactorUp*htDerivedweight*pileupWeight;

            weightISRUp *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*isrUp;
            weightISRDown *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*isrDown;
            weightFSRUp *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*fsrUp;
            weightFSRDown *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*fsrDown;

            weightPDFDown *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*pdfDown;
            weightPDFUp *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*pdfUp;
            weightPDFNom *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*pdfNom;

            weightScaleDown *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*scaleDown;
            weightScaleUp *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*scaleUp;
            weightScaleNom *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*scaleNom;

            weight *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
        }

        int NGenJets = 0, NGenJets_pt30 = 0;
        if(runtype == "MC")
        {
            //Get the numer of GenJets 
            const auto& GenJets = tr.getVec<TLorentzVector>("GenJets");           
            for(const auto& lv : GenJets)
            {
                NGenJets++;
                if(lv.Pt() > 30)
                {
                    NGenJets_pt30++;
                }
            }
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------
        bool pass_general    = passTriggerMC && passTrigger && passMadHT && passBlind && passMETFilters && passHEMVeto && correct2018Split;
        //bool pass_0l         = NGoodLeptons == 0;
        bool pass_1l         = NGoodLeptons == 1;
        bool pass_ht         = HT_trigger_pt30 > 300;
        bool pass_MBL        = (50 < Mbl && Mbl < 250);
        //bool pass_1e_1m      = (NGoodLeptons == 2) ? GoodLeptons[0].first != GoodLeptons[1].first : false;
        bool pass_njet_pt30  = NGoodJets_pt30 >= 7;
        bool pass_1btag_pt30 = NGoodBJets_pt30 >= 1;
        bool pass_2btag_pt30 = NGoodBJets_pt30 >= 2;
        bool pass_0btag_pt30 = NGoodBJets_pt30 == 0;
        bool pass_1e         = false;
        bool pass_1m         = false;
        bool pass_lBarrel    = false;
        if(pass_1l)
        { 
            if(GoodLeptons[0].first == "e") pass_1e = true;
            if(GoodLeptons[0].first == "m") pass_1m = true;
            pass_lBarrel = abs( GoodLeptons[0].second.Eta() ) <= 1.2;
        }
        //bool pass_5to6njet_pt30 = (NGoodJets_pt30 == 5 || NGoodJets_pt30 == 6);
        
        bool passBaseline1l_AllJets = passBaselineGoodOffline1l &&
                                      passTrigger               &&
                                      passTriggerMC             &&
                                      passBlind                 &&
            (
                ((runtype != "Data" || filetag.find("Data_SingleMuon")     != std::string::npos) && NGoodMuons == 1     && NGoodElectrons == 0)
                                                                                     ||
                ((runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) && NGoodElectrons == 1 && NGoodMuons == 0)
            );
        
        // ------------------------------------------------
        // --  Temporary Home of the W+Jets Control Region
        // ------------------------------------------------
        double mT = -1;
        bool pass_mT = false, pass_Mbl_all = false;
        if(pass_1l)
        {
            TLorentzVector metLV;
            metLV.SetPtEtaPhiE(MET,0,METPhi,MET);
            mT = utility::calcMT(GoodLeptons[0].second, metLV);
            pass_mT = mT > 50 && mT < 110;
        
            for(unsigned int i = 0; i < Jets.size(); ++i)
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

        const auto& NGoodBJets_pt30_loose = tr.getVar<int>("NGoodBJets_pt30_loose"+suffix);        
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
            //{""                                          , pass_general                                                                             },
            //{"_HT300"                                    , pass_general && pass_ht && JetID                                                         },
            //{"_HT300_0b"                                 , pass_general && pass_ht && JetID && pass_0btag_pt30                                      },
            //{"_HT300_1l"                                 , pass_general && pass_ht && JetID && pass_1l                                              },
            //{"_HT300_1l_0b"                              , pass_general && pass_ht && JetID && pass_1l && pass_0btag_pt30                           },
            //{"_HT300_1l_ge2b"                            , pass_general && pass_ht && JetID && pass_1l && pass_2btag_pt30                           },
            //{"_HT300_1l_ge1b_ge3j"                       , pass_general && pass_ht && JetID && pass_1l && pass_1btag_pt30 && NGoodJets_pt30 >= 3    },
            //{"_1l_HT300_ge7j"                            , pass_general && pass_1l && pass_ht && pass_njet_pt30  && JetID                           },
            //{"_1l_HT300_ge1b"                            , pass_general && pass_1l && pass_ht && pass_1btag_pt30 && JetID                           },
            //{"_1l_HT300_ge7j_ge1b"                       , pass_general && pass_1l && pass_ht && pass_njet_pt30  && pass_1btag_pt30 && JetID        },
            //{"_1l_HT300_ge1b_Mbl"                        , pass_general && passBaseline1l_AllJets                                                   },
            //{"_1l_HT300_ge3j_ge1b_Mbl"                   , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 3                            },
            //{"_1l_HT300_ge3j_ge1b_Mbl_noPuWeight"        , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 3                            },
            //{"_1l_HT300_ge3j_ge1b_Mbl_puWeightUp"        , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 3                            },
            //{"_1l_HT300_ge3j_ge1b_Mbl_puWeightDown"      , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 3                            },
            //{"_1e_HT300_ge4j_ge1b_Mbl"                   , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 4 && pass_1e                 },
            //{"_1m_HT300_ge4j_ge1b_Mbl"                   , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 4 && pass_1m                 },
            {"_1l_HT300_ge7j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good                                                      },                         
            {"_1e_HT300_ge7j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && pass_1e                                           },                         
            {"_1m_HT300_ge7j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && pass_1m                                           },   
            {"_1e_HT300_ge7j_ge1b_Mbl_noHTWeight"        , pass_general && passBaseline1l_Good && pass_1e                                           },                         
            {"_1m_HT300_ge7j_ge1b_Mbl_noHTWeight"        , pass_general && passBaseline1l_Good && pass_1m                                           },   
            //{"_1l_HT300_ge7j_ge1b_Mbl_posWeight"         , pass_general && passBaseline1l_Good && weight >  0.0                                     },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_negWeight"         , pass_general && passBaseline1l_Good && weight <= 0.0                                     },                         
            {"_1l_HT300_ge7j_ge1b_Mbl_noHTWeight"        , pass_general && passBaseline1l_Good                                                      }, 
            //{"_1l_HT300_ge7j_ge1b_Mbl_noLepWeight"       , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_noBTagWeight"      , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_noPuWeight"        , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_puWeightUp"        , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_puWeightDown"      , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_bTagWeightUp"        , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_bTagWeightDown"      , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_htWeightUp"        , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_htWeightDown"      , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_prefireWeightUp"   , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_prefireWeightDown" , pass_general && passBaseline1l_Good                                                      },                         

            //{"_1l_HT300_ge7j_ge1b_Mbl_pdfWeightUp"   , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_pdfWeightDown" , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_pdfWeightNom" , pass_general && passBaseline1l_Good                                                      },                         

            //{"_1l_HT300_ge7j_ge1b_Mbl_scaleWeightUp"   , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_scaleWeightDown" , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_scaleWeightNom" , pass_general && passBaseline1l_Good                                                      },                         

            //{"_1l_HT300_ge7j_ge1b_Mbl_fsrDown"           , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_fsrUp"             , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_isrDown"           , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_isrUp"             , pass_general && passBaseline1l_Good                                                      },                         
            //{"_1l_HT300_ge7j_ge2b_Mbl"                   , pass_general && passBaseline1l_Good && pass_2btag_pt30                                   },
            //{"_1l_HT300_ge7j_ge1b_Mbl_lBarrel"           , pass_general && passBaseline1l_Good && pass_lBarrel                                      },
            //{"_1e_HT300_ge7j_ge1b_Mbl_lBarrel"           , pass_general && passBaseline1l_Good && pass_lBarrel && pass_1e                           },
            //{"_1m_HT300_ge7j_ge1b_Mbl_lBarrel"           , pass_general && passBaseline1l_Good && pass_lBarrel && pass_1m                           },
            //{"_1l_HT300_ge7j_ge1b_Mbl_lEndCap"           , pass_general && passBaseline1l_Good && !pass_lBarrel                                     },
            //{"_1e_HT300_ge7j_ge1b_Mbl_lEndCap"           , pass_general && passBaseline1l_Good && !pass_lBarrel && pass_1e                          },
            //{"_1m_HT300_ge7j_ge1b_Mbl_lEndCap"           , pass_general && passBaseline1l_Good && !pass_lBarrel && pass_1m                          },
            //{"_1l_HT300_ge7j_ge1b_Mbl_d1_fsrDown"          , pass_general && passBaseline1l_Good && deepESM_bin1                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_d2_fsrUp"            , pass_general && passBaseline1l_Good && deepESM_bin2                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_d3_isrDown"          , pass_general && passBaseline1l_Good && deepESM_bin3                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_d4_isrUp"            , pass_general && passBaseline1l_Good && deepESM_bin4                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_d1"                  , pass_general && passBaseline1l_Good && deepESM_bin1                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_d2"                  , pass_general && passBaseline1l_Good && deepESM_bin2                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_d3"                  , pass_general && passBaseline1l_Good && deepESM_bin3                                      },                         
            //{"_1l_HT300_ge7j_ge1b_Mbl_d4"                  , pass_general && passBaseline1l_Good && deepESM_bin4                                      },                         

            //{"_1l_HT300_1j_ge1b_Mbl"                     , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 1                            },
            //{"_1l_HT300_2j_ge1b_Mbl"                     , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 2                            },
            //{"_1l_HT300_3j_ge1b_Mbl"                     , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 3                            },
            //{"_1l_HT300_4j_ge1b_Mbl"                     , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 4                            },
            //{"_1l_HT300_5j_ge1b_Mbl"                     , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 5                            },
            //{"_1l_HT300_6j_ge1b_Mbl"                     , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 6                            },
            //{"_1l_HT300_7j_ge1b_Mbl"                     , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 7                               },
            //{"_1l_HT300_8j_ge1b_Mbl"                     , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 8                               },
            //{"_1l_HT300_9j_ge1b_Mbl"                     , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 9                               },
            //{"_1l_HT300_10j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 10                              },
            //{"_1l_HT300_11j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 11                              },
            //{"_1l_HT300_12j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 12                              },
            //{"_1l_HT300_13j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 13                              },
            //{"_1l_HT300_14j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 14                              },
            //{"_1l_HT300_15j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 15                              },
            //{"_1l_HT300_n7j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 7        },
            //{"_1l_HT300_n8j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 8        },
            //{"_1l_HT300_n9j_ge1b_Mbl"                    , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 9        },
            //{"_1l_HT300_n10j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 10       },
            //{"_1l_HT300_n11j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 11       },
            //{"_1l_HT300_n12j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 12       },
            //{"_1l_HT300_n13j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 13       },
            //{"_1l_HT300_n14j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 14       },
            //{"_1l_HT300_n15j_ge1b_Mbl"                   , pass_general && passBaseline1l_Good && NGoodJets_pt30 != 15       },

            //{"_1l_HT300_5j_ge1b_Mbl_htCorr"              , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 5                            },
            //{"_1l_HT300_6j_ge1b_Mbl_htCorr"              , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 6                            },
            //{"_1l_HT300_7j_ge1b_Mbl_htCorr"              , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 7                               },
            //{"_1l_HT300_8j_ge1b_Mbl_htCorr"              , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 8                               },
            //{"_1l_HT300_ge7j_ge1b_Mbl_htCorr"            , pass_general && passBaseline1l_Good                                                      },
            //{"_1l_0b_ge300ht_50to110mt_ge30MET"          , pass_general && passBaseline1l_WCR                                                       },
            //{"_1l_0b_ge300ht_50to110mt_ge30MET_even"     , pass_general && passBaseline1l_WCR && evenEvent                                          },
            //{"_1l_0b_ge300ht_50to110mt_ge30MET_odd"      , pass_general && passBaseline1l_WCR && !evenEvent                                         },
            //{"_1e_1m_ge2b_le5j"                          , pass_general && passBaseline1e1m                                                         },
            {"_passQCDCR"                                , passBaseline1l_NonIsoMuon                                                                },
        };

        std::vector<TH1DInfo> histInfos = {
            {    "h_njets",               20,   0.0,   20.0},
            {"blind_njets",               20,   0.0,   20.0},
            {    "h_njetsQCDCR",          20,   0.0,   20.0},
            {    "h_ngjets",              20,   0.0,   20.0},
            {    "h_ngjets_pt30",         20,   0.0,   20.0},
            {    "h_ntops",               10,   0.0,   10.0},
            {"blind_ntops",               10,   0.0,   10.0},
            {    "h_nb",                  10,   0.0,   10.0},
            {"blind_nb",                  10,   0.0,   10.0},
            {    "h_deepESM",            200,   0.0,    1.0},
            {    "h_mva_all",            200,   0.0,    1.0},
            {    "h_mva_all_w2",         200,   0.0,    1.0},
            {    "h_mva_all_nEntries",   200,   0.0,    1.0},
            {    "h_mva_7j",             200,   0.0,    1.0},
            {    "h_mva_8j",             200,   0.0,    1.0},
            {    "h_mva_9j",             200,   0.0,    1.0},
            {    "h_mva_10j",            200,   0.0,    1.0},
            {    "h_mva_11j",            200,   0.0,    1.0},
            {    "h_mva_7j_w2",             200,   0.0,    1.0},
            {    "h_mva_8j_w2",             200,   0.0,    1.0},
            {    "h_mva_9j_w2",             200,   0.0,    1.0},
            {    "h_mva_10j_w2",            200,   0.0,    1.0},
            {    "h_mva_11j_w2",            200,   0.0,    1.0},
            {    "h_mva_7j_nEntries",             200,   0.0,    1.0},
            {    "h_mva_8j_nEntries",             200,   0.0,    1.0},
            {    "h_mva_9j_nEntries",             200,   0.0,    1.0},
            {    "h_mva_10j_nEntries",            200,   0.0,    1.0},
            {    "h_mva_11j_nEntries",            200,   0.0,    1.0},
            {"blind_deepESM",            200,   0.0,    1.0},
            {    "h_deepESMQCDCR",       200,   0.0,    1.0},
            {    "h_deepESMMerged",        4,   0.5,    4.5},
            {"blind_deepESMMerged",        4,   0.5,    4.5},
            {    "h_ht",                 500,   0.0, 5000.0},
            {    "h_htQCDCR",            500,   0.0, 5000.0},
            {"blind_ht",                 500,   0.0, 5000.0},
            {    "h_mbl",                300,   0.0,  300.0},
            {"blind_mbl",                300,   0.0,  300.0},
            {    "h_lPt",                200,   0.0, 2000.0},
            {"blind_lPt",                200,   0.0, 2000.0},
            {    "h_lEta",               200,  -6.0,    6.0},
            {"blind_lEta",               200,  -6.0,    6.0},
            {    "h_lPhi",               200,  -4.0,    4.0},
            {"blind_lPhi",               200,  -4.0,    4.0},
            {    "h_lMiniIso",          4000,   0.0,    5.0},
            {    "h_isomPt",             200,   0.0, 2000.0},
            {    "h_isomEta",            200,  -6.0,    6.0},
            {    "h_isomPhi",            200,  -4.0,    4.0},
            {    "h_jPt",                150,   0.0, 1500.0},
            {    "h_jPt_j1",             150,   0.0, 1500.0},
            {    "h_jPt_j2",             150,   0.0, 1500.0},
            {    "h_jPt_j3",             150,   0.0, 1500.0},
            {    "h_jPt_j4",             150,   0.0, 1500.0},
            {    "h_jPt_j5",             150,   0.0, 1500.0},
            {    "h_jPt_j6",             150,   0.0, 1500.0},
            {    "h_jPt_j7",             150,   0.0, 1500.0},
            {    "h_jPt_j8",             150,   0.0, 1500.0},
            {    "h_jPt_j9",             150,   0.0, 1500.0},
            {    "h_jPt_j10",             150,   0.0, 1500.0},
            {"blind_jPt",                150,   0.0, 1500.0},
            {    "h_jEta",               100,  -6.0,    6.0},
            {    "h_jEta_j1",               100,  -6.0,    6.0},
            {    "h_jEta_j2",               100,  -6.0,    6.0},
            {    "h_jEta_j3",               100,  -6.0,    6.0},
            {    "h_jEta_j4",               100,  -6.0,    6.0},
            {    "h_jEta_j5",               100,  -6.0,    6.0},
            {    "h_jEta_j6",               100,  -6.0,    6.0},
            {    "h_jEta_j7",               100,  -6.0,    6.0},
            {    "h_jEta_j8",               100,  -6.0,    6.0},
            {    "h_jEta_j9",               100,  -6.0,    6.0},
            {    "h_jEta_j10",               100,  -6.0,    6.0},
            {"blind_jEta",               100,  -6.0,    6.0},
            {    "h_jPhi",               80,  -4.0,    4.0},
            {    "h_jPhi_j1",               80,  -4.0,    4.0},
            {    "h_jPhi_j2",               80,  -4.0,    4.0},
            {    "h_jPhi_j3",               80,  -4.0,    4.0},
            {    "h_jPhi_j4",               80,  -4.0,    4.0},
            {    "h_jPhi_j5",               80,  -4.0,    4.0},
            {    "h_jPhi_j6",               80,  -4.0,    4.0},
            {    "h_jPhi_j7",               80,  -4.0,    4.0},
            {    "h_jPhi_j8",               80,  -4.0,    4.0},
            {    "h_jPhi_j9",               80,  -4.0,    4.0},
            {    "h_jPhi_j10",               80,  -4.0,    4.0},
            {"blind_jPhi",               80,  -4.0,    4.0},
            {    "h_jM",                  200,  0.0,    200.0},
            {    "h_jM_j1",               200,  0.0,    200.0},
            {    "h_jM_j2",               200,  0.0,    200.0},
            {    "h_jM_j3",               200,  0.0,    200.0},
            {    "h_jM_j4",               200,  0.0,    200.0},
            {    "h_jM_j5",               200,  0.0,    200.0},
            {    "h_jM_j6",               200,  0.0,    200.0},
            {    "h_jM_j7",               200,  0.0,    200.0},
            {    "h_jM_j8",               200,  0.0,    200.0},
            {    "h_jM_j9",               200,  0.0,    200.0},
            {    "h_jM_j10",               200,  0.0,    200.0},
            {"blind_jM",                200,  0.0,    200.0},
            {    "h_allMbl",             300,   0.0,  300.0},            
            {"blind_allMbl",             300,   0.0,  300.0},
            {"h_nvtx",                   101,  -0.5,  100.5},
            {"h_trueNumInteractions",    720,   0.0,  120.0},
            {"h_fixedGridRhoFastjetAll", 720,   0.0,  120.0},
            {"blind_fixedGridRhoFastjetAll", 720,   0.0,  120.0},
            {"h_weight",                 200,  -5.0,    5.0},
            {"h_leptonweight",           200,  -5.0,    5.0},
            {"h_pileupWeight",           200,  -5.0,   20.0},
            {"h_bTagWeight",             200,  -5.0,   20.0},
            {"h_htDerivedweight",        200,  -5.0,    5.0},
            {"h_prefiringScaleFactor",   200,  -5.0,    5.0},
            {"h_eventweight",            200,  -5.0,   50.0},            
        };

        std::vector<TH2DInfo> hist2DInfos = {
            {    "h_jPt_jM",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j1_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j2_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j3_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j4_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j5_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j6_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j7_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j8_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j9_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta00_05_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta05_10_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta10_15_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta15_20_p",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta20_Inf_p",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta00_05_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta05_10_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta10_15_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta15_20_n",       150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jPt_jM_j10_eta20_Inf_n",      150, 0.0, 1500.0, 200, 0.0, 200},
            {    "h_jEta_jPt_j1",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j2",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j3",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j4",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j5",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j6",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j7",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j8",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j9",       100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jEta_jPt_j10",      100, -6, 6,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j1",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j2",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j3",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j4",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j5",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j6",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j7",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j8",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j9",       80, -4, 4,  150,   0.0, 1500.0},
            {    "h_jPhi_jPt_j10",       80,-4, 4,  150,   0.0, 1500.0},
            {    "h_jEta_jM_j1",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j2",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j3",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j4",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j5",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j6",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j7",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j8",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j9",       100, -6, 6, 200, 0.0, 200},
            {    "h_jEta_jM_j10",      100, -6, 6, 200, 0.0, 200},
            {    "h_jPhi_jM_j1",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j2",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j3",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j4",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j5",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j6",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j7",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j8",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j9",       80, -4, 4,  200, 0.0, 200},
            {    "h_jPhi_jM_j10",       80, -4, 4,  200, 0.0, 200},
            {    "h_njets_deepESM", 15,    0,   15, 100,   0.0,   1.0},
            {"blind_njets_deepESM", 15,    0,   15, 100,   0.0,   1.0},
            {    "h_njets_mbl",     15,    0,   15, 100,   0.0, 300.0},
            {"blind_njets_mbl",     15,    0,   15, 100,   0.0, 300.0},
            {    "h_ht_deepESM",   100,    0, 3000, 100,   0.0,   1.0},
            {"blind_ht_deepESM",   100,    0, 3000, 100,   0.0,   1.0},
            {    "h_lEta_lPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {"blind_lEta_lPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {    "h_jEta_jPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {"blind_jEta_jPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {    "h_lEta_nb",      100, -6.0,  6.0,  10,   0.0,  10.0},
            {"blind_lEta_nb",      100, -6.0,  6.0,  10,   0.0,  10.0},
            {    "h_njets_rho",    100, 0.0, 100.0, 15,  -0.5,  15.5},
            {"blind_njets_rho",    100, 0.0, 100.0, 15,  -0.5,  15.5},
            {"h_lPt_lMiniIso",     200, 0.0,  2000.0, 4000, 0.0, 5.0},
            {"h_lEta_lMiniIso",    200, -6.0, 6.0,    4000, 0.0, 5.0},
            {"h_lPhi_lMiniIso",    200, -4.0, 4.0,    4000, 0.0, 5.0},

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

        my_histos["EventCounter"]->Fill(eventCounter);

        for(auto& kv : cut_map_1l)
        {
            if(kv.second)
            {
                double w = weight;
                if(kv.first.find("50to110mt")    != std::string::npos || 
                   kv.first.find("htCorr")       != std::string::npos || 
                   kv.first.find("noHTWeight")   != std::string::npos) w = weightNoHT;
                if(kv.first.find("passQCDCR")    != std::string::npos) w = weightQCDCR;
                if(kv.first.find("noLepWeight")  != std::string::npos) w = weightNoLepton;
                if(kv.first.find("noBTagWeight") != std::string::npos) w = weightNoBTag;
                if(kv.first.find("noPuWeight")   != std::string::npos) w = weightNoPU;
                if(kv.first.find("puWeightUp")   != std::string::npos) w = weightPUup;
                if(kv.first.find("puWeightDown") != std::string::npos) w = weightPUdown;
                if(kv.first.find("bTagWeightUp")   != std::string::npos) w = weightBtagUp;
                if(kv.first.find("bTagWeightDown") != std::string::npos) w = weightBtagDown;

                if(kv.first.find("noPrefireWeight")   != std::string::npos) w = weightNoPrefire;
                if(kv.first.find("prefireWeightUp")   != std::string::npos) w = weightPrefireUp;
                if(kv.first.find("prefireWeightDown") != std::string::npos) w = weightPrefireDown;
                if(kv.first.find("pdfWeightUp")   != std::string::npos) w = weightPDFUp;
                if(kv.first.find("pdfWeightDown") != std::string::npos) w = weightPDFDown;
                if(kv.first.find("htWeightUp")   != std::string::npos) w = weightHTUp;
                if(kv.first.find("htWeightDown") != std::string::npos) w = weightHTDown;
                if(kv.first.find("pdfWeightNom") != std::string::npos) w = weightPDFNom;

                if(kv.first.find("scaleWeightUp")   != std::string::npos) w = weightScaleUp;
                if(kv.first.find("scaleWeightDown") != std::string::npos) w = weightScaleDown;
                if(kv.first.find("scaleWeightNom") != std::string::npos) w = weightScaleNom;


                if(kv.first.find("isrUp")   != std::string::npos) w = weightISRUp;
                if(kv.first.find("isrDown") != std::string::npos) w = weightISRDown;
                if(kv.first.find("fsrUp")   != std::string::npos) w = weightFSRUp;
                if(kv.first.find("fsrDown") != std::string::npos) w = weightFSRDown;

                if (NGoodJets_pt30 >= 7) {
                    my_histos["h_mva_all"+kv.first]->Fill(deepESM_val, w);
                    my_histos["h_mva_all_w2"+kv.first]->Fill(deepESM_val, std::pow(w,2.0));
                    my_histos["h_mva_all_nEntries"+kv.first]->Fill(deepESM_val, 1);
                }

                if (NGoodJets_pt30 >= 7 && NGoodJets_pt30 < 11) {
                    my_histos["h_mva_"+std::to_string(NGoodJets_pt30)+"j"+kv.first]->Fill(deepESM_val, w);
                    my_histos["h_mva_"+std::to_string(NGoodJets_pt30)+"j_w2"+kv.first]->Fill(deepESM_val, std::pow(w,2.0));
                    my_histos["h_mva_"+std::to_string(NGoodJets_pt30)+"j_nEntries"+kv.first]->Fill(deepESM_val, 1);

                } else if (NGoodJets_pt30 >= 11) {
                    my_histos["h_mva_11j"+kv.first]->Fill(deepESM_val, w);
                    my_histos["h_mva_11j_w2"+kv.first]->Fill(deepESM_val, std::pow(w,2.0));
                    my_histos["h_mva_11j_nEntries"+kv.first]->Fill(deepESM_val, 1);
                }

                my_histos["h_njets"               +kv.first]->Fill(NGoodJets_pt30, w);
                my_histos["h_ngjets"              +kv.first]->Fill(NGenJets, eventweight);
                my_histos["h_ngjets_pt30"         +kv.first]->Fill(NGenJets_pt30, eventweight);
                my_histos["h_njetsQCDCR"          +kv.first]->Fill(NNonIsoMuonJets_pt30, w);
                my_histos["h_ntops"               +kv.first]->Fill(ntops, w);
                my_histos["h_nb"                  +kv.first]->Fill(NGoodBJets_pt30, w);
                my_histos["h_deepESM"             +kv.first]->Fill(deepESM_val, w);
                my_histos["h_deepESMQCDCR"        +kv.first]->Fill(deepESM_valNonIsoMuon, w);
                my_histos["h_deepESMMerged"       +kv.first]->Fill(deepESM_binNum, w);
                my_histos["h_ht"                  +kv.first]->Fill(HT_trigger_pt30, w);
                my_histos["h_htQCDCR"             +kv.first]->Fill(HT_NonIsoMuon_pt30, w);
                my_histos["h_mbl"                 +kv.first]->Fill(Mbl, w);
                my_histos["h_weight"              +kv.first]->Fill(weight, w);
                my_histos["h_leptonweight"        +kv.first]->Fill(leptonweight, w);
                my_histos["h_pileupWeight"        +kv.first]->Fill(pileupWeight, w);
                my_histos["h_bTagWeight"          +kv.first]->Fill(bTagWeight, w);
                my_histos["h_htDerivedweight"     +kv.first]->Fill(htDerivedweight, w);
                my_histos["h_prefiringScaleFactor"+kv.first]->Fill(prefiringScaleFactor, w);
                my_histos["h_eventweight"         +kv.first]->Fill(eventweight, w);
                my_histos["h_nvtx"                +kv.first]->Fill(NVtx, w);
                my_histos["h_fixedGridRhoFastjetAll"+kv.first]->Fill(fixedGridRhoFastjetAll, w);
                my_2d_histos["h_njets_rho"+kv.first]->Fill(fixedGridRhoFastjetAll, NGoodJets_pt30, w);
                for(unsigned int l = 0; l < GoodElectrons.size(); l++)
                {
                    if (!GoodElectrons[l]) continue;
                    my_histos["h_lMiniIso"+kv.first]->Fill(Electrons_MiniIso.at(l), w);
                    my_2d_histos["h_lPt_lMiniIso"+kv.first]->Fill(Electrons.at(l).Pt(), Electrons_MiniIso.at(l), w);
                    my_2d_histos["h_lPhi_lMiniIso"+kv.first]->Fill(Electrons.at(l).Phi(), Electrons_MiniIso.at(l), w);
                    my_2d_histos["h_lEta_lMiniIso"+kv.first]->Fill(Electrons.at(l).Eta(), Electrons_MiniIso.at(l), w);

                }
                for(unsigned int l = 0; l < GoodMuons.size(); l++)
                {
                    if (!GoodMuons[l]) continue;
                    my_histos["h_lMiniIso"+kv.first]->Fill(Muons_MiniIso.at(l), w);
                    my_2d_histos["h_lPt_lMiniIso"+kv.first]->Fill(Muons.at(l).Pt(), Muons_MiniIso.at(l), w);
                    my_2d_histos["h_lPhi_lMiniIso"+kv.first]->Fill(Muons.at(l).Phi(), Muons_MiniIso.at(l), w);
                    my_2d_histos["h_lEta_lMiniIso"+kv.first]->Fill(Muons.at(l).Eta(), Muons_MiniIso.at(l), w);

                }

                for(const auto& l : GoodLeptons)
                {
                    my_histos["h_lPt"+kv.first]->Fill(l.second.Pt(), w);
                    my_histos["h_lEta"+kv.first]->Fill(l.second.Eta(), w);
                    my_histos["h_lPhi"+kv.first]->Fill(l.second.Phi(), w);
                    my_2d_histos["h_lEta_lPhi"+kv.first]->Fill(l.second.Eta(), l.second.Phi(), w);
                    my_2d_histos["h_lEta_nb"+kv.first]->Fill(l.second.Eta(), NGoodBJets_pt30, w);
                }
                for(const auto& isoMuon : GoodNonIsoMuons)
                {
                    my_histos["h_isomPt"+kv.first]->Fill(isoMuon.second.Pt(), w);
                    my_histos["h_isomEta"+kv.first]->Fill(isoMuon.second.Eta(), w);
                    my_histos["h_isomPhi"+kv.first]->Fill(isoMuon.second.Phi(), w);
                }
                for(const auto& mbl : MblVec)
                {
                    my_histos["h_allMbl"+kv.first]->Fill(mbl, w);
                }

                // Pt ranked jets are for good jets only, keep track of good jet rank
                unsigned int goodJetRank = 1;
                double jeta = 0.0;
                for(unsigned int j = 0; j < Jets.size(); j++)
                {

                    if (!GoodJets_pt30[j]) continue;

                    if (goodJetRank < 10) {
                        jeta = Jets.at(j).Eta();
                        my_histos["h_jPt_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Pt(),w);
                        my_histos["h_jEta_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Eta(),w);
                        my_histos["h_jPhi_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Phi(),w);
                        my_histos["h_jM_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).M(),w);

                        my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Pt(),Jets.at(j).M(),w);
                        my_2d_histos["h_jEta_jPt_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Eta(),Jets.at(j).Pt(),w);
                        my_2d_histos["h_jEta_jM_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Eta(),Jets.at(j).M(),w);
                        my_2d_histos["h_jPhi_jPt_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Phi(),Jets.at(j).Pt(),w);
                        my_2d_histos["h_jPhi_jM_j"+std::to_string(goodJetRank)+kv.first]->Fill(Jets.at(j).Phi(),Jets.at(j).M(),w);

                        if      (jeta > 0.0 and jeta <= 0.5)   my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta00_05_p"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta > 0.5 and jeta <= 1.0)   my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta05_10_p"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta > 1.0 and jeta <= 1.5)   my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta10_15_p"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta > 1.5 and jeta <= 2.0)   my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta15_20_p"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta > 2.0)                   my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta20_Inf_p"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta <  0.0 and jeta >= -0.5) my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta00_05_n"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta < -0.5 and jeta >= -1.0) my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta05_10_n"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta < -1.0 and jeta >= -1.5) my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta10_15_n"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta < -1.5 and jeta >= -2.0) my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta15_20_n"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);
                        else if (jeta < -2.0)                  my_2d_histos["h_jPt_jM_j"+std::to_string(goodJetRank)+"_eta20_Inf_n"+kv.first]->Fill(Jets.at(j).Pt(), Jets.at(j).M(), w);

                    }

                    my_2d_histos["h_jPt_jM"+kv.first]->Fill(Jets.at(j).Pt(),Jets.at(j).M(),w);
                    my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                    my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                    my_histos["h_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                    my_2d_histos["h_jEta_jPhi"+kv.first]->Fill(Jets.at(j).Eta(), Jets.at(j).Phi(), w);
                    
                    goodJetRank++;
                }
                my_2d_histos["h_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, w);
                my_2d_histos["h_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, w);
                my_2d_histos["h_ht_deepESM"+kv.first]->Fill(HT_trigger_pt30, deepESM_val, w);
                my_2d_tp_histos["h_njets_deepESMMerged_preFireSF"+kv.first]->Fill(NJet, deepESM_binNum, prefiringScaleFactor, w);

                if( NGoodJets_pt30 <= 8 )
                {
                    my_histos["blind_njets"         +kv.first]->Fill(NGoodJets_pt30, w);
                    my_histos["blind_ntops"         +kv.first]->Fill(ntops, w);
                    my_histos["blind_nb"            +kv.first]->Fill(NGoodBJets_pt30, w);
                    my_histos["blind_deepESM"       +kv.first]->Fill(deepESM_val, w);
                    my_histos["blind_deepESMMerged" +kv.first]->Fill(deepESM_binNum, w);
                    my_histos["blind_ht"            +kv.first]->Fill(HT_trigger_pt30, w);
                    my_histos["blind_mbl"           +kv.first]->Fill(Mbl, w);
                    my_histos["blind_fixedGridRhoFastjetAll"+kv.first]->Fill(fixedGridRhoFastjetAll, w);
                    my_2d_histos["blind_njets_rho"+kv.first]->Fill(fixedGridRhoFastjetAll, NGoodJets_pt30, w);

                    for(const auto l : GoodLeptons)
                    {
                        my_histos["blind_lPt"+kv.first]->Fill(l.second.Pt(), w);
                        my_histos["blind_lEta"+kv.first]->Fill(l.second.Eta(), w);
                        my_histos["blind_lPhi"+kv.first]->Fill(l.second.Phi(), w);
                        my_2d_histos["blind_lEta_lPhi"+kv.first]->Fill(l.second.Eta(), l.second.Phi(), w);
                        my_2d_histos["blind_lEta_nb"+kv.first]->Fill(l.second.Eta(), NGoodBJets_pt30, w);
                    }
                    for(const auto& mbl : MblVec)
                    {
                        my_histos["blind_allMbl"+kv.first]->Fill(mbl, w);
                    }
                    for(unsigned int j = 0; j < Jets.size(); j++)
                    {
                        if(!GoodJets_pt30[j]) continue;
                        my_histos["blind_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                        my_histos["blind_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                        my_histos["blind_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                        my_2d_histos["blind_jEta_jPhi"+kv.first]->Fill(Jets.at(j).Eta(), Jets.at(j).Phi(), w);
                    }
                    my_2d_histos["blind_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, w);
                    my_2d_histos["blind_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, w);
                    my_2d_histos["blind_ht_deepESM"+kv.first]->Fill(HT_trigger_pt30, deepESM_val, w);
                }
            }
        }

        // Fill histos for deepESM training
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
            my_histos["GoodLeptons_pt_1_1l_ge7j_ge1b"]->Fill(GoodLeptons.at(0).second.Pt(), weight);
            my_histos["GoodLeptons_eta_1_1l_ge7j_ge1b"]->Fill(GoodLeptons.at(0).second.Eta(), weight);
            my_histos["GoodLeptons_phi_1_1l_ge7j_ge1b"]->Fill(GoodLeptons.at(0).second.Phi(), weight);
            my_histos["GoodLeptons_m_1_1l_ge7j_ge1b"]->Fill(GoodLeptons.at(0).second.M(), weight);

            for(unsigned int i = 0; i < Jets_cm_top6.size(); i++)
            {
                double jeta = Jets_cm_top6.at(i).Eta();
                if      (jeta > 0.0 and jeta <= 0.5)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta00_05_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta > 0.5 and jeta <= 1.0)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta05_10_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta > 1.0 and jeta <= 1.5)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta10_15_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta > 1.5 and jeta <= 2.0)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta15_20_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta > 2.0)                   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta20_Inf_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta <  0.0 and jeta >= -0.5) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta00_05_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta < -0.5 and jeta >= -1.0) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta05_10_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta < -1.0 and jeta >= -1.5) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta10_15_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta < -1.5 and jeta >= -2.0) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta15_20_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                else if (jeta < -2.0)                  my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_eta20_Inf_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);

                my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weight);
                my_2d_histos["Jet_cm_eta_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(Jets_cm_top6.at(i).Eta(), Jets_cm_top6.at(i).Pt(), weight);
                my_2d_histos["Jet_cm_phi_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(Jets_cm_top6.at(i).Phi(), Jets_cm_top6.at(i).Pt(), weight);
                my_2d_histos["Jet_cm_eta_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(Jets_cm_top6.at(i).Eta(), Jets_cm_top6.at(i).M(), weight);
                my_2d_histos["Jet_cm_phi_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(Jets_cm_top6.at(i).Phi(), Jets_cm_top6.at(i).M(), weight);

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

            my_histos["fwm2_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm2_top6, weightNoHT);
            my_histos["fwm3_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm3_top6, weightNoHT);
            my_histos["fwm4_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm4_top6, weightNoHT);
            my_histos["fwm5_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm5_top6, weightNoHT);
            my_histos["fwm6_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm6_top6, weightNoHT);
            my_histos["fwm7_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm7_top6, weightNoHT);
            my_histos["fwm8_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm8_top6, weightNoHT);
            my_histos["fwm9_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm9_top6, weightNoHT);
            my_histos["fwm10_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(fwm10_top6, weightNoHT);
            my_histos["jmt_ev0_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(jmt_ev0_top6, weightNoHT);
            my_histos["jmt_ev1_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(jmt_ev1_top6, weightNoHT);
            my_histos["jmt_ev2_top6_1l_ge7j_ge1b_NoHTweight"]->Fill(jmt_ev2_top6, weightNoHT);
            my_histos["GoodLeptons_pt_1_1l_ge7j_ge1b_NoHTweight"]->Fill(GoodLeptons.at(0).second.Pt(), weightNoHT);
            my_histos["GoodLeptons_eta_1_1l_ge7j_ge1b_NoHTweight"]->Fill(GoodLeptons.at(0).second.Eta(), weightNoHT);
            my_histos["GoodLeptons_phi_1_1l_ge7j_ge1b_NoHTweight"]->Fill(GoodLeptons.at(0).second.Phi(), weightNoHT);
            my_histos["GoodLeptons_m_1_1l_ge7j_ge1b_NoHTweight"]->Fill(GoodLeptons.at(0).second.M(), weightNoHT);

            for(unsigned int i = 0; i < Jets_cm_top6.size(); i++)
            {
                double jeta = Jets_cm_top6.at(i).Eta();
                if      (jeta > 0.0 and jeta <= 0.5)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta > 0.5 and jeta <= 1.0)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta > 1.0 and jeta <= 1.5)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta > 1.5 and jeta <= 2.0)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta > 2.0)                   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta <  0.0 and jeta >= -0.5) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta00_05_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta < -0.5 and jeta >= -1.0) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta05_10_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta < -1.0 and jeta >= -1.5) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta10_15_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta < -1.5 and jeta >= -2.0) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta15_20_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                else if (jeta < -2.0)                  my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight_eta20_Inf_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);

                my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightNoHT);
                my_2d_histos["Jet_cm_eta_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(Jets_cm_top6.at(i).Eta(), Jets_cm_top6.at(i).M(), weightNoHT);
                my_2d_histos["Jet_cm_eta_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(Jets_cm_top6.at(i).Eta(), Jets_cm_top6.at(i).Pt(), weightNoHT);
                my_2d_histos["Jet_cm_phi_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(Jets_cm_top6.at(i).Phi(), Jets_cm_top6.at(i).M(), weightNoHT);
                my_2d_histos["Jet_cm_phi_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(Jets_cm_top6.at(i).Phi(), Jets_cm_top6.at(i).Pt(), weightNoHT);

                my_histos["Jet_cm_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Pt()), weightNoHT);
                my_histos["Jet_cm_eta_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Eta()), weightNoHT);
                my_histos["Jet_cm_phi_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Phi()), weightNoHT);
                my_histos["Jet_cm_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b_NoHTweight"]->Fill(static_cast<double>(Jets_cm_top6.at(i).M()), weightNoHT);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b_NoHTweight"]->Fill(deepESM_val, Jets_cm_top6.at(i).Pt(), weightNoHT);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b_NoHTweight"]->Fill(deepESM_val, Jets_cm_top6.at(i).Eta(), weightNoHT);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b_NoHTweight"]->Fill(deepESM_val, Jets_cm_top6.at(i).Phi(), weightNoHT);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b_NoHTweight"]->Fill(deepESM_val, Jets_cm_top6.at(i).M(), weightNoHT);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b_NoHTweight"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Pt(), weightNoHT);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b_NoHTweight"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Eta(), weightNoHT);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b_NoHTweight"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Phi(), weightNoHT);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b_NoHTweight"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).M(), weightNoHT);
            }

        }

        // Fill histos for deepESM training
        if(passBaseline1l_NonIsoMuon)
        {
            my_histos["fwm2_top6_QCDCR"]->Fill(fwm2_top6, weightQCDCR);
            my_histos["fwm3_top6_QCDCR"]->Fill(fwm3_top6, weightQCDCR);
            my_histos["fwm4_top6_QCDCR"]->Fill(fwm4_top6, weightQCDCR);
            my_histos["fwm5_top6_QCDCR"]->Fill(fwm5_top6, weightQCDCR);
            my_histos["fwm6_top6_QCDCR"]->Fill(fwm6_top6, weightQCDCR);
            my_histos["fwm7_top6_QCDCR"]->Fill(fwm7_top6, weightQCDCR);
            my_histos["fwm8_top6_QCDCR"]->Fill(fwm8_top6, weightQCDCR);
            my_histos["fwm9_top6_QCDCR"]->Fill(fwm9_top6, weightQCDCR);
            my_histos["fwm10_top6_QCDCR"]->Fill(fwm10_top6, weightQCDCR);
            my_histos["jmt_ev0_top6_QCDCR"]->Fill(jmt_ev0_top6, weightQCDCR);
            my_histos["jmt_ev1_top6_QCDCR"]->Fill(jmt_ev1_top6, weightQCDCR);
            my_histos["jmt_ev2_top6_QCDCR"]->Fill(jmt_ev2_top6, weightQCDCR);

            for(unsigned int i = 0; i < Jets_cm_top6.size(); i++)
            {

                double jeta = Jets_cm_top6.at(i).Eta();
                if      (jeta > 0.0 and jeta <= 0.5)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta00_05_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta > 0.5 and jeta <= 1.0)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta05_10_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta > 1.0 and jeta <= 1.5)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta10_15_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta > 1.5 and jeta <= 2.0)   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta15_20_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta > 2.0)                   my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta20_Inf_p"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta <  0.0 and jeta >= -0.5) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta00_05_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta < -0.5 and jeta >= -1.0) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta05_10_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta < -1.0 and jeta >= -1.5) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta10_15_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta < -1.5 and jeta >= -2.0) my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta15_20_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                else if (jeta < -2.0)                  my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR_eta20_Inf_n"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);

                my_2d_histos["Jet_cm_pt_m_"+std::to_string(i+1)+"_QCDCR"]->Fill(Jets_cm_top6.at(i).Pt(), Jets_cm_top6.at(i).M(), weightQCDCR);
                my_2d_histos["Jet_cm_eta_m_"+std::to_string(i+1)+"_QCDCR"]->Fill(Jets_cm_top6.at(i).Eta(), Jets_cm_top6.at(i).M(), weightQCDCR);
                my_2d_histos["Jet_cm_eta_pt_"+std::to_string(i+1)+"_QCDCR"]->Fill(Jets_cm_top6.at(i).Eta(), Jets_cm_top6.at(i).Pt(), weightQCDCR);
                my_2d_histos["Jet_cm_phi_m_"+std::to_string(i+1)+"_QCDCR"]->Fill(Jets_cm_top6.at(i).Phi(), Jets_cm_top6.at(i).M(), weightQCDCR);
                my_2d_histos["Jet_cm_phi_pt_"+std::to_string(i+1)+"_QCDCR"]->Fill(Jets_cm_top6.at(i).Phi(), Jets_cm_top6.at(i).Pt(), weightQCDCR);

                my_histos["Jet_cm_pt_"+std::to_string(i+1)+"_QCDCR"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Pt()), weightQCDCR);
                my_histos["Jet_cm_eta_"+std::to_string(i+1)+"_QCDCR"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Eta()), weightQCDCR);
                my_histos["Jet_cm_phi_"+std::to_string(i+1)+"_QCDCR"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Phi()), weightQCDCR);
                my_histos["Jet_cm_m_"+std::to_string(i+1)+"_QCDCR"]->Fill(static_cast<double>(Jets_cm_top6.at(i).M()), weightQCDCR);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_deepESM_QCDCR"]->Fill(deepESM_val, Jets_cm_top6.at(i).Pt(), weightQCDCR);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_deepESM_QCDCR"]->Fill(deepESM_val, Jets_cm_top6.at(i).Eta(), weightQCDCR);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_deepESM_QCDCR"]->Fill(deepESM_val, Jets_cm_top6.at(i).Phi(), weightQCDCR);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_deepESM_QCDCR"]->Fill(deepESM_val, Jets_cm_top6.at(i).M(), weightQCDCR);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_njets_QCDCR"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Pt(), weightQCDCR);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_njets_QCDCR"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Eta(), weightQCDCR);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_njets_QCDCR"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Phi(), weightQCDCR);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_njets_QCDCR"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).M(), weightQCDCR);
            }

        }

        // ------------
        // -- Cut flow
        // ------------
        if(true) my_histos["h_cutFlow"]->Fill(0.5, weight);
        if(true && pass_general) my_histos["h_cutFlow"]->Fill(1.5, weight);  
        if(true && pass_general && pass_1l) my_histos["h_cutFlow"]->Fill(2.5, weight);
        if(true && pass_general && pass_1l && pass_ht) my_histos["h_cutFlow"]->Fill(3.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID) my_histos["h_cutFlow"]->Fill(4.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30) my_histos["h_cutFlow"]->Fill(5.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL) my_histos["h_cutFlow"]->Fill(6.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL && pass_njet_pt30) my_histos["h_cutFlow"]->Fill(7.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL && pass_njet_pt30 && passHEMVeto) my_histos["h_cutFlow"]->Fill(8.5, weight);   

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
