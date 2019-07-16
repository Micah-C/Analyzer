#define Semra_Analyzer_cxx
#include "Analyzer/Analyzer/include/Semra_Analyzer.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>

Semra_Analyzer::Semra_Analyzer() : inithisto(false) // SEMRA / define inithisto variable
{
}


/// Define histos
void Semra_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) // SEMRA / define variable map
{
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();

	my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;


    	for (const auto& cutVar : cutmap) {

		my_histos.emplace( "h_ntops_"+cutVar.first, std::make_shared<TH1D> ( ("h_ntops_"+cutVar.first).c_str(), ("h_ntops_"+cutVar.first).c_str(), 10, 0, 10 ) );
        	my_histos.emplace( "h_ht_"+cutVar.first, std::make_shared<TH1D> ( ("h_ht_"+cutVar.first).c_str(), ("h_ht_"+cutVar.first).c_str(), 60, 0, 3000 ) );
        	my_histos.emplace( "h_njets_"+cutVar.first, std::make_shared<TH1D> ( ("h_njets_"+cutVar.first).c_str(), ("h_njets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        	my_histos.emplace( "h_nbjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nbjets_"+cutVar.first).c_str(), ("h_nbjets_"+cutVar.first).c_str(), 15, 0, 15 ) );
		my_2d_histos.emplace( "h_njets_MVA_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_MVA_"+cutVar.first).c_str(), ("h_njets_MVA_"+cutVar.first).c_str(), 8, 7, 15, 50, 0, 1.0 ) );

	}


	//Define TEfficiencies if you are doing trigger studies (for proper error bars) or cut flow charts.
    	my_efficiencies.emplace("event_sel_weight", std::make_shared<TEfficiency>("event_sel_weight","event_sel_weight",9,0,9));
	

}


/// Put everything you want to do per event 
void Semra_Analyzer::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
 
   while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );
        
        const auto& runtype             = tr.getVar<std::string>("runtype");     
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC       = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt45     = tr.getVar<int>("NGoodBJets_pt45");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt45     = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodJets_pt45      = tr.getVar<int>("NGoodJets_pt45");
        
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline        = tr.getVar<bool>("passBaseline1l_Good");
       

	// -----------------------------------------------
	// -- SEMRA / Define Top Tag variables
	// -----------------------------------------------

	const auto& deepESM_val            = tr.getVar<double>("deepESM_val");
	const bool& pass_0l                = NGoodLeptons==0;
	
	const auto& ntops                  = tr.getVar<int>("ntops");
	const auto& ntops_1jet             = tr.getVar<int>("ntops_1jet");
	const auto& ntops_2jet             = tr.getVar<int>("ntops_2jet");
	const auto& ntops_3jet             = tr.getVar<int>("ntops_3jet");
	const auto& pass_HT_trigger        = HT_trigger_pt45 > 500;
	const auto& pass_NJets_pt45        = NGoodJets_pt45 >= 6;
	const auto& pass_NBJets_pt45_2b    = NGoodBJets_pt45 >= 2;

	

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double topPtScaleFactor     = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            topPtScaleFactor = tr.getVar<double>("topPtScaleFactor");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
        

        // SEMRA / Make cuts and fill histograms here & cutmap
        // cutmap
        const std::map<std::string, bool>& cutmap
	{
		{"", true}, // no cut
		{"0l", pass_0l},

		{"0l_HT500", pass_0l && pass_HT_trigger},
		{"0l_ge6j", pass_0l && pass_NJets_pt45},
		{"0l_ge1b", pass_0l && NGoodBJets_pt45 >= 1},
		{"0l_ge1b_HT500", pass_0l && NGoodBJets_pt45 >= 1 && pass_HT_trigger},  
		{"0l_ge2b", pass_0l && pass_NBJets_pt45_2b},
		{"0l_1t", pass_0l && ntops==1},
		{"0l_2t", pass_0l && ntops==2},
		{"0l_1t1j", pass_0l && ntops_1jet==1},
		{"0l_1t3j", pass_0l && ntops_3jet==1},
		{"0l_2t1j", pass_0l && ntops_1jet==2},
                {"0l_2t3j", pass_0l && ntops_3jet==2},
	

	};

	if (!inithisto) {
		InitHistos(cutmap);
		inithisto = true;
	}


	my_histos["EventCounter" ]->Fill( eventCounter );


	for (const auto& cutVar: cutmap) {

		if (cutVar.second) {
			my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight );
			my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt45, weight );		
			my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt45, weight );
			my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt45, weight );
			my_2d_histos["h_njets_MVA_"+cutVar.first]->Fill( NGoodJets_pt45, deepESM_val, weight );
		}
	}


        // Example Fill event selection efficiencies
        my_efficiencies["event_sel_weight"]->SetUseWeightedEvents();
        my_efficiencies["event_sel_weight"]->FillWeighted(true,eventweight,0);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID,eventweight,1);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1,eventweight,2);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC,eventweight,3);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt45 >= 1,eventweight,4);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt45 >= 1 && 50 < Mbl && Mbl < 250,eventweight,5);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt45 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt45 > 300,eventweight,6);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt45 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt45 > 300 && NGoodJets_pt45 >= 7,eventweight,7);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt45 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt45 > 300 && NGoodJets_pt45 >= 7,weight,8);
    } 
}


void Semra_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

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