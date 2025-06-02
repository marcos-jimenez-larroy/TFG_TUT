#define AnalyzeBackground_PC_cxx
#include "AnalyzeBackground_PC.h"
#include <TH2.h>
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <map>


//#include "DCB.C"
//#include "DCB_RF.C"
#include <TCanvas.h>



void AnalyzeBackground_PC::Loop()
{

    Double_t cuts[31] = {27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,41.5,42};

    double num_cuts = 31;

    

    
        //TH1* h0yymass_PC = new TH1D("h0yymass_PC","H0mass_PC",55,105.0,160.0); 
        TH1* h0yymass_PC_2 = new TH1D("h0yymass_PC_2","H0mass_PC_2",55,105.0,160.0); 
        TH1* h0yymass_RC_2 = new TH1D("h0yymass_RC_2","H0mass_RC_2",55,105.0,160.0); 

        //TH1* h0yymass = new TH1D("h0yymass","H0mass",55,105.0,160.0);
        //TH1* h0yymass_isPassed = new TH1D("h0yymass_isPassed","H0mass_isPassed",55,105.0,160.0);
        //TH1* h0yymass_isPassed_manual = new TH1D("h0yymass_isPassed_manual","H0mass_isPassed_manual",55,105.0,160.0);


        std::map<std::string, TH1D*> hMappT_PC;
        std::map<std::string, TH1D*> hMappT_RC;

        for (int i = 0; i < num_cuts; ++i) {
            std::string histName = "PC_" + std::to_string(cuts[i]);
            hMappT_PC[histName] = new TH1D(Form("PC_%d",i), Form("Cut of %f", cuts[i]),
                     55, 105., 160.);

            histName = "RC_" + std::to_string(cuts[i]);
            hMappT_RC[histName] = new TH1D(Form("RC_%d",i), Form("Relative Cut of %f", cuts[i]),
                    55, 105., 160.);
          }



        if (fChain == 0) return;

        Long64_t nentries = fChain->GetEntriesFast();

        Long64_t nbytes = 0, nb = 0;

        double SumOfWeights = 1.02238e15; //SumofWeights del ggH

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);   nbytes += nb;


            //Char_t pass_trigger = (passTrig_HLT_2g20_tight_icaloloose_L12EM15VHI || passTrig_HLT_2g22_tight_L12EM15VHI || passTrig_HLT_g120_loose_L1EM22VHI || passTrig_HLT_g140_loose_L1EM22VHI || passTrig_HLT_g35_medium_g25_medium_L12EM20VH);


            Char_t SetCuts = (isPassedPreselection && isPassedTriggerMatch && isPassedPID && isPassedIsolation && isPassedMassCut);

            //if(isPassed_PC){
                //h0yymass_PC->Fill(m_yy/1000.,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights); //This only for the 0.32 cut
            //}


            if (SetCuts && pT_y2 > 0.25*m_yy){
                double ww = weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights; 

                for (int j = 0; j < num_cuts; ++j) {

                   if (sqrt(pT_y1*pT_y2) > (cuts[j]/100.)*m_yy) {

                     std::string histName = "PC_" + std::to_string(cuts[j]);
                     hMappT_PC[histName]->Fill(m_yy/1000.,ww);
                   }

                   if (pT_y1 > (cuts[j]/100.)*m_yy) {

                    std::string histName = "RC_" + std::to_string(cuts[j]);
                    hMappT_RC[histName]->Fill(m_yy/1000.,ww);
                  }


                      }
                //std::cout<<" myy "<<m_yy/1000.<<std::endl;
             }


            /*
            if(SetCuts && pT_y2 > 0.25*m_yy){ //&& pT_y2 > 0.25*m_yy   // EN EL ISPASSED_PC FALTAN LOS DOS CORTES DE PTYY

                h0yymass->Fill(m_yy/1000.,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights);

                if(sqrt(pT_y1*pT_y2) > (cut/100.)*m_yy){
                    h0yymass_PC_2->Fill(m_yy/1000.,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights); //This only for the 0.32 cut
                }*/
                

                /*
                if(pT_y1 > 0.35*m_yy){
                    h0yymass_isPassed_manual->Fill(m_yy/1000.,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights);
                }
                
                if(sqrt(pT_y1*pT_y2) > 0.01*cuts[6]*m_yy){
                    h0yymass_PC->Fill(m_yy/1000.,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights);
                    pT2->Fill(pT_y2/1000.,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights);
                    pT1->Fill(pT_y1/1000.,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights);
                    pT->Fill(pT_yy/1000.0,weight*27.48*1000*crossSectionBRfilterEff/SumOfWeights);

                }
            } */

            
            
            }

            TFile* myfile = new TFile("/home/marcosj/HGamCore/New_Root/mc23a_background_ProductCuts_RelativeCuts.root", "RECREATE"); //Recreate significa que si no está lo cree y si está que lo reescriba.

            myfile->cd();
            for (auto& pair : hMappT_PC) {
                pair.second->Write();
                delete pair.second;  // Clean up memory
            }
            for (auto& pair : hMappT_RC) {
                pair.second->Write();
                delete pair.second;  // Clean up memory
            }
            myfile->Close();

}