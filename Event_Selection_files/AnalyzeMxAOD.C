#define AnalyzeMxAOD_cxx
#include "AnalyzeMxAOD.h"
#include <TH2.h>
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <map>


//#include "DCB.C"
//#include "DCB_RF.C"
#include <TCanvas.h>



void AnalyzeMxAOD::Loop()
{
//   In a ROOT session, you can do:                                                                                                            
//      root> .L AnalyzeMxAOD.C                                                                                                                   
//      root> AnalyzeMxAOD t                                                                                                                      
//      root> t.GetEntry(12); // Fill t data members with entry number 12                                                                      
//      root> t.Show();       // Show values of entry 12                                                                                       
//      root> t.Show(16);     // Read and show values of entry 16                                                                              
//      root> t.Loop();       // Loop on all entries                                                                                           
//                                                                                                                                             

//     This is the loop skeleton where:                                                                                                        
//    jentry is the global entry number in the chain                                                                                           
//    ientry is the entry number in the current Tree                                                                                           
//  Note that the argument to GetEntry must be:                                                                                                
//    jentry for TChain::GetEntry                                                                                                              
//    ientry for TTree::GetEntry and TBranch::GetEntry                                                                                         
//                                                                                                                                             
//       To read only selected branches, Insert statements like:                                                                               
// METHOD1:                                                                                                                                    
//    fChain->SetBranchStatus("*",0);  // disable all branches                                                                                 
//    fChain->SetBranchStatus("branchname",1);  // activate branchname                                                                         
// METHOD2: replace line                                                                                                                       
//    fChain->GetEntry(jentry);       //read all branches                                                                                      
//by  b_branchname->GetEntry(ientry); //read only this branch


    //TH1* h0yymass = new TH1D("h0yymass","H0mass",55,105,160); 

    //TH1* h0yypt_fid = new TH1D("h0yypt_fid","PTyy_fid",500,0,1500);


    Double_t xbins[30] = {0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0,36.0,39.0,42.0,45.0,48.0,51.0,54.0,57.0,60.0,69.0,78.0,87.0,96.0,132.0,186.0,258.0,618.0,1500.0};
    Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};

    double numBinspT = 16;

    TH1D* h_pTTruth = new TH1D("h_pTTruth", "h_pTTruth", numBinspT, new_xbins);


    std::map<TString, vector<int>> variations={
      {"ggH", {0,88,94,95,96,97,98,99,100,101,102,111 }}, // el 0 es NLO y el 88 NNLO nominal
      {"VBF", {0,2,3,6,7,8,9}},
      {"WmH", {0,2,3,4,5,7,9}},
      {"WpH", {0,2,3,4,5,7,9}},
      {"ZH", {0,2,3,4,5,7,9}},
      {"ttH", {0,2,3,6,7,8,9}}
    };


    //TH2* ptyy_migration = new TH2D("ptyy_migration","PTyy migration",16,new_xbins,16,new_xbins);


   if (fChain == 0) return;

   //Map con las variaciones (por ahora solo para ggH)
   std::map<std::string, TH1D*> hMappTTruth;


   TString prodMode = "ZH"; //Change when not using ggH
   for (int i =0; i< (variations[prodMode].size()); i++)
     {
       std::cout<< " i "<<i<<" total size "<<variations[prodMode].size()<<std::endl;
       std::string hname = "h_pT_truth_var_" + std::to_string(variations[prodMode].at(i));
       std::cout<<hname<<std::endl;
       hMappTTruth[hname] = new TH1D(hname.c_str(),hname.c_str(),numBinspT,new_xbins);
     }



   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   double ww = 1;


   //Sum of weights hechos a mano por ahora
   double w_ggH = 226630780.;
   double w_VBF = 7597700.7;
   double w_ttH = 441624.5;
   double w_WmH = 115366.78;
   double w_WpH = 184533.81;
   double w_ZH = 419690.66;

   //mcChannelNumber característico
   Int_t mcCN_ggH = 602421;
   Int_t mcCN_VBF = 601482;
   Int_t mcCN_ttH = 602422;
   Int_t mcCN_WmH = 601483;
   Int_t mcCN_WpH = 601484;
   Int_t mcCN_ZH = 601523;

   double SumOfWeights = 0.0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (mcChannelNumber == mcCN_ggH){
         SumOfWeights = w_ggH;
      } else if (mcChannelNumber == mcCN_VBF){
         SumOfWeights = w_VBF;
      } else if (mcChannelNumber == mcCN_ttH){
         SumOfWeights = w_ttH;
      } else if (mcChannelNumber == mcCN_WmH){
         SumOfWeights = w_WmH;
      } else if (mcChannelNumber == mcCN_WpH){
         SumOfWeights = w_WpH;
      } else if (mcChannelNumber == mcCN_ZH){
         SumOfWeights = w_ZH;
      }else{std::cout<<"Revisa esto"<<std::endl;}
   


      
      
      //ww = weight*30.0*1000*crossSectionBRfilterEff/SumOfWeights; //La L de 30 1/fb igual es 32/33 pero por ahora vale

      // if (Cut(ientry) < 0) continue;

      if (isFiducial){
         h_pTTruth->Fill(pT_yy_truth/1000.,weightInitial*30.0*1000*crossSectionBRfilterEff/SumOfWeights);
      
   
         if(prodMode == "ggH"){
           for (int i =0; i< (variations[prodMode].size()); i++)
             {
               std::string hname = "h_pT_truth_var_" + std::to_string(variations[prodMode].at(i));
               hMappTTruth[hname]->Fill(pT_yy_truth/1000., (weightInitial*MCweights->at(0).at(variations[prodMode].at(i))/MCweights->at(0).at(88))*30.0*1000*crossSectionBRfilterEff/SumOfWeights);
               //else hMappTTruth[hname]->Fill(pT_yy_truth/1000., (weightInitial*MCweights->at(0).at(variations[prodMode].at(i))/MCweights->at(0).at(0))*30.0*1000*crossSectionBRfilterEff/SumOfWeights);
               
             }
         }

      /*if (isPassed&&isFiducial){
         ptyy_migration->Fill(pT_yy/1000.0,pT_yy_truth/1000.0,ww);
      }*/
      }
   }
 




   //int dcb = DCB_RF(h0yymass,115.,135.,1.5e6);
   

   //Change name based on which production mode you are using!!
   TFile* myfile = new TFile("/home/marcosj/HGamCore/New_Root/ggH.root", "RECREATE"); //Recreate significa que si no está lo cree y si está que lo reescriba.

   //h0yymass->Write();
   h_pTTruth->Write();
   if(prodMode == "ggH"){
   myfile->cd();
   for (auto& pair : hMappTTruth) {
     pair.second->Write();
     delete pair.second;  // Clean up memory
   }
   }
   
   //dcb->Write();
   myfile->Close();

   /*TFile* migration_file = new TFile("/home/marcosj/HGamCore/New_Root/mc23a_ptyy_migration.root", "RECREATE");

   ptyy_migration->Write();
   migration_file->Close();*/
}