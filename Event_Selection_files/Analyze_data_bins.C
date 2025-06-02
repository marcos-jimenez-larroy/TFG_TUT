#define Analyze_data_bins_cxx
#include "Analyze_data_bins.h"
#include <TH2.h>
#include <TStyle.h>
//#include "DCB.C"
//#include "DCB_RF.C"
#include <TCanvas.h>
#include <iostream>
#include <map>



void Analyze_data_bins::Loop()
{
//   In a ROOT session, you can do:                                                                                                            
//      root> .L Analyze_data_bins.C                                                                                                                   
//      root> Analyze_data_bins t                                                                                                                      
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


Double_t xbins[30] = {0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0,36.0,39.0,42.0,45.0,48.0,51.0,54.0,57.0,60.0,69.0,78.0,87.0,96.0,132.0,186.0,258.0,618.0,1500.0};
Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};

   int numBins_pT = 16;  // Number of pT_yy bins  
   int numBins_myy = 55; // Number of bins for m_yy (En la inclusiva eran 80 Bines pero ahora habrá menos estadística así que congemos bines de 1 GeV)
   double m_yy_min = 105.0, m_yy_max = 160.0;
   
   // Map to store histograms for each pT_yy bin
   std::map<std::string, TH1D*> hMappT;

// Create histograms for each pT_yy bin
for (int i = 0; i < numBins_pT; ++i) {
   std::string histName = "h_myy_pT_" + std::to_string(i);
   hMappT[histName] = new TH1D(histName.c_str(), Form("m_yy in pT bin [%g, %g]", new_xbins[i], new_xbins[i+1]),
            numBins_myy, m_yy_min, m_yy_max);
 }
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;





   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      

      
      
      Char_t SetCuts = (isPassedPreselection && isPassedTriggerMatch && isPassedPID && isPassedIsolation && isPassedMassCut);
      // if (Cut(ientry) < 0) continue;

      if (SetCuts && pT_y2 > 0.25*m_yy && (pT_y1*pT_y2) > (0.32*0.32)*m_yy*m_yy){
         for (int j = 0; j < numBins_pT; ++j) {
            if (pT_yy/1000. >= new_xbins[j] && pT_yy/1000. < new_xbins[j+1]) {
              std::string histName = "h_myy_pT_" + std::to_string(j);
              hMappT[histName]->Fill(m_yy/1000.);
              break; // Stop after finding the correct bin
            }
               }
         //std::cout<<" myy "<<m_yy/1000.<<std::endl;
      }
      // if (Cut(ientry) < 0) continue;                                                                                                        
      //std::cout<<" myy "<<m_yy/1000.<<std::endl;
   }

   //int dcb = DCB_RF(h0yymass,115.,135.,1.5e6);

   //Change name based on which production mode you are using!!
   TFile* myfile = new TFile("/home/marcosj/HGamCore/New_Root/data22_23_bins_ProductCuts.root", "RECREATE"); //Recreate significa que si no está lo cree y si está que lo reescriba.

   myfile->cd();
   for (auto& pair : hMappT) {
     pair.second->Write();
     delete pair.second;  // Clean up memory
   }
   myfile->Close();
}