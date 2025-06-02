#define Analyze_data_cxx
#include "Analyze_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Analyze_data::Loop()
{
//   In a ROOT session, you can do:                                                                                                            
//      root> .L Analyze_data.C                                                                                                                   
//      root> Analyze_data t                                                                                                                      
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


    TH1* h0yymass = new TH1D("h0yymass","H0mass",55,105,160);

    TH1* h0yypt = new TH1D("h0yypt","PTyy",500,0.0,1500);
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;



      // if (Cut(ientry) < 0) continue;

      if (isPassed){
         h0yymass->Fill(m_yy/1000.0);
         h0yypt->Fill(pT_yy/1000.0);
      }
      // if (Cut(ientry) < 0) continue;                                                                                                        
      //std::cout<<" myy "<<m_yy/1000.<<std::endl;
   }

   TFile* myfile = new TFile("/home/marcosj/HGamCore/New_Root/data22_23_histos.root", "RECREATE"); //Recreate significa que si no está lo cree y si está que lo reescriba.

   h0yymass->Write();
   h0yypt->Write();
}