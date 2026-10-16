#ifndef Analyze_data_bins_h
#define Analyze_data_bins_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Analyze_data_bins {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain                                                                        
   Int_t           fCurrent; //!current Tree number in a TChain  

   Char_t          isPassed;
   Float_t         m_yy;
   Float_t         pT_yy;
   Float_t         weight;
   Float_t         crossSectionBRfilterEff;

   Char_t         isPassedPreselection;
   Char_t         isPassedTriggerMatch;
   Char_t         isPassedPID;
   Char_t         isPassedIsolation;
   Char_t         isPassedMassCut;
   Char_t         isPassed_PC;

   Float_t       pT_y1;
   Float_t       pT_y2;




   TBranch        *b_HGamEventInfoAuxDyn_isPassed; //! 
   TBranch        *b_HGamEventInfoAuxDyn_m_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_yy;  //!
   TBranch        *b_HGamEventInfoAuxDyn_weight;  //!
   TBranch        *b_HGamEventInfoAuxDyn_crossSectionBRfilterEff; //!


   TBranch        *b_HGamEventInfoAuxDyn_isPassedPreselection;
   TBranch        *b_HGamEventInfoAuxDyn_isPassedTriggerMatch;
   TBranch        *b_HGamEventInfoAuxDyn_isPassedPID;
   TBranch        *b_HGamEventInfoAuxDyn_isPassedIsolation;
   TBranch        *b_HGamEventInfoAuxDyn_isPassedMassCut;
   TBranch        *b_HGamEventInfoAuxDyn_isPassed_PC;
   TBranch        *b_HGamEventInfoAuxDyn_pT_y1;
   TBranch        *b_HGamEventInfoAuxDyn_pT_y2;



   Analyze_data_bins(TTree *tree=0);
   virtual ~Analyze_data_bins();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analyze_data_bins_cxx
Analyze_data_bins::Analyze_data_bins(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file                                                                              
// used to generate this class and read the Tree.                                                                                              
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/userdata/acueto/HGam_datas/data22/data22_13p6TeV.periodAllYear.physics_Main.MxAODDetailed.p6456.h031.001.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/userdata/acueto/HGam_datas/data22/data22_13p6TeV.periodAllYear.physics_Main.MxAODDetailed.p6456.h031.001.root");
      }
      f->GetObject("CollectionTree",tree);

   }
   Init(tree);
}

Analyze_data_bins::~Analyze_data_bins()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyze_data_bins::GetEntry(Long64_t entry)
{
// Read contents of entry.                                                                                                                     
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyze_data_bins::LoadTree(Long64_t entry)
{
// Set the environment to read one entry                                                                                                       
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
void Analyze_data_bins::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize                                                                      
   // a new tree or chain. Typically here the branch addresses and branch                                                                      
   // pointers of the tree will be set.                                                                                                        
   // It is normally not necessary to make changes to the generated                                                                            
   // code, but the routine can be extended by the user if needed.                                                                             
   // Init() will be called many times when running on PROOF                                                                                   
   // (once per file to be processed).  

    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassed",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.m_yy",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.pT_yy",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.weight",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.crossSectionBRfilterEff",1);

    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassedPreselection",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassedTriggerMatch",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassedPID",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassedIsolation",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassedMassCut",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassed_PC",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.pT_y1",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.pT_y2",1);




    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassed", &isPassed, &b_HGamEventInfoAuxDyn_isPassed);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy", &m_yy, &b_HGamEventInfoAuxDyn_m_yy);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_yy", &pT_yy, &b_HGamEventInfoAuxDyn_pT_yy);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.weight", &weight, &b_HGamEventInfoAuxDyn_weight);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.crossSectionBRfilterEff", &crossSectionBRfilterEff, &b_HGamEventInfoAuxDyn_crossSectionBRfilterEff);

    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedPreselection", &isPassedPreselection, &b_HGamEventInfoAuxDyn_isPassedPreselection);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedTriggerMatch", &isPassedTriggerMatch, &b_HGamEventInfoAuxDyn_isPassedTriggerMatch);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedPID", &isPassedPID, &b_HGamEventInfoAuxDyn_isPassedPID);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedIsolation", &isPassedIsolation, &b_HGamEventInfoAuxDyn_isPassedIsolation);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedMassCut", &isPassedMassCut, &b_HGamEventInfoAuxDyn_isPassedMassCut);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassed_PC", &isPassed_PC, &b_HGamEventInfoAuxDyn_isPassed_PC);


    fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_y1", &pT_y1, &b_HGamEventInfoAuxDyn_pT_y1);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_y2", &pT_y2, &b_HGamEventInfoAuxDyn_pT_y2);


    Notify();
}

bool Analyze_data_bins::Notify()
{
   // The Notify() function is called when a new file is opened. This                                                                          
   // can be either for a new TTree in a TChain or when when a new TTree                                                                       
   // is started when using PROOF. It is normally not necessary to make changes                                                                
   // to the generated code, but the routine can be extended by the                                                                            
   // user if needed. The return value is currently not used.                                                                                  

   return true;
}

void Analyze_data_bins::Show(Long64_t entry)
{
// Print contents of entry.                                                                                                                    
// If entry is not specified, print current entry                                                                                              
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyze_data_bins::Cut(Long64_t entry)
{
// This function may be called from Loop.                                                                                                      
// returns  1 if entry is accepted.                                                                                                            
// returns -1 otherwise.                                                                                                                       
   return 1;
}
#endif // #ifdef Analyze_data_bins_cxx  
