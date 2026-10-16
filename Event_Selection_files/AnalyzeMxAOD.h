#ifndef AnalyzeMxAOD_h
#define AnalyzeMxAOD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>


class AnalyzeMxAOD {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain                                                                        
   Int_t           fCurrent; //!current Tree number in a TChain  

   Char_t          isPassed;
   Char_t          isFiducial;
   Float_t         m_yy;
   Float_t         pT_yy;
   Float_t         pT_yy_truth;
   Float_t         weight;
   Float_t         weightInitial;
   Float_t         crossSectionBRfilterEff;
   Int_t           mcChannelNumber;
   vector<vector<float>> *MCweights;


   TBranch        *b_HGamEventInfoAuxDyn_isPassed; //! 
   TBranch        *b_HGamTruthEventInfoAuxDyn_isFiducial;
   TBranch        *b_HGamEventInfoAuxDyn_m_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_yy;  //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_pT_yy;  //!
   TBranch        *b_HGamEventInfoAuxDyn_weight;  //!
   TBranch        *b_HGamEventInfoAuxDyn_crossSectionBRfilterEff; //!
   TBranch        *b_HGamEventInfoAuxDyn_weightInitial;//!
   TBranch        *b_EventInfoAuxDyn_mcChannelNumber; //!
   TBranch        *b_TruthEventsAuxDyn_weights;//!

   AnalyzeMxAOD(TTree *tree=0);
   virtual ~AnalyzeMxAOD();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalyzeMxAOD_cxx
AnalyzeMxAOD::AnalyzeMxAOD(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file                                                                              
// used to generate this class and read the Tree.                                                                                              
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/userdata/acueto/HGamMxAODs/mc23a/mc23a.PhPy8_VBFH125.MxAODDetailedNoSkim.e8532_s4162_r14622_p6454_h031.root/mc23a.PhPy8_VBFH125.MxAODDetailedNoSkim.e8532_s4162_r14622_p6454_h031.001.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/userdata/acueto/HGamMxAODs/mc23a/mc23a.PhPy8_VBFH125.MxAODDetailedNoSkim.e8532_s4162_r14622_p6454_h031.root/mc23a.PhPy8_VBFH125.MxAODDetailedNoSkim.e8532_s4162_r14622_p6454_h031.001.root");
      }
      f->GetObject("CollectionTree",tree);

   }
   Init(tree);
}

AnalyzeMxAOD::~AnalyzeMxAOD()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeMxAOD::GetEntry(Long64_t entry)
{
// Read contents of entry.                                                                                                                     
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeMxAOD::LoadTree(Long64_t entry)
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
void AnalyzeMxAOD::Init(TTree *tree)
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

    MCweights = 0;

    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.isPassed",1);
    fChain->SetBranchStatus("HGamTruthEventInfoAuxDyn.isFiducial",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.m_yy",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.pT_yy",1);
    fChain->SetBranchStatus("HGamTruthEventInfoAuxDyn.pT_yy",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.weight",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.crossSectionBRfilterEff",1);
    fChain->SetBranchStatus("HGamEventInfoAuxDyn.weightInitial",1);
    fChain->SetBranchStatus("EventInfoAuxDyn.mcChannelNumber",1);
    fChain->SetBranchStatus("TruthEventsAuxDyn.weights",1);


    fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassed", &isPassed, &b_HGamEventInfoAuxDyn_isPassed);
    fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.isFiducial", &isFiducial, &b_HGamTruthEventInfoAuxDyn_isFiducial);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy", &m_yy, &b_HGamEventInfoAuxDyn_m_yy);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_yy", &pT_yy, &b_HGamEventInfoAuxDyn_pT_yy);
    fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.pT_yy", &pT_yy_truth, &b_HGamTruthEventInfoAuxDyn_pT_yy);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.weight", &weight, &b_HGamEventInfoAuxDyn_weight);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.crossSectionBRfilterEff", &crossSectionBRfilterEff, &b_HGamEventInfoAuxDyn_crossSectionBRfilterEff);
    fChain->SetBranchAddress("HGamEventInfoAuxDyn.weightInitial", &weightInitial, &b_HGamEventInfoAuxDyn_weightInitial);
    fChain->SetBranchAddress("EventInfoAuxDyn.mcChannelNumber", &mcChannelNumber, &b_EventInfoAuxDyn_mcChannelNumber);
    fChain->SetBranchAddress("TruthEventsAuxDyn.weights", &MCweights, &b_TruthEventsAuxDyn_weights);
    
    Notify();
}

bool AnalyzeMxAOD::Notify()
{
   // The Notify() function is called when a new file is opened. This                                                                          
   // can be either for a new TTree in a TChain or when when a new TTree                                                                       
   // is started when using PROOF. It is normally not necessary to make changes                                                                
   // to the generated code, but the routine can be extended by the                                                                            
   // user if needed. The return value is currently not used.                                                                                  

   return true;
}

void AnalyzeMxAOD::Show(Long64_t entry)
{
// Print contents of entry.                                                                                                                    
// If entry is not specified, print current entry                                                                                              
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeMxAOD::Cut(Long64_t entry)
{
// This function may be called from Loop.                                                                                                      
// returns  1 if entry is accepted.                                                                                                            
// returns -1 otherwise.                                                                                                                       
   return 1;
}
#endif // #ifdef AnalyzeMxAOD_cxx  
