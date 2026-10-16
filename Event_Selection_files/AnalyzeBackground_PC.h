#ifndef AnalyzeBackground_PC_h
#define AnalyzeBackground_PC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>


class AnalyzeBackground_PC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain                                                                        
   Int_t           fCurrent; //!current Tree number in a TChain  

   Char_t          isPassed;
   //Char_t          isFiducial;
   Float_t         m_yy;
   Float_t         pT_yy;
   //Float_t         pT_yy_truth;
   Float_t         weight;
   //Float_t         weightInitial;
   Float_t         crossSectionBRfilterEff;
   //Int_t           mcChannelNumber;
   //vector<vector<float>> *MCweights;

   //For Product Cut analysis

   Char_t         isPassedPreselection;
   //Char_t         passTrig_HLT_2g20_tight_icaloloose_L12EM15VHI;
   Char_t         passTrig_HLT_2g22_tight_L12EM15VHI;
   Char_t         passTrig_HLT_g120_loose_L1EM22VHI;
   Char_t         passTrig_HLT_g140_loose_L1EM22VHI;
   Char_t         passTrig_HLT_g35_medium_g25_medium_L12EM20VH;
   Char_t         isPassedTriggerMatch;
   Char_t         isPassedPID;
   Char_t         isPassedIsolation;
   Char_t         isPassedMassCut;
   Char_t         isPassed_PC;

   Float_t       pT_y1;
   Float_t       pT_y2;



   TBranch        *b_isPassed; //! 
   //TBranch        *b_isFiducial;
   TBranch        *b_m_yy;   //!
   TBranch        *b_pT_yy;  //!
   //TBranch        *b_pT_yy;  //!
   TBranch        *b_weight;  //!
   TBranch        *b_crossSectionBRfilterEff; //!
   //TBranch        *b_weightInitial;//!
   //TBranch        *b_mcChannelNumber; //!
   //TBranch        *b_weights;//!
   

   //For Product Cut Analysis

   TBranch        *b_isPassedPreselection;
   //TBranch        *b__passTrig_HLT_2g20_tight_icaloloose_L12EM15VHI;
   TBranch        *b_passTrig_HLT_2g22_tight_L12EM15VHI;
   TBranch        *b_passTrig_HLT_g120_loose_L1EM22VHI;
   TBranch        *b_passTrig_HLT_g140_loose_L1EM22VHI;
   TBranch        *b_passTrig_HLT_g35_medium_g25_medium_L12EM20VH;

   TBranch        *b_isPassedTriggerMatch;
   TBranch        *b_isPassedPID;
   TBranch        *b_isPassedIsolation;
   TBranch        *b_isPassedMassCut;
   TBranch        *b_isPassed_PC;
   TBranch        *b_pT_y1;
   TBranch        *b_pT_y2;

   AnalyzeBackground_PC(TTree *tree=0);
   virtual ~AnalyzeBackground_PC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalyzeBackground_PC_cxx
AnalyzeBackground_PC::AnalyzeBackground_PC(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file                                                                              
// used to generate this class and read the Tree.                                                                                              
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/userdata/acueto/HGamMxAODs/mc23a/BkgTest/hist-TestBkg.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/userdata/acueto/HGamMxAODs/mc23a/BkgTest/hist-TestBkg.root");
      }
      f->GetObject("BkgTree",tree);

   }
   Init(tree);
}

AnalyzeBackground_PC::~AnalyzeBackground_PC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeBackground_PC::GetEntry(Long64_t entry)
{
// Read contents of entry.                                                                                                                     
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeBackground_PC::LoadTree(Long64_t entry)
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
void AnalyzeBackground_PC::Init(TTree *tree)
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

    //MCweights = 0;

    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("isPassed",1);
    //fChain->SetBranchStatus("isFiducial",1);
    fChain->SetBranchStatus("m_yy",1);
    fChain->SetBranchStatus("pT_yy",1);
    //fChain->SetBranchStatus("pT_yy",1);
    fChain->SetBranchStatus("weight",1);
    fChain->SetBranchStatus("crossSectionBRfilterEff",1);
    //fChain->SetBranchStatus("weightInitial",1);
    //fChain->SetBranchStatus("mcChannelNumber",1);
    //fChain->SetBranchStatus("weights",1);
    //Product Cut Analysis

    fChain->SetBranchStatus("isPassedPreselection",1);
    fChain->SetBranchStatus("isPassed_PC",1);
    fChain->SetBranchStatus("isPassedTriggerMatch",1);
    fChain->SetBranchStatus("isPassedPID",1);
    fChain->SetBranchStatus("isPassedIsolation",1);
    fChain->SetBranchStatus("isPassedMassCut",1);
    //fChain->SetBranchStatus("passTrig_HLT_2g20_tight_icaloloose_L12EM15VHI",1);
    fChain->SetBranchStatus("passTrig_HLT_2g22_tight_L12EM15VHI",1);
    fChain->SetBranchStatus("passTrig_HLT_g120_loose_L1EM22VHI",1);
    fChain->SetBranchStatus("passTrig_HLT_g140_loose_L1EM22VHI",1);
    fChain->SetBranchStatus("passTrig_HLT_g35_medium_g25_medium_L12EM20VH",1);
    fChain->SetBranchStatus("pT_y1",1);
    fChain->SetBranchStatus("pT_y2",1);







    fChain->SetBranchAddress("isPassed", &isPassed, &b_isPassed);
    //fChain->SetBranchAddress("HGamTruth.isFiducial", &isFiducial, &b_HGamTruth_isFiducial);
    fChain->SetBranchAddress("m_yy", &m_yy, &b_m_yy);
    fChain->SetBranchAddress("pT_yy", &pT_yy, &b_pT_yy);
    //fChain->SetBranchAddress("pT_yy", &pT_yy_truth, &b_pT_yy);
    fChain->SetBranchAddress("weight", &weight, &b_weight);
    fChain->SetBranchAddress("crossSectionBRfilterEff", &crossSectionBRfilterEff, &b_crossSectionBRfilterEff);
    //fChain->SetBranchAddress("weightInitial", &weightInitial, &b_weightInitial);
    //fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
    //fChain->SetBranchAddress("weights", &MCweights, &b_weights);

    //Product Cut Analysis

    fChain->SetBranchAddress("isPassedPreselection", &isPassedPreselection, &b_isPassedPreselection);
    fChain->SetBranchAddress("isPassed_PC", &isPassed_PC, &b_isPassed_PC);

    fChain->SetBranchAddress("isPassedTriggerMatch", &isPassedTriggerMatch, &b_isPassedTriggerMatch);
    fChain->SetBranchAddress("isPassedPID", &isPassedPID, &b_isPassedPID);
    fChain->SetBranchAddress("isPassedIsolation", &isPassedIsolation, &b_isPassedIsolation);
    fChain->SetBranchAddress("isPassedMassCut", &isPassedMassCut, &b_isPassedMassCut);

    //fChain->SetBranchAddress(".passTrig_HLT_2g20_tight_icaloloose_L12EM15VHI", &passTrig_HLT_2g20_tight_icaloloose_L12EM15VHI, &b__passTrig_HLT_2g20_tight_icaloloose_L12EM15VHI);
    fChain->SetBranchAddress("passTrig_HLT_2g22_tight_L12EM15VHI", &passTrig_HLT_2g22_tight_L12EM15VHI, &b_passTrig_HLT_2g22_tight_L12EM15VHI);
    fChain->SetBranchAddress("passTrig_HLT_g120_loose_L1EM22VHI", &passTrig_HLT_g120_loose_L1EM22VHI, &b_passTrig_HLT_g120_loose_L1EM22VHI);
    fChain->SetBranchAddress("passTrig_HLT_g140_loose_L1EM22VHI", &passTrig_HLT_g140_loose_L1EM22VHI, &b_passTrig_HLT_g140_loose_L1EM22VHI);
    fChain->SetBranchAddress("passTrig_HLT_g35_medium_g25_medium_L12EM20VH", &passTrig_HLT_g35_medium_g25_medium_L12EM20VH, &b_passTrig_HLT_g35_medium_g25_medium_L12EM20VH);


    fChain->SetBranchAddress("pT_y1", &pT_y1, &b_pT_y1);
    fChain->SetBranchAddress("pT_y2", &pT_y2, &b_pT_y2);







    
    Notify();
}

bool AnalyzeBackground_PC::Notify()
{
   // The Notify() function is called when a new file is opened. This                                                                          
   // can be either for a new TTree in a TChain or when when a new TTree                                                                       
   // is started when using PROOF. It is normally not necessary to make changes                                                                
   // to the generated code, but the routine can be extended by the                                                                            
   // user if needed. The return value is currently not used.                                                                                  

   return true;
}

void AnalyzeBackground_PC::Show(Long64_t entry)
{
// Print contents of entry.                                                                                                                    
// If entry is not specified, print current entry                                                                                              
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeBackground_PC::Cut(Long64_t entry)
{
// This function may be called from Loop.                                                                                                      
// returns  1 if entry is accepted.                                                                                                            
// returns -1 otherwise.                                                                                                                       
   return 1;
}
#endif // #ifdef AnalyzeBackground_PC_cxx  
