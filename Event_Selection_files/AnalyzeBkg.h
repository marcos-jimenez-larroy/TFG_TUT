#ifndef AnalyzeBkg_h
#define AnalyzeBkg_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class AnalyzeBkg {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain                                                                        
   Int_t           fCurrent; //!current Tree number in a TChain  

   Char_t          isPassed;
   Float_t         m_yy;
   Float_t         pT_yy;
   Float_t         weight;
   Float_t         crossSectionBRfilterEff;
   

   TBranch        *b_isPassed; //! 
   TBranch        *b_m_yy;   //!
   TBranch        *b_pT_yy;  //!
   TBranch        *b_weight;  //!
   TBranch        *b_crossSectionBRfilterEff; //!
  

   AnalyzeBkg(TTree *tree=0);
   virtual ~AnalyzeBkg();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalyzeBkg_cxx
AnalyzeBkg::AnalyzeBkg(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file                                                                              
// used to generate this class and read the Tree.                                                                                              
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("file_path");
      if (!f || !f->IsOpen()) {
         f = new TFile("file_path");
      }
      f->GetObject("BkgTree",tree);

   }
   Init(tree);
}

AnalyzeBkg::~AnalyzeBkg()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeBkg::GetEntry(Long64_t entry)
{
// Read contents of entry.                                                                                                                     
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeBkg::LoadTree(Long64_t entry)
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
void AnalyzeBkg::Init(TTree *tree)
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
    fChain->SetBranchStatus("isPassed",1);
    fChain->SetBranchStatus("m_yy",1);
    fChain->SetBranchStatus("pT_yy",1);
    fChain->SetBranchStatus("weight",1);
    fChain->SetBranchStatus("crossSectionBRfilterEff",1);



    fChain->SetBranchAddress("isPassed", &isPassed, &b_isPassed);
    fChain->SetBranchAddress("m_yy", &m_yy, &b_m_yy);
    fChain->SetBranchAddress("pT_yy", &pT_yy, &b_pT_yy);
    fChain->SetBranchAddress("weight", &weight, &b_weight);
    fChain->SetBranchAddress("crossSectionBRfilterEff", &crossSectionBRfilterEff, &b_crossSectionBRfilterEff);

    
    Notify();
}

bool AnalyzeBkg::Notify()
{
   // The Notify() function is called when a new file is opened. This                                                                          
   // can be either for a new TTree in a TChain or when when a new TTree                                                                       
   // is started when using PROOF. It is normally not necessary to make changes                                                                
   // to the generated code, but the routine can be extended by the                                                                            
   // user if needed. The return value is currently not used.                                                                                  

   return true;
}

void AnalyzeBkg::Show(Long64_t entry)
{
// Print contents of entry.                                                                                                                    
// If entry is not specified, print current entry                                                                                              
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeBkg::Cut(Long64_t entry)
{
// This function may be called from Loop.                                                                                                      
// returns  1 if entry is accepted.                                                                                                            
// returns -1 otherwise.                                                                                                                       
   return 1;
}
#endif // #ifdef AnalyzeBkg_cxx  
