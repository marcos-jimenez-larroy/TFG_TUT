// .L Analyze_data.C
// .x runtest_data22_23.C
void runtest_data22_23()
{
  TChain *p1 = new TChain("CollectionTree","chain");
  p1->Add("/dir/data22/data22_13p6TeV.periodAllYear.physics_Main.MxAODDetailed.p6456.h031.001.root");
  p1->Add("/dir/data22/data22_13p6TeV.periodAllYear.physics_Main.MxAODDetailed.p6456.h031.002.root");
  p1->Add("/dir/data22/data22_13p6TeV.periodAllYear.physics_Main.MxAODDetailed.p6456.h031.003.root");
  p1->Add("/dir/data22/data22_13p6TeV.periodAllYear.physics_Main.MxAODDetailed.p6456.h031.004.root");
  p1->Add("/dir/data22/data22_13p6TeV.periodAllYear.physics_Main.MxAODDetailed.p6456.h031.005.root");
  

  //
  // run over the files
  TTree* mytree = (TTree*)p1;
  //Analyze_data t(mytree); //Depending on the purpose!
  Analyze_data_bins t(mytree);
  //
  t.Loop();
}