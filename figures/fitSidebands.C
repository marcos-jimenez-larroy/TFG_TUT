#include "setup_categories_XGBoost_ttH.h"

void fitSidebands(TString datafile="files/skimmed_dataall.root"){
  gROOT->SetBatch(kTRUE);

  // Categories to fit
  std::vector<TString> categories = categories_XGBoost_ttH::catLabels;
  std::vector<TString> function = categories_XGBoost_ttH::bkgFunctions;
  std::vector<double> N_SS = categories_XGBoost_ttH::N_SS;

  // Functions to consider
  std::vector<TString> listOfFunctions = {"expPoly2","exponential","power","linear"};
  //std::vector<TString> listOfFunctions = {"expPoly2"};
  std::map<TString, Color_t> color;
  color["expPoly2"] = kOrange;
  color["exponential"] = kBlue;
  color["power"] = kRed;
  color["linear"] = kGreen+1;

  // Input tree containing data
  TFile infile(datafile);
  TTree *tree = (TTree*)infile.Get("tree");

  // Observable (same name than in tree)
  RooRealVar m_yy("m_yy","m_yy",105,160);
  m_yy.setBins(55);
  m_yy.setRange("TIlow",105,120);
  m_yy.setRange("TIhigh",130,160);
  m_yy.setRange("Full",105,160);
  m_yy.setRange("CBFitRange",115,135);


  std::map<TString,double> minbounds;
  std::map<TString,double> maxbounds;


  ifstream minfile("minBound.txt");
  if(minfile.good()){
    cout<< " minBound.txt exists "<<endl;

    std::string line;
    int index = -1;
    while( !minfile.eof() ){
      TString a;
      double b;
      minfile >> a >> b;
      minbounds[a]=b;

    }

  }
  ifstream maxfile("maxBound.txt");
  if(maxfile.good()){
    cout<< " maxBound.txt exists "<<endl;

    std::string line;
    int index = -1;
    while( !maxfile.eof() ){
      TString a;
      double b;
      maxfile >> a >> b;
      maxbounds[a]=b;

    }

  }

  for(TString cat: categories){
    cout<< "category "<<cat<<" minbound "<<minbounds[cat]<<std::endl;


  }

  // Control plots to check validity of fits
  TCanvas *c = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  c->Print("plot_bkgFit.pdf[");

  
  // Loop over categories and make fit on data sidebands
  for(int i=0; i<32; i++){

    
    int cat = categories_XGBoost_ttH::MxAODtoHfitter[i];
    if(cat==-999) continue;
    cout<<"mxaod "<<i <<" then "<<cat<<endl;
    m_yy.setRange("S90",minbounds[categories[cat]],maxbounds[categories[cat]]);
  
    cout<<"min "<<minbounds[categories[cat]]<<" max "<<maxbounds[categories[cat]]<<endl;
    
    std::cout << "\n\nMake fit on category: " << categories[cat] << std::endl;
    std::vector<double> yieldVec;
    std::vector<double> chi2Vec;
    std::vector<double> norm;
    std::vector<RooDataHist*> histVec;

    // Make histograms with TI sidebands only (blind signal region)

    ostringstream histo;
    histo.str(""); histo<<"hmyy_"<<i;
		     
    // TH1F *h = new TH1F(histo.str().c_str(),histo.str().c_str(),55,105,260);
    tree->Draw("m_yy/1000.>>hist",Form("(m_yy/1000.<120 || m_yy/1000.>130) && catCoup_idx==%i",i));
    TH1F *h = (TH1F*)gDirectory->Get("hist");
    cout<<"Integral "<<h->Integral()<<" catCoup  "<<i<<endl;

    if(h->Integral()==0.0){
      continue;
    }
    RooDataHist data("data","",RooArgList(m_yy),h);

    // Nominal frame for initial fits
    RooPlot *frame = m_yy.frame(RooFit::Title(categories[cat]));
    frame->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    frame->GetYaxis()->SetTitle("Events / GeV");
    data.plotOn(frame, RooFit::Range("TIlow,TIhigh"), RooFit::Name("data"));

    // Try all functions
    for(TString fct : listOfFunctions){

      // Only try expPoly2 if chosen by SS test
      if(fct=="expPoly2" && function[i]!="expPoly2"){
	yieldVec.push_back(0);
	chi2Vec.push_back(999);
	histVec.push_back(0);
	continue;
      }

      RooRealVar c0("c0","",-1.,-10000.,10000.);
      RooRealVar c1("c1","",-1.,-100.,10.);
      RooRealVar c2("c2","",-1.,-10.,10000.);
      //RooRealVar N("N","",1.0,0.,100000.);

      TString expression;
      RooArgList arguments(m_yy);
      if(fct=="exponential"){
	expression = "exp((m_yy-100.)/100.*c1)";
	arguments.add(c1);
      }
      else if(fct=="power"){
	expression = "pow(m_yy/100.,c1)";
	arguments.add(c1);
      }
      else if(fct=="linear"){
	expression = "c1*m_yy+c2";
	arguments.add(c1);
	c2.setRange(0.,10000.);
	arguments.add(c2);
      }
      else if(fct=="expPoly2"){
	expression = "exp((m_yy-100.)/100.*(c1+c2*(m_yy-100.)/100.))";
	//expression = "exp(c0+c1*m_yy+c2*m_yy*m_yy)"
	  //	f = ROOT.TF1(f_name, "TMath::Exp([0]+[1]*x+[2]*x*x)",xmin,xmax)
	arguments.add(c1);
	arguments.add(c2);
      }
      else{
	std::cout << "Not valid functional form: " << fct << std::endl;
	continue;
      }

      // Build PDF; normalised and fitted to sidebands only
      RooGenericPdf pdf(Form("pdf_%s",fct.Data()), expression, arguments);
      //pdf.setNormRange("TIlow,TIhigh");
      //  pdf.fitTo(data, RooFit::Range("TIlow,TIhigh"), RooFit::Minimizer("Minuit2"));
      RooFitResult *fr  = pdf.fitTo(data, RooFit::Range("TIlow,TIhigh"), RooFit::Minimizer("Minuit"), RooFit::Save());
      fr->Print();
      // Calculate total background normalisation, extrapolated from integral inside TI sidebands and write to datacard
      pdf.setNormRange("massRange");
      double normTIlow = pdf.createIntegral(m_yy,m_yy,"TIlow")->getVal();
      double normTIhigh = pdf.createIntegral(m_yy,m_yy,"TIhigh")->getVal();
      double normTot =  h->Integral()/(normTIlow+normTIhigh);
      // double normTot = ((RooRealVar*)fr->floatParsFinal().find("N"))->getVal();
      cout<<" value of norm "<<normTot<<endl;
      double normSignalRegion = pdf.createIntegral(m_yy,m_yy,"S90")->getVal();
      //double normSignalRegion = normTot - normTIlow - normTIhigh;

      
      
      pdf.plotOn(frame, RooFit::Range("TIlow,TIhigh"), RooFit::LineColor(color[fct]), RooFit::Name(Form("pdf_%s",fct.Data())));
      pdf.plotOn(frame, RooFit::Range(105,160), RooFit::LineStyle(kDashed), RooFit::LineColor(color[fct]));
      
      yieldVec.push_back(normSignalRegion);
      chi2Vec.push_back(frame->chiSquare(Form("pdf_%s",fct.Data()),"data"));
      histVec.push_back(pdf.generateBinned(m_yy,normTot,RooFit::ExpectedData()));
      norm.push_back(normTot);
    }

    frame->Draw();

    // Quantities for chosen function, based on smalles chi2
    TString selFct = "exponential";
    double bestChi2 = 1000;
    double yieldBest = 0;
    double normBest = 0;
    TH1 *selHist;
    for(int j=0; j<listOfFunctions.size(); j++){
      if(yieldVec[j]==0) continue;
      if(chi2Vec[j]<bestChi2){
	bestChi2 = chi2Vec[j];
	selFct = listOfFunctions[j];
	yieldBest = yieldVec[j];
	selHist = histVec[j]->createHistogram("m_yy");
	normBest = norm[j];
      }
    }

    TLatex l; l.SetTextSize(0.03); l.SetTextFont(42);
    for(int j=0; j<listOfFunctions.size(); j++){
      if(yieldVec[j]==0) continue;
      l.DrawLatexNDC(0.4,0.85-0.04*j,Form("#color[%i]{%s}",color[listOfFunctions[j]],listOfFunctions[j].Data()));
    }

    //selHist->Scale(h->Integral()/selHist->Integral());
    //selHist->Draw("same");
    
    l.SetTextSize(0.04);
    //l.DrawLatexNDC(0.65,0.85,Form("SS test: %s",function[i].Data()));
    //    histo.str(""); histo<<"Category: "<<categories.cat
    l.DrawLatexNDC(0.65,0.85,Form("Category: %s",categories[cat].Data()));
    //l.DrawLatexNDC(0.65,0.80,Form("From #chi^{2}: %s",selFct.Data()));
    l.DrawLatexNDC(0.65,0.8,Form("(#chi^{2}=%.2f)",bestChi2));
    l.DrawLatexNDC(0.65,0.75,Form("(Expected B90=%.2f)",normBest*yieldBest));


    c->Print("plot_bkgFit.pdf");
    
    /*
    // Compute difference between chosen function and others
    RooPlot *frameDiff = m_yy.frame(RooFit::Title(categories[i]));
    frameDiff->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    frameDiff->GetYaxis()->SetTitle("#DeltaN");


    // Define background PDF for best function
    RooRealVar c1("c1","",-1.,-100.,10.);
    RooRealVar c2("c2","",-1.,-10.,10000.);

    TString expression;
    RooArgList arguments(m_yy);
    if(selFct=="exponential"){
      expression = "exp((m_yy-100.)/100.*c1)";
      arguments.add(c1);
    }
    else if(selFct=="power"){
      expression = "pow(m_yy/100.,c1)";
      arguments.add(c1);
    }
    else if(selFct=="linear"){
      expression = "c1*m_yy+c2";
      arguments.add(c1);
      c2.setRange(0.,10000.);
      arguments.add(c2);
    }
    else if(selFct=="expPoly2"){
      expression = "exp((m_yy-100.)/100.*(c1+c2*(m_yy-100.)/100.))";
      arguments.add(c1);
      arguments.add(c2);
    }
    else{
      std::cout << "Not valid functional form: " << selFct << std::endl;
      continue;
    }

    // Build PDF; normalised and fitted to sidebands only
    RooGenericPdf bestBkgPDF("bestBkgPDF", expression, arguments);

    // compute uncertainty from fit and counting
    TString opt = "hist";
    double maxN_SS = 0;
    double maxYieldDiff = 0;
    for(int j=0; j<listOfFunctions.size(); j++){
      if(listOfFunctions[j]==selFct) continue;
      if(yieldVec[j]==0 || chi2Vec[j]>5) continue;
      TH1 *h = histVec[j]->createHistogram("m_yy");
      RooDataHist hist("histo_"+listOfFunctions[j],"",RooArgList(m_yy),h);
      hist.plotOn(frameDiff,RooFit::LineColor(color[listOfFunctions[j]]),RooFit::MarkerColor(color[listOfFunctions[j]]), RooFit::DataError(RooAbsData::None), RooFit::MarkerSize(0.2));
      // Signal PDF
      RooRealVar mu("mu","",125090);
      RooRealVar sigma("sigma","",1.7);
      RooRealVar alpha("alpha","",1.5);
      RooRealVar n("n","",3);
      RooCBShape CB("CB","",m_yy,mu,sigma,alpha,n);
      // Signal + background PDF
      RooRealVar N_SS("N_SS","",10,-1000,1000);
      RooRealVar N_bkg("N_bkg","",100,-1000,100000000);
      RooAddPdf addPDF("addPDF","",RooArgList(CB,bestBkgPDF),RooArgList(N_SS,N_bkg));
      addPDF.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::SumW2Error(kFALSE));
      addPDF.plotOn(frameDiff, RooFit::LineColor(color[listOfFunctions[j]]), RooFit::Name(Form("CBext_%s",listOfFunctions[j].Data())));
      maxN_SS = std::max(maxN_SS, N_SS.getValV());
      maxYieldDiff = std::max(maxYieldDiff, fabs(yieldBest-yieldVec[j]));
    }

    frameDiff->Draw();
    l.DrawLatexNDC(0.50,0.85,Form("max(N_{sig}) = %.2f",maxN_SS));
    l.DrawLatexNDC(0.50,0.78,Form("max#DeltaN([120,130]GeV) = %.2f",maxYieldDiff));
    l.DrawLatexNDC(0.50,0.70,Form("N_SS (current) = %.2f",N_SS[i]));
    TLine line;
    line.DrawLine(105,0,160,0);
    c->Print("plot_bkgFit.pdf");

    delete h;*/
  }


  // Close all opened files
  infile.Close();
  c->Print("plot_bkgFit.pdf]");
  delete c;

}
