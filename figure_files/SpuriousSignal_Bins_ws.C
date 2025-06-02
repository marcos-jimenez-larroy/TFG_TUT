#include "TLorentzVector.h"
#include <vector>
//#include "MyEvent.h"

//#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
//#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
//#include "RooGenericPdf.h"
#include "RooTFnBinding.h"

//Script to fit a function+signal parametrization to background MC only! (RooAddPdf)
//Look which function makes the signal the smallest

//Conditions:

//Nsp (Spurious Signal) < 10% of expected events (350 con RooFit en el PDF 'check.pdf')
//Nsp < 20% total statistical error del MC



using namespace RooFit;





RooWorkspace * BackgroundFit(TH1* h0yymass,double xmin, double xmax, TString func, double *p_value,RooWorkspace* w_in, double *chi2a)
{


    RooRealVar hMass("myy","myy [GeV]",xmin,xmax);
    RooDataHist hist("hist","hist",hMass,Import(*h0yymass));


    int nBins = 55;

    RooWorkspace * w_out = new RooWorkspace("w_out") ; 

    w_out->import(hist); //So I do not save it five times
    

    w_in->loadSnapshot("postfit");

    RooFitResult *result = (RooFitResult*)w_in->obj("fit_result");

    RooArgList SignalParameterList = result->floatParsFinal();

//Function definitions:

//Double-Sided-Crystal-Ball
double par[6]; //These parameters result from the previous fit of the signal-only MC


    RooRealVar* var_mu = (RooRealVar*) SignalParameterList.find("mu");
    RooRealVar* var_width = (RooRealVar*) SignalParameterList.find("width");
    RooRealVar* var_a1 = (RooRealVar*) SignalParameterList.find("a1");
    RooRealVar* var_p1 = (RooRealVar*) SignalParameterList.find("p1");
    RooRealVar* var_a2 = (RooRealVar*) SignalParameterList.find("a2");
    RooRealVar* var_p2 = (RooRealVar*) SignalParameterList.find("p2");

    par[0] = var_mu->getVal();
    par[1] = var_width->getVal();
    par[2] = var_a1->getVal();
    par[3] = var_p1->getVal();
    par[4] = var_a2->getVal();
    par[5] = var_p2->getVal();







//DCB parameters (Not allowed to change, only normalization)
RooRealVar mu("mu","mu",par[0]);
RooRealVar width("width","width",par[1]);
RooRealVar a1("a1","a1",par[2]);
RooRealVar p1("p1","p1",par[3]);
RooRealVar a2("a2","a2",par[4]);
RooRealVar p2("p2","p2",par[5]);

RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);

//Exponential

RooRealVar c("c","c",-1.0,-10.0,0.0);

RooExponential ExpPol1("ExpPol1","ExpPol1",hMass,c);

//Bernstein Polynomial

double coef[4]; //These parameters result from the previous fit of the signal-only MC
coef[0] = 1.0;
coef[1] = 1.0;
coef[2] = 1.0;
coef[3] = 1.0;

//Bernstein parameters

RooRealVar c0("c0","c0",coef[0],-5.0,5.0);
RooRealVar c1("c1","c1",coef[1],-5.0,5.0);
RooRealVar c2("c2","c2",coef[2],-5.0,5.0);
RooRealVar c3("c3","c3",coef[3],-5.0,5.0);


RooBernstein BernsteinPol("BernsteinPol","BernsteinPol",hMass,RooArgList(c0,c1,c2,c3));

//Power Function (Simple polynomial)

RooRealVar powerExponent1("powerExponent1","powerExponent1",-1.0,-30.0,10.0);



RooGenericPdf PowerFunction("PowerFunction","TMath::Power(x[0]/100.,x[1])",RooArgList(hMass, powerExponent1));

//Second Order Polynomial in Exponential (ExpPol2)

RooRealVar b1("b1","b1",-1.0,-10.0,0.0);
RooRealVar b2("b2","b2",-1.0,-10.0,0.0);

RooGenericPdf ExpPol2("ExpPol2","TMath::Exp(x[1]*TMath::Power(x[0]-100.,1)/100. + x[2]*TMath::Power(x[0]-100.,2)/10000.)",RooArgList(hMass,b1,b2));

//Third Order Polynomial in Exponential (ExpPol3)

RooRealVar d1("d1","d1",-1.0,-50.0,50.0);
RooRealVar d2("d2","d2",-1.0,-50.0,50.0);
RooRealVar d3("d3","d3",-1.0,-50.0,50.0);

RooGenericPdf ExpPol3("ExpPol3","TMath::Exp(x[1]*TMath::Power(x[0]-100.,1)/100. + x[2]*TMath::Power(x[0]-100.,2)/10000. + x[3]*TMath::Power(x[0]-100.,3)/1000000.)",RooArgList(hMass,d1,d2,d3));


//Fit the functions to our histograms



  //Normalization coefficients

RooRealVar N_SS("N_SS","N_SS",10,-1000,1000);
RooRealVar N_bkg("N_bkg","N_bkg",100,-1000,100000000);


if (func == "ExpPol1")
{


    RooAddPdf FitFunctionExpPol1("FitFunctionExpPol1","FitFunctionExpPol1",RooArgList(dcbPdf,ExpPol1),RooArgList(N_SS,N_bkg));

    w_out->import(FitFunctionExpPol1) ; 
    //object that includes all parameters of the pdf (and their errors)
    RooArgSet * fitParams = w_out->pdf("FitFunctionExpPol1")->getParameters(*w_out->data("hist")) ;
    //snapshot of the parameters before the fit
    w_out->saveSnapshot("prefit" , *fitParams ) ; 
    //replaced auto with RooFitResult
    RooFitResult * result = w_out->pdf("FitFunctionExpPol1")->fitTo(*w_out->data("hist"), RooFit::Save(true));
    result->Print();
    w_out->import(*result,"ExpPol1_fit_result") ;
    //snapshot of the parameters after the fit
    w_out->saveSnapshot("postfit" , *fitParams ) ; 

    


    RooAbsReal* chi2 = w_out->pdf("FitFunctionExpPol1")->createChi2(hist, Range("fullRange"),Extended(true), DataError(RooAbsData::Poisson));
    
    chi2a[0] = chi2->getVal();
    p_value[0] = TMath::Prob(chi2->getVal(),nBins-1.0);//TMath::Prob(chi2,NDF)



    auto pl = hMass.frame(RooFit::Name("ExpPol1"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 

    hist.plotOn(pl,RooFit::Name("hist"));


    w_out->pdf("FitFunctionExpPol1")->plotOn(pl,RooFit::Name("FitFunctionExpPol1"),RooFit::LineColor(kRed));

    RooArgList SignalParameterList1 = result->floatParsFinal();
    RooRealVar* N_Sig = (RooRealVar*) SignalParameterList1.find("N_SS");
    RooRealVar* N_Bkg = (RooRealVar*) SignalParameterList1.find("N_bkg");
    pl->GetYaxis()->SetTitle("Events/GeV");
    pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");
    pl->GetXaxis()->SetTitleSize(0.06);
    pl->GetXaxis()->SetTitleOffset(1);
    pl->Draw();
    TLatex l; l.SetTextSize(0.07); l.SetTextFont(42);
    l.DrawLatexNDC(0.42,0.86,"ExpPol1");

    l.DrawLatexNDC(0.20,0.46,Form("N_{SS} = %.2f",N_Sig->getVal()));
    l.DrawLatexNDC(0.20,0.38,Form("N_{bkg} = %.2f",N_Bkg->getVal()));
    l.DrawLatexNDC(0.20,0.30,Form("Bkg uncertainty = %.2f",N_Sig->getError()));
    l.DrawLatexNDC(0.20,0.22,Form("#chi^{2} = %.2f",chi2a[0]));


    
    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol1_bin9.pdf");

    
    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol1.pdf");
    return w_out;
}

else if(func == "Bernstein")
{

    RooAddPdf FitFunctionBernstein("FitFunctionBernstein","FitFunctionBernstein",RooArgList(dcbPdf,BernsteinPol),RooArgList(N_SS,N_bkg));

    w_out->import(FitFunctionBernstein); //Problem here

    //object that includes all parameters of the pdf (and their errors)
    RooArgSet * fitParams = w_out->pdf("FitFunctionBernstein")->getParameters(*w_out->data("hist")) ;
    //snapshot of the parameters before the fit
    w_out->saveSnapshot("prefit" , *fitParams ) ; 
    //replaced auto with RooFitResult
    RooFitResult * result = w_out->pdf("FitFunctionBernstein")->fitTo(*w_out->data("hist"), RooFit::Save(true));
    result->Print();
    w_out->import(*result,"Bernstein_fit_result") ;
    //snapshot of the parameters after the fit
    w_out->saveSnapshot("postfit" , *fitParams ) ; 


    auto pl = hMass.frame(RooFit::Name("Bernstein Polynomial"));
   
    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl,RooFit::Name("hist"));


    w_out->pdf("FitFunctionBernstein")->plotOn(pl,RooFit::Name("FitFunctionBernstein"),RooFit::LineColor(kRed));

    RooAbsReal* chi2 = w_out->pdf("FitFunctionBernstein")->createChi2(hist, Range("fullRange"),Extended(true), DataError(RooAbsData::Poisson));
    chi2a[1] = chi2->getVal();
    p_value[1] = TMath::Prob(chi2->getVal(),nBins-4.0);//TMath::Prob(chi2,NDF)


    RooArgList SignalParameterList1 = result->floatParsFinal();
    RooRealVar* N_Sig = (RooRealVar*) SignalParameterList1.find("N_SS");
    RooRealVar* N_Bkg = (RooRealVar*) SignalParameterList1.find("N_bkg");
    pl->GetYaxis()->SetTitle("Events/GeV");
    pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");
    pl->GetXaxis()->SetTitleSize(0.06);
    pl->GetXaxis()->SetTitleOffset(1);
    pl->Draw();
    TLatex l; l.SetTextSize(0.07); l.SetTextFont(42);
    l.DrawLatexNDC(0.42,0.86,"Bernstein");

    l.DrawLatexNDC(0.20,0.46,Form("N_{SS} = %.2f",N_Sig->getVal()));
    l.DrawLatexNDC(0.20,0.38,Form("N_{bkg} = %.2f",N_Bkg->getVal()));
    l.DrawLatexNDC(0.20,0.30,Form("Bkg uncertainty = %.2f",N_Sig->getError()));
    l.DrawLatexNDC(0.20,0.22,Form("#chi^{2} = %.2f",chi2a[1]));



    

    

    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/Bernstein_bin9.pdf");
    
    
    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/BernsteinPol.pdf");

    return w_out;
}
else if(func=="Power")
{
    RooAddPdf FitFunctionPower("FitFunctionPower","FitFunctionPower",RooArgList(dcbPdf,PowerFunction),RooArgList(N_SS,N_bkg));

    w_out->import(FitFunctionPower) ; 
    //object that includes all parameters of the pdf (and their errors)
    RooArgSet * fitParams = w_out->pdf("FitFunctionPower")->getParameters(*w_out->data("hist")) ;
    //snapshot of the parameters before the fit
    w_out->saveSnapshot("prefit" , *fitParams ) ; 
    //replaced auto with RooFitResult


   // RooBinSamplingPdf binSampler("PowerBinSampler", "PBS", hMass, FitFunctionPower);

    //RooFitResult * result =binSampler.fitTo(*w_out->data("hist"), RooFit::Save(true));


    RooFitResult * result = w_out->pdf("FitFunctionPower")->fitTo(*w_out->data("hist"), RooFit::Save(true));
    result->Print();

    RooHistFunc histFunc("histFunc", "interpolated data", hMass, hist);

    RooArgList SignalParameterList1 = result->floatParsFinal();
    RooRealVar* N_Sig = (RooRealVar*) SignalParameterList1.find("N_SS");
    RooRealVar* N_Bkg = (RooRealVar*) SignalParameterList1.find("N_bkg");
    
    /*double chi2_mano = 0.0;
    for(int i=1;i<=55;i++){
      hMass.setVal(104.5+i);
      hMass.setRange(Form("bin_%i",i),104 +i,105+i);
      double ex_val = histFunc.getVal(hMass);
      double ob_val = N_Bkg->getVal()*w_out->pdf("FitFunctionPower")->createIntegral(hMass,hMass,Form("bin_%i",i))->getVal();
      chi2_mano += (ex_val-ob_val)*(ex_val-ob_val)/ex_val;

      std::cout<<ob_val<<std::endl;


    }
    std::cout<<"Chi 2 a mano: "<<chi2_mano<<std::endl;*/


    w_out->import(*result,"Power_fit_result") ;
    //snapshot of the parameters after the fit
    w_out->saveSnapshot("postfit" , *fitParams ) ; 





    





    auto pl = hMass.frame(RooFit::Name("Power Function"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl,RooFit::Name("hist"));

    RooAbsReal* chi2 = w_out->pdf("FitFunctionPower")->createChi2(hist, Range("fullRange"),Extended(true), DataError(RooAbsData::Poisson));
    chi2a[2] = chi2->getVal();
    p_value[2] = TMath::Prob(chi2->getVal(),nBins-1.0);//TMath::Prob(chi2,NDF)

    


    w_out->pdf("FitFunctionPower")->plotOn(pl,RooFit::Name("FitFunctionPower"),RooFit::LineColor(kRed));
    
    pl->GetYaxis()->SetTitle("Events/GeV");
    pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");
    pl->GetXaxis()->SetTitleSize(0.06);
    pl->GetXaxis()->SetTitleOffset(1);
    pl->Draw();
    TLatex l; l.SetTextSize(0.07); l.SetTextFont(42);
    l.DrawLatexNDC(0.42,0.86,"Power");

    l.DrawLatexNDC(0.20,0.46,Form("N_{SS} = %.2f",N_Sig->getVal()));
    l.DrawLatexNDC(0.20,0.38,Form("N_{bkg} = %.2f",N_Bkg->getVal()));
    l.DrawLatexNDC(0.20,0.30,Form("Bkg uncertainty = %.2f",N_Sig->getError()));
    l.DrawLatexNDC(0.20,0.22,Form("#chi^{2} = %.2f",chi2a[2]));

    

    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/Power_bin9.pdf");


    


    
    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/PowerFunction.pdf");
    
    return w_out;

}
else if(func=="ExpPol2")
{
    RooAddPdf FitFunctionExpPol2("FitFunctionExpPol2","FitFunctionExpPol2",RooArgList(dcbPdf,ExpPol2),RooArgList(N_SS,N_bkg));

    w_out->import(FitFunctionExpPol2) ; 
    //object that includes all parameters of the pdf (and their errors)
    RooArgSet * fitParams = w_out->pdf("FitFunctionExpPol2")->getParameters(*w_out->data("hist")) ;
    //snapshot of the parameters before the fit
    w_out->saveSnapshot("prefit" , *fitParams ) ; 
    //replaced auto with RooFitResult
    RooFitResult * result = w_out->pdf("FitFunctionExpPol2")->fitTo(*w_out->data("hist"), RooFit::Save(true));
    result->Print();
    w_out->import(*result,"ExpPol2_fit_result") ;
    //snapshot of the parameters after the fit
    w_out->saveSnapshot("postfit" , *fitParams ) ; 


    auto pl = hMass.frame(RooFit::Name("ExpPol2"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl,RooFit::Name("hist"));

    RooAbsReal* chi2 = w_out->pdf("FitFunctionExpPol2")->createChi2(hist, Range("fullRange"),Extended(true), DataError(RooAbsData::Poisson));
    chi2a[3] = chi2->getVal();
    p_value[3] = TMath::Prob(chi2->getVal(),nBins-2.0);//TMath::Prob(chi2,NDF)


    w_out->pdf("FitFunctionExpPol2")->plotOn(pl,RooFit::Name("FitFunctionExpPol2"),RooFit::LineColor(kRed));
    RooArgList SignalParameterList1 = result->floatParsFinal();
    RooRealVar* N_Sig = (RooRealVar*) SignalParameterList1.find("N_SS");
    RooRealVar* N_Bkg = (RooRealVar*) SignalParameterList1.find("N_bkg");
    pl->GetYaxis()->SetTitle("Events/GeV");
    pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");
    pl->GetXaxis()->SetTitleSize(0.06);
    pl->GetXaxis()->SetTitleOffset(1);
    pl->Draw();
    TLatex l; l.SetTextSize(0.07); l.SetTextFont(42);
    l.DrawLatexNDC(0.42,0.86,"ExpPol2");

    l.DrawLatexNDC(0.20,0.46,Form("N_{SS} = %.2f",N_Sig->getVal()));
    l.DrawLatexNDC(0.20,0.38,Form("N_{bkg} = %.2f",N_Bkg->getVal()));
    l.DrawLatexNDC(0.20,0.30,Form("Bkg uncertainty = %.2f",N_Sig->getError()));
    l.DrawLatexNDC(0.20,0.22,Form("#chi^{2} = %.2f",chi2a[3]));

    

    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol2_bin9.pdf");
    
    
    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol2.pdf");
    

    return w_out;
}
else if(func=="ExpPol3")
{
    RooAddPdf FitFunctionExpPol3("FitFunctionExpPol3","FitFunctionExpPol3",RooArgList(dcbPdf,ExpPol3),RooArgList(N_SS,N_bkg));

    w_out->import(FitFunctionExpPol3) ; 
    //object that includes all parameters of the pdf (and their errors)
    RooArgSet * fitParams = w_out->pdf("FitFunctionExpPol3")->getParameters(*w_out->data("hist")) ;
    //snapshot of the parameters before the fit
    w_out->saveSnapshot("prefit" , *fitParams ) ; 
    //replaced auto with RooFitResult
    RooFitResult * result = w_out->pdf("FitFunctionExpPol3")->fitTo(*w_out->data("hist"), RooFit::Save(true));
    result->Print();
    w_out->import(*result,"ExpPol3_fit_result") ;
    //snapshot of the parameters after the fit
    w_out->saveSnapshot("postfit" , *fitParams ) ; 



    auto pl = hMass.frame(RooFit::Name("ExpPol3"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 

    hist.plotOn(pl,RooFit::Name("hist"));

    RooAbsReal* chi2 = w_out->pdf("FitFunctionExpPol3")->createChi2(hist, Range("fullRange"),Extended(true), DataError(RooAbsData::Poisson));
    chi2a[4] = chi2->getVal();
    p_value[4] = TMath::Prob(chi2->getVal(),nBins-3.0);//TMath::Prob(chi2,NDF)


    w_out->pdf("FitFunctionExpPol3")->plotOn(pl,RooFit::Name("FitFunctionExpPol3"),RooFit::LineColor(kRed));
    RooArgList SignalParameterList1 = result->floatParsFinal();
    RooRealVar* N_Sig = (RooRealVar*) SignalParameterList1.find("N_SS");
    RooRealVar* N_Bkg = (RooRealVar*) SignalParameterList1.find("N_bkg");
    pl->GetYaxis()->SetTitle("Events/GeV");
    pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");
    pl->GetXaxis()->SetTitleSize(0.06);
    pl->GetXaxis()->SetTitleOffset(1);
    pl->Draw();
    TLatex l; l.SetTextSize(0.07); l.SetTextFont(42);
    l.DrawLatexNDC(0.42,0.86,"ExpPol3");

    l.DrawLatexNDC(0.20,0.46,Form("N_{SS} = %.2f",N_Sig->getVal()));
    l.DrawLatexNDC(0.20,0.38,Form("N_{bkg} = %.2f",N_Bkg->getVal()));
    l.DrawLatexNDC(0.20,0.30,Form("Bkg uncertainty = %.2f",N_Sig->getError()));
    l.DrawLatexNDC(0.20,0.22,Form("#chi^{2} = %.2f",chi2a[4]));

   

    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol3_bin9.pdf");
    
    

    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol3.pdf");
    
    return w_out;

}

else{
  cout<<"Function Not Defined"<<endl;
  return nullptr;
}

}




Double_t SpuriousSignal(TString histname) //Fit all at the same time
{
  TString thePath = "/home/mjlarroy/TFG_TUT/Root_Results/";

  TFile file(thePath + "mc23a_background_bins.root","READ");
  TFile file_signal(thePath + "mc23a_signal_bins.root","READ");

  TString OutputName= "/outFile_Signal/outFileSignal_" + histname + ".root";
  TString WorkspaceName= "w_" + histname;

  TFile file_ws(thePath + OutputName,"READ");
  RooWorkspace *w = (RooWorkspace*)file_ws.Get(WorkspaceName);


  std::vector<TString> funcList = {"ExpPol1","Bernstein","Power","ExpPol2","ExpPol3"};
  
  TH1D *hmyy_signal = (TH1D*)file_signal.Get(histname);
  TH1D *hmyy = (TH1D*)file.Get(histname);


  //N_expected from signal (120. and 130. hard-coded because we used 55 bins)

  double jmin = (120. - 105.0)/hmyy_signal->GetBinWidth(1); 
  double jmax = (130.-105.0)/hmyy_signal->GetBinWidth(1);

  double N_expected = hmyy_signal->Integral(jmin,jmax); //Must coincide with the N_bkg in every fit


  double p_value [5];
  double chi2[5];
  double N_signal[5]; double N_err[5];
  
  TString outTag = histname;

  for(int i=0;i<funcList.size();i++){

  
  TFile * outfile = new TFile( thePath +"outFile_Background/" +"outFileBackground_"+ funcList[i] + "_" + outTag + ".root" , "update" ) ;

  RooWorkspace * w_out = BackgroundFit(hmyy, 105., 160., funcList[i], p_value, w, chi2);

  w_out->loadSnapshot("postfit");
  RooFitResult *result = (RooFitResult*)w_out->obj(funcList[i] + "_fit_result");

  RooArgList SignalParameterList = result->floatParsFinal();
  RooRealVar* N_SS = (RooRealVar*) SignalParameterList.find("N_SS");

  N_signal[i] = N_SS->getVal();
  N_err[i] = N_SS->getError();



  //writes workspace to file, overwriting previous one with identical name
  w_out->Write("w_out_" +funcList[i] + "_" + histname , TObject::kOverwrite ) ;
  outfile->Close() ; 
  }


  //Imprimir las cosas en pantalla para ver parámetros de cada uno pero comentado cuando aparecen todfos los histogramas.
  
  cout<<endl;
  cout<<"Ratio N_SS/N_expected and N_SS/N_SS_error for each function (in pairs)"<<endl;
  for(int i = 0; i<5;i++){
    cout<<funcList[i]<<endl;
    cout<<"NSS/N_expected ="<<abs(N_signal[i])/N_expected<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"NSS/NSS_error ="<<abs(N_signal[i])/N_err[i]<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"Chi2 ="<<chi2[i]<<endl;
    cout<<endl;
  }

  //First of all, if the exponential (0) or the power fucntion (2) are enoguh we will use them

  Double_t best_func = 0.0;

  

    //Posibilidad para relajar el constraint


  if ((std::abs(N_signal[0])/N_expected <0.14 || std::abs(N_signal[0])/N_err[0] < 0.2) && p_value[0]>0.01){
    if((std::abs(N_signal[2])/N_expected < 0.14 || std::abs(N_signal[2])/N_err[2] < 0.2) && p_value[2]>0.01){
        if (std::abs(N_signal[0])<std::abs(N_signal[2])) {
            best_func = N_err[0];
        }
        else{
            best_func = N_err[2];
        }
    }
    else{
      best_func = N_err[0];
    }
  }
  else if ((std::abs(N_signal[2])/N_expected < 0.14 || std::abs(N_signal[2])/N_err[2] < 0.2) && p_value[2]>0.01){
    best_func =N_err[2];
  }
  else if ((std::abs(N_signal[3])/N_expected < 0.14 || std::abs(N_signal[3])/N_err[3] < 0.2) && p_value[3]>0.01){
    best_func = N_err[3];
  }
  else if ((std::abs(N_signal[4])/N_expected < 0.14 || std::abs(N_signal[4])/N_err[4] < 0.2) && p_value[4]>0.01){
    best_func = N_err[4];
  }
  else if ((std::abs(N_signal[1])/N_expected < 0.14 || std::abs(N_signal[1])/N_err[1] < 0.2) && p_value[1]>0.01){
    best_func = N_err[1];
    
  }

  if(histname == "h_myy_pT_0"){
    for (int i=9;i<10;i++){
      std::cout<<"NSS "<<N_signal[i]<<std::endl;
      std::cout<<"NSS/Nexpected: "<<N_signal[i]/N_expected<<std::endl;
      std::cout<<"Nerr: "<<N_err[i]<<std::endl;
      std::cout<<"Chi2: "<<chi2[i]<<std::endl;
      std::cout<<"p-value: "<<p_value[i]<<std::endl;


    }
  }

  return best_func;
  

  





}


void SpuriousSignal_Bins_ws()
{

    Double_t BestFunctions[16];

    for (int i=9;i<10;i++){
        std::string histName= "h_myy_pT_" + std::to_string(i);

        BestFunctions[i] = SpuriousSignal(histName);

    }
    for(int i=0;i<16;i++){
        cout<<"h_myy_pT_"<<i<<": "<<BestFunctions[i]<<endl;

    }
}