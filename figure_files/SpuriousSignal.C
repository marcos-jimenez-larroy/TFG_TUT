#include "TLorentzVector.h"
#include <vector>
//#include "MyEvent.h"

//#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
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

int BackgroundFit(TH1* h0yymass,double xmin, double xmax, TString func, double *N_signal, double *N_err,double *chi2)
{
    RooRealVar hMass("myy","myy [GeV]",xmin,xmax);
    RooDataHist hist("hist","hist",hMass,Import(*h0yymass));


//Function definitions:

//Double-Sided-Crystal-Ball
double par[6]; //These parameters result from the previous fit of the signal-only MC
par[0] = 125.49;
par[1] = 1.77;
par[2] = 1.72;
par[3] = 4.6;
par[4] = 1.54;
par[5] = 20.2;

//Import data



//DCB parameters (Not allowed to change, only normalization)
RooRealVar mu("mu","mu",par[0]);
RooRealVar width("width","width",par[1]);
RooRealVar a1("a1","a1",par[2]);
RooRealVar p1("p1","p1",par[3]);
RooRealVar a2("a2","a2",par[4]);
RooRealVar p2("p2","p2",par[5]);

RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);

//Normalization coefficients

RooRealVar N_SS("N_SS","",10,-1000,1000);
RooRealVar N_bkg("N_bkg","",100,-1000,100000000);


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



if (func == "ExpPol1")
{


    RooAddPdf FitFunction("FitFunction","FitFunction",RooArgList(dcbPdf,ExpPol1),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"),RooFit::Save(true));
    result->Print();

    N_signal[0] = N_SS.getValV();
    N_err[0] = N_SS.getError();

    auto pl = hMass.frame(RooFit::Name("ExpPol1"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 

    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[0] = pl->chiSquare("FitFunction","hist");

    pl->Draw();

    TLatex l; l.SetTextSize(0.04); l.SetTextFont(42);
    l.DrawLatexNDC(0.55,0.8,"ExpPol1");
    l.DrawLatexNDC(0.55,0.7,Form("#chi^{2}/NDF = %.2f",chi2[0]/1.0));


    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol1_bin0.pdf");

    return 0;

}

else if(func == "Bernstein")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,BernsteinPol),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"),RooFit::Save(true));
    result->Print();

    N_signal[1] = N_SS.getValV();
    N_err[1] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("Bernstein"));
   
    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[1] = pl->chiSquare("FitFunction","hist");
    pl->Draw();

    TLatex l; l.SetTextSize(0.04); l.SetTextFont(42);
    l.DrawLatexNDC(0.55,0.8,"Bernstein");
    l.DrawLatexNDC(0.55,0.7,Form("#chi^{2}/NDF = %.4f",chi2[1]/4.0));

    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/BernsteinPol_bin0.pdf");

    return 0;
}
else if(func=="Power")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,PowerFunction),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    result->Print();

    N_signal[2] = N_SS.getValV();
    N_err[2] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("Power"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[2] = pl->chiSquare("FitFunction","hist");
    pl->Draw();

    TLatex l; l.SetTextSize(0.04); l.SetTextFont(42);
    l.DrawLatexNDC(0.55,0.8,"Power");
    l.DrawLatexNDC(0.55,0.7,Form("#chi^{2}/NDF = %.2f",chi2[2]/1.0));

    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/PowerFunction_bin0.pdf");
    
    return 0;
}
else if(func=="ExpPol2")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,ExpPol2),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    result->Print();

    N_signal[3] = N_SS.getValV();
    N_err[3] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("ExpPol2"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[3] = pl->chiSquare("FitFunction","hist");
    pl->Draw();

    TLatex l; l.SetTextSize(0.04); l.SetTextFont(42);
    l.DrawLatexNDC(0.55,0.8,"ExpPol2");
    l.DrawLatexNDC(0.55,0.7,Form("#chi^{2}/NDF = %.2f",chi2[3]/2.0));


    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol2_bin0.pdf");
    
    return 0;
}
else if(func=="ExpPol3")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,ExpPol3),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    result->Print();

    N_signal[4] = N_bkg.getValV(); //Ojo cambiar. Es un check rápido
    N_err[4] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("ExpPol3"));

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 

    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[4] = pl->chiSquare("FitFunction","hist");

    pl->Draw();

    TLatex l; l.SetTextSize(0.04); l.SetTextFont(42);
    l.DrawLatexNDC(0.55,0.8,"ExpPol3");
    l.DrawLatexNDC(0.55,0.7,Form("#chi^{2}/NDF = %.2f",chi2[4]/3.0));

    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol3_bin0.pdf");

    /*TFile* file = new TFile("/home/mjlarroy/TFG_TUT/figures/fitresultexp3.root","RECREATE");
    result->Write("NAME");
    file->Close();*/ 
    
    //Para salvar el resultado del fit, aunque ahora preguntaré qué manera es la mejor para luego trabajarlo
    
    return 0;
}
cout<<"Function incorrectly specified";
return 0;


}

void SpuriousSignal() //Fit all at the same time
{
  TFile file_bins("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_background_bins.root","READ");
  TFile file("/home/mjlarroy/TFG_TUT/Root_Results/MC_reweighted_background_histos.root","READ");
  TFile file_signal("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_histos.root","READ");
  
  TH1D *hmyy_signal = (TH1D*)file_signal.Get("h0yymass");
  TH1D *hmyy = (TH1D*)file_bins.Get("h_myy_pT_0");
  //TH1D *hmyy = (TH1D*)file.Get("histmc_bkg_fit");



  //N_expected from signal (120.125 and 129.75 hard-coded because we used 80 bins)

  double jmin = (120.125 - 105.0)/hmyy_signal->GetBinWidth(1); 
  double jmax = (129.75-105.0)/hmyy_signal->GetBinWidth(1);

  double N_exp = hmyy_signal->Integral(jmin,jmax); //Must coincide with the N_bkg in every fit


  std::vector<TString> funcList = {"ExpPol1","Bernstein","Power","ExpPol2","ExpPol3"};

  double N_signal[5];
  double chi2 [5];
  double N_err[5];

  for(int i=0;i<5;i++){
    BackgroundFit(hmyy, 105., 160., funcList[i],N_signal,N_err,chi2);
  }


  cout<<endl;
  cout<<"Ratio N_SS/N_expected and N_SS/N_SS_error for each function (in pairs)"<<endl;
  for(int i = 0; i<5;i++){
    cout<<funcList[i]<<endl;
    cout<<"NSS/N_expected ="<<abs(N_signal[i])/N_exp<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"NSS_error/NSS ="<<N_err[i]/abs(N_signal[i])<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"Chi2 ="<<chi2[i]<<endl;
    cout<<N_signal[i]<<endl;
    cout<<endl;
  }


}