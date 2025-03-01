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

int BackgroundFit(TH1* h0yymass,double xmin, double xmax, TString func, double *N_signal, double *chi2)
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
RooRealVar mu("mu","mu",par[0],par[0]-0.001,par[0]+0.001);
RooRealVar width("width","width",par[1],par[1]-0.001,par[1]+0.001);
RooRealVar a1("a1","a1",par[2],par[2]-0.001,par[2]+0.001);
RooRealVar p1("p1","p1",par[3],par[3]-0.001,par[3]+0.001);
RooRealVar a2("a2","a2",par[4],par[4]-0.001,par[4]+0.001);
RooRealVar p2("p2","p2",par[5],par[5]-0.001,par[5]+0.001);

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

RooRealVar powerExponent1("powerExponent1","powerExponent1",-1.0,-3.0,0.0);
RooRealVar powerExponent2("powerExponent2","powerExponent2",-1.0,-3.0,0.0);

RooRealVar constant("constant","constant",1.0,0.0,1000.0);

RooRealVar A1("A1","A1",1.0,0.0,100.0);
RooRealVar A2("A2","A2",1.0,0.0,100.0);


RooGenericPdf PowerFunction("PowerFunction","x[3]*TMath::Power(x[0],x[1])+x[4]*TMath::Power(x[0],x[2]) + x[5]",RooArgList(hMass, powerExponent1,powerExponent2,A1,A2,constant));

//Second Order Polynomial in Exponential (ExpPol2)

RooRealVar b1("b1","b1",-1.0,-10.0,0.0);
RooRealVar b2("b2","b2",-1.0,-10.0,0.0);

RooGenericPdf ExpPol2("ExpPol2","TMath::Exp(x[1]*TMath::Power(x[0],1) + x[2]*TMath::Power(x[0],2))",RooArgList(hMass,b1,b2));

//Third Order Polynomial in Exponential (ExpPol3)

RooRealVar d1("d1","d1",-1.0,-5.0,0.0);
RooRealVar d2("d2","d2",-1.0,-5.0,0.0);
RooRealVar d3("d3","d3",-1.0,-5.0,0.0);

RooGenericPdf ExpPol3("ExpPol3","TMath::Exp(x[1]*TMath::Power(x[0],1) + x[2]*TMath::Power(x[0],2) + x[3]*TMath::Power(x[0],3))",RooArgList(hMass,d1,d2,d3));


//Fit the functions to our histograms



if (func == "ExpPol1")
{


    RooAddPdf FitFunction("FitFunction","FitFunction",RooArgList(dcbPdf,ExpPol1),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"),RooFit::Save(true));
    result->Print();

    N_signal[0] = N_SS.getValV();

    auto pl = hMass.frame();
    chi2[0] = pl->chiSquare(Form("pdf_%s",func.Data()),"data");

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl);
    FitFunction.plotOn(pl);
    pl->Draw();
    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol1.pdf");

    return 0;

}

else if(func == "Bernstein")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,BernsteinPol),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"),RooFit::Save(true));
    result->Print();

    N_signal[1] = N_SS.getValV();
    

    auto pl = hMass.frame();
    chi2[1] = pl->chiSquare();

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl);
    FitFunction.plotOn(pl);
    pl->Draw();
    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/BernsteinPol.pdf");

    return 0;
}
else if(func=="Power")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,PowerFunction),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    result->Print();

    N_signal[2] = N_SS.getValV();

    auto pl = hMass.frame();
    chi2[2] = pl->chiSquare();

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl);
    FitFunction.plotOn(pl);
    pl->Draw();
    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/PowerFunction.pdf");
    
    return 0;
}
else if(func=="ExpPol2")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,ExpPol2),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    result->Print();

    N_signal[3] = N_SS.getValV();

    auto pl = hMass.frame();
    chi2[3] = pl->chiSquare();

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl);
    FitFunction.plotOn(pl);
    pl->Draw();
    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol2.pdf");
    
    return 0;
}
else if(func=="ExpPol3")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,ExpPol3),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    result->Print();

    N_signal[4] = N_SS.getValV();

    auto pl = hMass.frame();
    chi2[4] = pl->chiSquare();

    double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); 
    hist.plotOn(pl);
    FitFunction.plotOn(pl);
    pl->Draw();
    canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol3.pdf");
    
    return 0;
}
cout<<"Function incorrectly specified";
return 0;


}

void SpuriousSignal() //Fit all at the same time
{

  TFile file("/home/mjlarroy/TFG_TUT/Root_Results/MC_reweighted_background_histos.root","READ");
  TFile file_signal("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_histos.root","READ");
  
  TH1D *hmyy = (TH1D*)file.Get("histmc_bkg_fit");
  TH1D *hmyy_signal = (TH1D*)file_signal.Get("h0yymass");

  double Bkg_error = 0;
  int nb = hmyy->GetNbinsX();

  for(int i=1;i<=nb;i++){
    Bkg_error += pow(hmyy->GetBinError(i),2);
  }

  Bkg_error = sqrt(Bkg_error);

  //N_expected from signal (120.125 and 129.75 hard-coded because we used 80 bins)

  double jmin = (120.125 - 105.0)/hmyy_signal->GetBinWidth(1); 
  double jmax = (129.75-105.0)/hmyy_signal->GetBinWidth(1);

  double N_exp = hmyy_signal->Integral(jmin,jmax);


  std::vector<TString> funcList = {"ExpPol1","Bernstein","Power","ExpPol2","ExpPol3"};

  double N_signal[5];
  double chi2 [5];

  for(int i=0;i<5;i++){
    BackgroundFit(hmyy, 105., 160., funcList[i],N_signal,chi2);
  }
  cout<<endl;
  cout<<"Ratio N_SS/N_expected and N_SS/Bkg_error for each function (in pairs)"<<endl;
  for(int i = 0; i<5;i++){
    cout<<funcList[i]<<endl;
    cout<<"NSS/N_expected ="<<abs(N_signal[i])/N_exp<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"NSS/Bkg_error ="<<abs(N_signal[i])/Bkg_error<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"Chi2 ="<<chi2[i]<<endl;
    cout<<endl;

  }


}