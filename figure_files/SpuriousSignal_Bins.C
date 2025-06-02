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



int BackgroundFit(TH1* h0yymass,double xmin, double xmax, TString func, double *N_signal, double *N_err,double *chi2, double *signal_fit_parameters)
{
    RooRealVar hMass("myy","myy [GeV]",xmin,xmax);
    RooDataHist hist("hist","hist",hMass,Import(*h0yymass));


//Function definitions:

//Double-Sided-Crystal-Ball
double par[6]; //These parameters result from the previous fit of the signal-only MC
par[0] = signal_fit_parameters[0];
par[1] = signal_fit_parameters[1];
par[2] = signal_fit_parameters[2];
par[3] = signal_fit_parameters[3];
par[4] = signal_fit_parameters[4];
par[5] = signal_fit_parameters[5];




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
    //result->Print();

    N_signal[0] = N_SS.getValV();
    N_err[0] = N_SS.getError();

    auto pl = hMass.frame(RooFit::Name("ExpPol1"));

    /*double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); */

    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[0] = pl->chiSquare("FitFunction","hist");

    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol1.pdf");

    return 0;

}

else if(func == "Bernstein")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,BernsteinPol),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"),RooFit::Save(true));
    //result->Print();

    N_signal[1] = N_SS.getValV();
    N_err[1] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("Bernstein Polynomial"));
   
    /*double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); */
    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[1] = pl->chiSquare("FitFunction","hist");
    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/BernsteinPol.pdf");

    return 0;
}
else if(func=="Power")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,PowerFunction),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    //result->Print();

    N_signal[2] = N_SS.getValV();
    N_err[2] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("Power Function"));

    /*double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); */
    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[2] = pl->chiSquare("FitFunction","hist");
    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/PowerFunction.pdf");
    
    return 0;
}
else if(func=="ExpPol2")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,ExpPol2),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    //result->Print();

    N_signal[3] = N_SS.getValV();
    N_err[3] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("ExpPol2"));

    /*double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); */
    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[3] = pl->chiSquare("FitFunction","hist");
    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol2.pdf");
    
    return 0;
}
else if(func=="ExpPol3")
{
    RooAddPdf FitFunction("Total","Total",RooArgList(dcbPdf,ExpPol3),RooArgList(N_SS,N_bkg));

    auto result = FitFunction.fitTo(hist, RooFit::Minimizer("Minuit2"), RooFit::Save(true));
    //result->Print();

    N_signal[4] = N_SS.getValV();
    N_err[4] = N_SS.getError();


    auto pl = hMass.frame(RooFit::Name("ExpPol3"));

    /*double w = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
    canvasn->cd(); */

    hist.plotOn(pl,RooFit::Name("hist"));


    FitFunction.plotOn(pl,RooFit::Name("FitFunction"));
    
    chi2[4] = pl->chiSquare("FitFunction","hist");

    //pl->Draw();
    //canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ExpPol3.pdf");
    
    return 0;
}
cout<<"Function incorrectly specified";
return 0;


}

TString SpuriousSignal(TString histname, double *fit_parameters) //Fit all at the same time
{
  TFile file("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_background_bins.root","READ");
  TFile file_signal("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_bins.root","READ");
  
  TH1D *hmyy_signal = (TH1D*)file_signal.Get(histname);
  TH1D *hmyy = (TH1D*)file.Get(histname);


  //N_expected from signal (120.125 and 129.75 hard-coded because we used 80 bins)

  double jmin = (120.125 - 105.0)/hmyy_signal->GetBinWidth(1); 
  double jmax = (129.75-105.0)/hmyy_signal->GetBinWidth(1);

  double N_expected = hmyy_signal->Integral(jmin,jmax); //Must coincide with the N_bkg in every fit


  std::vector<TString> funcList = {"ExpPol1","Bernstein","Power","ExpPol2","ExpPol3"};

  double N_signal[5];
  double chi2 [5];
  double N_err[5];

  for(int i=0;i<5;i++){
    BackgroundFit(hmyy, 105., 160., funcList[i],N_signal,N_err,chi2, fit_parameters);
  }
  //Imprimir las cosas en pantalla para ver parámetros de cada uno pero comentado cuando aparecen todfos los histogramas.
  /*
  cout<<endl;
  cout<<"Ratio N_SS/N_expected and N_SS/N_SS_error for each function (in pairs)"<<endl;
  for(int i = 0; i<5;i++){
    cout<<funcList[i]<<endl;
    cout<<"NSS/N_expected ="<<abs(N_signal[i])/N_exp<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"NSS/NSS_error ="<<abs(N_signal[i])/N_err[i]<<endl; //Not very nice but it is only a 5 elemtn loop
    cout<<"Chi2 ="<<chi2[i]<<endl;
    cout<<endl;
  }*/

  //First of all, if the exponential (0) or the power fucntion (2) are enoguh we will use them

  std::string best_func = "none";
  if ((std::abs(N_signal[0])/N_expected < 0.1 || std::abs(N_signal[0])/N_err[0] < 0.2) && chi2[0]>0.01){
    if((std::abs(N_signal[2])/N_expected < 0.1 || std::abs(N_signal[2])/N_err[2] < 0.2) && chi2[2]>0.01){
        if (std::abs(N_signal[0])<std::abs(N_signal[2])) {
            best_func = funcList[0];
        }
        else{
            best_func = funcList[2];
        }
    }
    else{
      best_func = funcList[0];
    }
  }
  else if ((std::abs(N_signal[2])/N_expected < 0.1 || std::abs(N_signal[2])/N_err[2] < 0.2) && chi2[2]>0.01){
    best_func = funcList [2];
  }
  else if ((std::abs(N_signal[3])/N_expected < 0.1 || std::abs(N_signal[3])/N_err[3] < 0.2) && chi2[3]>0.01){
    best_func = funcList[3];
  }
  else if ((std::abs(N_signal[4])/N_expected < 0.1 || std::abs(N_signal[4])/N_err[4] < 0.2) && chi2[4]>0.01){
    best_func = funcList[4];
  }
  else if ((std::abs(N_signal[1])/N_expected < 0.1 || std::abs(N_signal[1])/N_err[1] < 0.2) && chi2[1]>0.01){
    best_func = funcList[1];
  }

  return best_func;





}


void SpuriousSignal_Bins()
{


  double signal_fit_parameters[6]; double signal_fit_errors[6];


    std::string BestFunctions[29];

    for (int i=12;i<13;i++){
        std::string histName= "h_myy_pT_" + std::to_string(i);

        MakeFit_Bins(histName, signal_fit_parameters, signal_fit_errors);

        BestFunctions[i] = SpuriousSignal(histName, signal_fit_parameters);

    }
    for(int i=12;i<13;i++){
        cout<<"h_myy_pT_"<<i<<": "<<BestFunctions[i]<<endl;

    }
}