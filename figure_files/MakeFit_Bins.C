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

using namespace RooFit;


/*
Double_t DoubleSidedCB2_RF(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  double u   = (x-mu)/width;
  double A1  = TMath::Power(p1/TMath::Abs(a1),-p1)*TMath::Exp(-a1*a1/2.);
  double A2  = TMath::Power(p2/TMath::Abs(a2),-p2)*TMath::Exp(-a2*a2/2.);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2.);
  else            result *= A2*TMath::Power(B2-u,-p2);
  return result;
}


double DoubleSidedCB_RF(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2_RF(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}
*/

//Roofit

int DCB_RF(TH1 *h0yymass, double xmin, double xmax, double N, double *fit_parameters, double *fit_errors){

double par[6];

par[0]=125.0;
par[1]=0.8;
par[2]=1.;
par[3]=1.;
par[4]=1.;
par[5]=1.;

RooRealVar hMass("myy","myy [GeV]",xmin,xmax);
RooDataHist hist("hist","hist",hMass,Import(*h0yymass));

//DCB parameters
RooRealVar mu("mu","mu",par[0],xmin,xmax);
RooRealVar width("width","width",par[1],0.,2.);
RooRealVar a1("a1","a1",par[2],0.,100.);
RooRealVar p1("p1","p1",par[3],0.,100.);
RooRealVar a2("a2","a2",par[4],0.,100.);
RooRealVar p2("p2","p2",par[5],0.,100.);

RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);


auto result = dcbPdf.fitTo(hist, RooFit::Save(true));
result->Print();

fit_parameters[0] = mu.getValV();
fit_errors[0] = mu.getError();
fit_parameters[1] = width.getValV();
fit_errors[1] = width.getError();
fit_parameters[2] = a1.getValV();
fit_errors[2] = a1.getError();
fit_parameters[3] = p1.getValV();
fit_errors[3] = p1.getError();
fit_parameters[4] = a2.getValV();
fit_errors[4] = a2.getError();
fit_parameters[5] = p2.getValV();
fit_errors[5] = p2.getError();


return 0;
}

void MakeFit_Bins(TString histname, double *fit_parameters, double *fit_errors)
{

  TFile file("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_bins.root","READ");
  
  TH1D *hmyy_signal = (TH1D*)file.Get(histname);

  DCB_RF(hmyy_signal, 105., 160., 10000.,fit_parameters, fit_errors);


}