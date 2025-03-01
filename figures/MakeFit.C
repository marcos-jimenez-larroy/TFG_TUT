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

int DCB_RF(TH1 *h0yymass, double xmin, double xmax, double N){

double par[7];
par[0]=N; //Right now this is not doing anything
par[1]=125.0;
par[2]=0.8;
par[3]=1.;
par[4]=1.;
par[5]=1.;
par[6]=1.;

RooRealVar hMass("myy","myy [GeV]",xmin,xmax);
RooDataHist hist("hist","hist",hMass,Import(*h0yymass));

//DCB parameters
RooRealVar mu("mu","mu",par[1],xmin,xmax);
RooRealVar width("width","width",par[2],0.,2.);
RooRealVar a1("a1","a1",par[3],0.,100.);
RooRealVar p1("p1","p1",par[4],0.,100.);
RooRealVar a2("a2","a2",par[5],0.,100.);
RooRealVar p2("p2","p2",par[6],0.,100.);

RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);


auto result = dcbPdf.fitTo(hist, RooFit::Save(true));
result->Print();

auto pl = hMass.frame();

double w = 650;
double he = 450;
TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
 canvasn->cd(); 
hist.plotOn(pl);
dcbPdf.plotOn(pl);
pl->Draw();
 canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/check.pdf");

return 0;
}

void MakeFit()
{

  TFile file("mc23a_signal_histos.root","READ");
  
  TH1D *hmyy = (TH1D*)file.Get("h0yymass");

  DCB_RF(hmyy, 105., 160., 10000.);


}
