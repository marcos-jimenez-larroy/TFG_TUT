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

RooRealVar hMass("myy","m_{#gamma#gamma} [GeV]",xmin,xmax);
RooDataHist hist("hist","hist",hMass,Import(*h0yymass));

//DCB parameters
RooRealVar mu("mu","mu",par[1],xmin,xmax);
RooRealVar width("width","width",par[2],0.,2.);
RooRealVar a1("a1","a1",par[3],0.,100.);
RooRealVar p1("p1","p1",par[4],0.,100.);
RooRealVar a2("a2","a2",par[5],0.,100.);
RooRealVar p2("p2","p2",par[6],0.,100.);

RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);

//Fit utilizando la pdf "dcbPdf" a los datos guardados en "hist"
auto result = dcbPdf.fitTo(hist, RooFit::Save(true));
result->Print();

//Para acceder a los valores del fit (aquí un ejemplo para la media "mu" pero el resto son iguales)

//Esto es una forma de hacerlo un poco más rollo pero por si quieres guardarlas aparte. Si no, puedes hacerlo como hace más abajo (líneas 80-85) para imprimirlos en el canvas solo con el mu.getVal
RooArgList ParameterList = result->floatParsFinal();
RooRealVar* var_mu = (RooRealVar*) ParameterList.find("mu");

double fitted_mu = var_mu->getVal();
double fitted_mu_error = var_mu->getError();




RooAbsReal* chi2 = dcbPdf.createChi2(hist,Extended(true), DataError(RooAbsData::Poisson));

double Xndf = chi2->getVal()/(55.0-6.0); //number of degrees of freedom = N# bins - N# free parameters in the fit

auto pl = hMass.frame();

double w = 650;
double he = 450;
TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,w,he);
 canvasn->cd(); 

 pl->GetYaxis()->SetTitle("Events/GeV");

hist.plotOn(pl);
dcbPdf.plotOn(pl,RooFit::Name("DSCB Signal Fit"));

pl->Draw();

TLatex l; l.SetTextSize(0.05); l.SetTextFont(42);
l.DrawLatexNDC(0.52,0.86,Form("#chi^{2}/ndf = %.2f",Xndf));
l.DrawLatexNDC(0.52,0.8,Form("m_{H} = %.2f #pm %.2f GeV",mu.getValV(), mu.getError()));
l.DrawLatexNDC(0.52,0.73,Form("#sigma_{H} = %.2f #pm %.2f GeV",width.getValV(), width.getError()));
l.DrawLatexNDC(0.52,0.66,Form("#alpha_{low} = %.2f #pm %.2f GeV",a1.getValV(), a1.getError()));
l.DrawLatexNDC(0.52,0.60,Form("#alpha_{high} = %.2f #pm %.2f GeV",a2.getValV(), a2.getError()));
l.DrawLatexNDC(0.52,0.54,Form("n_{low} = %.2f #pm %.2f GeV",p1.getValV(), p1.getError()));
l.DrawLatexNDC(0.52,0.48,Form("n_{high} = %.2f #pm %.2f GeV",p2.getValV(), p2.getError()));

l.SetTextFont(72);
  l.SetTextColor(1);
  l.SetTextSize(0.05);
  l.DrawLatexNDC(0.20,0.85,"ATLAS");
  l.SetTextFont(42);
  l.SetTextSize(0.05);

  l.DrawLatexNDC(0.20,0.78,"#sqrt{s} = 13,6 TeV");
//double norm = dcbPdf.getNorm(RooArgSet(hMass));

//std::cout<<norm<<std::endl;
//l.DrawLatexNDC(0.55,0.66,Form("N=%.2f",norm));

//RooAbsReal* normIntegral = dcbPdf.createIntegral(RooArgSet(hMass));

// Then propagate the error
//double normValue = normIntegral->getVal();
//double normError = normIntegral->getPropagatedError(*result);


//std::cout<<normValue<<std::endl;
//std::cout<<normError<<std::endl;

 canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/check.pdf");

return 0;
}

void MakeFit()
{

  TFile file("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_histos.root","READ");
  
  TH1D *hmyy = (TH1D*)file.Get("h0yymass");

  //std::cout<<hmyy->Integral()<<std::endl;

  DCB_RF(hmyy, 105., 160., 10000.);


}
