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


RooWorkspace * DCB_RF_ws(
  TH1 *h0yymass, 
  double xmin, 
  double xmax, 
  double N, 
  double *fit_parameters, 
  double *fit_errors, 
  bool doPlots ){

  double par[6];

  par[0]=125.0;
  par[1]=0.8;
  par[2]=1.;
  par[3]=1.;
  par[4]=1.;
  par[5]=1.;

  RooWorkspace * w = new RooWorkspace("w") ; 
  RooRealVar hMass("myy","myy [GeV]",xmin,xmax);
  RooDataHist hist("hist","hist",hMass,Import(*h0yymass));

  w->import(hist)  ;

  //DCB parameters
  RooRealVar mu("mu","mu",par[0],xmin,xmax);
  RooRealVar width("width","width",par[1],0.,10.);
  RooRealVar a1("a1","a1",par[2],0.,100.);
  RooRealVar p1("p1","p1",par[3],0.,100.);
  RooRealVar a2("a2","a2",par[4],0.,100.);
  RooRealVar p2("p2","p2",par[5],0.,100.);

  RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);
  w->import(dcbPdf) ; 
  //object that includes all parameters of the pdf (and their errors)
  RooArgSet * fitParams = w->pdf("dcbPdf")->getParameters(*w->data("hist")) ;
  //snapshot of the parameters before the fit
  w->saveSnapshot("prefit" , *fitParams ) ; 
  //replaced auto with RooFitResult
  RooFitResult * result = w->pdf("dcbPdf")->fitTo(*w->data("hist"), RooFit::Save(true));
  result->Print();
  w->import(*result,"fit_result") ;
  //snapshot of the parameters after the fit
  w->saveSnapshot("postfit" , *fitParams ) ; 


  return w;
}

        //Ahora mismo está puesto para que use el binado 'new', que es un poco más pequeño que el original!


void MakeFit_Bins_ws(
    TString histname = "h_myy_pT_0", 
    double * fit_parameters = nullptr , 
    double * fit_errors = nullptr )
{
  TString thePath = "/home/mjlarroy/TFG_TUT/Root_Results/" ;
  TFile file(thePath + "mc23a_signal_bins.root","READ");

  TH1D *hmyy_signal = (TH1D*)file.Get(histname);

  //in the future, you can use the outTag to identify prod mechanism, eta region, conv/unconv...
  //with this logic, you can run over all histograms and save the workspaces in the same file
  TString outTag = histname ; 
  TFile * outfile = new TFile( thePath + "outFileSignal_"+ outTag + ".root" , "update" ) ;

  RooWorkspace * w = DCB_RF_ws(hmyy_signal, 105., 160., 10000.,fit_parameters, fit_errors, 1);

  //writes workspace to file, overwriting previous one with identical name
  w->Write("w_" + histname , TObject::kOverwrite ) ;
  outfile->Close() ; 
}
