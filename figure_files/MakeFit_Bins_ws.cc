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


RooWorkspace * DCB_RF_ws(
  TH1 *h0yymass, 
  double xmin, 
  double xmax, 
  double N, 
  double *fit_parameters, 
  double *fit_errors, 
  bool doPlots,
  TString number){

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

  RooArgList SignalParameterList = result->floatParsFinal();


    RooRealVar* var_mu = (RooRealVar*) SignalParameterList.find("mu");
    RooRealVar* var_width = (RooRealVar*) SignalParameterList.find("width");
    RooRealVar* var_a1 = (RooRealVar*) SignalParameterList.find("a1");
    RooRealVar* var_p1 = (RooRealVar*) SignalParameterList.find("p1");
    RooRealVar* var_a2 = (RooRealVar*) SignalParameterList.find("a2");
    RooRealVar* var_p2 = (RooRealVar*) SignalParameterList.find("p2");





  auto pl = hMass.frame();
   
    double wi = 650;
    double he = 450;
    TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,wi,he);
    canvasn->cd(); 
    hist.plotOn(pl,RooFit::Name("hist"));


    w->pdf("dcbPdf")->plotOn(pl,RooFit::Name("dcbPdf"),RooFit::LineColor(kBlue));

    RooAbsReal* chi2 = w->pdf("dcbPdf")->createChi2(hist,Extended(true), DataError(RooAbsData::Poisson));

    double Xndf = chi2->getVal()/(55.0-6.0);

    pl->GetYaxis()->SetTitle("Events/GeV");
    pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");


    pl->Draw();

  TLatex l; l.SetTextSize(0.05); l.SetTextFont(42);
  l.DrawLatexNDC(0.52,0.86,Form("#chi^{2}/ndf = %.2f",Xndf));
  l.DrawLatexNDC(0.52,0.8,Form("m_{H} = %.2f #pm %.2f GeV",var_mu->getValV(), var_mu->getError()));
  l.DrawLatexNDC(0.52,0.73,Form("#sigma_{H} = %.2f #pm %.2f GeV",var_width->getValV(), var_width->getError()));
  l.DrawLatexNDC(0.52,0.66,Form("#alpha_{low} = %.2f #pm %.2f GeV",var_a1->getValV(), var_a1->getError()));
  l.DrawLatexNDC(0.52,0.60,Form("#alpha_{high} = %.2f #pm %.2f GeV",var_a2->getValV(), var_a2->getError()));
  l.DrawLatexNDC(0.52,0.54,Form("n_{low} = %.2f #pm %.2f GeV",var_p1->getValV(), var_p1->getError()));
  l.DrawLatexNDC(0.52,0.48,Form("n_{high} = %.2f #pm %.2f GeV",var_p2->getValV(), var_p2->getError()));

  l.SetTextFont(72);
    l.SetTextColor(1);
    l.SetTextSize(0.05);
    l.DrawLatexNDC(0.20,0.85,"ATLAS");
    l.SetTextFont(42);
    l.SetTextSize(0.05);

  l.DrawLatexNDC(0.20,0.78,"#sqrt{s} = 13,6 TeV");

  canvasn->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/Signal_fit_" + number + ".pdf");







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

  TString Number = histname;

  //in the future, you can use the outTag to identify prod mechanism, eta region, conv/unconv...
  //with this logic, you can run over all histograms and save the workspaces in the same file
  TString outTag = histname ; 
  TFile * outfile = new TFile( thePath + "outFile_Signal/outFileSignal_"+ outTag + ".root" , "update" ) ;

  RooWorkspace * w = DCB_RF_ws(hmyy_signal, 105., 160., 10000.,fit_parameters, fit_errors, 1, Number);

  //writes workspace to file, overwriting previous one with identical name
  w->Write("w_" + histname , TObject::kOverwrite ) ;
  outfile->Close() ; 
}
