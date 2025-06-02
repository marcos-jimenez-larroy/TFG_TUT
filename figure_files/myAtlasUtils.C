#include <iostream>
#include <cmath>

#include "AtlasUtils.h"

#include "TLatex.h"
#include "TMarker.h"
#include "TLine.h"
#include "TPave.h"
#include "TH1.h"

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
#include <map>

using namespace RooFit;
//
TColor *color1 = new TColor(131, 0.90, 0.00, 0.20); // reddd
TColor *color2 = new TColor(132, 0.14, 0.42, 0.56); // steelblue

void myhhdate(TH1 *hh,const char *opt,float ms,int color,int st,float wd,int ls,int lc,int fills,int color2)
{
  hh->SetMarkerSize(ms);
  hh->SetMarkerColor(color);
  hh->SetMarkerStyle(st);
  hh->SetLineWidth(wd);
  hh->SetLineStyle(ls);
  hh->SetLineColor(lc);
  hh->SetFillStyle(fills);
  hh->SetFillColor(color2);
  hh->Draw(opt);
}


void myepsfile(TCanvas *canvas,const char *file,const char *epp)
{
  canvas->Print(file,epp);
}

void myText(Double_t x,Double_t y,Color_t color,const char *text,Float_t mytsize) {

  //Double_t tsize=0.05;
  TLatex l; //l.SetTextAlign(12);
  l.SetTextSize(mytsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}


void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,Int_t mcolor_line,Int_t mstyle,const char *text,Float_t mytsize,Double_t tsize) 
{

  //Double_t tsize=0.06;

  TLatex l; l.SetTextAlign(12); 
  l.SetTextSize(mytsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  //printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

  //TPave *mbox= new TPave(x1,y1,x2,y2,0,"NDC");
  TPave *mbox= new TPave(x1,y1,x2,y2,1,"NDC");

  mbox->SetFillColor(mcolor);
  mbox->SetLineColor(mcolor_line);
  //mbox->SetFillStyle(1001);
  mbox->SetFillStyle(mstyle);
  mbox->SetLineWidth(2.);
  mbox->Draw();

  TLine mline;
  mline.SetLineWidth(4);
  mline.SetLineColor(1);
  mline.SetLineStyle(1);
  Double_t y_new=(y1+y2)/2.;
  //mline.DrawLineNDC(x1,y_new,x2,y_new);

}

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,const char *text,Float_t msize,Float_t mytsize,Float_t tsize) 
{
  //Double_t tsize=0.06;
  TMarker *marker = new TMarker(x-(0.4*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);
  marker->Draw();

  TLatex l; l.SetTextAlign(12); 
  l.SetTextSize(mytsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);
}



void myLine(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,Int_t mstyle,const char *text,Float_t mytsize,Double_t tsize) 
{

  TLatex l; l.SetTextAlign(12); 
  l.SetTextSize(mytsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  //printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

  TLine mline;
  mline.SetLineWidth(4);
  mline.SetLineColor(mcolor);
  mline.SetLineStyle(mstyle);
  Double_t y_new=(y1+y2)/2.;
  mline.DrawLineNDC(x1,y_new,x2,y_new);

}

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par){
if (reject && x[0] > 122 && x[0]<128){
    TF1::RejectPoint();
    return 0;
}
return par[0] + (x[0]-100.0)*par[1] + (x[0]-100.)*(x[0]-100.)*par[2] + (x[0]-100.)*(x[0]-100.)*(x[0]-100.)*par[3];
}


TF1 *fitExclude(TH1 *h0yymass,double xmin, double xmax) {

TF1 *ajuste = new TF1("ajuste",fline,xmin,xmax,4); //4 parameters
ajuste->SetParameters(100000.0,-300,300,300);

reject = true;
h0yymass->Fit(ajuste,"0");
reject = false;

TF1 *f1 = new TF1("f1",fline,xmin,xmax,4);
f1->SetParameters(ajuste->GetParameters());
f1->SetBit(TF1::kNotDraw);
h0yymass->GetListOfFunctions()->Add(f1);
gROOT->GetListOfFunctions()->Remove(f1);

return ajuste;

}
/*
RooFitResult *fitExclude_RooFit(TH1 *h0yymass, double xmin, double xmax, RooGenericPdf *pdf){



  RooFitResult result = pdf->fitTo(h0yymass, Minimizer("Minuit2"), Save(true), Range(105.,120.125), Range(129.75,160.));
  return result;
}*/

