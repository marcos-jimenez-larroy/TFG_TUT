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

void Inclusive_Representation()
{
    ostringstream histo;
  //
  // titulo del eje x
  string xt;
  xt = "m_{#gamma#gamma} [GeV]"; 
  // titulo eje y
  string yt;
  yt = "Events/GeV";

  TFile *mc;  // root de todos los MC del que voy a leer
  TF1 *ajuste_bkg; // histograma de todos los MC que voy a mostrar
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/MC_reweighted_background_histos_22_23.root"; // nombre del root de los MC
  mc = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "ajuste"; // nombre del histograma que quiero leer
  ajuste_bkg = (TF1*)mc->Get(histo.str().c_str()); // leer histograma MC

  TFile *data;  // root de todos los MC del que voy a leer
  TH1D *hist_data; // histograma de todos los MC que voy a mostrar
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/data22_23_histos.root"; // nombre del root de los MC
  data = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "h0yymass"; // nombre del histograma que quiero leer
  hist_data = (TH1D*)data->Get(histo.str().c_str()); // leer histograma MC

  int ni = 1; int nb = hist_data->GetNbinsX(); // numero de bines en el eje X









  //Resultados del DSCB de la región inclusiva
  /*Double_t a1 = 1.7735;
  Double_t a2 = 1.7375;
  Double_t mu = 125.50;
  Double_t width = 1.816;
  Double_t p1 = 4.4322;
  Double_t p2 = 10.352;*/

  double L1 = 27.58;
  double Lt = 31.40+27.58;
  //Double_t N = 1475.54*Lt/L1;

  TH1D* hist_dif;
  hist_dif = (TH1D*)hist_data->Clone(histo.str().c_str()); 
  hist_dif->Add(ajuste_bkg,-1);

  TH1D* hist_ajuste;
  hist_ajuste = (TH1D*)hist_data->Clone(histo.str().c_str()); 
  hist_ajuste->Reset();
  double value = 0.0;

  TH1D* hist_sum;
  hist_sum = (TH1D*)hist_data->Clone(histo.str().c_str()); 
  hist_sum->Reset();

  TH1D* hist_signal;
  hist_signal = (TH1D*)hist_data->Clone(histo.str().c_str()); 
  hist_signal->Reset();




  RooRealVar hMass("myy","myy [GeV]",105.,160.);
  RooDataHist hist("hist","hist",hMass,Import(*hist_data));

  RooRealVar mu("mu","mu",125.50,125.5,125.5);
  RooRealVar width("width","width",1.816,1.816,1.816);
  RooRealVar a1("a1","a1",1.7735,1.7735,1.7735);
  RooRealVar p1("p1","p1",4.4322,4.4322,4.4322);
  RooRealVar a2("a2","a2",1.7375,1.7375,1.7375);
  RooRealVar p2("p2","p2",10.352,10.352,10.352);

  RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);

  RooRealVar b1("b1","b1",300,-1000.0,400.0);
  RooRealVar b2("b2","b2",300,-1000.0,400.0);
  RooRealVar b3("b3","b3",300,-1000.0,400.0);
  RooRealVar c("c","c",20000,0.0,100000.0);


  RooGenericPdf ThirdPol("ThirdPol","x[1] + (x[0]-100.0)*x[2] + (x[0]-100.)*(x[0]-100.)*x[3] + (x[0]-100.)*(x[0]-100.)*(x[0]-100.)*x[4]",RooArgList(hMass,c,b1,b2,b3));


  RooRealVar N_SS("N_SS","N_SS",1000,-10000,10000);
  RooRealVar N_bkg("N_bkg","N_bkg",100,-1000,100000000);

  RooAddPdf FitFunction("FitFunction","FitFunction",RooArgList(dcbPdf,ThirdPol),RooArgList(N_SS,N_bkg));

  


  RooFitResult * result = FitFunction.fitTo(hist, RooFit::Save(true));

  std::cout<<"Hi"<<std::endl;

  RooArgList SignalParameterList1 = result->floatParsFinal();
  RooRealVar* N_Sig = (RooRealVar*) SignalParameterList1.find("N_SS");



  Double_t N = N_Sig->getVal();
  std::cout<<N<<std::endl;
  Double_t N_err = N_Sig->getError();
  std::cout<<N_err<<std::endl;



  RooRandom::randomGenerator()->SetSeed(123);  // Set a random seed for reproducibility

  RooDataSet* signal_data = dcbPdf.generate(RooArgSet(hMass), 10000);
    
  for (int i = 1; i <= 55; ++i) {
    double xValue = hist_signal->GetXaxis()->GetBinCenter(i); // Get bin center
    hMass.setVal(xValue); // Set the value of x
    double pdfValue = dcbPdf.getVal(RooArgSet(hMass)); // Evaluate the PDF

    // Fill the histogram with the PDF value
    hist_signal->Fill(xValue, N*pdfValue);

}
double ndat = 0.0;

for(int i=1;i<=nb;i++){
  value = ajuste_bkg->Eval(105+i);
  hist_ajuste->SetBinContent(i,value);

  double dd = hist_data->GetBinContent(i);
  ndat += dd;
  hist_data->SetBinError(i,sqrt(dd));
  hist_dif->SetBinError(i,sqrt(dd));

  hist_signal->SetBinError(i,0.0);
}

    hist_sum->Add(hist_ajuste,1.);
    hist_sum->Add(hist_sum,hist_signal);

for(int i=1;i<=nb;i++){
  hist_sum->SetBinError(i,0.0);

}

  // PARTE ESTÉTICA DE LA REPRESENTACIÓN


  // anchura y altura del canvas
  double w = 650;
  double h = 650;
  double xx = 0.; // posicion x de la esquina superior izquierda del canvas
  double yy = 0.; // posicion y de la esquina superior izquierda del canvas
  // crea canvas
  TCanvas *canvasp = new TCanvas("canvasp","canvasp",xx,yy,w,h);
  canvasp->SetWindowSize(w+(w-canvasp->GetWw()),h+(h-canvasp->GetWh()));
  canvasp->SetFillStyle(4000); // para hacer transparente
  canvasp->cd();
  //
  // crear pads (como unos subplots de toda la vida)
  TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.0,0.29,0.9,0.99);
  TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0.0,0.10,0.9,0.41);
  //
  // pintar pads
  canvas_1->SetFillStyle(4000); //transparentes
  canvas_1->Draw();
  canvas_2->SetFillStyle(4000);
  canvas_2->Draw();
  //
  canvas_1->cd();
  //canvas_1->SetLogy(); // escala logaritmica en eje Y
  //
  // crear frame
  double bl = hist_data->GetBinCenter(ni)-hist_data->GetBinWidth(ni)/2.; // principio del eje X
  double bu = hist_data->GetBinCenter(nb)+hist_data->GetBinWidth(nb)/2.; // final del eje X
  //double yl = histdat->GetMinimum();  // minimo eje Y
  double yl = 0.0001; // minimo eje Y. el 0.5 por logaithmic scale
  double yu = hist_data->GetMaximum();  // maximo eje Y
  yu = yu*1.4; // maximo eje Y
  TH1F *h1 = canvas_1->DrawFrame(bl,yl,bu,yu); // crea el frame
  //
  TLine ll(bl,0.,bu,0.);
  ll.SetLineWidth(2.);
  ll.SetLineStyle(1); //Crea la línea que vemos luego debajo
  ll.SetLineColor(0);
  //

  //Te renta hacer el Frame para poner todos los settings de una vez y no sobre cada histograma y estas cosas.


  // titulo eje y
  histo.str(""); histo << yt;
  h1->SetYTitle(histo.str().c_str());
  // settings eje y
  h1->GetYaxis()->SetTitleOffset(1.6);
  h1->GetYaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetNdivisions(510);
  //
  // settings eje x
  h1->GetXaxis()->SetLabelOffset(2000); //Esto lo hace para no poner título ni uds en el eje de la foto de arriba
  h1->GetXaxis()->SetTitleSize(0.0);
  h1->GetXaxis()->SetLabelSize(0.0);
  h1->GetXaxis()->SetNdivisions(606);
  h1->GetXaxis()->SetTickLength(0.04);
  // titulo eje x
  histo.str(""); histo << xt;
  h1->SetXTitle(histo.str().c_str());

  h1->SetTitle("Invariant mass of the two most energetic photons"); //Ask tmb

  //
  // pintar el frame
  h1->Draw();
  //
  // pintar los histogramas (función en el archivo myAtlasUtils.C)
  //Pintas del más grande al más pequeño

  myhhdate(hist_data,"esamex0",1.0,1,20,2.,1,1,1,1); //Falta el detalle de los ejes! (Preguntar a Ana)
  //myhhdate(histmc,"same][",0.0,1,0,2.,1,1,1000,0);
  //
  //
  myhhdate(hist_signal,"same][L",0.0,6,0,2.,1,6,1000,0);


  
  
  myhhdate(hist_ajuste,"same][L",0.0,4,0,2.,7,4,1000,0);
  myhhdate(hist_sum,"same][L",0.0,2,0,2.,1,2,1000,0);

  



  //
  //myhhdate(hist_dif,"esamex0",1.0,1,20,2.,1,1,1,1);

  // labels
  double x1 = 0.31;double x2 = 0.50;
  double yy1 = 0.89; 
  double y2 = yy1 - 0.04;
  
  double y3 = y2 - 0.04;
  
  double y4 = y3 - 0.04;

  histo.str(""); histo<<"Data (entries = ";histo << int(ndat);histo<<")";
  myMarkerText(x1,yy1,1,20,histo.str().c_str(),1.0,0.04,0.04);
  myLine(x1,y2,0.05,2,20,"Signal + Background",0.04,0.04);
  myLine(x1,y3,0.05,4,7,"Background",0.04,0.04);
  myLine(x1,y4,0.05,6,1000,"Signal",0.04,0.04);
  /*
  double y5 = y4 - 0.04;
  double y6 = y5 - 0.04;
  double y7 = y6 - 0.04;*/
  

  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(1);
  l.DrawLatex(x1+0.37,y2,"ATLAS");
  l.SetTextFont(42);

  l.DrawLatex(x1+0.37,y4,"#sqrt{s} = 13,6 TeV");
  l.DrawLatex(x1+0.37,y4-0.08,"L = 58.98 fb^{-1}");

  //
  // re-dibujar los ejes
  gPad->RedrawAxis();
  //
  canvas_2->cd();
  canvas_2->SetBottomMargin(0.30);
  //
  //double mx = 0.29;
  double mx = 0.09;
  TH1F *h2 = canvas_2->DrawFrame(bl,-400,bu,800);
  h2->SetYTitle("Events - fitted bkg");
  //
  histo.str(""); histo << xt;
  h2->SetXTitle(histo.str().c_str());
  //
  h2->GetXaxis()->SetTitleSize(0.11);
  h2->GetXaxis()->SetTitleOffset(1.1);
  //
  h2->GetXaxis()->SetLabelSize(0.09);
  h2->GetXaxis()->SetNdivisions(505);
  h2->GetXaxis()->SetMoreLogLabels(1);
  h2->GetXaxis()->SetNoExponent(1);
  //
  h2->GetYaxis()->SetTitleSize(0.10);
  h2->GetYaxis()->SetTitleOffset(0.7);
  h2->GetYaxis()->SetLabelSize(0.09);
  //h2->GetYaxis()->SetMoreLogLabels(1);


  //
  TAxis *axis_x2 = h2->GetXaxis();
  axis_x2->SetNdivisions(606);
  axis_x2->SetMoreLogLabels(1);
  axis_x2->SetNoExponent(1);
  axis_x2->SetTickLength(0.08);
  //
  TAxis *axis_y2 = h2->GetYaxis();
  axis_y2->SetNdivisions(507);
  //
  h2->Draw();
  //
  ll.Draw();
  //ll.DrawLine(bl,200.0,bu,200.0);
  //ll.DrawLine(bl,-200.0,bu,-200.0);
  //ll.DrawLine(bl,1.0-0.05,bu,1.0-0.05);





  myhhdate(hist_dif,"esamex0",1.0,1,20,2.,1,1,1,1);
  myhhdate(hist_signal,"same][L",0.0,6,0,2.,1,6,1000,0);
  //
  gPad->RedrawAxis();
  //
  canvasp->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/Prueba_InclusiveRegion_22_23.pdf");

}