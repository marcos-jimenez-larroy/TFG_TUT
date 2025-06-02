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


ostringstream histo;

string xt;
xt = "m_{#gamma#gamma} [GeV]"; 
// titulo eje y
string yt;
yt = "Events";


void DrawRatio(TString histName)
{



    TFile *data_RC;  

    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/data22_23_bins.root"; 
    data_RC = new TFile(histo.str().c_str(),"read"); 

    TFile *data_PC;

    histo.str(""); histo << "/home/mjlarroy/TFG_TUT/Root_Results/data22_23_bins_ProductCuts.root";
    data_PC = new TFile(histo.str().c_str(),"read");

    double w = 650;
    double h = 650;
    double xx = 0.; // posicion x de la esquina superior izquierda del canvas
    double yy = 0.; 


        TH1D *data_RC_hist; // histograma de todos los MC que voy a mostrar
        data_RC_hist = (TH1D*)data_RC->Get(histName); // leer histograma MC

        TH1D *data_PC_hist;
        data_PC_hist = (TH1D*)data_PC->Get(histName);

        TH1* ratio = new TH1D("ratio","ratio",55,105.0,160.0);

        for(int j=1;j<=55;j++){
            ratio->SetBinContent(j,data_PC_hist->GetBinContent(j)/data_RC_hist->GetBinContent(j));
            ratio->SetBinError(j,sqrt((data_PC_hist->GetBinError(j)/data_RC_hist->GetBinContent(j))*(data_PC_hist->GetBinError(j)/data_RC_hist->GetBinContent(j)) + (ratio->GetBinContent(j)*data_RC_hist->GetBinError(j)/data_RC_hist->GetBinContent(j))*(ratio->GetBinContent(j)*data_RC_hist->GetBinError(j)/data_RC_hist->GetBinContent(j))));
        }

        TCanvas *canvasp = new TCanvas("canvasp","canvasp",xx,yy,w,h);
        canvasp->SetWindowSize(w+(w-canvasp->GetWw()),h+(h-canvasp->GetWh()));
        canvasp->SetFillStyle(4000); // para hacer transparente
        canvasp->cd();

        TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.0,0.29,0.9,0.99);
        TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0.0,0.10,0.9,0.41);
        

        // pintar pads
        canvas_1->SetFillStyle(4000); //transparentes
        canvas_1->Draw();
        canvas_2->SetFillStyle(4000);
        canvas_2->Draw();
        

        canvas_1->cd();
        //canvas_1->SetLogy();
        //canvas_1->SetLogx();
        // crear frame
        double bl = 105; // principio del eje X
        double bu = data_RC_hist->GetBinCenter(55)+data_RC_hist->GetBinWidth(55)/2.;; // final del eje X

        //double yl = histdat->GetMinimum();  // minimo eje Y
        double yl = 0.0; // minimo eje Y. el 0.5 por logaithmic scale
        double yu = data_RC_hist->GetMaximum();  // maximo eje Y

        yu = yu*1.2; yl = yl*0.9;

        TH1F *h1 = canvas_1->DrawFrame(bl,yl,bu,yu); // crea el frame
        //

        //
        //Te renta hacer el Frame para poner todos los settings de una vez y no sobre cada histograma y estas cosas.

        // titulo eje y
        histo.str(""); histo << yt;
        h1->SetYTitle(histo.str().c_str());
        // settings eje y
        h1->GetYaxis()->SetTitleSize(0.05);
        h1->GetYaxis()->SetTitleOffset(1.4);
        h1->GetYaxis()->SetLabelSize(0.05);
        h1->GetYaxis()->SetNdivisions(606);

        h1->GetXaxis()->SetNdivisions(705);
        h1->GetXaxis()->SetMoreLogLabels(1);
        h1->GetXaxis()->SetNoExponent(1);
        //
        // settings eje x
        histo.str(""); histo << xt;
        h1->SetXTitle(histo.str().c_str());
        //
        
        
        h1->GetXaxis()->SetTitleSize(0.0);
        h1->GetXaxis()->SetTitleOffset(1.1);

        //
        h1->GetXaxis()->SetLabelSize(0.00);
        h1->GetXaxis()->SetLabelOffset(2000);

        //h1->GetYaxis()->SetRangeUser(0.0011,0.2);

        gPad->RedrawAxis();
         

        // pintar el frame
        h1->Draw();
        // titulo eje x

        myhhdate(data_PC_hist,"same][",0.0,kBlue,0,2.,1,kBlue,0,0);
        myhhdate(data_RC_hist,"same][",0.0,kRed,0,2.,5,kRed,0,0);

        myLine(0.68,0.75,0.05,kBlue,20,"Product Cut 0.32",0.04,0.04);
        myLine(0.68,0.7,0.05,kRed,5,"Relative Cut 0.35",0.04,0.04);

        TLatex l1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
        l1.SetNDC();
        l1.SetTextFont(72);
        l1.SetTextColor(1);
        l1.SetTextSize(0.04);
        l1.DrawLatex(0.50,0.86,"ATLAS");
        l1.SetTextFont(42);

        l1.DrawLatex(0.50,0.81,"#sqrt{s} = 13,6 TeV");
        l1.DrawLatex(0.70,0.81,"L = 58.98 fb^{-1}");



        canvas_2->cd();
    //canvas_2->SetLogx();

        canvas_2->SetBottomMargin(0.30);
        //
        //double mx = 0.29;

        TH1F *h2;

        if (histName == "h_myy_pT_0"){
            h2 = canvas_2->DrawFrame(bl,0.9,bu,1.15);
        }

        else {
            h2 = canvas_2->DrawFrame(bl,0.95,bu,1.05);
        }

        h2->SetYTitle("Ratio PC/RC");
        //
        histo.str(""); histo << xt;
        h2->SetXTitle(histo.str().c_str());
        //
        h2->GetXaxis()->SetTitleSize(0.11);
        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetXaxis()->SetLabelSize(0.1);
        //
        h2->GetXaxis()->SetNdivisions(705);
        h2->GetXaxis()->SetMoreLogLabels(1);
        h2->GetXaxis()->SetNoExponent(1);
        //
        h2->GetYaxis()->SetTitleSize(0.10);
        h2->GetYaxis()->SetTitleOffset(0.7);
        h2->GetYaxis()->SetLabelSize(0.1);
        h2->GetYaxis()->SetNdivisions(606);

        myhhdate(ratio,"esamex0",1.0,1,20,2.,1,1,1,1);

        TString FileName = "/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ratio_PC_RC_" + histName + ".pdf";
        canvasp->Print(FileName);
    
        return;
}

void PC_histograms()

{


    for(int i=0;i<16;i++){
        TString histName = "h_myy_pT_" + std::to_string(i);
        DrawRatio(histName);
    }

}