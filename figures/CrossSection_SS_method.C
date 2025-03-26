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
using namespace RooFit;


    // clase que convierte cualquier tipo de variable a variable string
    ostringstream histo;
    //
    // titulo del eje x
    string xt;
    xt = "pT_{#gamma#gamma} [GeV]"; 
    // titulo eje y
    string yt;
    yt = "#frac{d#sigma}{d#Omega}";

    std::vector<TString> funcList = {
        "ExpPol3",
        "ExpPol3",
        "Power",
        "Power",
        "ExpPol3",
        "ExpPol3",
        "ExpPol3",
        "ExpPol3",
        "ExpPol3",
        "ExpPol3",
        "ExpPol3", 
        "ExpPol3",
        "Power", 
        "ExpPol3",
        "ExpPol3",
        "ExpPol1",
        "ExpPol1",
        "ExpPol3",
        "ExpPol1",
        "ExpPol1",
        "Bernstein",
        "ExpPol2",
        "ExpPol2",
        "ExpPol2",
        "Power",
        "ExpPol1",
        "ExpPol1",
        "Power",
        "ExpPol1"};  
    //10,17,20,21 -> With the relaxed constraint becomes ExpPol1
    //Needed to relax  in bin 12 the constraintof the NSS/NExpected condition form <10% to <14%.

    std::vector<TString> new_funcList = {
        "Power",
        "Power",
        "Power",
        "Power",
        "Power",
        "Power",
        "ExpPol3",
        "ExpPol3",
        "ExpPol1",
        "Power",
        "ExpPol1", 
        "ExpPol3",
        "Power", 
        "ExpPol3",
        "ExpPol3",
        "ExpPol1"};  
    

        //Ahora mismo está puesto para que use el binado normal. Las funciones 'new' se refieren al binado nuevo, más pequeño.

    void CrossSection_SS_method(){

        Double_t xbins[30] = {0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0,36.0,39.0,42.0,45.0,48.0,51.0,54.0,57.0,60.0,69.0,78.0,87.0,96.0,132.0,186.0,258.0,618.0,1500.0};
        Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};
        
        TH1* h0yypT = new TH1D("h0yypT","PTyy",16,new_xbins);

        TString thePath = "/home/mjlarroy/TFG_TUT/Root_Results/outFile_Background/";
        TString thePathData = "/home/mjlarroy/TFG_TUT/Root_Results/data22_bins.root";
        TString thePathSignal = "/home/mjlarroy/TFG_TUT/Root_Results/";

        TFile *data;
        data = new TFile(thePathData,"read");

        for(int i=0;i<29;i++){


            //Traer el Workspace correspondiente del ajuste al bcg
            TString FileName= "outFileBackground_" + funcList[i] + "_h_myy_pT_" + std::to_string(i) + ".root";

            TFile *file_ws;
            file_ws = new TFile(thePath + FileName,"READ");

            TString WorkspaceName = "w_out_" + funcList[i] + "_h_myy_pT_" + std::to_string(i);

            RooWorkspace *w = (RooWorkspace*)file_ws->Get(WorkspaceName);


            //Datos
            TH1D *hist_data;
            TString HistName = "h_myy_pT_" + std::to_string(i);
            hist_data = (TH1D*)data->Get(HistName); // leer histograma de datos

            RooRealVar hMass("myy","myy [GeV]",105.0,160.0);


            RooDataHist hist("hist","hist",hMass,Import(*hist_data));


            hMass.setRange("Left",105.,120.125);
            hMass.setRange("Right",129.75,160.);
            hMass.setRange("Complete",105.,160.);


            //Sacar los resultados y la PDF del Workspace
            RooFitResult *result_bkg = (RooFitResult*)w->obj(funcList[i] + "_fit_result");

            RooArgList BackgroundParameterList = result_bkg->floatParsFinal();
            RooGenericPdf *pdf;


            double jmin = (120.125 - 105.0)/hist_data->GetBinWidth(1); 
            double jmax = (129.75-105.0)/hist_data->GetBinWidth(1);

            double L = 30.0;

            //Obtaining fit results from the signal

            TString SignalFileName= "outFileSignal_h_myy_pT_" + std::to_string(i) + ".root";
            TFile *file_signal;
            file_signal = new TFile(thePathSignal + SignalFileName,"READ");

            TString SignalWorkspaceName = "w_h_myy_pT_" + std::to_string(i);

            RooWorkspace *w_signal = (RooWorkspace*)file_signal->Get(SignalWorkspaceName);

            RooFitResult *result_signal = (RooFitResult*)w_signal->obj("fit_result");

            RooArgList SignalParameterList = result_signal->floatParsFinal();

            //Double-Sided-Crystal-Ball
            double par[6]; //These parameters result from the previous fit of the signal-only MC


            RooRealVar* var_mu = (RooRealVar*) SignalParameterList.find("mu");
            RooRealVar* var_width = (RooRealVar*) SignalParameterList.find("width");
            RooRealVar* var_a1 = (RooRealVar*) SignalParameterList.find("a1");
            RooRealVar* var_p1 = (RooRealVar*) SignalParameterList.find("p1");
            RooRealVar* var_a2 = (RooRealVar*) SignalParameterList.find("a2");
            RooRealVar* var_p2 = (RooRealVar*) SignalParameterList.find("p2");

            par[0] = var_mu->getVal();
            par[1] = var_width->getVal();
            par[2] = var_a1->getVal();
            par[3] = var_p1->getVal();
            par[4] = var_a2->getVal();
            par[5] = var_p2->getVal();

            //DCB parameters (Not allowed to change, only normalization)
            RooRealVar mu("mu","mu",par[0],123.,127.);
            RooRealVar width("width","width",par[1]);
            RooRealVar a1("a1","a1",par[2]);
            RooRealVar p1("p1","p1",par[3]);
            RooRealVar a2("a2","a2",par[4]);
            RooRealVar p2("p2","p2",par[5]);

            RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);


            //Variable definitions:

            RooRealVar c("c","c",1.0,-10.0,10.0);
            RooRealVar powerExponent1("powerExponent1","powerExponent1",0.,-5.0,5.0);
            RooRealVar c0("c0","c0",1.0,-5.0,5.0);
            RooRealVar c1("c1","c1",1.0,-5.0,5.0);
            RooRealVar c2("c2","c2",1.0,-5.0,5.0);
            RooRealVar c3("c3","c3",1.0,-5.0,5.0);
            RooRealVar b1("b1","b1",1.0,-10.0,2.0);
            RooRealVar b2("b2","b2",1.0,-10.0,2.0);
            RooRealVar d1("d1","d1",1.0,-50.0,50.0);
            RooRealVar d2("d2","d2",1.0,-50.0,50.0);
            RooRealVar d3("d3","d3",1.0,-50.0,50.0);


            RooRealVar N_SS("N_SS","N_SS",1500,0,2000);
            RooRealVar N_bkg("N_bkg","N_bkg",6000,0,25000);

            RooAddPdf *total_pdf;

            //Cargarlos según la función necesaria
            if (funcList[i] == "Power"){
                

                RooRealVar* var_powerExponent1 = (RooRealVar*) BackgroundParameterList.find("powerExponent1");
                powerExponent1.setVal(var_powerExponent1->getVal());
                pdf = new RooGenericPdf("Power","TMath::Power(x[0]/100.,x[1])",RooArgList(hMass,powerExponent1));

                total_pdf = new RooAddPdf("Total","Sum of PDFs",RooArgList(dcbPdf,*pdf),RooArgList(N_SS,N_bkg));

            }
            else if(funcList[i] == "Bernstein"){

                RooRealVar* var_c0 = (RooRealVar*) BackgroundParameterList.find("c0");
                RooRealVar* var_c1 = (RooRealVar*) BackgroundParameterList.find("c1");
                RooRealVar* var_c2 = (RooRealVar*) BackgroundParameterList.find("c2");
                RooRealVar* var_c3 = (RooRealVar*) BackgroundParameterList.find("c3");

                c0.setVal(var_c0->getVal());
                c1.setVal(var_c1->getVal());
                c2.setVal(var_c2->getVal());
                c3.setVal(var_c3->getVal());   

                pdf = new RooGenericPdf("Bernstein","x[1]*TMath::Power(1-x[0],3) + x[2]*3*x[0]*(1-x[0])*(1-x[0]) + x[3]*3*x[0]*x[0]*(1-x[0]) + x[4]*TMath::Power(x[0],3)",RooArgList(hMass,c0,c1,c2,c3));

                total_pdf = new RooAddPdf("Total","Total",RooArgList(dcbPdf,*pdf),RooArgList(N_SS,N_bkg));

            }
            else if(funcList[i] == "ExpPol1"){

                RooRealVar* var_c = (RooRealVar*) BackgroundParameterList.find("c");

                c.setVal(var_c->getVal());                

                pdf = new RooGenericPdf("ExpPol1","TMath::Exp(x[0]*x[1])",RooArgList(hMass,c));

                total_pdf = new RooAddPdf("Total","Total",RooArgList(dcbPdf,*pdf),RooArgList(N_SS,N_bkg));
                
            }
            else if(funcList[i] == "ExpPol2"){

                RooRealVar* var_b1 = (RooRealVar*) BackgroundParameterList.find("b2");
                RooRealVar* var_b2 = (RooRealVar*) BackgroundParameterList.find("b2");

                b1.setVal(var_b1->getVal());
                b2.setVal(var_b2->getVal());

                pdf = new RooGenericPdf("ExpPol2","TMath::Exp(x[1]*TMath::Power(x[0]-100.,1)/100. + x[2]*TMath::Power(x[0]-100.,2)/10000.)",RooArgList(hMass,b1,b2));

                total_pdf = new RooAddPdf("Total","Total",RooArgList(dcbPdf,*pdf),RooArgList(N_SS,N_bkg));

            }
            else if(funcList[i] == "ExpPol3"){

                RooRealVar* var_d1 = (RooRealVar*) BackgroundParameterList.find("d1");
                RooRealVar* var_d2 = (RooRealVar*) BackgroundParameterList.find("d2");
                RooRealVar* var_d3 = (RooRealVar*) BackgroundParameterList.find("d3");

                d1.setVal(var_d1->getVal());
                d2.setVal(var_d2->getVal());
                d3.setVal(var_d3->getVal());

                pdf = new RooGenericPdf("ExpPol3","TMath::Exp(x[1]*TMath::Power(x[0]-100.,1)/100. + x[2]*TMath::Power(x[0]-100.,2)/10000. + x[3]*TMath::Power(x[0]-100.,3)/1000000.)",RooArgList(hMass,d1,d2,d3));
                
                total_pdf = new RooAddPdf("Total","Total",RooArgList(dcbPdf,*pdf),RooArgList(N_SS,N_bkg));

            }
            else{
                return;
            }


            RooFitResult *result = total_pdf->fitTo(hist, RooFit::Range("Left,Right"), RooFit::Minimizer("Minuit2"), RooFit::Save(true) );

            result->Print();


            RooArgList ParameterList = result->floatParsFinal();
            RooRealVar* var_NSignal = (RooRealVar*) ParameterList.find("N_SS");
            double NSignal = var_NSignal->getVal();
    


            double wid = 650;
            double he = 450;
            TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,wid,he);
            canvasn->cd(); 


            auto pl = hMass.frame(RooFit::Name("frame"));
            

            //FitHist->plotOn(pl,RooFit::Name("FitHist"));
            hist.plotOn(pl,RooFit::Name("hist"));
            total_pdf->plotOn(pl,RooFit::Name("pdf"),RooFit::Range("Complete")); //Normaliza a los eventos o esto vale directamente?
            

            pl->Draw();

            TString OutFileName = "/home/mjlarroy/TFG_TUT/figures/pdf_graphs/XS/CrossSection_SS_method_bin" + std::to_string(i) + ".pdf";


            canvasn->Print(OutFileName);

            //Parte de ver cuántos eventos de señal encontramos. Entre 120.125 y 129.75 en principio. Lo haremos restando los histogramas e integrando entre jmin y jmax

            h0yypT->SetBinContent(i,NSignal/((h0yypT->GetXaxis()->GetBinWidth(i))*L));


            
            delete pdf;
            
        }
        TCanvas *canvas_cs = new TCanvas("canvas_cs","canvas_cs",0,0,650,450);
        canvas_cs->cd();
        canvas_cs->SetLogx();

        h0yypT->Draw("P");

        canvas_cs->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/XS/CrossSection_SS_method_complete.pdf");


        
        return;
        



    }

