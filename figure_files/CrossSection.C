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
    

    void CrossSection(){

        Double_t xbins[30] = {0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0,36.0,39.0,42.0,45.0,48.0,51.0,54.0,57.0,60.0,69.0,78.0,87.0,96.0,132.0,186.0,258.0,618.0,1500.0};
        Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};

        TH1* h0yypT = new TH1D("h0yypT","PTyy",16,new_xbins);

        TString thePath = "/home/mjlarroy/TFG_TUT/Root_Results/outFile_Background/";
        TString thePathData = "/home/mjlarroy/TFG_TUT/Root_Results/new_data22_bins.root";

        TFile *data;
        data = new TFile(thePathData,"read");

        for(int i=2;i<3;i++){


            //Traer el Workspace correspondiente del ajuste al bcg
            TString FileName= "outFileBackground_" + new_funcList[i] + "_new_h_myy_pT_" + std::to_string(i) + ".root";

            TFile *file_ws;
            file_ws = new TFile(thePath + FileName,"READ");

            TString WorkspaceName = "w_out_" + new_funcList[i] + "_new_h_myy_pT_" + std::to_string(i);

            RooWorkspace *w = (RooWorkspace*)file_ws->Get(WorkspaceName);


            //Datos
            TH1D *hist_data;
            TString HistName = "new_h_myy_pT_" + std::to_string(i);
            hist_data = (TH1D*)data->Get(HistName); // leer histograma de datos

            RooRealVar hMass("myy","myy [GeV]",105.0,160.0);


            RooDataHist hist("hist","hist",hMass,Import(*hist_data));


            hMass.setRange("Left",105.,120.125);
            hMass.setRange("Right",129.75,160.);
            hMass.setRange("Complete",105.,160.);
            hMass.setRange("Reduced",120.125,129.75);



            //Sacar los resultados y la PDF del Workspace
            RooFitResult *result_bkg = (RooFitResult*)w->obj(new_funcList[i] + "_fit_result");

            RooArgList SignalParameterList = result_bkg->floatParsFinal();
            RooGenericPdf *pdf;

            RooArgList *ArgList;

            double jmin = (120.125 - 105.0)/hist_data->GetBinWidth(1); 
            double jmax = (129.75-105.0)/hist_data->GetBinWidth(1);

            double L = 30.0; //fb^-1


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
            RooRealVar N("N","N",10000.,2000.,1000000.);


            //Cargarlos según la función necesaria
            if (new_funcList[i] == "Power"){
                

                RooRealVar* var_powerExponent1 = (RooRealVar*) SignalParameterList.find("powerExponent1");
                powerExponent1.setVal(var_powerExponent1->getVal());
                pdf = new RooGenericPdf("Power","TMath::Power(x[0]/100.,x[1])",RooArgList(hMass,powerExponent1));

            }
            else if(new_funcList[i] == "Bernstein"){

                RooRealVar* var_c0 = (RooRealVar*) SignalParameterList.find("c0");
                RooRealVar* var_c1 = (RooRealVar*) SignalParameterList.find("c1");
                RooRealVar* var_c2 = (RooRealVar*) SignalParameterList.find("c2");
                RooRealVar* var_c3 = (RooRealVar*) SignalParameterList.find("c3");

                c0.setVal(var_c0->getVal());
                c1.setVal(var_c1->getVal());
                c2.setVal(var_c2->getVal());
                c3.setVal(var_c3->getVal());   

                pdf = new RooGenericPdf("Bernstein","x[1]*TMath::Power(1-x[0],3) + x[2]*3*x[0]*(1-x[0])*(1-x[0]) + x[3]*3*x[0]*x[0]*(1-x[0]) + x[4]*TMath::Power(x[0],3)",RooArgList(hMass,c0,c1,c2,c3));

            }
            else if(new_funcList[i] == "ExpPol1"){

                RooRealVar* var_c = (RooRealVar*) SignalParameterList.find("c");

                c.setVal(var_c->getVal());                

                pdf = new RooGenericPdf("ExpPol1","TMath::Exp(x[0]*x[1])",RooArgList(hMass,c));                
            }
            else if(new_funcList[i] == "ExpPol2"){

                RooRealVar* var_b1 = (RooRealVar*) SignalParameterList.find("b2");
                RooRealVar* var_b2 = (RooRealVar*) SignalParameterList.find("b2");

                b1.setVal(var_b1->getVal());
                b2.setVal(var_b2->getVal());

                pdf = new RooGenericPdf("ExpPol2","TMath::Exp(x[1]*TMath::Power(x[0]-100.,1)/100. + x[2]*TMath::Power(x[0]-100.,2)/10000.)",RooArgList(hMass,b1,b2));

            }
            else if(new_funcList[i] == "ExpPol3"){

                RooRealVar* var_d1 = (RooRealVar*) SignalParameterList.find("d1");
                RooRealVar* var_d2 = (RooRealVar*) SignalParameterList.find("d2");
                RooRealVar* var_d3 = (RooRealVar*) SignalParameterList.find("d3");

                d1.setVal(var_d1->getVal());
                d2.setVal(var_d2->getVal());
                d3.setVal(var_d3->getVal());

                pdf = new RooGenericPdf("ExpPol3","TMath::Exp(x[1]*TMath::Power(x[0]-100.,1)/100. + x[2]*TMath::Power(x[0]-100.,2)/10000. + x[3]*TMath::Power(x[0]-100.,3)/1000000.)",RooArgList(hMass,d1,d2,d3));
                
            }
            else{
                return;
            }


            RooFitResult *result = pdf->fitTo(hist, RooFit::Range("Left,Right"), RooFit::Minimizer("Minuit2"), RooFit::Save(true) );

            result->Print();

            double norm = pdf->getNorm(RooArgSet(hMass));
    
           
            // Print the normalization
            std::cout << "Normalization: " << norm << std::endl;


            double NEvents = hist_data->Integral();
            cout<<NEvents<<endl;

            auto FitHist = pdf->createHistogram("FitHist", hMass, RooFit::Binning(80, 105., 160.));

            //FitHist->Scale(NEvents); //Número de eventos de los datos
            //FitHist->Scale(norm_param->getVal()); //Parámetro de normalización del Fit
            cout<<FitHist->GetBinContent(1)<<endl;
            FitHist->Scale(NEvents);
            cout<<FitHist->GetBinContent(1)<<endl;


            double wi = 650;
            double he = 450;
            TCanvas *canvasn = new TCanvas("canvasn","canvasn",0,0,wi,he);
            canvasn->cd(); 


            auto pl = hMass.frame(RooFit::Name("frame"));
            

            //FitHist->plotOn(pl,RooFit::Name("FitHist"));
            hist.plotOn(pl,RooFit::Name("hist"));
            pdf->plotOn(pl,RooFit::Name("pdf"),RooFit::Range("Complete")); //Normaliza a los eventos o esto vale directamente?
            

            pl->Draw();

            TString OutFileName = "/home/mjlarroy/TFG_TUT/figures/pdf_graphs/XS/new_CrossSection_bin" + std::to_string(i) + ".pdf";


            canvasn->Print(OutFileName);

            //Parte de ver cuántos eventos de señal encontramos. Entre 120.125 y 129.75 en principio. Lo haremos restando los histogramas e integrando entre jmin y jmax

            TH1D* hist_dif;
            hist_dif = (TH1D*)hist_data->Clone(HistName);
            hist_dif->SetName("hist_dif");

            hist_dif->Add(FitHist,-1);


            //
            //
            /////////////////////////////////////////////////////////////////////////////////////////////////////
            //PRUEBA DE FITAR EL DSCB A LA RESTA:

            //Obtaining fit results from the signal
            TString thePathSignal = "/home/mjlarroy/TFG_TUT/Root_Results/";
            TString SignalFileName= "outFileSignal_h_myy_pT_" + std::to_string(i) + ".root";
            TFile *file_signal;
            file_signal = new TFile(thePathSignal + SignalFileName,"READ");

            TString SignalWorkspaceName = "w_h_myy_pT_" + std::to_string(i);

            RooWorkspace *w_signal = (RooWorkspace*)file_signal->Get(SignalWorkspaceName);

            RooFitResult *result_signal = (RooFitResult*)w_signal->obj("fit_result");

            RooArgList SignalParameterList2 = result_signal->floatParsFinal();

            //Double-Sided-Crystal-Ball
            double par[6]; //These parameters result from the previous fit of the signal-only MC


            RooRealVar* var_mu = (RooRealVar*) SignalParameterList2.find("mu");
            RooRealVar* var_width = (RooRealVar*) SignalParameterList2.find("width");
            RooRealVar* var_a1 = (RooRealVar*) SignalParameterList2.find("a1");
            RooRealVar* var_p1 = (RooRealVar*) SignalParameterList2.find("p1");
            RooRealVar* var_a2 = (RooRealVar*) SignalParameterList2.find("a2");
            RooRealVar* var_p2 = (RooRealVar*) SignalParameterList2.find("p2");

            par[0] = var_mu->getVal();
            par[1] = var_width->getVal();
            par[2] = var_a1->getVal();
            par[3] = var_p1->getVal();
            par[4] = var_a2->getVal();
            par[5] = var_p2->getVal();

            

            //DCB parameters (Not allowed to change, only normalization)
            RooRealVar mu("mu","mu",par[0],122.,128.);
            RooRealVar width("width","width",par[1]);
            RooRealVar a1("a1","a1",par[2]);
            RooRealVar p1("p1","p1",par[3]);
            RooRealVar a2("a2","a2",par[4]);
            RooRealVar p2("p2","p2",par[5]);

            RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB_RF",hMass,mu,width,a1,p1,a2,p2);



            RooDataHist hist2("hist","hist",hMass,Import(*hist_dif));


            RooFitResult *result_signal2 = dcbPdf.fitTo(hist2,  RooFit::Range("Reduced"),RooFit::Minimizer("Minuit2"), RooFit::Save(true) );

            result_signal2->Print();

            norm = dcbPdf.getNorm(RooArgSet(hMass));
    
           
            // Print the normalization
            std::cout << "Normalization 2: " << norm << std::endl;



            TCanvas *canvasp = new TCanvas("canvasp","canvasp",0,0,wi,he);
            canvasp->cd();
            auto pl2 = hMass.frame(RooFit::Name("frame"));
            

            //FitHist->plotOn(pl,RooFit::Name("FitHist"));
            hist2.plotOn(pl2,RooFit::Name("hist2"));
            dcbPdf.plotOn(pl2,RooFit::Name("dcbpdf"),RooFit::Range("Complete"));

            pl2->Draw();


            //Aquí acabe la prueba cuidado.
            /////////////////////////////////////////////////////////////////////////////////////
/*
            TCanvas *canvasp = new TCanvas("canvasp","canvasp",0,0,width,he);
            canvasp->cd();
            
            hist_data->Draw();
            FitHist->Draw("same");
            hist_dif->Draw("same");

            TString OutFileName2 = "/home/mjlarroy/TFG_TUT/figures/pdf_graphs/XS/new_CrossSection_TH1_bin" + std::to_string(i) + ".pdf";

            canvasp->Print(OutFileName2);

            h0yypT->SetBinContent(i,(hist_dif->Integral(jmin,jmax))/((h0yypT->GetXaxis()->GetBinWidth(i))*L));


            */
            delete pdf;
            
        }
        /*
        TCanvas *canvas_cs = new TCanvas("canvas_cs","canvas_cs",0,0,650,450);
        canvas_cs->cd();
        canvas_cs->SetLogx();

        h0yypT->Draw("P");

        canvas_cs->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/XS/new_CrossSection_complete.pdf");*/


        
        return;
        



    }

