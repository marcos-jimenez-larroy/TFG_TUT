/Para correrlo hacemos root -l reco.C

{
    // clase que convierte cualquier tipo de variable a variable string
    ostringstream histo;
    //
    // titulo del eje x
    string xt;
    xt = "Product Cut Purity"; 
    // titulo eje y
    string yt;
    yt = "Product Cut Efficiency";

    double cuts[12] = {20,22,24,26,28,30,32,34,36,38,40,42};
    double efficiency[12] = {};
    double purity[12] = {};


    for(int i=0; i++; i<12){


    

    TFile *mc_signal;
    TH1D *hist_signal;
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_ProductCuts_" << std::to_string(i) << ".root";
    mc_signal = new TFile(histo.str().c_str(),"read");
  
    histo.str("");  histo << "h0yymass";
    hist = (TH1D*)mc_signal->Get(histo.str().c_str());

    TH1D *hist_signal_PC;
    histo.str("");  histo << "h0yymass_PC";
    hist_signal_PC = (TH1D*)mc_signal->Get(histo.str().c_str());


    TFile *mc_background;
    TH1D *hist_background_PC;

    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/mc23a_background_ProductCuts_32.root";
    mc_background = new TFile(histo.str().c_str(),"read");
    histo.str("");  histo << "h0yymass_PC";
    hist_background_PC = (TH1D*)mc_background->Get(histo.str().c_str());


    efficiency[i] = hist_signal_PC->Integral()/hist_signal->Integral();
    purity[i] = hist_signal_PC->Integral()/hist_background_PC->Integral();


    }






}