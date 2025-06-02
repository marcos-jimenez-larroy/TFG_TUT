{
    ostringstream histo;
  
    string xt;
    xt = "Purity"; 
    // titulo eje y
    string yt;
    yt = "Efficiency";

    Double_t cuts[31] = {27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,41.5,42}; //In percentage

    Double_t Efficiency_PC[31] = {};
    Double_t Purity_PC[31] = {};

    Double_t Efficiency_RC[31] = {};
    Double_t Purity_RC[31] = {};

    double num_cuts = 31;



    TFile *signal;  // root de todos los MC del que voy a leer
    TH1D *setcuts_hist; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_ProductCuts_RelativeCuts_VBF.root"; // nombre del root de los MC
    signal = new TFile(histo.str().c_str(),"read"); // abrir fichero

    histo.str("");  histo << "h0yymass"; // nombre del histograma que quiero leer
    setcuts_hist = (TH1D*)signal->Get(histo.str().c_str()); // leer histograma MC

    TFile *background;
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/mc23a_background_ProductCuts_RelativeCuts.root"; // nombre del root de los MC
    background = new TFile(histo.str().c_str(),"read"); // abrir fichero


    for(int i=0; i<num_cuts;i++){

        TH1D *PC_signal_hist; 
        TH1D *PC_background_hist;
        TH1D *RC_signal_hist; 
        TH1D *RC_background_hist;


        histo.str("");  histo << "PC_" << std::to_string(i); // nombre del histograma que quiero leer
        PC_signal_hist = (TH1D*)signal->Get(histo.str().c_str()); // leer histograma MC
        PC_background_hist = (TH1D*)background->Get(histo.str().c_str());

        Efficiency_PC[i] = PC_signal_hist->Integral()/setcuts_hist->Integral();
        Purity_PC[i] = PC_signal_hist->Integral()/PC_background_hist->Integral();

        histo.str("");  histo << "RC_" << std::to_string(i); // nombre del histograma que quiero leer
        RC_signal_hist = (TH1D*)signal->Get(histo.str().c_str()); // leer histograma MC
        RC_background_hist = (TH1D*)background->Get(histo.str().c_str());

        Efficiency_RC[i] = RC_signal_hist->Integral()/setcuts_hist->Integral();
        Purity_RC[i] = RC_signal_hist->Integral()/RC_background_hist->Integral();
    }

    for(int i=0; i<num_cuts; i++){
        std::cout<<Form("Efficiency for PC of %.3f = %f", cuts[i]/100., Efficiency_PC[i])<<std::endl;
    }
    for(int i=0; i<num_cuts; i++){
        std::cout<<Form("Purity for PC of %.3f =  %.8f", cuts[i]/100., Purity_PC[i])<<std::endl;
    }
    for(int i=0; i<num_cuts; i++){
        std::cout<<Form("Efficiency for RC of %.3f = %f", cuts[i]/100., Efficiency_RC[i])<<std::endl;
    }
    for(int i=0; i<num_cuts; i++){
        std::cout<<Form("Purity for RC of %.3f =  %.8f", cuts[i]/100., Purity_RC[i])<<std::endl;
    }
/*
    TGraph *g_PC = new TGraph(num_cuts,Purity_PC,Efficiency_PC);
    g_PC->Draw("ap");
    TGraph *g_RC = new TGraph(num_cuts,Purity_RC,Efficiency_RC);
    g_RC->Draw("ap");*/






















}