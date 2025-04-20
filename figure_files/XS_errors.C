//Para correrlo hacemos root -l macro_figure_ptyy.C

{
    // clase que convierte cualquier tipo de variable a variable string
    ostringstream histo;
    //
    // titulo del eje x
    string xt;
    xt = "p^{T}_{#gamma#gamma} [GeV]"; 
    // titulo eje y
    string yt;
    yt = "d#sigma/dp^{T}_{#gamma#gamma} [fb/GeV]";

    int ggH_variations[11] = {88,94,95,96,97,98,99,100,101,102,111}; //0 es NLO y 88 es NNLO que en caso de ggH es el nominal. El 0 lo omito ahora

    //LEER HISTOGRAMA CON EL PASOD E TRUTH A RECO (COMPARAMOS CON RECO LUEGO LA SECCIÑÓN EFICAZ)

    TFile *reco;  // root de todos los MC del que voy a leer
    TH1D *hist_reco; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/reco.root"; // nombre del root de los MC
    reco = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h"; // nombre del histograma que quiero leer
    hist_reco = (TH1D*)reco->Get(histo.str().c_str()); // leer histograma MC

    //READ MC HISTOGRAMS FOR EVERY PRODUCTION MODE
  
    TFile *VBF;  // root de todos los MC del que voy a leer
    TH1D *hist_VBF; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/VBF.root"; // nombre del root de los MC
    VBF = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h_pTTruth"; // nombre del histograma que quiero leer
    hist_VBF = (TH1D*)VBF->Get(histo.str().c_str()); // leer histograma MC

    TFile *ttH;  // root de todos los MC del que voy a leer
    TH1D *hist_ttH; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/ttH.root"; // nombre del root de los MC
    ttH = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h_pTTruth"; // nombre del histograma que quiero leer
    hist_ttH = (TH1D*)ttH->Get(histo.str().c_str()); // leer histograma MC

    TFile *ZH;  // root de todos los MC del que voy a leer
    TH1D *hist_ZH; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/ZH.root"; // nombre del root de los MC
    ZH = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h_pTTruth"; // nombre del histograma que quiero leer
    hist_ZH = (TH1D*)ZH->Get(histo.str().c_str()); // leer histograma MC

    TFile *WpH;  // root de todos los MC del que voy a leer
    TH1D *hist_WpH; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/WpH.root"; // nombre del root de los MC
    WpH = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h_pTTruth"; // nombre del histograma que quiero leer
    hist_WpH = (TH1D*)WpH->Get(histo.str().c_str()); // leer histograma MC

    TFile *WmH;  // root de todos los MC del que voy a leer
    TH1D *hist_WmH; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/WmH.root"; // nombre del root de los MC
    WmH = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h_pTTruth"; // nombre del histograma que quiero leer
    hist_WmH = (TH1D*)WmH->Get(histo.str().c_str()); // leer histograma MC

    


    Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};
    double numBinspT = 16;

    TH1D* hist_total = new TH1D("h_pTTruth_total", "h_pTTruth_total", numBinspT, new_xbins);
    hist_total->Add(hist_VBF);
    hist_total->Add(hist_ttH);
    hist_total->Add(hist_WmH);
    hist_total->Add(hist_WpH);
    hist_total->Add(hist_ZH);

    TH1D* hist_total_superior = new TH1D("h_pTTruth_superior", "h_pTTruth_superior", numBinspT, new_xbins);
    TH1D* hist_total_inferior = new TH1D("h_pTTruth_inferior", "h_pTTruth_inferior", numBinspT, new_xbins);

    std::map<std::string, TH1D*> hMappT;
    std::map<std::string, TH1D*> hMappT_ggH;

    //88,94,95,96,97,98,99,100,101,102,111




    for(int i=0;i<11;i++){


        std::string histName = "h_pT_truth_var_" + std::to_string(ggH_variations[i]);


        hMappT[histName] = new TH1D(histName.c_str(), histName.c_str(), numBinspT, new_xbins);
        
        hMappT[histName]->Add(hist_VBF);
        hMappT[histName]->Add(hist_ttH);
        hMappT[histName]->Add(hist_WmH);
        hMappT[histName]->Add(hist_WpH);
        hMappT[histName]->Add(hist_ZH);

        TFile *ggH;  // root de todos los MC del que voy a leer
        histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/ggH.root"; // nombre del root de los MC
        ggH = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
        TH1D *hist_ggH; // histograma de todos los MC que voy a mostrar


        hist_ggH = (TH1D*)ggH->Get(histName.c_str()); // leer histograma MC
        

        hMappT[histName]->Add(hist_ggH);

    }

    double value = 0.0;
    double anchura = 0.0;
    double L = 27.58; //Luiminosidad integrada
    double min_value = 0.0;
    double max_value = 0.0;


    for(int j=1;j<=numBinspT;j++){

        min_value = 10000.0;
        max_value = -10.0;

        
        for(int i=0;i<11;i++){
            std::string histName = "h_pT_truth_var_" + std::to_string(ggH_variations[i]);

            value = hMappT[histName]->GetBinContent(j);
            anchura = hMappT[histName]->GetBinWidth(j);

            hMappT[histName]->SetBinContent(j,value/(anchura*L)); //Los pongo a nivel particle/truth

            hMappT[histName]->SetBinError(j,0.0); //POR RAZONES ESTÉTICAS Y ADEMÁS ES EL PROPIO VALOR EL QUE USAREMOS COMO ERROR LUEGO


            min_value = std::min(min_value,hMappT[histName]->GetBinContent(j));
            max_value = std::max(max_value,hMappT[histName]->GetBinContent(j));


        }


        hist_total_inferior->SetBinContent(j,min_value);
        hist_total_superior->SetBinContent(j,max_value);
    }
    

    int ni = 1; int nb = numBinspT; // numero de bines en el eje X
    cout<<"n bins = "<<nb<<endl;
  
  
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
  
    TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.0,0.0,0.9,0.99);
  
    // pintar pads
    canvas_1->SetFillStyle(4000); //transparentes
    canvas_1->Draw();

  
    canvas_1->cd();
    //canvas_1->SetLogy();
    canvas_1->SetLogx();
    // crear frame
    double bl = 0.; // principio del eje X
    double bu = 1600.; // final del eje X
    //double yl = histdat->GetMinimum();  // minimo eje Y
    double yl = 0.0; // minimo eje Y. el 0.5 por logaithmic scale
    double yu = 1.8;
  
    TH1F *h1 = canvas_1->DrawFrame(bl,yl,bu,yu); // crea el frame
  
    //Te renta hacer el Frame para poner todos los settings de una vez y no sobre cada histograma y estas cosas.
  
  
    // titulo eje y
    histo.str(""); histo << yt;
    h1->SetYTitle(histo.str().c_str());
    // settings eje y
    h1->GetYaxis()->SetTitleOffset(1.5);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetYaxis()->SetNdivisions(505);
    //
    // settings eje x
    histo.str(""); histo << xt;
    h1->SetXTitle(histo.str().c_str());
    //
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleOffset(1.1);
    //
    h1->GetXaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetNdivisions(505);
    h1->GetXaxis()->SetMoreLogLabels(1);
    h1->GetXaxis()->SetNoExponent(1);
    //
    // titulo eje x
  
  
  
    //
    // pintar el frame
    h1->Draw();
  
    //Plotear

    
/*
    for(int i=0;i<11;i++){

        std::string histName = "h_pT_truth_var_" + std::to_string(ggH_variations[i]);

        if (i>9){

        myhhdate(hMappT[histName],"SAME][",0.0,i+4,0,2.,1,i+4,1000,0);
        continue;
        }

        myhhdate(hMappT[histName],"SAME][",0.0,i+1,0,2.,1,i+1,1000,0);

        
    }
    //88,94,95,96,97,98,99,100,101,102,111

    myLine(x1+0.35,0.85,0.05,1,20,"88",0.03,0.03);
    myLine(x1+0.35,0.82,0.05,2,20,"94",0.03,0.03);
    myLine(x1+0.35,0.79,0.05,3,20,"95",0.03,0.03);
    myLine(x1+0.35,0.76,0.05,4,20,"96",0.03,0.03);
    myLine(x1+0.35,0.73,0.05,5,20,"97",0.03,0.03);
    myLine(x1+0.35,0.70,0.05,6,20,"98",0.03,0.03);
    myLine(x1+0.35,0.67,0.05,7,20,"99",0.03,0.03);
    myLine(x1+0.35,0.64,0.05,8,20,"100",0.03,0.03);
    myLine(x1+0.35,0.61,0.05,9,20,"101",0.03,0.03);
    myLine(x1+0.35,0.58,0.05,13,20,"102",0.03,0.03);
    myLine(x1+0.35,0.55,0.05,14,20,"111",0.03,0.03);*/

    myhhdate(hist_total_superior,"same][",0.0,2,0,2.,1,2,1000,0); //definido en myAtlasUtils.C
    myhhdate(hist_total_inferior,"same][",0.0,4,0,2.,1,4,1000,0); //definido en myAtlasUtils.C


  
    double x1 = 0.35;double x2 = 0.50;
    double yy1 = 0.90; 
    double y2 = yy1 - 0.04;
  
    TLatex l1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
    l1.SetNDC();
    l1.SetTextFont(72);
    l1.SetTextColor(1);
    l1.SetTextSize(0.035);
    l1.DrawLatex(x1+0.30,y2-0.3,"ATLAS");
    l1.SetTextFont(42);
  
    l1.DrawLatex(x1+0.30,y2-0.35,"#sqrt{s} = 13,6 TeV");
    l1.DrawLatex(x1+0.30,y2-0.4,"L = 27.58 fb^{-1}");
  
  
    gPad->RedrawAxis();
  
    myLine(x1+0.36,0.8,0.05,2,1000,"Superior",0.03,0.03);
    myLine(x1+0.36,0.85,0.05,4,20,"Inferior",0.03,0.03);
    
  
    canvasp->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/XS_errors.pdf");

    TFile* myfile = new TFile("/home/mjlarroy/TFG_TUT/Root_Results/XS_errors.root", "RECREATE"); //Recreate significa que si no está lo cree y si está que lo reescriba.

   myfile->cd();
   for (auto& pair : hMappT){
     delete pair.second;  // Clean up memory
   }
   hist_total_inferior->Write();
   hist_total_superior->Write();
   
   myfile->Close();

  }