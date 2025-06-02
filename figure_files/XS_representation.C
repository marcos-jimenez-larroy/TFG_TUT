{
  ostringstream histo;

  string xt;
  xt = "p^{T}_{#gamma#gamma} [GeV]"; 
  // titulo eje y
  string yt;
  yt = "d#sigma/dp^{T}_{#gamma#gamma} [fb/GeV]";

  TFile *mc_xs;  // root de todos los MC del que voy a leer
  TH1D *mc_xs_hist; // histograma de todos los MC que voy a mostrar
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/MC_ptyy_fiducial_cross_section.root"; // nombre del root de los MC
  mc_xs = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "hist_rebinned"; // nombre del histograma que quiero leer
  mc_xs_hist = (TH1D*)mc_xs->Get(histo.str().c_str()); // leer histograma MC


  TFile *exp_xs;  // root de todos los MC del que voy a leer
  TH1D *exp_xs_hist; // histograma de todos los MC que voy a mostrar
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/C_Experimental_CrossSection_22_23.root"; // nombre del root de los MC
  exp_xs = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "h0yypT"; // nombre del histograma que quiero leer
  exp_xs_hist = (TH1D*)exp_xs->Get(histo.str().c_str()); // leer histograma MC

  TFile *errors;  // root de todos los MC del que voy a leer
  TH1D *superior; // histograma de todos los MC que voy a mostrar
  TH1D *inferior;
  TH1D *NLO;
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/XS_errors.root"; // nombre del root de los MC
  errors = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "h_pTTruth_superior"; // nombre del histograma que quiero leer
  superior = (TH1D*)errors->Get(histo.str().c_str()); // leer histograma MC

  histo.str("");  histo << "h_pTTruth_inferior"; // nombre del histograma que quiero leer
  inferior = (TH1D*)errors->Get(histo.str().c_str()); // leer histograma MC

  histo.str("");  histo << "h_pT_truth_var_0"; // nombre del histograma que quiero leer
  NLO = (TH1D*)errors->Get(histo.str().c_str()); // leer histograma MC





  int ni = 1; int nb = superior->GetNbinsX(); // numero de bines en el eje X

  TFile *reco;  // root de todos los MC del que voy a leer
    TH1D *hist_reco; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/reco.root"; // nombre del root de los MC
    reco = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h"; // nombre del histograma que quiero leer
    hist_reco = (TH1D*)reco->Get(histo.str().c_str()); // leer histograma MC


  //Como los errores son asimétricos usaremos un TGraphAsymmErrors. Resultados del SS test tmb

  Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};

  Double_t sys_error[16] = {22.235,13.6853,12.5925,6.81755,3.15433,2.8733,13.3374,11.9845,12.4012,11.7149,2.04462,2.08365,0.643965,0.372638,0.375956,0.0266472}; //Con Power al frente->Valores del NSS

  Double_t N_expected[16] = {83.4173,167.805,168.445,147.158,123.663,102.871,123.312,94.5311,94.9008,121.134,69.751,43.8252,52.8856,34.3136,8.42659,0.366741};

  Double_t Bkg_uncert[16] = {50.9912,77.6312,74.9757,70.3685,64.7363,59.3115,64.6205,49.5746,54.9216,57.2969,36.5536,22.2011,19.6753,11.7251,4.12488,0.682029}; //Con Power al frente


  TH1D *his = new TH1D("his","his",16,new_xbins);




  Double_t x_positions[16] = {};
  Double_t y_positions[16] = {};
  Double_t y_positions_exp[16] = {};
  Double_t ratio[16] = {};
  Double_t ratio_NLO[16] = {};

  Double_t x_errors[16] = {};
  Double_t y_sup_errors[16] = {};
  Double_t y_sup_rel_errors[16] = {};
  Double_t y_inf_errors[16] = {};
  Double_t y_inf_rel_errors[16] = {};
  Double_t rel_error[16] = {};
  Double_t val = 0.0;
  Double_t exp_stat_error[16] = {};
  Double_t exp_total_error[16] = {};
  Double_t NLO_values[16] = {};
  Double_t zeros[16] = {};
  Double_t ratio_NLO_MC[16] = {};

  Double_t L = 58.98;

  Double_t total_xs = 0.0;
  Double_t total_err_xs = 0.0;

    for(int i=0;i<nb;i++){

      NLO_values[i] = NLO->GetBinContent(i+1)/NLO->Integral(); //Normalizada
      zeros[i] = 0;
      val = exp_xs_hist->GetBinContent(i+1);
      exp_xs_hist->SetBinContent(i+1,val*1.0/hist_reco->GetBinContent(i+1));
      
      his->SetBinContent(i+1,exp_xs_hist->GetBinContent(i+1));

      ratio[i] = mc_xs_hist->GetBinContent(i+1)/exp_xs_hist->GetBinContent(i+1);
      ratio_NLO[i] = NLO_values[i]/exp_xs_hist->GetBinContent(i+1);

      x_positions[i] = mc_xs_hist->GetBinCenter(i+1);
      y_positions[i] = mc_xs_hist->GetBinContent(i+1)/mc_xs_hist->Integral(); //Normalizada
      y_positions_exp[i] = exp_xs_hist->GetBinContent(i+1);

      ratio_NLO_MC[i] = NLO_values[i]/y_positions[i];


      std::cout<<exp_xs_hist->GetBinContent(i+1)<<std::endl;

      x_errors[i] = mc_xs_hist->GetBinWidth(i+1)/2.0;

      y_sup_errors[i] = (superior->GetBinContent(i+1) - mc_xs_hist->GetBinContent(i+1))/mc_xs_hist->Integral();
      y_inf_errors[i] = (-1.0*(inferior->GetBinContent(i+1) - mc_xs_hist->GetBinContent(i+1)))/mc_xs_hist->Integral(); //It is defined as positive for the TGraph
    

      y_sup_rel_errors[i] = y_sup_errors[i]/mc_xs_hist->GetBinContent(i+1);
      y_inf_rel_errors[i] = y_inf_errors[i]/mc_xs_hist->GetBinContent(i+1);

      rel_error[i] = (superior->GetBinContent(i+1)-inferior->GetBinContent(i+1))/mc_xs_hist->GetBinContent(i+1);


      exp_stat_error[i] = exp_xs_hist->GetBinError(i+1);

      exp_total_error[i] = sqrt(exp_stat_error[i]*exp_stat_error[i] + (sys_error[i]*sys_error[i])/((L*exp_xs_hist->GetBinWidth(i+1))*(L*exp_xs_hist->GetBinWidth(i+1))));
      

      total_xs += val*exp_xs_hist->GetBinWidth(i+1);
      total_err_xs += (exp_stat_error[i]*exp_xs_hist->GetBinWidth(i+1))*(exp_stat_error[i]*exp_xs_hist->GetBinWidth(i+1));

    }

    total_err_xs = sqrt(total_err_xs);

    std::cout<<"Experimental total cross section: "<<total_xs<<std::endl;
    std::cout<<"Error: "<<total_err_xs<<std::endl;

    TGraphAsymmErrors *mc_xs_complete = new TGraphAsymmErrors(16,x_positions,y_positions,x_errors,x_errors,y_inf_errors,y_sup_errors);
    TGraphAsymmErrors *ratio_xs_complete = new TGraphAsymmErrors(16,x_positions,ratio,x_errors,x_errors,y_inf_rel_errors,y_sup_rel_errors);
    TGraphAsymmErrors *NLO_xs_complete = new TGraphAsymmErrors(16,x_positions,NLO_values,x_errors,x_errors,zeros,zeros);
    TGraphAsymmErrors *ratio_NLO_xs_complete = new TGraphAsymmErrors(16,x_positions,ratio_NLO,x_errors,x_errors,zeros,zeros);
    TGraphAsymmErrors *ratio_NLO_MC_complete = new TGraphAsymmErrors(16,x_positions,ratio_NLO_MC,x_errors,x_errors,zeros,zeros);

    TGraphAsymmErrors *exp_xs_complete = new TGraphAsymmErrors(16,x_positions,y_positions_exp,x_errors,x_errors,exp_total_error,exp_total_error);

    //Por si lo quieres dibujar aparte el TGraph. Nosotros crearemos dos histogramas (uno sup y otro inf)




    //  PARTE ESTÉTICA DE LA REPRESENTACIÓN
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

  TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.0,0.29,0.9,0.99);
  TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0.0,0.10,0.9,0.41);
  

  // pintar pads
  canvas_1->SetFillStyle(4000); //transparentes
  canvas_1->Draw();
  canvas_2->SetFillStyle(4000);
  canvas_2->Draw();
  

  canvas_1->cd();
  canvas_1->SetLogy();
  canvas_1->SetLogx();
  // crear frame
  double bl = 1.0; // principio del eje X
  double bu = exp_xs_hist->GetBinCenter(16)+exp_xs_hist->GetBinWidth(16)/2.;; // final del eje X
  //double yl = histdat->GetMinimum();  // minimo eje Y
  double yl = 0.00; // minimo eje Y. el 0.5 por logaithmic scale
  double yu = 2.6;  // maximo eje Y

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

  h1->GetXaxis()->SetRangeUser(1,1500);
  h1->GetYaxis()->SetRangeUser(0.0011,0.2);

  gPad->RedrawAxis();
  gPad->Modified(); gPad->Update(); 

  // pintar el frame
  h1->Draw();
  // titulo eje x


  mc_xs_complete->SetMarkerColor(4);
  mc_xs_complete->SetMarkerStyle(1);
  mc_xs_complete->SetFillColor(kRed);
  mc_xs_complete->SetFillStyle(3001);
  mc_xs_complete->Draw("same][2");
  mc_xs_complete->Draw("same][p");


  NLO_xs_complete->SetMarkerColor(4);
  NLO_xs_complete->SetMarkerStyle(20);
  NLO_xs_complete->Draw("same][p");
  //



  
  


  /*

  exp_xs_complete->SetMarkerColor(kGreen);
  exp_xs_complete->SetFillColor(kGreen);
  exp_xs_complete->Draw("same][E2");*/

  /*for (int i=0; i<16;i++){
    his->GetXaxis()->SetBinLabel(i+1,Form("%.0f - %.0f", new_xbins[i],new_xbins[i+1]));
  }*/

  //Plotear
  //myhhdate(mc_xs_hist,"same][",0.0,kRed,0,2.,1,kRed,1000,0); //definido en myAtlasUtils.C
  //myhhdate(superior,"same][",0.0,kRed,0,2.,1,kRed,1001,kRed); //definido en myAtlasUtils.C
  //myhhdate(inferior,"same][",0.0,10,0,2.,1,10,1001,10); //definido en myAtlasUtils.C
  
  //myhhdate(exp_xs_hist,"same][",0.0,kBlack,0,2.,1,kBlack,1000,0); //definido en myAtlasUtils.C


 


  double x1 = 0.35;double x2 = 0.50;
  double yy1 = 0.90; 
  double y2 = yy1 - 0.04;

  TLatex l1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l1.SetNDC();
  l1.SetTextFont(72);
  l1.SetTextColor(1);
  l1.SetTextSize(0.04);
  l1.DrawLatex(x1+0.37,y2-0.0,"ATLAS");
  l1.SetTextFont(42);

  l1.DrawLatex(x1+0.37,y2-0.06,"#sqrt{s} = 13,6 TeV");
  l1.DrawLatex(x1+0.37,y2-0.11,"L = 58.98 fb^{-1}");

  
  myBoxText(x1-0.1,0.33,0.05,kRed,0,3001,"Theoretical NNLO Result",0.04,0.04);
  //myLine(x1-0.1,0.28,0.05,kBlack,1000,"Experimental Result",0.04,0.03);
  myMarkerText(x1-0.1,0.38,4,20,"Theoretical NLO Result",1.0,0.04,0.04);
  myMarkerText(x1-0.1,0.28,1,20,"Ratio NLO/NNLO",1.0,0.04,0.04);

  /*
  myBoxText(x1+0.26,0.83,0.05,kRed,0,3001,"Theoretical NNLO Result",0.04,0.04);
  //myLine(x1-0.1,0.28,0.05,kBlack,1000,"Experimental Result",0.04,0.03);
  myMarkerText(x1+0.26,0.88,4,20,"Theoretical NLO Result",1.0,0.04,0.04);
  myMarkerText(x1+0.26,0.78,1,20,"Ratio NLO/NNLO",1.0,0.04,0.04);*/



  canvas_2->cd();
  canvas_2->SetLogx();

  canvas_2->SetBottomMargin(0.30);
  //
  //double mx = 0.29;
  TH1F *h2 = canvas_2->DrawFrame(bl,0.5,bu,2.0);

  h2->SetYTitle("Ratio NLO/NNLO");
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

  gPad->RedrawAxis();

  h2->Draw();

  TLine l(bl,1.,bu,1.);
  l.SetLineStyle(7);
  l.DrawLine(bl,1,bu,1);
/*
  ratio_xs_complete->SetMarkerColor(4);
  ratio_xs_complete->SetMarkerStyle(1);
  ratio_xs_complete->SetFillColor(kRed);
  ratio_xs_complete->SetFillStyle(3001);
  ratio_xs_complete->Draw("2");
  ratio_xs_complete->Draw("p");

  ratio_NLO_xs_complete->SetMarkerColor(4);
  ratio_NLO_xs_complete->SetMarkerStyle(20);
  ratio_NLO_xs_complete->Draw("same][p");*/

  ratio_NLO_MC_complete->SetMarkerColor(1);
  ratio_NLO_MC_complete->SetMarkerStyle(20);
  ratio_NLO_MC_complete->Draw("same][p");

  
  

  
  



  canvasp->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/TWO_GRAPH_NLO_LOGY_NORMALIZED_EXP_cross_section_22_23.pdf");


  
    
}