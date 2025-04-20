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
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/Experimental_CrossSection_22_23.root"; // nombre del root de los MC
  exp_xs = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "h0yypT"; // nombre del histograma que quiero leer
  exp_xs_hist = (TH1D*)exp_xs->Get(histo.str().c_str()); // leer histograma MC

  TFile *errors;  // root de todos los MC del que voy a leer
  TH1D *superior; // histograma de todos los MC que voy a mostrar
  TH1D *inferior;
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/XS_errors.root"; // nombre del root de los MC
  errors = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "h_pTTruth_superior"; // nombre del histograma que quiero leer
  superior = (TH1D*)errors->Get(histo.str().c_str()); // leer histograma MC

  histo.str("");  histo << "h_pTTruth_inferior"; // nombre del histograma que quiero leer
  inferior = (TH1D*)errors->Get(histo.str().c_str()); // leer histograma MC

  int ni = 1; int nb = superior->GetNbinsX(); // numero de bines en el eje X

  TFile *reco;  // root de todos los MC del que voy a leer
    TH1D *hist_reco; // histograma de todos los MC que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/reco.root"; // nombre del root de los MC
    reco = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h"; // nombre del histograma que quiero leer
    hist_reco = (TH1D*)reco->Get(histo.str().c_str()); // leer histograma MC


  //Como los errores son asimétricos usaremos un TGraphAsymmErrors

  Double_t sys_error[16] = {2.32365,13.6853,12.5925,6.81755,3.15433,2.8733,13.3374,11.9845,12.4012,11.7149,2.04462,2.08365,0.643965,0.372638,0.375956,0.0266472};

  Double_t x_positions[16] = {};
  Double_t y_positions[16] = {};
  Double_t y_positions_exp[16] = {};

  Double_t x_errors[16] = {};
  Double_t y_sup_errors[16] = {};
  Double_t y_inf_errors[16] = {};
  Double_t rel_error[16] = {};
  Double_t val = 0.0;
  Double_t exp_stat_error[16] = {};
  Double_t exp_total_error[16] = {};

  Double_t L = 58.98;

    for(int i=0;i<nb;i++){

      val = exp_xs_hist->GetBinContent(i+1);
      exp_xs_hist->SetBinContent(i+1,val*1.0/hist_reco->GetBinContent(i+1));

      x_positions[i] = mc_xs_hist->GetBinCenter(i+1);
      y_positions[i] = mc_xs_hist->GetBinContent(i+1);
      y_positions_exp[i] = exp_xs_hist->GetBinContent(i+1);

      x_errors[i] = mc_xs_hist->GetBinWidth(i+1)/2.0;

      y_sup_errors[i] = superior->GetBinContent(i+1) - mc_xs_hist->GetBinContent(i+1);
      y_inf_errors[i] = -1.0*(inferior->GetBinContent(i+1) - mc_xs_hist->GetBinContent(i+1)); //It is defined as positive for the TGraph
    
      rel_error[i] = (superior->GetBinContent(i+1)-inferior->GetBinContent(i+1))/mc_xs_hist->GetBinContent(i+1);


      exp_stat_error[i] = exp_xs_hist->GetBinError(i+1);

      exp_total_error[i] = sqrt(exp_stat_error[i]*exp_stat_error[i] + (sys_error[i]*sys_error[i])/((L*exp_xs_hist->GetBinWidth(i+1))*(L*exp_xs_hist->GetBinWidth(i+1))));

    }

    TGraphAsymmErrors *mc_xs_complete = new TGraphAsymmErrors(16,x_positions,y_positions,x_errors,x_errors,y_inf_errors,y_sup_errors);


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

  TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.0,0.01,0.9,0.99);
  //TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0.0,0.10,0.9,0.41);

  // pintar pads
  canvas_1->SetFillStyle(4000); //transparentes
  canvas_1->Draw();
  

  canvas_1->cd();
  //canvas_1->SetLogy();
  canvas_1->SetLogx();
  // crear frame
  double bl = 1.0; // principio del eje X
  double bu = exp_xs_hist->GetBinCenter(16)-exp_xs_hist->GetBinWidth(16)/2.;; // final del eje X
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
  h1->GetYaxis()->SetTitleOffset(1.4);
  h1->GetYaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetNdivisions(606);

  h1->GetXaxis()->SetNdivisions(505);
  h1->GetXaxis()->SetMoreLogLabels(1);
  h1->GetXaxis()->SetNoExponent(1);
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
  // titulo eje x



  //
  // pintar el frame
  h1->Draw();

  mc_xs_complete->SetMarkerColor(kRed);
  mc_xs_complete->SetFillColor(kRed);
  mc_xs_complete->Draw("same][E2");

  /*

  exp_xs_complete->SetMarkerColor(kGreen);
  exp_xs_complete->SetFillColor(kGreen);
  exp_xs_complete->Draw("same][E2");*/

  //Plotear
  //myhhdate(mc_xs_hist,"same][",0.0,kRed,0,2.,1,kRed,1000,0); //definido en myAtlasUtils.C
  //myhhdate(superior,"same][",0.0,kRed,0,2.,1,kRed,1001,kRed); //definido en myAtlasUtils.C
  //myhhdate(inferior,"same][",0.0,10,0,2.,1,10,1001,10); //definido en myAtlasUtils.C
  myhhdate(exp_xs_hist,"same][",0.0,kBlue,0,2.,1,kBlue,1000,0); //definido en myAtlasUtils.C


 


  double x1 = 0.35;double x2 = 0.50;
  double yy1 = 0.90; 
  double y2 = yy1 - 0.04;

  TLatex l1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l1.SetNDC();
  l1.SetTextFont(72);
  l1.SetTextColor(1);
  l1.SetTextSize(0.04);
  l1.DrawLatex(x1+0.30,y2-0.20,"ATLAS");
  l1.SetTextFont(42);

  l1.DrawLatex(x1+0.30,y2-0.26,"#sqrt{s} = 13,6 TeV");
  l1.DrawLatex(x1+0.30,y2-0.31,"L = 58.98 fb^{-1}");

  myLine(x1+0.36,0.85,0.05,kRed,20,"MC",0.04,0.03);
  myLine(x1+0.36,0.8,0.05,kBlue,1000,"Exp",0.04,0.03);

  canvasp->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/EXP_cross_section_22_23.pdf");


  gPad->RedrawAxis();
    
}