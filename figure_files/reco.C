//Para correrlo hacemos root -l reco.C

{
    // clase que convierte cualquier tipo de variable a variable string
    ostringstream histo;
    //
    // titulo del eje x
    string xt;
    xt = "p^{T}_{#gamma#gamma} [GeV]"; 
    // titulo eje y
    string yt;
    yt = "Reco Efficiency";


    //Detector cross section measurements
    TFile *mc;
    TH1D *hist;
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/MC_ptyy_cross_section.root";
    mc = new TFile(histo.str().c_str(),"read");
  
    histo.str("");  histo << "hist_rebinned";
    hist = (TH1D*)mc->Get(histo.str().c_str());

    //Real cross section (reco)
    TFile *mc_fiducial;
    TH1D *hist_fiducial;
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/MC_ptyy_fiducial_cross_section.root";
    mc_fiducial = new TFile(histo.str().c_str(),"read");
    
    histo.str("");  histo << "hist_rebinned";
    hist_fiducial = (TH1D*)mc_fiducial->Get(histo.str().c_str());


    TH1D *histratio;
    histratio = (TH1D*)hist_fiducial->Clone(histo.str().c_str());
    histratio->Reset("ICESM");
    histratio->SetName("histratio");

    int ni = 1; int nb = hist_rebinned->GetNbinsX();

    double var = 0.0;
    double var_fiducial = 0.0;
    double var_error = 0.0;
    double var_fiducial_error = 0.0;

    Double_t xbins[30] = {0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0,36.0,39.0,42.0,45.0,48.0,51.0,54.0,57.0,60.0,69.0,78.0,87.0,96.0,132.0,186.0,258.0,618.0,1500.0};
    Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};

    TH1D *h = new TH1D("h","pT_migrations",16,0.0,16.0);

  double media = 0.0;
  double errors = 0.0;

    for(int j=ni; j<=nb; j++){

        var = hist->GetBinContent(j);
        var_fiducial = hist_fiducial->GetBinContent(j);

        var_error = hist->GetBinError(j);
        var_fiducial_error = hist_fiducial->GetBinError(j);

        media += var/var_fiducial;

        histratio->SetBinContent(j,var/var_fiducial);

        histratio->SetBinError(j, sqrt(1.0/(var_fiducial*var_fiducial) * ((1-2*var/var_fiducial)*var_error*var_error + var*var/(var_fiducial*var_fiducial)*var_fiducial_error*var_fiducial_error)));//ESTE ERROR HAY QUE CAMBIARLO

        errors += histratio->GetBinError(j);

        h->SetBinContent(j,histratio->GetBinContent(j));
        h->SetBinError(j,histratio->GetBinError(j));
    }
    media = media/16.0; //SALE 0,74
    errors = sqrt(errors)/16.0; //SALE 0,07
    std::cout<<errors<<std::endl;



    //REPRESENTATIONS

    // anchura y altura del canvas
  double w = 650;
  double he = 650;
  double xx = 0.; // posicion x de la esquina superior izquierda del canvas
  double yy = 0.; // posicion y de la esquina superior izquierda del canvas
  // crea canvas
  TCanvas *canvasp = new TCanvas("canvasp","canvasp",xx,yy,w,he);
  canvasp->SetWindowSize(w+(w-canvasp->GetWw()),he+(he-canvasp->GetWh()));
  canvasp->SetFillStyle(4000); // para hacer transparente
  canvasp->cd();


  TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.0,0.0,0.9,0.99);
  canvas_1->SetFillStyle(4000); //transparentes
  canvas_1->Draw();
  canvas_1->cd();

  gPad->SetBottomMargin(0.3);
  TH1F *h1 = canvas_1->DrawFrame(0,0,16,1);


    //h->GetXaxis()->SetRangeUser(0,1500);//Rango de eje X
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.04);

    h->GetXaxis()->SetTitleOffset(3.5);
    h->GetYaxis()->SetRangeUser(0,1);//Rango de eje X

    h->GetXaxis()->SetTitle("p^{T}_{#gamma#gamma truth} [GeV]");//Título del hist
    h->GetYaxis()->SetTitle("Unfolding Coefficient");//Sale lejos xd

  //
  h1->Draw();
  h->Draw("EX0");
  gPad->RedrawAxis();




  for (int i=0; i<16;i++){
    h->GetXaxis()->SetBinLabel(i+1,Form("%.0f - %.0f", new_xbins[i],new_xbins[i+1]));
    h->GetXaxis()->ChangeLabel(i+1,90,-1,-1,-1,-1);
  }

  h->GetXaxis()->LabelsOption("v");

  canvas_1->Update();

  canvasp->Update();

  canvasp->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/reco_efficiency.pdf");

  TFile* myfile = new TFile("/home/mjlarroy/TFG_TUT/Root_Results/reco.root", "RECREATE"); //Recreate significa que si no está lo cree y si está que lo reescriba.

  h->Write();

  myfile->Close();



  gPad->RedrawAxis();



}