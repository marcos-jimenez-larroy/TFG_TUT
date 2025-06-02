{

  
    // clase que convierte cualquier tipo de variable a variable string
    ostringstream histo;
    //
    // titulo del eje x
    string xt;
    xt = "p^{T}_{#gamma#gamma reco} (GeV)"; //En el de verdad habría que sumar también la aportación de los muones ojo
    //
    // titulo eje y
    string yt;
    yt = "p^{T}_{#gamma#gamma truth} (GeV)";

    TFile *dat;  // root de datos del que voy a leer
    TH2D *hist; // histograma de datos que voy a mostrar
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/mc23a_ptyy_migration.root"; // nombre del root de los datos
    dat = new TFile(histo.str().c_str(),"read"); // abrir fichero datos
  
    histo.str("");  histo << "ptyy_migration"; // nombre del histograma que quiero leer
    hist = (TH2D*)dat->Get(histo.str().c_str()); // leer histograma datos
    hist->SetName("p^{T}_{#gamma#gamma} migrations");

    TFile *fid_bins;
    fid_bins = new TFile("/home/mjlarroy/TFG_TUT/Root_Results/mc23a_signal_fiducial_histos.root","READ");

    TH1D *fid_ptyy;
    fid_ptyy = (TH1D*)fid_bins->Get("h0yypt");

    fid_ptyy->Rebin(16,"fid_ptyy_rebinned",new_xbins);

    //REPRESENTATIONS


    Double_t xbins[30] = {0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0,36.0,39.0,42.0,45.0,48.0,51.0,54.0,57.0,60.0,69.0,78.0,87.0,96.0,132.0,186.0,258.0,618.0,1500.0};
    Double_t new_xbins [17] = {0.0,6.0,12.0,18.0,24.0,30.0,36.0,45.0,54.0,66.0,90.0,114.0,138.0,192.0,312.0,624.0,1500.0};

    TH2D *h = new TH2D("h","p^{T}_{#gamma#gamma} migrations",16,0.0,16.0,16,0.0,16.0);

    

    for(int i=1;i<=16;i++){

      //double width_i = xbins[i]-xbins[i-1];
      //double width_j = 0.0;
      double sum = fid_ptyy_rebinned->GetBinContent(i);

        for(int j=1;j<=16;j++){
          //width_j = xbins[i]-xbins[i-1];
        h->SetBinContent(i,j,int(hist->GetBinContent(i,j)/sum*100));
        }
    }


    // anchura y altura del canvas
  double w = 700;
  double ha = 700;
  double xx = 0.; // posicion x de la esquina superior izquierda del canvas
  double yy = 0.; // posicion y de la esquina superior izquierda del canvas
  // crea canvas
  TCanvas *canvasp = new TCanvas("canvasp","canvasp",xx,yy,w,ha);
  gStyle->SetOptTitle(0);
  canvasp->SetWindowSize(w+(w-canvasp->GetWw()),ha+(ha-canvasp->GetWh()));
  canvasp->SetFillStyle(4000); // para hacer transparente
  //auto *axis = new TGaxis(0,1.0,0,1.0,0.0,100,530,"S");
  

  /*for(int i=1;i<=29;i++){
    axis->ChangeLabel(i,-1,-1,-1,3,-1,xbins[i]);
  }*/

  canvasp->cd();


  canvasp->SetGrid();
  //canvasp->SetLogx();
  //canvasp->SetLogy();

    //h->GetXaxis()->SetRangeUser(0,30);//Rango de eje X
    //h->GetYaxis()->SetRangeUser(0,30);//Rango de eje X
    gPad->SetBottomMargin(0.3);
    gPad->SetLeftMargin(0.3);

    h->GetXaxis()->SetLabelSize(0.06);    
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);


    h->GetXaxis()->SetTickLength(0.03);

    h->GetXaxis()->SetTitle("p^{T}_{#gamma#gamma reco} [GeV]");//Título del hist
    h->GetYaxis()->SetTitle("p^{T}_{#gamma#gamma truth} [GeV]");//Sale lejos xd
    h->GetYaxis()->SetTitleOffset(3.4);
    h->GetXaxis()->SetTitleOffset(3.5);


    h->Draw("TEXT COLZ");

    
    


 

  for (int i=0; i<16;i++){
    h->GetXaxis()->SetBinLabel(i+1,Form("%.0f - %.0f", new_xbins[i],new_xbins[i+1]));
    h->GetYaxis()->SetBinLabel(i+1,Form("%.0f - %.0f", new_xbins[i],new_xbins[i+1]));
  }

  h->GetXaxis()->LabelsOption("v");

  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextSize(0.04);
    l.SetTextColor(1);
    l.DrawLatex(0.2,0.96,"ATLAS");
    l.SetTextFont(42);

    l.DrawLatex(0.37,0.96,"#sqrt{s} = 13,6 TeV");


  canvasp->Update();
  canvasp->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/ptyy_migration.pdf");


 

}