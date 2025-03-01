//Para correrlo hacemos root -l macro_figures_ptyy.C

{
  // clase que convierte cualquier tipo de variable a variable string
  ostringstream histo;
  //
  // titulo del eje x
  string xt;
  xt = "p^{T}_{#gamma#gamma} [GeV]"; 
  // titulo eje y
  string yt;
  yt = "diff_cross_section";

  TFile *mc;  // root de todos los MC del que voy a leer
  TH1D *histmc_all; // histograma de todos los MC que voy a mostrar
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/mc23a_signal_histos.root"; // nombre del root de los MC
  mc = new TFile(histo.str().c_str(),"read"); // abrir fichero MC

  histo.str("");  histo << "h0yypt"; // nombre del histograma que quiero leer
  histmc_all = (TH1D*)mc->Get(histo.str().c_str()); // leer histograma MC


  //Originally bins were of 6.0 GeV



  int len = 30;
  Double_t xbins[30] = {0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0,36.0,39.0,42.0,45.0,48.0,51.0,54.0,57.0,60.0,69.0,78.0,87.0,96.0,132.0,186.0,258.0,618.0,1500.0};
  
  /*int n = 0; double content = 0.0; xbins[0] = 0.0;

  for (int i=1;i<=len;i++){
    if (n<25){
      xbins[i] = xbins[i-1] + 3.0;
      n ++;
    }
    else if (n<82){
      xbins[i] = xbins[i-1] + 6.0;
      n++;
    }
    else if  (n<273){
      xbins[i] = xbins[i-1] + 9.0;
      n++;
    }
    else{
      xbins[i] = xbins[i] +12.0;
      n++;
    }
  } //No es un mal binado este.*/




  histmc_all->Rebin(29,"hist_rebinned",xbins);
  /*TH1D* hist_rebinned;
  hist_rebinned = (TH1D*)histmc_all->Clone(histo.str().c_str());*/
  
  TH1D *histratio;
  histratio = (TH1D*)hist_rebinned->Clone(histo.str().c_str());

  /*TFile *data;  // root de todos los MC del que voy a leer
  TH1D *histdata; // histograma de todos los MC que voy a mostrar
  histo.str("");   histo << "/home/mjlarroy/TFG_TUT/data22_histos.root"; // nombre del root de los MC
  data = new TFile(histo.str().c_str(),"read"); // abrir fichero datos

  
  histo.str("");  histo << "h0yypt"; // nombre del histograma que quiero leer
  histdata = (TH1D*)data->Get(histo.str().c_str()); // leer histograma MC*/


  int ni = 1; int nb = hist_rebinned->GetNbinsX(); // numero de bines en el eje X
  cout<<"n bins = "<<nb<<endl;

  double nmc = 0.0;
  double ndat = 0.0;
  double value = 0.0;
  double error = 0.0;
  double anchura = 0.0;

  for(int j=ni; j<=nb; j++){

    anchura = (hist_rebinned->GetXaxis()->GetBinWidth(j))/3.0; //6.0 xq es la anchura más pequeña

    value = hist_rebinned->GetBinContent(j); //Valor en el bin j
    error = sqrt(value);

    hist_rebinned->SetBinContent(j,value/anchura);

    //double dat = histdata->GetBinContent(j);

    //ndat += dat;

    nmc += value; //Actualizar el número de datos

    hist_rebinned->SetBinError(j,error/anchura);

    histratio->SetBinContent(j,error/value);

    
      cout<<histratio->GetBinContent(j)<<endl;
    

  }

  // normalizacion del histograma de MC al histograma de los datos.
  /*
  for(int j=ni;j<=nb;j++)
    {
      double xmc = histmc_all->GetBinContent(j); // signal del MC
      xmc = xmc*ndat/nmc; // aqui se realiza la normalizacion del MC al número de datos
      histmc_all->SetBinContent(j,xmc);
      histmc_all->SetBinError(j,0.); // para que salga solido el histograma (cosa técnica de root)
    }*/
  //

  cout<<"n events = "<<nmc<<endl;


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

  TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.0,0.29,0.9,0.99);
  TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0.0,0.10,0.9,0.41);

  // pintar pads
  canvas_1->SetFillStyle(4000); //transparentes
  canvas_1->Draw();
  canvas_2->SetFillStyle(4000);
  canvas_2->Draw();

  canvas_1->cd();
  //canvas_1->SetLogy();
  canvas_1->SetLogx();
  // crear frame
  double bl = hist_rebinned->GetBinCenter(ni)-hist_rebinned->GetBinWidth(ni)/2.; // principio del eje X
  double bu = hist_rebinned->GetBinCenter(nb)+hist_rebinned->GetBinWidth(nb)/2.+100.0; // final del eje X
  //double yl = histdat->GetMinimum();  // minimo eje Y
  double yl = 0.0015; // minimo eje Y. el 0.5 por logaithmic scale
  double yu = hist_rebinned->GetMaximum();  // maximo eje Y
  yu = yu*1.3; // maximo eje Y

  TH1F *h1 = canvas_1->DrawFrame(bl,yl,bu,yu); // crea el frame

  
  //
  TLine l(bl,1.,bu,1.);
  l.SetLineWidth(1.);
  l.SetLineStyle(2); //Crea la línea que vemos luego debajo
  //
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
  // titulo eje x


  h1->SetTitle("pT distribution for the signal events (MC)"); //Ask tmb

  //
  // pintar el frame
  h1->Draw();

  //Plotear
  myhhdate(hist_rebinned,"same][",0.0,132,0,2.,1,132,1000,132); //definido en myAtlasUtils.C

  double x1 = 0.35;double x2 = 0.50;
  double yy1 = 0.90; 
  double y2 = yy1 - 0.04;

  TLatex l1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l1.SetNDC();
  l1.SetTextFont(72);
  l1.SetTextColor(1);
  l1.SetTextSize(0.045);
  l1.DrawLatex(x1+0.35,y2-0.25,"ATLAS");
  l1.SetTextFont(42);

  l1.DrawLatex(x1+0.35,y2-0.31,"#sqrt{s} = 13,6 TeV");
  l1.DrawLatex(x1+0.35,y2-0.36,"L = 30 fb^{-1}");


  gPad->RedrawAxis();

  //Segundo Pad

  canvas_2->cd();
  canvas_2->SetLogx();

  canvas_2->SetBottomMargin(0.30);
  //
  //double mx = 0.29;
  double mx = 0.09;
  TH1F *h2 = canvas_2->DrawFrame(bl,0.0,bu,1.45);
  h2->SetYTitle("Relative error");
  //
  histo.str(""); histo << xt;
  h2->SetXTitle(histo.str().c_str());
  //
  h2->GetXaxis()->SetTitleSize(0.11);
  h2->GetXaxis()->SetTitleOffset(1.1);
  //
  h2->GetXaxis()->SetLabelSize(0.09);
  h2->GetXaxis()->SetNdivisions(505);
  h2->GetXaxis()->SetMoreLogLabels(1);
  h2->GetXaxis()->SetNoExponent(1);
  //
  h2->GetYaxis()->SetTitleSize(0.10);
  h2->GetYaxis()->SetTitleOffset(0.7);
  h2->GetYaxis()->SetLabelSize(0.09);
  //
  TAxis *axis_x2 = h2->GetXaxis();
  axis_x2->SetNdivisions(505);
  axis_x2->SetMoreLogLabels(1);
  axis_x2->SetNoExponent(1);
  axis_x2->SetTickLength(0.08);
  //
  TAxis *axis_y2 = h2->GetYaxis();
  axis_y2->SetNdivisions(505);
  //
  h2->Draw();


  myhhdate(histratio,"esamex0",1.0,1,20,2.,1,1,1,1); //1.0,1,20,2.,1,1,1,1
  l.DrawLine(bl,0.2,bu,0.2);

  

  gPad->RedrawAxis();
  cout << " end " << endl;
}