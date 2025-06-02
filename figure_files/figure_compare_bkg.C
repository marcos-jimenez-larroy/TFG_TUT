//Para correrlo hacemos root -l figure_compare_bkg.C

{
    // clase que convierte cualquier tipo de variable a variable string
    ostringstream histo;
    //
    // titulo del eje x
    string xt;
    xt = "m_{#gamma#gamma} [GeV]"; 
    // titulo eje y
    string yt;
    yt = "Events";
  
    TFile *mc;  // root del BKG que voy a leer
    TH1D *histmc_bkg; // histograma del BKG que voy a leer
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/mc23a_bkg_histos.root"; // nombre del root del MC
    mc = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h0yymass"; // nombre del histograma que quiero leer
    histmc_bkg = (TH1D*)mc->Get(histo.str().c_str()); // leer histograma MC
  
  
    // bins were of 1 GeV

    TFile *data;  // root de datos que voy a leer
    TH1D *histdata; // histograma de datos que voy a leer
    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/Root_Results/data22_23_histos.root"; // nombre del root de los datos
    data = new TFile(histo.str().c_str(),"read"); // abrir fichero MC
  
    histo.str("");  histo << "h0yymass"; // nombre del histograma que quiero leer
    histdata = (TH1D*)data->Get(histo.str().c_str()); // leer histograma datos
  
    //COMPARE THE MC FOR BACKGROUND AND THE DATA IN THE SIDEBAND REGION TO SEE HOW RELEVANT THE MISSING PHOTON-JET CONTRIBUTION IS
  
    //We have already seen that the photon-jet events are relevant. We will fit the data and reweight the MC in order to 'count' this effect
    
    TF1 *ajuste = fitExclude(histdata,105,160); //3rd degree polinomial
  
  

    
    TH1D *histratio;
    histratio = (TH1D*)histdata->Clone(histo.str().c_str());
    histratio->SetName("histratio");

    TH1D *histratio_ajuste;
    histratio_ajuste = (TH1D*)histdata->Clone(histo.str().c_str());
    histratio_ajuste->SetName("histratio_ajuste");

    TH1D *histmc_bkg_fit;
    histmc_bkg_fit = (TH1D*)histmc_bkg->Clone(histo.str().c_str());
    histmc_bkg_fit->SetName("histmc_bkg_fit");
  

  
  
    int ni = 1; int nb = histdata->GetNbinsX(); // numero de bines en el eje X
    cout<<"n bins = "<<nb<<endl;
  
    double nmc = 0.0;
    double ndat = 0.0;
    double value_mc = 0.0;
    double value_data = 0.0;
    double value_fit = 0.0;
    double error = 0.0;

    double L1 = 27.58;
    double Lt = 31.40+27.58;
    for(int j=ni; j<=nb; j++){
  
  
      value_mc = (histmc_bkg->GetBinContent(j))*Lt/L1; //Valor del MC en el bin j 
      value_data = histdata->GetBinContent(j); //Valor del dato en el bin j
      value_fit = ajuste->Eval(105+j); //El histograma de 105 a 160 en 55 bines.
      
      histmc_bkg_fit->SetBinContent(j,value_fit);
  
      nmc += value_mc; //Actualizar el número de eventos en el MC
      ndat += value_data; //Actualizar el número de eventos en los datos (ojo que estamos contando el número de eventos de señal pero en comparación son pocos)
      
  
    }
  
    // normalizacion del histograma de MC al histograma de los datos.
    
    for(int j=ni;j<=nb;j++)
      {
        double xmc = (histmc_bkg->GetBinContent(j))*Lt/L1; // signal del MC
        double xdat = histdata->GetBinContent(j);
        xmc = xmc*ndat/nmc; // aqui se realiza la normalizacion del MC al número de datos
        histmc_bkg->SetBinContent(j,xmc);
        histmc_bkg->SetBinError(j,0.); // para que salga solido el histograma (cosa técnica de root)

        error = sqrt(xdat)/xmc + xdat/(xmc*xmc) * sqrt(xmc);

        histratio->SetBinContent(j,xdat/xmc);
        histratio->SetBinError(j,error);

        value_fit = ajuste->Eval(105+j);

        histratio_ajuste->SetBinContent(j,value_fit/xmc);
        histratio_ajuste->SetBinError(j,value_fit/(xmc*xmc)*sqrt(xmc));
      }
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
    //canvas_1->SetLogx();
    // crear frame
    double bl = histdata->GetBinCenter(ni)-histdata->GetBinWidth(ni)/2.; // principio del eje X
    double bu = histdata->GetBinCenter(nb)+histdata->GetBinWidth(nb)/2.; // final del eje X
    //double yl = histdat->GetMinimum();  // minimo eje Y
    double yl = 0.0015; // minimo eje Y. el 0.5 por logaithmic scale
    double yu = histdata->GetMaximum();  // maximo eje Y
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
    h1->GetYaxis()->SetTitleOffset(1.55);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetYaxis()->SetNdivisions(606);
    h1->GetYaxis()->SetMaxDigits(4);
    //
    // settings eje x
    histo.str(""); histo << xt;
    h1->SetXTitle(histo.str().c_str());
    //
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleOffset(1.1);
    //
    h1->GetXaxis()->SetLabelSize(0.01);
    h1->GetXaxis()->SetNdivisions(505);
    // titulo eje x
  
  
    h1->SetTitle("m_{#gamma#gamma}"); //Ask tmb
  
    //
    // pintar el frame
    h1->Draw();
  
    //Plotear
    myhhdate(histdata,"same][",0.0,132,0,2.,1,132,1000,0); //definido en myAtlasUtils.C
    myhhdate(histmc_bkg,"same][",0.0,2,0,2.,1,2,1000,0); //definido en myAtlasUtils.C
    myhhdate(histmc_bkg_fit,"same][",0.0,3,0,2.,1,3,1000,0); //definido en myAtlasUtils.C



  
    double x1 = 0.35;double x2 = 0.50;
    double yy1 = 0.90; 
    double y2 = yy1 - 0.045;
    double y3 = y2 - 0.045;
  
    myBoxText(x1,yy1,0.05,0,132,1000,"Data",0.04,0.03);
    myBoxText(x1,y2,0.05,0,2,1000,"MC #gamma#gamma Background",0.04,0.03);
    myBoxText(x1,y3,0.05,0,3,1000,"MC reweighted Background",0.04,0.03);
    


    TLatex l1; //l.SetTextAlign(12); l.SetTextSize(tsize); 
    l1.SetNDC();
    l1.SetTextFont(72);
    l1.SetTextColor(1);
    l1.SetTextSize(0.045);
    l1.DrawLatex(x1+0.35,y2-0.25,"ATLAS");
    l1.SetTextFont(42);
  
    l1.DrawLatex(x1+0.35,y2-0.31,"#sqrt{s} = 13,6 TeV");
    l1.DrawLatex(x1+0.35,y2-0.36,"L = 58,98 fb^{-1}");
  
  
    gPad->RedrawAxis();
  
    //Segundo Pad
  
    canvas_2->cd();
    //canvas_2->SetLogx();
  
    canvas_2->SetBottomMargin(0.30);
    //
    //double mx = 0.29;
    double mx = 0.09;
    TH1F *h2 = canvas_2->DrawFrame(bl,0.90,bu,1.05);
    h2->SetYTitle("Ratio Over MC#gamma#gamma");
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
  
  
    myhhdate(histratio_ajuste,"esamex0",1.0,1,20,2.,1,1,1,1); //1.0,1,20,2.,1,1,1,1
    //myhhdate(histratio,"esamex0",1.0,1,20,2.,1,1,1,1); //1.0,1,20,2.,1,1,1,1
    l.SetLineStyle(7);
    l.DrawLine(bl,1,bu,1);


    //myLine(x1+0.4,y2,0.05,1,1000,"Fit/MC#gamma#gamma",0.08,0.03);
    //myLine(x1+0.4,y2-0.08,0.05,1,1000,"Data/MC#gamma#gamma",0.08,0.03);
  
    
  
    gPad->RedrawAxis();

    histo.str("");   histo << "/home/mjlarroy/TFG_TUT/figures/pdf_graphs/MC_background_22_23.pdf";
    myepsfile(canvasp,histo.str().c_str(),"pdf");

    TFile* myfile = new TFile("/home/mjlarroy/TFG_TUT/Root_Results/MC_reweighted_background_histos_22_23.root","RECREATE");

    histmc_bkg_fit->Write();
    ajuste->Write();
    histratio_ajuste->Write();


    cout << " end " << endl;
  }