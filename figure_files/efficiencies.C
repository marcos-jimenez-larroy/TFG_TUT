{


    TH1F *h = new TH1F("h","Production Mode Efficiency",6,0.,6.);
    TH1F *h2 = new TH1F("h2","Production Mode Efficiency",6,0.,6.);


    vector<TString> ProdMode = {
        "ggH",
        "VBF",
        "ttH",
        "WmH",
        "WpH",
        "ZH"}; 

        Double_t Efficiencies[6] = {36.88,36.70,37.54,33.63,29.92,32.32};
        Double_t ATLAS_Efficiencies[6] = {35.47,35.83,36.35,32.98,29.66,31.53};


        for(int i=0;i<6;i++){
            h->SetBinContent(i+1,Efficiencies[i]);
            h2->SetBinContent(i+1,ATLAS_Efficiencies[i]);

        }


    double w = 650;
    double he = 650;
    double xx = 0.; // posicion x de la esquina superior izquierda del canvas
    double yy = 0.; // posicion y de la esquina superior izquierda del canvas
    // crea canvas
    TCanvas *canvasp = new TCanvas("canvasp","canvasp",xx,yy,w,he);
    canvasp->SetGrid();
    canvasp->SetWindowSize(w+(w-canvasp->GetWw()),he+(he-canvasp->GetWh()));
    canvasp->SetFillStyle(4000); // para hacer transparente
    canvasp->cd();

    h->GetXaxis()->SetRangeUser(0,6);//Rango de eje X
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->SetTitle("Production Mode Efficiency");

    h->SetFillColor(38);
    h->SetBarWidth(0.40);
    h->SetBarOffset(0.1);


    h->GetXaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetRangeUser(0,50);//Rango de eje X

    h->GetXaxis()->SetTitle("Production Mode");//Título del hist
    h->GetYaxis()->SetTitle("Efficiency (%)");//Sale lejos xd

    h->Draw("BARS");

    


        for (int i=0; i<6;i++){
            h->GetXaxis()->SetBinLabel(i+1,ProdMode[i]);
          }

    TLatex *latex = new TLatex(0.5, 0.97, "Production Mode Efficiency");
    latex->SetNDC();
    latex->SetTextAlign(22); // Center text
    latex->SetTextSize(0.05); // Size of the title text
    latex->SetTextFont(42);   // Font (Helvetica)
    latex->Draw();

    h2->SetFillColor(4);
    h2->SetBarWidth(0.40);
    h2->SetBarOffset(0.5);
    h2->Draw("BARS SAME");

    myBoxText(0.7,0.9,0.05,38,38,1000,"Used MC",0.03,0.03);
    myBoxText(0.7,0.85,0.05,4,4,1000,"ATLAS Result",0.03,0.03);


    // Update the canvas
    canvasp->Update();
    canvasp->Print("/home/mjlarroy/TFG_TUT/figures/pdf_graphs/Production_Modes_Efficiencies.pdf");

}