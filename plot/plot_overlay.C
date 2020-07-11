void normalizeHist(TH1F*Hist){
  int     nbins    = Hist->GetNbinsX();
  double  nentries = Hist->GetEntries();
  printf("nentries = %2.3f\n", nentries);

  for (int i=0; i<nbins; ++i){
    double content = Hist->GetBinContent(i+1)/nentries;
    if ( content > 0) Hist->SetBinContent(i+1, content);
  }
}


//-----------------------------------------------------------------------------
void plot_trkseed_tprdem_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_tkrseed_tprdem_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("trkseed_115/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("trkseed_120/%s", Variable));

  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);
  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");


  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}


//-----------------------------------------------------------------------------
void plot_helix_tprdem_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_helix_tprdem_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("helix_0/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("helix_10/%s", Variable));
  
  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");

  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}


//-----------------------------------------------------------------------------
void plot_trkseed_tprdep_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_trkseed_tprdep_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("trkseed_315/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("trkseed_320/%s", Variable));

  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);
  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");


  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}


//-----------------------------------------------------------------------------
void plot_helix_tprdep_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_helix_tprdep_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("helix_300/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("helix_310/%s", Variable));
  
  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");

  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//                                    CalPatRec
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void plot_trkseed_cprdem_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_trkseed_cprdem_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("trkseed_15/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("trkseed_20/%s", Variable));

  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);
  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");


  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}


//-----------------------------------------------------------------------------
void plot_helix_cprdem_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_helix_cprdem_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("helix_100/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("helix_110/%s", Variable));
  
  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");

  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}


//-----------------------------------------------------------------------------
void plot_trkseed_cprdep_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_trkseed_cprdep_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("trkseed_215/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("trkseed_220/%s", Variable));

  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);
  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");


  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}


//-----------------------------------------------------------------------------
void plot_helix_cprdep_overlay( const char* Fn0, const char* Fn1, const char* Variable) {

  TCanvas* c1 = new TCanvas(Form("c_helix_cprdep_%s", Variable),Form("c_%s", Variable),900,600);
  
  c1->cd();

  TH1F* he0_trig    = gh1(Fn0   ,"TriggerAna", Form("helix_200/%s", Variable));
  TH1F* he1_trig    = gh1(Fn1   ,"TriggerAna", Form("helix_210/%s", Variable));
  
  he1_trig->SetFillColor(623);
  he1_trig->SetFillStyle(3003);
  he1_trig->SetLineColor(623);

  //  he0_trig->GetYaxis()->SetRangeUser(0.8, 1.0);

  TLegend *leg = new TLegend(0.5957, 0.7474, 0.8363, 0.8676, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(he0_trig, "Background only", "L");
  leg->AddEntry(he1_trig, "CE with final track", "L");

  normalizeHist(he0_trig);
  normalizeHist(he1_trig);

  double    he0_max = he0_trig->GetMaximum();
  double    he1_max = he1_trig->GetMaximum();
  //  he0_trig->GetXaxis()->SetTitle("E_{cluster} [MeV]");
  if (he0_max > he1_max){
    he0_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he0_trig->Draw();
    he1_trig->Draw("sames"); he0_trig->Draw("same");
  }else {
    he1_trig->GetYaxis()->SetTitle(Form("Entries / (%2.3f)", he0_trig->GetBinWidth(1)));
    he1_trig->Draw();
    he0_trig->Draw("sames"); he1_trig->Draw("same");  
  }
  leg->Draw();

  c1->Update();
  TPaveStats *stats1 = (TPaveStats*)he0_trig ->GetListOfFunctions()->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats*)he1_trig ->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue); 
  stats2->SetTextColor(kRed); 
  stats2->SetX1NDC(0.53); stats2->SetX2NDC(0.70);
  stats2->SetY1NDC(0.46); stats2->SetY2NDC(0.70);
  stats1->SetX1NDC(0.72); stats1->SetX2NDC(0.89);
  stats1->SetY1NDC(0.46); stats1->SetY2NDC(0.70);
  c1->Modified();
}
