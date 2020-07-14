void    plot_stntuple(const char*FnStn){

  TString       names[2]  = {"TrkPatRec","CalPatRec"};//{"","MuMinus","MuPlus","IPAMu","Proton"};
  EColor        colors[5] = {kBlue,kGreen,kMagenta,kRed,kOrange};
  
  TH1F*         hp[5] = {NULL};
  TH1F*         hdp[5] = {NULL};
  TH1F*         hpt[5] = {NULL};
  TH1F*         hpzMC[5] = {NULL};
  TH1F*         hchi2xy[5] = {NULL};
  TH1F*         hchi2pz[5] = {NULL};
  TH1F*         hnsh[5] = {NULL};
  TH1F*         hlambda[5] = {NULL};
  TH1F*         hd0[5] = {NULL};
  TH1F*         hgenz[5] = {NULL};

  TH1F*         hSeedP    [2] = {NULL};
  TH1F*         hSeedDp   [2] = {NULL};
  TH1F*         hSeedNSh  [2] = {NULL};
  TH1F*         hSeedD0   [2] = {NULL};
  TH1F*         hSeedFrNSh[2] = {NULL};
  TH1F*         hSeedChi2 [2] = {NULL};
  
  
  int  nhist(2);
  int  stnTrkSeedDir[2] = {100, 0};
  for (int i=0; i<nhist; ++i){
    hSeedP   [i]  = gh1(FnStn   ,"TriggerAna", Form("trkseed_%i/p"      , stnTrkSeedDir[i]));
    hSeedP   [i]->GetXaxis()->SetRangeUser(80, 120);
    hSeedDp  [i]  = gh1(FnStn   ,"TriggerAna", Form("trkseed_%i/deltaP1", stnTrkSeedDir[i]));
    hSeedDp  [i]->GetXaxis()->SetRangeUser(-30, 30);
    hSeedNSh [i]  = gh1(FnStn   ,"TriggerAna", Form("trkseed_%i/nhits"  , stnTrkSeedDir[i]));
    hSeedNSh [i]->GetXaxis()->SetRangeUser(0, 100);
    hSeedNSh [i]->GetXaxis()->SetTitle("nStrawHits");
    hSeedD0  [i]  = gh1(FnStn   ,"TriggerAna", Form("trkseed_%i/d0"  , stnTrkSeedDir[i]));
    //hSeedD0->GetXaxis()->SetRangeUser(0, 100);
    hSeedChi2 [i]  = gh1(FnStn   ,"TriggerAna", Form("trkseed_%i/chi2"   , stnTrkSeedDir[i]));
    hSeedFrNSh[i] = gh1(FnStn   ,"TriggerAna", Form("trkseed_%i/FrNHits1"   , stnTrkSeedDir[i]));
    hSeedFrNSh[i]->GetXaxis()->SetRangeUser(0, 1.1);
    hSeedFrNSh[i]->GetXaxis()->SetTitle("nSh_{MC-matched}/nSh");
  }

  int  stnHelixDir[2] = {0, 100};
  for (int i=0; i<nhist; ++i){
    hp[i]  = gh1(FnStn  ,"TriggerAna", Form("helix_%i/p"        , stnHelixDir[i]));
    hdp[i] = gh1(FnStn  ,"TriggerAna", Form("helix_%i/deltaP1"  , stnHelixDir[i]));
    hdp[i]->GetXaxis()->SetRangeUser(-40,40);

    hpt[i] = gh1(FnStn   ,"TriggerAna", Form("helix_%i/pT"      , stnHelixDir[i]));
    //hpt[i]->GetXaxis()->SetRangeUser(-40,40);
    hpzMC[i] = gh1(FnStn   ,"TriggerAna", Form("helix_%i/pzMC1" , stnHelixDir[i]));
    //    hpzMC[i]->GetXaxis()->SetRangeUser(-20,20);

    hnsh[i]  = gh1(FnStn   ,"TriggerAna", Form("helix_%i/nhits" , stnHelixDir[i]));
    hchi2xy[i] = gh1(FnStn   ,"TriggerAna", Form("helix_%i/chi2XY" , stnHelixDir[i]));
    hchi2xy[i]->GetXaxis()->SetRangeUser(0,10);
    hchi2pz[i] = gh1(FnStn   ,"TriggerAna", Form("helix_%i/chi2PhiZ" , stnHelixDir[i]));
    hchi2pz[i]->GetXaxis()->SetRangeUser(0,10);

    hlambda[i] = gh1(FnStn   ,"TriggerAna", Form("helix_%i/lambda" , stnHelixDir[i]));
    // hlambda[i]->Rebin(4);
    //    hlambda[i]->SetLineColor(colors[i]);
    hd0[i] = gh1(FnStn   ,"TriggerAna", Form("helix_%i/d0" , stnHelixDir[i]));
    //    hd0[i]->Rebin(4);
    //    hd0[i]->SetLineColor(colors[i]);
  }
  
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(111);
  TPaveText*tp = new TPaveText(0.19, 0.899, 0.27, 0.939, "NDC");//(0.1, 0.2, 0.96, 0.99, "NDC");
  tp->SetBorderSize(0);
  tp->SetFillStyle(0);
  tp->SetLineStyle(0);
  TText* mu2eLabel = tp->AddText("#font[72]{Mu2e} #font[42]{Simulation preliminary}");
  mu2eLabel->SetTextSize(0.030);
  mu2eLabel->SetTextAlign(13);
    
  TCanvas*c0[2];
  for (int i=0; i<nhist; ++i){
    c0[i] = new TCanvas(Form("c0_%i",i),Form("c0_p_%i", i), 900, 900);
    c0[i]->SetLeftMargin(0.15);
    c0[i]->SetRightMargin(0.05);
    c0[i]->SetTopMargin(0.05);
    c0[i]->SetBottomMargin(0.1);

    //  c0->cd()->SetGrid();
    hp[i]->GetXaxis()->SetRangeUser(60, 140);
    hp[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f MeV/c)", hp[i]->GetBinWidth(1)));
    double max = hp[i]->GetMaximum()*1.2;
    hp[i]->SetMaximum(max);
    hp[i]->Draw();
    gPad->Update();

    tp   ->Draw();
    
    c0[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c0[i]->GetTitle()));
  }


  TCanvas*c1[2], *c2[2], *c3[2], *c4[2], *c5[2], *c6[2], *c7[2], *c8[2];
  TF1*f;
  for (int i=0; i<nhist; ++i){ 
    c1[i] = new TCanvas(Form("c1_%i",i),Form("c1_dp_%i", i), 900, 900);
    c1[i]->SetLeftMargin(0.15);
    c1[i]->SetRightMargin(0.05);
    c1[i]->SetTopMargin(0.05);
    c1[i]->SetBottomMargin(0.1);

    //  c1->cd()->SetGrid();
    //    hdp[i]->GetXaxis()->SetRangeUser(60, 140);
    hdp[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f MeV/c)", hdp[i]->GetBinWidth(1)));
    double max = hdp[i]->GetMaximum()*1.2;
    hdp[i]->SetMaximum(max);
    hdp[i]->Draw();
    hdp[i]->Fit("gaus");
    f = (TF1*)hdp[i]->GetListOfFunctions()->FindObject("gaus");
    f->SetNpx(1000);
    tp   ->Draw();

    c1[i]->Update();
    TPaveStats *stats1 = (TPaveStats*)hdp[i]->GetListOfFunctions()->FindObject("stats");
    stats1->SetTextColor(kBlack); 
    stats1->SetX1NDC(0.67); stats1->SetX2NDC(0.913);
    stats1->SetY1NDC(0.62); stats1->SetY2NDC(0.91);

    c1[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c1[i]->GetTitle()));
    
    //--------------------------------------------------------------------------------
    c2[i] = new TCanvas(Form("c2_%i",i),Form("c2_dpt_%i", i), 900, 900);
    c2[i]->SetLeftMargin(0.15);
    c2[i]->SetRightMargin(0.05);
    c2[i]->SetTopMargin(0.05);
    c2[i]->SetBottomMargin(0.1);

    //  c2->cd()->SetGrid();
    //    hpt[i]->GetXaxis()->SetRangeUser(60, 140);
    hpt[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f MeV/c)", hpt[i]->GetBinWidth(1)));
    max = hpt[i]->GetMaximum()*1.2;
    hpt[i]->SetMaximum(max);
    hpt[i]->Draw();
    hpt[i]->Fit("gaus");
    f = (TF1*)hpt[i]->GetListOfFunctions()->FindObject("gaus");
    f->SetNpx(1000);    
    tp   ->Draw();
    gPad->Update();
    c2[i]->Update();
    TPaveStats *stats2 = (TPaveStats*)hpt[i]->GetListOfFunctions()->FindObject("stats");
    stats2->SetTextColor(kBlack); 
    stats2->SetX1NDC(0.67); stats2->SetX2NDC(0.913);
    stats2->SetY1NDC(0.62); stats2->SetY2NDC(0.91);
    
    c2[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c2[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    c3[i] = new TCanvas(Form("c3_%i",i),Form("c3_dpz_%i", i), 900, 900);
    c3[i]->SetLeftMargin(0.15);
    c3[i]->SetRightMargin(0.05);
    c3[i]->SetTopMargin(0.05);
    c3[i]->SetBottomMargin(0.1);

    //  c3->cd()->SetGrid();
    //    hpzMC[i]->GetXaxis()->SetRangeUser(60, 140);
    hpzMC[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f MeV/c)", hpzMC[i]->GetBinWidth(1)));
    max = hpzMC[i]->GetMaximum()*1.2;
    hpzMC[i]->SetMaximum(max);
    hpzMC[i]->Draw();
    hpzMC[i]->Fit("gaus");
    f = (TF1*)hpzMC[i]->GetListOfFunctions()->FindObject("gaus");
    f->SetNpx(1000);    
    
    tp   ->Draw();
    gPad->Update();

    c3[i]->Update();
    TPaveStats *stats3 = (TPaveStats*)hpzMC[i]->GetListOfFunctions()->FindObject("stats");
    stats3->SetTextColor(kBlack); 
    stats3->SetX1NDC(0.67); stats3->SetX2NDC(0.913);
    stats3->SetY1NDC(0.62); stats3->SetY2NDC(0.91);
    
    c3[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c3[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    c4[i] = new TCanvas(Form("c4_%i",i),Form("c4_nsh_%i", i), 900, 900);
    c4[i]->SetLeftMargin(0.15);
    c4[i]->SetRightMargin(0.05);
    c4[i]->SetTopMargin(0.05);
    c4[i]->SetBottomMargin(0.1);

    //  c4->cd()->SetGrid();
    //    hnsh[i]->GetXaxis()->SetRangeUser(60, 140);
    hnsh[i]->GetYaxis()->SetTitle("Entries");
    hnsh[i]->Rebin();
    max = hnsh[i]->GetMaximum()*1.2;
    hnsh[i]->SetMaximum(max);
    hnsh[i]->Draw();
    tp   ->Draw();
    c4[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c4[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    c5[i] = new TCanvas(Form("c5_%i",i),Form("c5_d0_%i", i), 900, 900);
    c5[i]->SetLeftMargin(0.15);
    c5[i]->SetRightMargin(0.05);
    c5[i]->SetTopMargin(0.05);
    c5[i]->SetBottomMargin(0.1);

    //  c5->cd()->SetGrid();
    //    hd0[i]->GetXaxis()->SetRangeUser(60, 140);
    hd0[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f mm)", hd0[i]->GetBinWidth(1)));
    //hd0[i]->Rebin();
    max = hd0[i]->GetMaximum()*1.2;
    hd0[i]->SetMaximum(max);
    hd0[i]->Draw();
  
    tp   ->Draw();
    c5[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c5[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    c6[i] = new TCanvas(Form("c6_%i",i),Form("c6_lambda_%i", i), 900, 900);
    c6[i]->SetLeftMargin(0.15);
    c6[i]->SetRightMargin(0.05);
    c6[i]->SetTopMargin(0.05);
    c6[i]->SetBottomMargin(0.1);

    //  c6->cd()->SetGrid();
    //    hlambda[i]->GetXaxis()->SetRangeUser(60, 140);
    hlambda[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f mm/rad)", hlambda[i]->GetBinWidth(1)));
    //hlambda[i]->Rebin();
    max = hlambda[i]->GetMaximum()*1.2;
    hlambda[i]->SetMaximum(max);
    hlambda[i]->Draw();
  
    tp   ->Draw();
    c6[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c6[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    c7[i] = new TCanvas(Form("c7_%i",i),Form("c7_chi2xy_%i", i), 900, 900);
    c7[i]->SetLeftMargin(0.15);
    c7[i]->SetRightMargin(0.05);
    c7[i]->SetTopMargin(0.05);
    c7[i]->SetBottomMargin(0.1);

    //  c7->cd()->SetGrid();
    //    hchi2xy[i]->GetXaxis()->SetRangeUser(60, 140);
    hchi2xy[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f mm/rad)", hchi2xy[i]->GetBinWidth(1)));
    //hchi2xy[i]->Rebin();
    max = hchi2xy[i]->GetMaximum()*1.2;
    hchi2xy[i]->SetMaximum(max);
    hchi2xy[i]->Draw();

    tp   ->Draw();
    c7[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c7[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    c8[i] = new TCanvas(Form("c8_%i",i),Form("c8_chi2pz_%i", i), 900, 900);
    c8[i]->SetLeftMargin(0.15);
    c8[i]->SetRightMargin(0.05);
    c8[i]->SetTopMargin(0.05);
    c8[i]->SetBottomMargin(0.1);

    //  c8->cd()->SetGrid();
    //    hchi2pz[i]->GetXaxis()->SetRangeUser(60, 140);
    hchi2pz[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f mm/rad)", hchi2pz[i]->GetBinWidth(1)));
    //hchi2pz[i]->Rebin();
    max = hchi2pz[i]->GetMaximum()*1.2;
    hchi2pz[i]->SetMaximum(max);
    hchi2pz[i]->Draw();
  
    tp   ->Draw();
    c8[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c8[i]->GetTitle()));

  }


  //--------------------------------------------------------------------------------
  TCanvas*cc1[2], *cc2[2], *cc3[2], *cc4[2], *cc5[2], *cc6[2];
  for (int i=0; i<nhist; ++i){ 
    cc1[i] = new TCanvas(Form("cc1_%i",i),Form("cc1_SeedDp_%i", i), 900, 900);
    cc1[i]->SetLeftMargin(0.15);
    cc1[i]->SetRightMargin(0.05);
    cc1[i]->SetTopMargin(0.05);
    cc1[i]->SetBottomMargin(0.1);

    //  cc1->cd()->SetGrid();
    //    hSeedDp[i]->GetXaxis()->SetRangeUser(60, 140);
    hSeedDp[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f MeV/c)", hSeedDp[i]->GetBinWidth(1)));
    double max = hSeedDp[i]->GetMaximum()*1.2;
    hSeedDp[i]->SetMaximum(max);
    hSeedDp[i]->Draw();
    hSeedDp[i]->Fit("gaus");
    TF1*f = (TF1*)hSeedDp[i]->GetListOfFunctions()->FindObject("gaus");
    f->SetNpx(1000);
    tp   ->Draw();

    cc1[i]->Update();
    TPaveStats *stats1 = (TPaveStats*)hSeedDp[i]->GetListOfFunctions()->FindObject("stats");
    stats1->SetTextColor(kBlack); 
    stats1->SetX1NDC(0.67); stats1->SetX2NDC(0.913);
    stats1->SetY1NDC(0.62); stats1->SetY2NDC(0.91);
    
    cc1[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), cc1[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    cc2[i] = new TCanvas(Form("cc2_%i",i),Form("cc2_SeedP_%i", i), 900, 900);
    cc2[i]->SetLeftMargin(0.15);
    cc2[i]->SetRightMargin(0.05);
    cc2[i]->SetTopMargin(0.05);
    cc2[i]->SetBottomMargin(0.1);

    //  cc2->ccd()->SetGrid();
    //    hSeedP[i]->GetXaxis()->SetRangeUser(60, 140);
    hSeedP[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f MeV/c)", hSeedP[i]->GetBinWidth(1)));
    max = hSeedP[i]->GetMaximum()*1.2;
    hSeedP[i]->SetMaximum(max);
    hSeedP[i]->Draw();
    tp   ->Draw();
    
    cc2[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), cc2[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    cc3[i] = new TCanvas(Form("cc3_%i",i),Form("cc3_seedFrNSh_%i", i), 900, 900);
    cc3[i]->SetLeftMargin(0.15);
    cc3[i]->SetRightMargin(0.05);
    cc3[i]->SetTopMargin(0.05);
    cc3[i]->SetBottomMargin(0.1);

    //  cc3->cd()->SetGrid();
    //    hSeedFrNSh[i]->GetXaxis()->SetRangeUser(60, 140);
    hSeedFrNSh[i]->GetYaxis()->SetTitle("Entries");
    max = hSeedFrNSh[i]->GetMaximum()*1.2;
    hSeedFrNSh[i]->SetMaximum(max);
    hSeedFrNSh[i]->Draw();
    tp   ->Draw();
    
    cc3[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), cc3[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    cc4[i] = new TCanvas(Form("cc4_%i",i),Form("cc4_SeedNSh_%i", i), 900, 900);
    cc4[i]->SetLeftMargin(0.15);
    cc4[i]->SetRightMargin(0.05);
    cc4[i]->SetTopMargin(0.05);
    cc4[i]->SetBottomMargin(0.1);

    //  cc4->cd()->SetGrid();
    //    hSeedNSh[i]->GetXaxis()->SetRangeUser(60, 140);
    hSeedNSh[i]->GetYaxis()->SetTitle("Entries");
    hSeedNSh[i]->Rebin();
    max = hSeedNSh[i]->GetMaximum()*1.2;
    hSeedNSh[i]->SetMaximum(max);
    hSeedNSh[i]->Draw();
    tp   ->Draw();
    cc4[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), cc4[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    cc5[i] = new TCanvas(Form("cc5_%i",i),Form("cc5_SeedD0_%i", i), 900, 900);
    cc5[i]->SetLeftMargin(0.15);
    cc5[i]->SetRightMargin(0.05);
    cc5[i]->SetTopMargin(0.05);
    cc5[i]->SetBottomMargin(0.1);

    //  cc5->cd()->SetGrid();
    //    hSeedD0[i]->GetXaxis()->SetRangeUser(60, 140);
    hSeedD0[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f mm)", hSeedD0[i]->GetBinWidth(1)));
    //hSeedD0[i]->Rebin();
    max = hSeedD0[i]->GetMaximum()*1.2;
    hSeedD0[i]->SetMaximum(max);
    hSeedD0[i]->Draw();
  
    tp   ->Draw();
    cc5[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), cc5[i]->GetTitle()));

    //--------------------------------------------------------------------------------
    cc6[i] = new TCanvas(Form("cc6_%i",i),Form("cc6_SeedChi2_%i", i), 900, 900);
    cc6[i]->SetLeftMargin(0.15);
    cc6[i]->SetRightMargin(0.05);
    cc6[i]->SetTopMargin(0.05);
    cc6[i]->SetBottomMargin(0.1);

    //  cc6->cd()->SetGrid();
    //    hSeedChi2[i]->GetXaxis()->SetRangeUser(60, 140);
    hSeedChi2[i]->GetYaxis()->SetTitle(Form("Entries / (%2.2f mm/rad)", hSeedChi2[i]->GetBinWidth(1)));
    //hSeedChi2[i]->Rebin();
    max = hSeedChi2[i]->GetMaximum()*1.2;
    hSeedChi2[i]->SetMaximum(max);
    hSeedChi2[i]->Draw();
  
    tp   ->Draw();
    cc6[i]->Print(Form("Figures/%s_%s.pdf", names[i].Data(), cc6[i]->GetTitle()));

  }



  // TCanvas*c1 = new TCanvas("c1","CRY_d0", 900, 600);
  // hd0  [0]->Draw();
  // for (int i=1; i<nhist; ++i) hd0[i]->Draw("same");
  // leg->Draw();
  // c1->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c1->GetTitle()));
 
  // TCanvas*c2 = new TCanvas("c2","CRY_lambda", 900, 600);
  // hlambda  [0]->Draw();
  // for (int i=1; i<nhist; ++i) hlambda[i]->Draw("same");
  // leg->Draw();
  // c2->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c2->GetTitle()));
 
  // TCanvas*c3 = new TCanvas("c2","CRY_dp", 900, 600);
  // hdp  [0]->Draw();
  // for (int i=1; i<nhist; ++i) hdp[i]->Draw("same");
  // leg->Draw();
  // c3->Print(Form("Figures/%s_%s.pdf", names[i].Data(), c3->GetTitle()));
 
}
