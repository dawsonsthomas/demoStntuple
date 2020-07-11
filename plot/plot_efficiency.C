void plot_rejectionVsEff(const char* signalFile,  const char* bkgFile){
  TCanvas* c1 = new TCanvas("rejectionVsEff","rejectionVsEff",900,600);
  //TCanvas* c2 = new TCanvas("Cal_Thresh_Eff_Dem","Cal_Thresh_Eff_Dem",900,600);
  //TCanvas* c3 = new TCanvas("no_tracker_only_Dem","no_tracker_only__Dem",900,600);
  c1->cd();

  TFile* signal = TFile::Open(signalFile);
  TFile* bkg    = TFile::Open(bkgFile);
  TH1F*  fbkgLhdCorTot       = (TH1F*)bkg->Get("RPCPhotonAna/likelihood/lh_tot");
  TH1F*  fsignalLhdCorTot    = (TH1F*)signal->Get("RPCPhotonAna/likelihood/lh_tot");

  TH1F*  fbkgintegral        = (TH1F*)bkg->Get("RPCPhotonAna/likelihood/lh_tot");
  TH1F*  frpcintegral        = (TH1F*)signal->Get("RPCPhotonAna/likelihood/lh_tot");

  int nbins = fbkgLhdCorTot->GetNbinsX();
  
  Double_t bkgnorm = fbkgLhdCorTot->GetEntries();
  Double_t signorm = fsignalLhdCorTot->GetEntries();

  fbkgLhdCorTot->Scale(1/bkgnorm);
  fsignalLhdCorTot->Scale(1/signorm);

  TGraphErrors*gr = new TGraphErrors();
  gr->SetNameTitle("rejectionVsEff","rejectionVsEff");
  
  for (int i=1; i<=nbins; ++i){
    double rejection = 30000;
    double nbkg = fbkgLhdCorTot->Integral(i,nbins);
    if (nbkg > 0) rejection = 1/nbkg;
    double eff       = fsignalLhdCorTot->Integral(i,nbins);
    gr->SetPoint(i-1, eff, rejection);
    fbkgintegral->SetBinContent(i, nbkg);
    frpcintegral->SetBinContent(i, eff);
  }

  gr->Draw("AP");
  //c2 -> cd();
  //fbkgintegral->Draw();
  //c3 -> cd();
  //frpcintegral->Draw();
}

void plot_timingVsEff(const char* signalFile){
  TCanvas* c1 = new TCanvas("timingVsEff","timingVsEff",900,600);
  TCanvas* c2 = new TCanvas("aghgfd","asdfgh",900,600);
  c1->cd();

  TFile* signal = TFile::Open(signalFile);
  TH1F*  fsignalgent0    = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_time1");
  fsignalgent0->Rebin(4);
  TH1F*  fgent0integral   = (TH1F*)fsignalgent0->Clone();
  fgent0integral -> SetName("clAccumulatTime");

  int nbins = fsignalgent0 ->GetNbinsX();
  double intscale = 1.571;
  Double_t signorm = fsignalgent0->Integral(1,nbins);
  signorm = intscale*signorm;
  fsignalgent0->Scale(1/signorm);
  fsignalgent0->Draw();

  for (int i=1; i<=nbins; ++i){
    double eff  = fsignalgent0->Integral(i,nbins);
    fgent0integral->SetBinContent(i, eff);
    double err(0);
    if (eff>0) err = 0; //sqrt(eff*(1-eff)/(signorm));
    
    fgent0integral->SetBinError(i,err);
  }
  c2->cd();
  fgent0integral->Draw();
 
}

void plot_timmin_v_thresh_dem(){

  TCanvas* c1 = new TCanvas("Cal_Thresh_Eff","Cal_Thresh_Eff",900,600);
  c1->cd();

  TH1F* tim = new TH1F ("Timing", "Timing v. Threshhold",11,40,50);
  double mean[11] = {0.451,0.3952,0.3474,0.3032,0.27,0.2354,0.2073,0.1827,0.1592,0.1444,0.1292};
  double rms[11] = {0.5576,0.5312,0.4698,0.4494,0.4229,0.3852,0.3568,0.3315,0.2959,0.2718,0.2592};

  for(int i=0; i<11; i++){
    tim->SetBinContent(i+1,mean[i]);
  }

  tim->Draw();
  tim->SetTitle("Timing v. Threshhold Electrons (CalTimePeakFinder)");
  tim->GetXaxis()->SetTitle("Cal. E_{Thresh} MeV");
  tim->GetYaxis()->SetTitle("Mean Analysis Time (ms)");
}

void plot_timmin_v_thresh_dep(){

  TCanvas* c1 = new TCanvas("Cal_Thresh_Eff","Cal_Thresh_Eff",900,600);
  c1->cd();

  TH1F* tim = new TH1F ("Timing", "timing v. Threshhold",11,40,50);
  double mean[11] = {0.4415,0.3867,0.3421,0.2958,0.2615,0.2295,0.2014,0.1764,0.1522,0.1389,0.1235};
  double rms[11] = {0.5596,0.5278,0.4775,0.4432,0.4159,0.3885,0.3564,0.3312,0.2934,0.2732,0.2552};

  for(int i=0; i<11; i++){
    tim->SetBinContent(i+1,mean[i]);
  }

  tim->Draw();
  tim->SetTitle("Timing v. Threshhold Positrons (CalTimePeakFinder)");
  tim->GetXaxis()->SetTitle("Cal. E_{Thresh} MeV");
  tim->GetYaxis()->SetTitle("Mean Analysis Time (ms)");
}

//--------------------------------------------------------------------------------
void plot_calo_thresh_eff_dem( const char* Fn0) {

  TCanvas* c1 = new TCanvas("Cal_Thresh_Eff_Dem","Cal_Thresh_Eff_Dem",900,600);
  TCanvas* c2 = new TCanvas("no_tracker_only_Dem","no_tracker_only__Dem",900,600);
  
//-----------------------------------------------------------------------------
// differential efficiency plot for cut set C seletion
//-----------------------------------------------------------------------------
  c1->cd();
  TH1F*  fThreshEffDem    = gh1(Fn0   ,"TriggerAna", "evt_0/ThreshEffDem");//
  TH1F*  fTPRnoneDem    = gh1(Fn0   ,"TriggerAna", "evt_0/TPRnoneDem");//
  TH1F*  fcalThreshNormDem    = gh1(Fn0   ,"TriggerAna", "evt_0/calThreshNormDem");//

  fThreshEffDem->SetLineColor(kGreen);
  fThreshEffDem->SetMarkerColor(kGreen);
  fThreshEffDem->SetMarkerStyle(24);
  fThreshEffDem->SetMarkerSize(1);

  fTPRnoneDem->SetLineColor(kBlue);
  fTPRnoneDem->SetMarkerColor(kBlue);
  fTPRnoneDem->SetMarkerStyle(24);
  fTPRnoneDem->SetMarkerSize(1);

  double   NEvents    =  fThreshEffDem->GetBinContent(1); 
  int      nbins      =  fThreshEffDem->GetNbinsX();
  double eff = 0;
  double err = 0;

  for (int i=0; i<nbins; ++i){
    eff     = fThreshEffDem->GetBinContent(i+1)/NEvents;
    err     = sqrt(eff*(1.-eff)/NEvents);
    fThreshEffDem->SetBinError  (i+1, err    ); 
    fThreshEffDem->SetBinContent(i+1, eff    );

    eff     = fTPRnoneDem->GetBinContent(i+1)/NEvents;
    err     = sqrt(eff*(1.-eff)/NEvents);
    fTPRnoneDem->SetBinError  (i+1, err    ); 
    fTPRnoneDem->SetBinContent(i+1, eff    );
  }

  //TLegend *leg = new TLegend(0.34, 0.62, 0.8363, 0.8676, NULL, "brNDC");
  //leg->SetFillStyle(1001);
  //leg->SetFillColor(kWhite);
  //leg->SetBorderSize(1);
  //leg->SetShadowColor(0);
  //leg->AddEntry(fThreshEff, "Calo-seeded");;
  //leg->AddEntry(fTPRnone, "No Tracker-only");;

  fThreshEffDem->SetTitle("Cal-seeded");
  fThreshEffDem->GetXaxis()->SetTitle("E_{Thresh} For Cal Cluster");
  fThreshEffDem->GetYaxis()->SetTitle("Reconstructed e- tracks (normalized)");
  fTPRnoneDem->SetTitle("Cal-seeded w/ No Tracker-only");
  fTPRnoneDem->GetXaxis()->SetTitle("E_{Thresh} For Cal Cluster");
  fTPRnoneDem->GetYaxis()->SetTitle("Reconstructed e- tracks (normalized)");


  fThreshEffDem->Draw();
  c2 -> cd();
  fTPRnoneDem->Draw();
}


//--------------------------------------------------------------------------------
void plot_calo_thresh_eff_dep( const char* Fn0) {

  TCanvas* c1 = new TCanvas("Cal_Thresh_Eff_Dep","Cal_Thresh_Eff_Dep",900,600);
  TCanvas* c2 = new TCanvas("no_tracker_only_Dep","no_tracker_only_Dep",900,600);
  
//-----------------------------------------------------------------------------
// differential efficiency plot for cut set C seletion
//-----------------------------------------------------------------------------
  c1->cd();
  TH1F*  fThreshEffDep    = gh1(Fn0   ,"TriggerAna", "evt_0/ThreshEffDep");//
  TH1F*  fTPRnoneDep    = gh1(Fn0   ,"TriggerAna", "evt_0/TPRnoneDep");//
  TH1F*  fcalThreshNormDep    = gh1(Fn0   ,"TriggerAna", "evt_0/calThreshNormDep");//

  fThreshEffDep->SetLineColor(kGreen);
  fThreshEffDep->SetMarkerColor(kGreen);
  fThreshEffDep->SetMarkerStyle(24);
  fThreshEffDep->SetMarkerSize(1);

  fTPRnoneDep->SetLineColor(kBlue);
  fTPRnoneDep->SetMarkerColor(kBlue);
  fTPRnoneDep->SetMarkerStyle(24);
  fTPRnoneDep->SetMarkerSize(1);

  double   NEvents    =  fThreshEffDep->GetBinContent(1); 
  int      nbins      =  fThreshEffDep->GetNbinsX();
  double eff = 0;
  double err = 0;

  for (int i=0; i<nbins; ++i){
    eff     = fThreshEffDep->GetBinContent(i+1)/NEvents;
    err     = sqrt(eff*(1.-eff)/NEvents);
    fThreshEffDep->SetBinError  (i+1, err    ); 
    fThreshEffDep->SetBinContent(i+1, eff    );

    eff     = fTPRnoneDep->GetBinContent(i+1)/NEvents;
    err     = sqrt(eff*(1.-eff)/NEvents);
    fTPRnoneDep->SetBinError  (i+1, err    ); 
    fTPRnoneDep->SetBinContent(i+1, eff    );
  }

  //TLegend *leg = new TLegend(0.34, 0.62, 0.8363, 0.8676, NULL, "brNDC");
  //leg->SetFillStyle(1001);
  //leg->SetFillColor(kWhite);
  // leg->SetBorderSize(1);
  // leg->SetShadowColor(0);
  // leg->AddEntry(fThreshEff, "Calo-seeded");;
  // leg->AddEntry(fTPRnone, "No Tracker-only");;

  fThreshEffDep->SetTitle("Cal-seeded");
  fThreshEffDep->GetXaxis()->SetTitle("E_{Thresh} For Cal Cluster");
  fThreshEffDep->GetYaxis()->SetTitle("Reconstructed e+ tracks (normalized)");
  fTPRnoneDep->SetTitle("Cal-seeded w/ No Tracker-only");
  fTPRnoneDep->GetXaxis()->SetTitle("E_{Thresh} For Cal Cluster");
  fTPRnoneDep->GetYaxis()->SetTitle("Reconstructed e+ tracks (normalized)");

  fThreshEffDep->Draw();
  c2 -> cd();
  fTPRnoneDep->Draw();
}
