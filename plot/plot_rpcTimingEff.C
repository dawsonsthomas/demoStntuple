void plot_events_v_t0(){

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
void plot_events_v_t0( const char* file) {

  TCanvas* c1 = new TCanvas("RPCEvents_After_T0","RPCEvents_After_T0",900,600);
  
//-----------------------------------------------------------------------------
// Differential RPC photon production plot for event window t0 seletion
//-----------------------------------------------------------------------------
  c1->cd();
  TFile* f = new TFile(file);
  TH1F*  fGenTime0 = (TH1F*)f->Get("RPCPhotonAna/gen/gen_time0");

  fGenTime0->SetLineColor(kBlue);
  fGenTime0->SetMarkerColor(kBlue);
  fGenTime0->SetMarkerStyle(24);
  fGenTime0->SetMarkerSize(1);

  int      nbins      =  fGenTime0->GetNbinsX();
  double   NEvents    =  fGenTime0->Integral(1,nbins+1);
  double   binwidth   = fGenTime0->GetBinWidth(1);
  double   init       = fGenTime0->GetBinLowEdge(1);
  double   eff    = 0;
  double   err    = 0;

  for (int i=0; i<nbins; i++){
    eff     = fGenTime0->Integral(i+1,nbins+1)/NEvents; 
    err     = sqrt(eff*(1.-eff)/NEvents);
    fGenTime0->SetBinContent(i+1, eff);
    //fGenTime0->SetBinError  (i+1, err);
    fGenTime0->SetBinError (i+1, 0);
  }

  fGenTime0->SetTitle("Proportion RPC photons gen after t0");
  fGenTime0->GetXaxis()->SetTitle("t0");
  fGenTime0->GetYaxis()->SetTitle("Proportion RPC photons gen after t0");


  fGenTime0->Draw();
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
