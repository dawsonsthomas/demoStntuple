
//--------------------------------------------------------------------------------
void plot_rpc_photon_likelihood( const char* bkg, const char* signal) {

  TFile *bkg = TFile::Open(bkg);
  TFile *signal = TFile::Open(signal);

  TCanvas* c1 = new TCanvas("Cal_Thresh_Eff_Dem","Cal_Thresh_Eff_Dem",900,600);
  TCanvas* c2 = new TCanvas("no_tracker_only_Dem","no_tracker_only__Dem",900,600);
  
//import 1d and 2d histograms from signal and bkg root files
  c1->cd();
  TH1F*  bEnergy    = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_E0");
  TH1F*  bTime      = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_time0");
  TH1F*  bnCr       = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_nCr0");
  TH1F*  brDist     = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_rDist0");
  TH1F*  be1eRatio  = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e1eRatio0");
  TH1F*  be2eRatio  = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e2eRatio0");
  TH1F*  be3eRatio  = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e2eRatio0");

  TH2F*  beRDist    = (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/clDisk%i_eRDist%i");
  TH2F*  be1RDist   = (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/clDisk%i_e1RDist%i");
  TH2F*  be2RDist   = (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/clDisk%i_e2RDist%i");
  TH2F*  be3RDist   = (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/clDisk%i_e3RDist%i");
  TH2F*  bnCrRDist  = (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/clDisk%i_nCrRDist%i");

  TH1F*  sEnergy    = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_E1");
  TH1F*  sTime      = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_time1");
  TH1F*  snCr       = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_nCr1");
  TH1F*  srDist     = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_rDist1");
  TH1F*  se1eRatio  = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e1eRatio1");
  TH1F*  se2eRatio  = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e2eRatio1");
  TH1F*  se3eRatio  = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e2eRatio1");

  TH2F*  seRDist    = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/clDisk%i_eRDist%i");
  TH2F*  se1RDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/clDisk%i_e1RDist%i");
  TH2F*  se2RDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/clDisk%i_e2RDist%i");
  TH2F*  se3RDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/clDisk%i_e3RDist%i");
  TH2F*  snCrRDist  = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/clDisk%i_nCrRDist%i");



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
