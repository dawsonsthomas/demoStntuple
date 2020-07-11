//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with DpF > 1 MeV : obviously, misreconstructed ones
//  5  : events with N(tracks) > 1
//  6  : events trk_41 with 0.8< E/P < 1.1 - tracks missed by CalPatRec
//  7  : events (muo) with LogLHRCal >   20
//  8  : events (ele) with LogLHRCal < - 20
//  9  : events (muo) with 0.42 < E/P < 0.46
// 10  : events (muo) with Set C track with ECL > 80 MeV
// 28  : Set C DEM tracks with E/P > 1.1
// 29  : TRK_19 (Set C DEM tracks with a cluster) and LLHR(cal) < 0
// 31  : EVT_6 events with ce_costh > 0.8 
// 32  : TRK_1 events with chi2tcm > 100. 
// 33  : DU < -80mm - study edge effects
// 34  : EVT_7: events with E_CL > 60 and no tracks (makes sense only for single CE events)
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TTrigAna001Module.hh"

ClassImp(TTrigAna001Module)
//-----------------------------------------------------------------------------
TTrigAna001Module::TTrigAna001Module(const char* name, const char* title):
TStnModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);

  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;
  fFillDioHist     = 1;


  fMinTrigHelixNSh = 15;

  //fMinT0 = 0; // do not cut on time by default

  fTrackID = new TStnTrackID();
  fLogLH   = new TEmuLogLH();
  //-----------------------------------------------------------------------------
  // MC truth: define which MC particle to consider as signal
  //-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode = 2;			// conversionGun, 28:StoppedParticleReactionGun
  
  //----------------------------------------------------------------------
  // Get the Czarnecki spectrum 
  //----------------------------------------------------------------------
  fspec        =  TFile::Open("/grid/fermiapp/mu2e/personal/gianipez/Spectra/Czarnecki.root");
  fDIOspectrum = (TH1F*)fspec->Get("Spectrum");

  fEventCounter = 0;
}

//-----------------------------------------------------------------------------
TTrigAna001Module::~TTrigAna001Module() {
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookTimeClusterHistograms   (TimeClusterHist_t*   Hist, const char* Folder){
  
  HBook1F(Hist->fNHits         ,"nhits"      ,Form("%s: # of straw hits"              ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fNComboHits    ,"ncombohits" ,Form("%s: # of combo hits"              ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fClusterEnergy ,"clusterE"   ,Form("%s: cluster energy; E [MeV]      ",Folder), 400,   0,  200,Folder);
  HBook1F(Hist->fT0            ,"t0"         ,Form("%s: t0; t0[ns]"                   ,Folder), 800, 400,  1700,Folder);
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookHelixHistograms   (HelixHist_t*   Hist, const char* Folder){
  
  HBook1F(Hist->fNHits         ,"nhits"      ,Form("%s: # of straw hits"              ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fClusterTime   ,"clusterTime",Form("%s: cluster time; t_{cluster}[ns]",Folder), 800, 400,  1700,Folder);
  HBook1F(Hist->fClusterEnergy ,"clusterE"   ,Form("%s: cluster energy; E [MeV]      ",Folder), 400,   0,  200,Folder);
  HBook1F(Hist->fRadius        ,"radius"     ,Form("%s: curvature radius; r [mm]"     ,Folder), 500,   0,   500,Folder);
  HBook1F(Hist->fMom           ,"p"          ,Form("%s: momentum; p [MeV/c]"          ,Folder), 300,   50,   200,Folder);
  HBook1F(Hist->fPt            ,"pT"         ,Form("%s: pT; pT [MeV/c]"               ,Folder), 600,   0,   150,Folder);
  HBook1F(Hist->fLambda        ,"lambda"     ,Form("%s: abs(#lambda); |#lambda|"      ,Folder), 800,   0,   400,Folder);
  HBook1F(Hist->fT0            ,"t0"         ,Form("%s: t0; t0[ns]"                   ,Folder), 800, 400,  1700,Folder);
  HBook1F(Hist->fT0Err         ,"t0err"      ,Form("%s: t0err; t0err [ns]"            ,Folder), 100,   0,    10,Folder);
  HBook1F(Hist->fD0            ,"d0"         ,Form("%s: D0; d0 [mm]"                  ,Folder), 1600,   -400,    400,Folder);
  HBook1F(Hist->fPDG1          ,"pdg1"       ,Form("%s: pdg most frequent particle"   ,Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fPDGMother1    ,"pdgm1"      ,Form("%s: pdg most freq particle mother",Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fFrNhits1      ,"FrNHits1"   ,Form("%s: Fract CmbHits from most freq prtcl" ,Folder), 200,   0,   2, Folder);
  HBook1F(Hist->fPDG2          ,"pdg2"       ,Form("%s: pdg 2nd most freq particle"   ,Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fPDGMother2    ,"pdgm2"      ,Form("%s: pdg 2nd most freq part mother",Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fFrNhits2      ,"FrNHits2"   ,Form("%s: Frac CmbHits 2nd most freq prtcl"  ,Folder), 200,   0,   2,Folder);
  HBook1F(Hist->fz1Origin      ,"z1Origin"   ,Form("%s: Z Coordinate of Origin1"      ,Folder), 300,   0,   15000,Folder);
  HBook1F(Hist->fr1Origin      ,"r1Origin"   ,Form("%s: Radial distance to Origin1"   ,Folder), 500,   0,   5000,Folder);
  HBook1F(Hist->fz2Origin      ,"z2Origin"   ,Form("%s: Z Coordinate of Origin2"      ,Folder), 300,   0,   15000,Folder);
  HBook1F(Hist->fr2Origin      ,"r2Origin"   ,Form("%s: Radial distance to Origin2"   ,Folder), 500,   0,   5000,Folder);
  HBook1F(Hist->fp1MC          ,"p1MC"       ,Form("%s: Momentum Particle1 MC"        ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fpT1MC         ,"pT1MC"      ,Form("%s: TMomentum Particle1 MC"       ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fp2MC          ,"p2MC"       ,Form("%s: Momentum Particle2 MC"        ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fpT2MC         ,"pT2MC"      ,Form("%s: TMomentum Particle2 MC"       ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fdeltaP1       ,"deltaP1"    ,Form("%s: p minus p1MC"                 ,Folder), 170, -20, 150,Folder);
  HBook1F(Hist->fdeltaP2       ,"deltaP2"    ,Form("%s: p minus p2MC"                 ,Folder), 170, -20, 150,Folder);
  HBook1F(Hist->fpzMC1         ,"pzMC1"      ,Form("%s: Z Momentum"                   ,Folder), 400,-400, 400,Folder);
  HBook1F(Hist->fAlg           ,"alg"        ,Form("%s: algorithm mask"               ,Folder), 10,  0, 10,Folder);
  HBook1F(Hist->fBestAlg       ,"bestAlg"    ,Form("%s: best algorithm "              ,Folder), 10,  0, 10,Folder);
  HBook1F(Hist->fChi2XY        ,"chi2XY"     ,Form("%s: #chi^{2} XY; #chi^{2}_{XY}/ndof;Entries"        ,Folder), 400,  0, 20,Folder);
  HBook1F(Hist->fPhiZ          ,"chi2PhiZ"   ,Form("%s: #chi^{2} #phi-Z;#chi^{2}_{#phi-Z}/ndof;Entries" ,Folder),  400,  0, 20,Folder);
  HBook2F(Hist->fNHitsVsChi2PhiZ,"nHitsVsChi2PhiZ"   ,Form("%s: nhits vs #chi^{2}_{#phi-Z};nHits; #chi^{2}_{#phi-Z}/ndof" ,Folder), 150,   0,   150, 2000,  0, 100,Folder);
  HBook1F(Hist->fDT0           ,"dt0Track"   ,Form("%s: t0_{helix} - t0_{track}; #Delta t0=t0_{helix} - t0_{track}[ns]", Folder), 400, -100,  100,Folder);
  HBook1F(Hist->fDE            ,"dETrack"   ,Form("%s: E_{helix} - E_{track}; #Delta E=E_{helix} - E_{track}[MeV]", Folder), 400, -100,  100,Folder);

  
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookTrackSeedHistograms   (TrackSeedHist_t*   Hist, const char* Folder){
  
  HBook1F(Hist->fNHits       ,"nhits"      ,Form("%s: # of straw hits"              ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fClusterTime ,"clusterTime",Form("%s: cluster time; t_{cluster} [ns]",Folder), 800, 400,  1700,Folder);
  HBook1F(Hist->fClusterEnergy ,"clusterE"   ,Form("%s: cluster energy; E [MeV]      ",Folder), 400,   0,  200,Folder);
  HBook1F(Hist->fRadius      ,"radius"     ,Form("%s: curvature radius; r [mm]"     ,Folder), 500,   0,   500,Folder);
  HBook1F(Hist->fMom         ,"p"          ,Form("%s: momentum; p [MeV/c]"          ,Folder), 300,   50,   200,Folder);
  HBook1F(Hist->fPt         ,"pT"          ,Form("%s: pT; pT [MeV/c]"               ,Folder), 300,   0,   150,Folder);
  HBook1F(Hist->fTanDip      ,"tanDip"     ,Form("%s: abs(tanDip); |tanDip|"               ,Folder), 300,   0,     3,Folder);
  HBook1F(Hist->fChi2        ,"chi2"       ,Form("%s: #chi^{2}; #chi^{2}/ndof"      ,Folder), 100,   0,    10,Folder);
  HBook1F(Hist->fFitCons     ,"fitcons"   ,Form("%s: fit-consistency;Fit-con"      ,Folder), 100,   0,    1,Folder);
  HBook1F(Hist->fD0          ,"d0"         ,Form("%s: D0; d0 [mm]"            ,Folder), 400,   -400,    400,Folder);
  HBook1F(Hist->fPDG1          ,"pdg1"       ,Form("%s: pdg most frequent particle"   ,Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fPDGMother1    ,"pdgm1"      ,Form("%s: pdg most freq particle mother",Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fFrNhits1      ,"FrNHits1"   ,Form("%s: Fract CmbHits from most freq prtcl" ,Folder), 200,   0,   2, Folder);
  HBook1F(Hist->fPDG2          ,"pdg2"       ,Form("%s: pdg 2nd most freq particle"   ,Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fPDGMother2    ,"pdgm2"      ,Form("%s: pdg 2nd most freq part mother",Folder), 2252,   -30,   2222,Folder);
  HBook1F(Hist->fFrNhits2      ,"FrNHits2"   ,Form("%s: Frac CmbHits 2nd most freq prtcl"  ,Folder), 200,   0,   2,Folder);
  HBook1F(Hist->fz1Origin      ,"z1Origin"   ,Form("%s: Z Coordinate of Origin1"      ,Folder), 300,   0,   15000,Folder);
  HBook1F(Hist->fr1Origin      ,"r1Origin"   ,Form("%s: Radial distance to Origin1"   ,Folder), 500,   0,   5000,Folder);
  HBook1F(Hist->fz2Origin      ,"z2Origin"   ,Form("%s: Z Coordinate of Origin2"      ,Folder), 300,   0,   15000,Folder);
  HBook1F(Hist->fr2Origin      ,"r2Origin"   ,Form("%s: Radial distance to Origin2"   ,Folder), 500,   0,   5000,Folder);
  HBook1F(Hist->fp1MC          ,"p1MC"       ,Form("%s: Momentum Particle1 MC"        ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fpT1MC         ,"pT1MC"      ,Form("%s: TMomentum Particle1 MC"       ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fp2MC          ,"p2MC"       ,Form("%s: Momentum Particle2 MC"        ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fpT2MC         ,"pT2MC"      ,Form("%s: TMomentum Particle2 MC"       ,Folder), 150,   0,   150,Folder);
  HBook1F(Hist->fdeltaP1       ,"deltaP1"    ,Form("%s: p minus p1MC"                 ,Folder), 170, -20, 150,Folder);
  HBook1F(Hist->fdeltaP2       ,"deltaP2"    ,Form("%s: p minus p2MC"                 ,Folder), 170, -20, 150,Folder);
  HBook1F(Hist->fpzMC1         ,"pzMC1"      ,Form("%s: Z Momentum"                   ,Folder), 400,-400, 400,Folder);
 
  
}


//-----------------------------------------------------------------------------
void TTrigAna001Module::BookCaloHistograms(CaloHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------
  HBook1F(Hist->fVaneID ,"vane_id",Form("%s: Vane ID"       ,Folder), 10, 0,  10,Folder);

  for (int i=0; i<4; i++) {
    HBook1F(Hist->fEnergy  [i],Form("energy_%i",i),Form("%s: Hit Energy[%i]",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fTime    [i],Form("time_%i"  ,i),Form("%s: Hit time  [%i]",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fNHits   [i],Form("nhits_%i" ,i),Form("%s: NHits     [%i]",Folder,i), 50, 0,  50,Folder);
    HBook1F(Hist->fRadius  [i],Form("r_%i"     ,i),Form("%s: Radius    [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRadiusWE[i],Form("rwe_%i"   ,i),Form("%s: RadiusWE  [%i]",Folder,i),100, 0,1000,Folder);

    HBook1F(Hist->fE700    [i],Form("e700_%i",i),Form("%s: Hit Energy[%i] (T > 700ns)",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fT700    [i],Form("t700_%i",i),Form("%s: Hit time  [%i] (T > 700ns)",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fN700    [i],Form("n700_%i",i),Form("%s: NHits     [%i] (T > 700ns)",Folder,i), 50, 0,  50,Folder);

    HBook1F(Hist->fR700  [i],Form("r700_%i"  ,i),Form("%s: Radius (T>700) [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRWE700[i],Form("rwe700_%i",i),Form("%s: Radius*E(T>700)[%i]",Folder,i),100, 0,1000,Folder);
  }
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
  //   char name [200];
  //   char title[200];

  HBook1F(Hist->fVaneID ,"vane_id",Form("%s: Vane ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(Hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),150, 0, 300,Folder);
  HBook1F(Hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
  HBook1F(Hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
  HBook1F(Hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
  HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
  HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
  HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
  HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
  HBook1F(Hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
  HBook1F(Hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
  //   char name [200];
  //   char title[200];

  HBook1F(Hist->fP      ,"p"       ,Form("%s: Momentum"     ,Folder), 20000,     0, 2000,Folder);
  HBook1F(Hist->fPNorm   ,"pNorm"       ,Form("%s: Momentum (Normalized)"     ,Folder),5000,     0, 500,Folder);
  HBook1F(Hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
  HBook1F(Hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),5000, -2500, 2500,Folder);
  HBook1F(Hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(Hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 6000,     0, 6000,Folder);
  HBook1F(Hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
  HBook1F(Hist->fCosThNorm,"cos_thNorm"  ,Form("%s: Cos(Theta) (Norm weight)"   ,Folder), 200,   -1.,   1.,Folder);
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
  //   char name [200];
  //   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 300,  90  ,120. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 1200, 100.,106.,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),400,   0  ,200. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),400,   0  ,200. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),400,   0  ,200. ,Folder);
  HBook1D(Hist->fPNorm       ,"pNorm"     ,Form("%s: Track P(Norm WT)"   ,Folder), 10000,  0  ,1000. ,Folder);
  Hist->fPNorm->Sumw2(kTRUE);

  HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(Hist->fFitMomErrNorm,"momerrNorm",Form("%s: Track FitMomError(Normalized)" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
  HBook2F(Hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

  HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 500, 50,100,Folder);
  HBook1F(Hist->fPtNorm      ,"ptNorm"    ,Form("%s: Track Pt(Normalized)"  ,Folder), 500, 50,100,Folder);
  HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
  HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Norm    ,"chi2Norm"  ,Form("%s: Track chi2 total(Normalized)"  ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fNDof       ,"ndof"     ,Form("%s: Number of DOF"     ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);
  HBook1F(Hist->fChi2DofC   ,"chi2dc"   ,Form("%s: track chi2/N calc" ,Folder), 500, 0, 10,Folder);
  HBook1F(Hist->fTrkQual[0] ,"trkQual0" ,Form("%s: trkQual for TrkPatRec" ,Folder), 400, -2, 2,Folder);
  HBook1F(Hist->fTrkQual[1] ,"trkQual1" ,Form("%s: trkQual for CalPatRec" ,Folder), 400, -2, 2,Folder);
  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Norm      ,"t0Norm"    ,Form("%s: track T0(Normalized)"  ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 1000, 0,  10,Folder);
  HBook1F(Hist->fT0ErrNorm   ,"t0errNorm" ,Form("%s: track T0Err(Normalized)",Folder), 1000, 0,  10,Folder);
  HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
  HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
  HBook1F(Hist->fFitConsNorm[0],"fconNorm"     ,Form("%s: track fit cons [0](Normalized)",Folder), 200, 0,   1,Folder);
  HBook1F(Hist->fFitConsNorm[1],"fcon1Norm"    ,Form("%s: track fit cons [1](Normalized)",Folder), 1000, 0,   0.1,Folder);
  HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 1600,-400, 400,Folder);
  HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(Hist->fTanDipNorm  ,"tdipNorm"  ,Form("%s: track tan(dip)(Normalized)",Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);
  HBook1F(Hist->fBestAlg    ,"bestAlg"  ,Form("%s: best algorithm "    ,Folder),  10,  0, 10,Folder);


  HBook1F(Hist->fvDt         ,"vdt"      ,Form("%s: track virtual dT"  ,Folder), 2400, -200, 200 ,Folder);
  HBook1F(Hist->fvDt0        ,"vdt0"     ,Form("%s: track virtual dT0" ,Folder), 2400, -200, 200 ,Folder);
  HBook1F(Hist->fT0Pull      ,"t0Pull"   ,Form("%s: track T0 pull"     ,Folder), 400, -20, 20 ,Folder);
  HBook1F(Hist->fvtof        ,"vtof"     ,Form("%s: track tof"  ,Folder), 400,0  ,20 ,Folder);
  HBook1F(Hist->fvDp         ,"vdp"      ,Form("%s: track delta(P)"    ,Folder), 200,-5   ,5 ,Folder);
  HBook1F(Hist->fvDx         ,"vdx"      ,Form("%s: track virtual dX"  ,Folder), 500,-250 ,250,Folder);
  HBook1F(Hist->fvDy         ,"vdy"      ,Form("%s: track virtual dY"  ,Folder), 500,-250 ,250,Folder);
  HBook1F(Hist->fvDz         ,"vdz"      ,Form("%s: track virtual dZ"  ,Folder), 500,-250 ,250,Folder);
  HBook1F(Hist->fvDr         ,"vdr"      ,Form("%s: track virtual dR"  ,Folder), 500,0 ,500,Folder);
  HBook2F(Hist->fvDtVsDt0    ,"vdtVsDt0" ,Form("%s: Track dT Vs dt0"   ,Folder), 200, -20.,20,200,-20.,20,Folder);
  HBook2F(Hist->fvDtVsT0err  ,"vdtVsT0err",Form("%s: Track dT Vs T0err" ,Folder), 200, -20.,20,200,0.,1.,Folder);
  HBook2F(Hist->fvDtVsDr     ,"vdtVsDr"  ,Form("%s: Track dT Vs dR"    ,Folder), 200, -20.,20,250,0.,250,Folder);
  HBook2F(Hist->fvDtVsDp     ,"vdtVsDp"  ,Form("%s: Track dT Vs Dp"    ,Folder), 200, -20.,20,200,-5.,5,Folder);
  HBook2F(Hist->fvtofVsDp     ,"vtofVsDp"  ,Form("%s: Track tof Vs Dp"    ,Folder), 400, 0.,20,200,-5.,5,Folder);

  HBook1F(Hist->fvDtNorm      ,"vdtNorm"      ,Form("%s: track virtual dT(Normalized)"  ,Folder), 2400, -200, 200 ,Folder);
  HBook1F(Hist->fvDt0Norm     ,"vdt0Norm"     ,Form("%s: track virtual dT0(Normalized)" ,Folder), 2400, -200, 200 ,Folder);
  HBook1F(Hist->fvtofNorm    ,"vtofNorm"    ,Form("%s: track tof(Normalized)"  ,Folder), 400, 0  ,20 ,Folder);
  HBook1F(Hist->fvDpNorm      ,"vdpNorm"      ,Form("%s: track delta(P)(Normalized)"    ,Folder), 200,-5   ,5 ,Folder);
  HBook1F(Hist->fvDxNorm      ,"vdxNorm"      ,Form("%s: track virtual dX(Normalized)"  ,Folder), 500,-250 ,250,Folder);
  HBook1F(Hist->fvDyNorm      ,"vdyNorm"      ,Form("%s: track virtual dY(Normalized)"  ,Folder), 500,-250 ,250,Folder);
  HBook1F(Hist->fvDzNorm      ,"vdzNorm"      ,Form("%s: track virtual dZ(Normalized)"  ,Folder), 500,-250 ,250,Folder);
  HBook1F(Hist->fvDrNorm      ,"vdrNorm"      ,Form("%s: track virtual dR(Normalized)"  ,Folder), 500,0 ,500,Folder);
  HBook2F(Hist->fvDtVsDt0Norm ,"vdtVsDt0Norm" ,Form("%s: Track dT Vs dt0(Normalized)"   ,Folder), 200, -20.,20,200,-20.,20,Folder);
  HBook2F(Hist->fvDtVsT0errNorm  ,"vdtVsT0errNorm",Form("%s: Track dT Vs T0err(Normalized)" ,Folder), 200, -20.,20,200,0.,1.,Folder);
  HBook2F(Hist->fvDtVsDrNorm     ,"vdtVsDrNorm"  ,Form("%s: Track dT Vs dR(Normalized)"    ,Folder), 200, -20.,20,250,0.,250,Folder);
  HBook2F(Hist->fvDtVsDpNorm     ,"vdtVsDpNorm"  ,Form("%s: Track dT Vs Dp(Normalized)"    ,Folder), 200, -20.,20,200,-5.,5,Folder);
  HBook2F(Hist->fvtofVsDpNorm     ,"vtofVsDpNorm"  ,Form("%s: Track tof Vs Dp(Normalized)"    ,Folder), 200, -20.,20,200,-5.,5,Folder);

  HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 800,-80  ,80 ,Folder);
  HBook1F(Hist->fChi2Match  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);

  HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-250 ,250,Folder);
  HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-250 ,250,Folder);
  HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-100 ,100,Folder);
  HBook2F(Hist->fDvVsDu     ,"dv_vs_du" ,Form("%s: Track Dv Vs Du"    ,Folder), 100, -250,250,100,-100.,100,Folder);
  HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);
  HBook2F(Hist->fDuVsPath   ,"du_vs_path",Form("%s: Track Du Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDucVsPath  ,"duc_vs_path",Form("%s: T-C Duc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsPath   ,"dv_vs_path",Form("%s: T-C  Dv Vs Path"  ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvcVsPath  ,"dvc_vs_path",Form("%s: T-C Dvc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDtVsPath   ,"dt_vs_path",Form("%s: T-C DT Vs Path"   ,Folder),  50,   0 ,500,100,  -5.,  5.,Folder);
  HBook2F(Hist->fDuVsTDip   ,"du_vs_tdip",Form("%s: Track Du Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsTDip   ,"dv_vs_tdip",Form("%s: Track Dv Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);

  HBook1F(Hist->fZ1         ,"z1"       ,Form("%s: track Z1      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
  HBook1F(Hist->fVaneID     ,"vid"      ,Form("%s: track vane ID"     ,Folder),  10,-5   ,  5,Folder);
  HBook1F(Hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
  HBook1F(Hist->fXTrk       ,"xtrk"     ,Form("%s: track XTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fRTrk       ,"rtrk"     ,Form("%s: track RTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
  HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 200, 0   ,2,Folder);
  HBook1F(Hist->fEp[0]      ,"ep_0"     ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
  HBook1F(Hist->fEp[1]      ,"ep_1"     ,Form("%s: track E/P"         ,Folder), 3000, 0   ,15.,Folder);
  HBook2F(Hist->fEpVsPath   ,"ep_vs_path",Form("%s: E/P Vs Path"      ,Folder),  50,   0 ,500,150,  0.,  1.5,Folder);
  HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
  HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

  HBook1F(Hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

  HBook2F(Hist->fEpVsDt     ,"ep_vs_dt" ,Form("%s: E/P vs Dt"         ,Folder), 200, -10, 10,150,0.,1.5,Folder);
  HBook1F(Hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
  HBook1F(Hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) XSlope" ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHR     ,"llhr"     ,Form("%s: LogLH(e/m)"        ,Folder), 200,-100 ,100,Folder);

  HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 1000,-500,500,Folder);
  HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);

  HBook2F(Hist->fNEPlVsNHPl ,"nep_vs_nhp",Form("%s: Track NEXP vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fNDPlVsNHPl ,"ndp_vs_nhp",Form("%s: Track NDIF vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fChi2dVsNDPl,"chi2d_vs_ndp",Form("%s: Track Chi2/Dof vs NDP",Folder), 30, 0,30,100,0.,10,Folder);
  HBook2F(Hist->fDpFVsNDPl  ,"dpf_vs_ndp"  ,Form("%s: Track DpF vs NDP",Folder)     , 30, 0,30,100,-5,5,Folder);

  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  HBook1F(Hist->fcalThreshNormDem  ,"calThreshNormDem" ,Form("%s: Num of Events With Tracks from Cal-seeded or Tracker-only; E_{Calo}"  ,Folder),60,40,100,Folder);
  HBook1F(Hist->fTPRnoneDem  ,"TPRnoneDem" ,Form("%s: Num of Events With no Tracker-only Tracks; E_{Calo}"  ,Folder),60,40,100,Folder);
  HBook1F(Hist->fThreshEffDem  ,"ThreshEffDem" ,Form("%s: Efficiancy of Cal-seeded recon vs energy threshhold; E_{Calo}"  ,Folder),60,40,100,Folder);
  HBook1F(Hist->fcalThreshNormDep  ,"calThreshNormDep" ,Form("%s: Num of Events With Tracks from Cal-seeded or Tracker-only; E_{Calo}"  ,Folder),60,40,100,Folder);
  HBook1F(Hist->fTPRnoneDep  ,"TPRnoneDep" ,Form("%s: Num of Events With no Tracker-only Tracks; E_{Calo}"  ,Folder),60,40,100,Folder);
  HBook1F(Hist->fThreshEffDep  ,"ThreshEffDep" ,Form("%s: Efficiancy of Cal-seeded recon vs energy threshhold; E_{Calo}"  ,Folder),60,40,100,Folder);
  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1D(Hist->fNormMom    ,"dio_mom" ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(Hist->fInstLum    ,"instLum" ,Form("%s: instantaneous p beam intensity"  ,Folder), 1000, 1e6, 1e9,Folder);
  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  HBook1F(Hist->fNClusters[0] ,"ncl0"      ,Form("%s: Number of Reconstructed Clusters",Folder),101,-0.5,100.5,Folder);
  HBook1F(Hist->fNClusters[1] ,"ncl1"      ,Form("%s: Number of Reconstructed Clusters with E>50 MeV",Folder),101,-0.5,100.5,Folder);
  HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
  HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);
  HBook1F(Hist->fDNHits     ,"dNHits"    ,Form("%s: nHits_{trkseed} - nHits_{helix}; nHits_{trkseed} - nHits_{helix}",Folder), 201,-100.5,100.5,Folder);

  HBook1F(Hist->fNTClHits[0]   ,"nTClHits0",Form("%s: N_{comboHits} in the TimeCluster; n_{ComboHits}",Folder), 101, -0.5, 100.5,Folder);
  HBook1F(Hist->fNTClHits[1]   ,"nTClHits1",Form("%s: N_{comboHits} in the TimeCluster - track with #chi^{2}/ndof<5, nhits#geq20; n_{ComboHits}",Folder), 101, -0.5, 100.5,Folder);
  HBook1F(Hist->fNTClHits[2]   ,"nTClHits2",Form("%s: N_{comboHits} in the TimeCluster - track with #chi^{2}/ndof<5, nhits#geq25; n_{ComboHits}",Folder), 101, -0.5, 100.5,Folder);

  HBook1F(Hist->fCprDemTrigger[0]   ,"cprdemtrigger0",Form("%s: efficiency vs nhits-cut"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[1]   ,"cprdemtrigger1",Form("%s:3D-fit efficiency vs nhits-cut and #chi^{2} < 4",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[2]   ,"cprdemtrigger2",Form("%s 3D-fit: efficiency vs nhits-cut",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[3]   ,"cprdemtrigger3",Form("%s 3D-fit: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[4]   ,"cprdemtrigger4",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[5]   ,"cprdemtrigger5",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[6]   ,"cprdemtrigger6",Form("%s 3D-fit: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[8]   ,"cprdemtrigger8",Form("%s 3D-fit: ce efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[7]   ,"cprdemtrigger7",Form("%s 3D-fit: efficiency vs cluster energy, requiring 15 straw hits in the time peak",Folder), 8, 40, 80,Folder);
  HBook1F(Hist->fCprDemTrigger[9]   ,"cprdemtrigger9",Form("%s 3D-fit-ce: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[10]  ,"cprdemtrigger10",Form("%s 3D-fit-ce: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[11]  ,"cprdemtrigger11",Form("%s 3D-fit-ce: efficiency vs nhits-cut and p>80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[12]  ,"cprdemtrigger12",Form("%s: 3D-fit ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[13]  ,"cprdemtrigger13",Form("%s: pat-rec  efficiency vs nhits-cut"       ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[14]  ,"cprdemtrigger14",Form("%s: pat-rec ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDemTrigger[15]  ,"cprdemtrigger15",Form("%s: pat-rec ce  efficiency vs nhits-cut and #Chi^{2}<4 (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  
  HBook1F(Hist->fTprDemTrigger[0]   ,"tprdemtrigger0",Form("%s: efficiency vs nhits-cut"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[1]   ,"tprdemtrigger1",Form("%s:3D-fit efficiency vs nhits-cut and #chi^{2} < 4",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[2]   ,"tprdemtrigger2",Form("%s 3D-fit: efficiency vs nhits-cut",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[3]   ,"tprdemtrigger3",Form("%s 3D-fit: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[4]   ,"tprdemtrigger4",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[5]   ,"tprdemtrigger5",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[6]   ,"tprdemtrigger6",Form("%s 3D-fit: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[8]   ,"tprdemtrigger8",Form("%s 3D-fit: ce efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[7]   ,"tprdemtrigger7",Form("%s 3D-fit: efficiency vs cluster energy, requiring 15 straw hits in the time peak",Folder), 8, 40, 80,Folder);
  HBook1F(Hist->fTprDemTrigger[9]   ,"tprdemtrigger9",Form("%s 3D-fit-ce: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[10]  ,"tprdemtrigger10",Form("%s 3D-fit-ce: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[11]  ,"tprdemtrigger11",Form("%s 3D-fit-ce: efficiency vs nhits-cut and p>80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[12]  ,"tprdemtrigger12",Form("%s: 3D-fit ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[13]  ,"tprdemtrigger13",Form("%s: pat-rec  efficiency vs nhits-cut"       ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[14]  ,"tprdemtrigger14",Form("%s: pat-rec ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDemTrigger[15]  ,"tprdemtrigger15",Form("%s: pat-rec ce  efficiency vs nhits-cut and #Chi^{2}<4 (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  
  HBook1F(Hist->fCprDepTrigger[0]   ,"cprdeptrigger0",Form("%s: efficiency vs nhits-cut"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[1]   ,"cprdeptrigger1",Form("%s:3D-fit efficiency vs nhits-cut and #chi^{2} < 4",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[2]   ,"cprdeptrigger2",Form("%s 3D-fit: efficiency vs nhits-cut",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[3]   ,"cprdeptrigger3",Form("%s 3D-fit: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[4]   ,"cprdeptrigger4",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[5]   ,"cprdeptrigger5",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[6]   ,"cprdeptrigger6",Form("%s 3D-fit: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[8]   ,"cprdeptrigger8",Form("%s 3D-fit: ce efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[7]   ,"cprdeptrigger7",Form("%s 3D-fit: efficiency vs cluster energy, requiring 15 straw hits in the time peak",Folder), 8, 40, 80,Folder);
  HBook1F(Hist->fCprDepTrigger[9]   ,"cprdeptrigger9",Form("%s 3D-fit-ce: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[10]  ,"cprdeptrigger10",Form("%s 3D-fit-ce: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[11]  ,"cprdeptrigger11",Form("%s 3D-fit-ce: efficiency vs nhits-cut and p>80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[12]  ,"cprdeptrigger12",Form("%s: 3D-fit ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[13]  ,"cprdeptrigger13",Form("%s: pat-rec  efficiency vs nhits-cut"       ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[14]  ,"cprdeptrigger14",Form("%s: pat-rec ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fCprDepTrigger[15]  ,"cprdeptrigger15",Form("%s: pat-rec ce  efficiency vs nhits-cut and #Chi^{2}<4 (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  
  HBook1F(Hist->fTprDepTrigger[0]   ,"tprdeptrigger0",Form("%s: efficiency vs nhits-cut"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[1]   ,"tprdeptrigger1",Form("%s:3D-fit efficiency vs nhits-cut and #chi^{2} < 4",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[2]   ,"tprdeptrigger2",Form("%s 3D-fit: efficiency vs nhits-cut",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[3]   ,"tprdeptrigger3",Form("%s 3D-fit: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[4]   ,"tprdeptrigger4",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[5]   ,"tprdeptrigger5",Form("%s 3D-fit: efficiency vs nhits-cut and p > 80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[6]   ,"tprdeptrigger6",Form("%s 3D-fit: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[8]   ,"tprdeptrigger8",Form("%s 3D-fit: ce efficiency vs nhits-cut and p > 80 MeV/c",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[7]   ,"tprdeptrigger7",Form("%s 3D-fit: efficiency vs cluster energy, requiring 15 straw hits in the time peak",Folder), 8, 40, 80,Folder);
  HBook1F(Hist->fTprDepTrigger[9]   ,"tprdeptrigger9",Form("%s 3D-fit-ce: efficiency vs nhits-cut and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[10]  ,"tprdeptrigger10",Form("%s 3D-fit-ce: efficiency vs nhits-cut and #chi^{2}/ndof < 3 and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[11]  ,"tprdeptrigger11",Form("%s 3D-fit-ce: efficiency vs nhits-cut and p>80 MeV/c and d0#in[-200,200] mm",Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[12]  ,"tprdeptrigger12",Form("%s: 3D-fit ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[13]  ,"tprdeptrigger13",Form("%s: pat-rec  efficiency vs nhits-cut"       ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[14]  ,"tprdeptrigger14",Form("%s: pat-rec ce  efficiency vs nhits-cut (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  HBook1F(Hist->fTprDepTrigger[15]  ,"tprdeptrigger15",Form("%s: pat-rec ce  efficiency vs nhits-cut and #Chi^{2}<4 (normalized to N-good-tracks)"        ,Folder), 11, 9.5, 20.5,Folder);
  
  //  char  name[200];
  // for (int i=0; i<2; i++) {
  //   sprintf(name,"ncch_%i",i);
  //   HBook1F(Hist->fNCaloCrystalHits[i],name,Form("%s: N(calo crystal hits) [%i]",Folder,i),500,0,1000,Folder);
  //   sprintf(name,"ncch_vs_vane_%i",i);
  //   HBook2F(Hist->fNCaloHitsVsVane[i],name,Form("%s: N(calo crystal hits) vs vane[%i]",Folder,i),4,0,4,200,0,200,Folder);
  //   sprintf(name,"ncch_vs_row_%i",i);
  //   HBook2F(Hist->fNCaloHitsVsRow[i],name,Form("%s: N(calo crystal hits) vs row [%i]",Folder,i),20,0,20,200,0,200,Folder);
  //   sprintf(name,"ncch_vs_col_%i",i);
  //   HBook2F(Hist->fNCaloHitsVsCol[i],name,Form("%s: N(calo crystal hits) vs col [%i]",Folder,i),50,0,50,200,0,200,Folder);
  // }

  // for (int i=0; i<2; i++) {
  //   HBook1F(Hist->fETot        [i],Form("etot_%i"    ,i),Form("%s: Etot[%i]",Folder,i), 300, 0,150,Folder);
  //   HBook2F(Hist->fECrVsR      [i],Form("ecr_vs_r_%i",i),Form("%s: E Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,500,0,100,Folder);
  //   HBook2F(Hist->fNCrVsR      [i],Form("ncr_vs_r_%i",i),Form("%s: N Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,100,0,100,Folder);

  //   HBook2F(Hist->fNCrystalHitsVsR[i],Form("ncrh_vs_r_%i",i),Form("%s: N Crystal Hits[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
  //   HBook2F(Hist->fNHitCrystalsVsR[i],Form("nhcr_vs_r_%i",i),Form("%s: N Hit Crystals[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
  // }

  // HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot",Form("%s: NHit Crystals Tot",Folder), 100, 0,100,Folder);
  // HBook1F(Hist->fECal           ,"ecal",Form("%s: E(cal), sum over both disks",Folder), 500, 0,250,Folder);
  // HBook1F(Hist->fECalOverEKin   ,"ec_over_ek",Form("%s: E(cal)/E(kin)",Folder), 200, 0,2,Folder);
  HBook1F(Hist->fDNHitsHelix    ,"dNHits_hel",Form("%s: #Delta N_{hits}; N_{cpr} - N_{tpr}; Entries",Folder), 101, -50.5,50.5,Folder);
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode   ,"pdg"         ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fNStrawHits,"nsth"        ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
}
//_____________________________________________________________________________
void TTrigAna001Module::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

  //-----------------------------------------------------------------------------
  // book crystal histograms
  //-----------------------------------------------------------------------------
  HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");

  //--------------------------------------------------------------------------------
  // book timecluster histograms
  //--------------------------------------------------------------------------------
  int book_timecluster_histset[kNTimeClusterHistSets];
  for (int i=0; i<kNTimeClusterHistSets; ++i)  book_timecluster_histset[i] = 0;

  book_timecluster_histset[0] = 1;   // all events
  book_timecluster_histset[1] = 1;   // timeclusters with NHits>10
  book_timecluster_histset[2] = 1;   // timeclusters with NHits>15

  for (int i=0; i<kNTimeClusterHistSets; i++) {
    if (book_timecluster_histset[i] != 0) {
      sprintf(folder_name,"timecluster_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTimeCluster[i] = new TimeClusterHist_t;
      BookTimeClusterHistograms(fHist.fTimeCluster[i],Form("Hist/%s",folder_name));//FIXME!
    }
  }
  

  //--------------------------------------------------------------------------------
  // book trackSeed histograms
  //--------------------------------------------------------------------------------
  int book_trackSeed_histset[kNTrackSeedHistSets];
  for (int i=0; i<kNTrackSeedHistSets; ++i)  book_trackSeed_histset[i] = 0;

  //CalSeedFit Dem
  book_trackSeed_histset[0] = 1;   // events with at least one trackSeed
  book_trackSeed_histset[1] = 1;   // events with at least one trackSeed with p > 80 MeV/c
  book_trackSeed_histset[2] = 1;   // events with at least one trackSeed with p > 90 MeV/c
  book_trackSeed_histset[3] = 1;   // events with at least one trackSeed with p > 100 MeV/c
  book_trackSeed_histset[4] = 1;   // events with at least one trackSeed with 10 < nhits < 15
  book_trackSeed_histset[5] = 1;   // events with at least one trackSeed with nhits >= 15
  book_trackSeed_histset[6] = 1;   // events with one trck passing trkQual
  book_trackSeed_histset[7] = 1;   // events with at least one trackSeed with nhits >= 15 and chi2XY(ZPhi)<4
  book_trackSeed_histset[8] = 1;   // events with one trck passing trkQual and at least one trackSeed with nhits >= 15 and chi2XY(ZPhi)<4
  book_trackSeed_histset[9] = 1;   // events with one matched track from CalTrkFit
  book_trackSeed_histset[10] = 1;  // events with one trck passing trkQual and d0 in [-20, 20] cm
  book_trackSeed_histset[15] = 1;  // tracks in the events that passed the HelixTrigger cut
  book_trackSeed_histset[20] = 1;  // track from signal-generated particle
  
  //TrkSeedFit Dem
  book_trackSeed_histset[100] = 1;   // all trackSeeds
  book_trackSeed_histset[115] = 1;  // tracks in the events that passed the HelixTrigger cut
  book_trackSeed_histset[120] = 1;   // trackSeed from signal-generated particle

  //CalSeedFit Dep
  book_trackSeed_histset[200] = 1;   // events with at least one trackSeed
  book_trackSeed_histset[215] = 1;  // tracks in the events that passed the HelixTrigger cut
  book_trackSeed_histset[220] = 1;   // trackSeed from signal-generated particle

  //TrkSeedFit Dep
  book_trackSeed_histset[300] = 1;   // events with at least one trackSeed
  book_trackSeed_histset[315] = 1;  // tracks in the events that passed the HelixTrigger cut
  book_trackSeed_histset[320] = 1;   // trackSeed from signal-generated particle



  for (int i=0; i<kNTrackSeedHistSets; i++) {
    if (book_trackSeed_histset[i] != 0) {
      sprintf(folder_name,"trkseed_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrackSeed[i] = new TrackSeedHist_t;
      BookTrackSeedHistograms(fHist.fTrackSeed[i],Form("Hist/%s",folder_name));
    }
  }
  

  //--------------------------------------------------------------------------------
  // book helix histograms
  //--------------------------------------------------------------------------------
  int book_helix_histset[kNHelixHistSets];
  for (int i=0; i<kNHelixHistSets; ++i)  book_helix_histset[i] = 0;

  //Helices from  TrkPatRec downstrem e-
  book_helix_histset[0] = 1;   // events with at least one helix
  book_helix_histset[1] = 1;   // events with at least one helix with p > 80 MeV/c
  book_helix_histset[2] = 1;   // events with at least one helix with p > 90 MeV/c
  book_helix_histset[3] = 1;   // events with at least one helix with p > 100 MeV/c
  book_helix_histset[4] = 1;   // events with at least one helix with 10 < nhits < 15
  book_helix_histset[5] = 1;   // events with at least one helix with nhits >= 15
  book_helix_histset[6] = 1;   // events with at least one helix with nhits >= 15 and chi2XY(ZPhi)<=5
  
  book_helix_histset[8] = 1;   // helices associated with a reconstructed track
  book_helix_histset[9] = 1;   // helices associated with a reconstructed track from signal
  book_helix_histset[10] = 1;   // events from a particle generated in the Al target with 100 < p < 105

  //Helices from  TrkPatRec: positron downstream
  book_helix_histset[300] = 1;   // events with at least one helix
  book_helix_histset[301] = 1;   // events with at least one helix with p > 80 MeV/c
  book_helix_histset[302] = 1;   // events with at least one helix with p > 90 MeV/c
  book_helix_histset[303] = 1;   // events with at least one helix with p > 100 MeV/c
  book_helix_histset[304] = 1;   // events with at least one helix with 10 < nhits < 15
  book_helix_histset[305] = 1;   // events with at least one helix with nhits >= 15
  book_helix_histset[306] = 1;   // events with at least one helix with nhits >= 15 and chi2XY(ZPhi)<=5
  book_helix_histset[308] = 1;   // helices associated with a reconstructed track
  book_helix_histset[309] = 1;   // helices associated with a reconstructed track form signal
  book_helix_histset[310] = 1;   // events from a particle generated in the Al target with 90 < p < 100

  
  //Helices from CalPatRec Dem
  book_helix_histset[100] = 1;   // events with at least one helix
  book_helix_histset[101] = 1;   // events with at least one helix with p > 80 MeV/c
  book_helix_histset[102] = 1;   // events with at least one helix with p > 90 MeV/c
  book_helix_histset[103] = 1;   // events with at least one helix with p > 100 MeV/c
  book_helix_histset[104] = 1;   // events with at least one helix with 10 < nhits < 15
  book_helix_histset[105] = 1;   // events with at least one helix with nhits >= 15
  book_helix_histset[106] = 1;   // events with at least one helix with nhits >= 15 and chi2XY(ZPhi)<=5
  book_helix_histset[107] = 1;   // events with at least one track passing quality cuts
  book_helix_histset[108] = 1;   // helices associated with a reconstructed track
  book_helix_histset[109] = 1;   // helices associated with a reconstructed track from signal
  book_helix_histset[110] = 1;   // events from a particle generated in the Al target with 100 < p < 105

  book_helix_histset[120] = 1;   // helices with chi2_zphi > 5

  //Helices from CalPatRec Dep
  book_helix_histset[200] = 1;   // events with at least one helix
  book_helix_histset[201] = 1;   // events with at least one helix with p > 80 MeV/c
  book_helix_histset[202] = 1;   // events with at least one helix with p > 90 MeV/c
  book_helix_histset[203] = 1;   // events with at least one helix with p > 100 MeV/c
  book_helix_histset[204] = 1;   // events with at least one helix with 10 < nhits < 15
  book_helix_histset[205] = 1;   // events with at least one helix with nhits >= 15
  book_helix_histset[206] = 1;   // events with at least one helix with nhits >= 15 and chi2XY(ZPhi)<=5
  book_helix_histset[208] = 1;   // helices associated with a reconstructed track
  book_helix_histset[209] = 1;   // helices associated with a reconstructed track form signal
  book_helix_histset[210] = 1;   // events from a particle generated in the Al target with 90 < p < 100

  for (int i=0; i<kNHelixHistSets; i++) {
    if (book_helix_histset[i] != 0) {
      sprintf(folder_name,"helix_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fHelix[i] = new HelixHist_t;
      BookHelixHistograms(fHist.fHelix[i],Form("Hist/%s",folder_name));
    }
  }
  

  //-----------------------------------------------------------------------------
  // book event histograms
  //-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with a reconstructed track
  book_event_histset[ 2] = 1;	        // events without reconstructed tracks
  book_event_histset[ 3] = 1;	        // events with a reconstructed cluster
  book_event_histset[ 4] = 1;	        // events without reconstructed clusters
  book_event_histset[ 5] = 1;	        // events w/o reconstructed tracks, |costh|<0.4
  book_event_histset[ 6] = 1;	        // events with tracks passing "Set C" cuts
  book_event_histset[ 7] = 1;	        // events with E(cluster) > 60 MeV
  book_event_histset[ 8] = 1;	        // events with the highest energy cluster on the 1st disk
  book_event_histset[ 9] = 1;	        // events with the highest energy cluster on the 2nd disk
  book_event_histset[10] = 0;	        // 
  book_event_histset[11] = 1;	        // selection cuts
  book_event_histset[12] = 1;	        // 
  book_event_histset[13] = 1;	        // 
  book_event_histset[14] = 1;	        // 
  book_event_histset[15] = 1;	        // 
  book_event_histset[16] = 1;	        // 
  book_event_histset[17] = 1;	        // 
  book_event_histset[18] = 1;	        // 
  book_event_histset[19] = 0;	        // 
  book_event_histset[20] = 0;	        // 
  book_event_histset[21] = 0;	        // 
  book_event_histset[22] = 0;	        // 
  book_event_histset[23] = 0;	        // 
					// TrkPatRec tracks
  book_event_histset[24] = 1;	        // events with at least one reco track
  book_event_histset[25] = 1;	        // 
  book_event_histset[26] = 1;	        // 
  book_event_histset[27] = 1;	        // 
  book_event_histset[28] = 1;	        // 

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
      fBookedEvent[i] = 1;
    }else {
      fBookedEvent[i] = 0;
    }
  }//-----------------------------------------------------------------------------
  // book simp histograms
  //-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
  //-----------------------------------------------------------------------------
  // book track histograms
  //-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all tracks e-
  book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts 
  book_track_histset[  2] = 1;		// all tracks e- passing trkQual cut
  book_track_histset[  3] = 1;		// "Set C" tracks with 100 <= P < 110 
  book_track_histset[  4] = 1;		// "trkQual" tracks with 100 <= P < 110 

  book_track_histset[  5] = 1;		// "trkQual" tracks with 100 <= P < 110, and CalTrkFit e- with chi2 < 3 & d0 in [-20, 20] cm &  nhits >=15
  book_track_histset[  6] = 1;		// "trkQual" tracks with 100 <= P < 110, and CalTrkFit e- with chi2 < 4 & d0 in [-20, 20] cm &  nhits >=15

  book_track_histset[  7] = 1;          // "trkQual" tracks with 100 <= P < 110 found by CalPatRec
  book_track_histset[  8] = 1;          // "trkQual" tracks with 100 <= P < 110 found by CalPatRec with 100 <= P < 110,
  book_track_histset[  9] = 1;          // "trkQual" tracks with 100 <= P < 110 found by CalPatRec with 100 <= P < 110, and CalTrkFit e- with chi2 < 3 & d0 in [-20, 20] cm &  nhits >=15
  book_track_histset[ 10] = 1;          // "trkQual" tracks with 100 <= P < 110 found by CalPatRec with 100 <= P < 110, and CalTrkFit e- with chi2 < 4 & d0 in [-20, 20] cm &  nhits >=15


  book_track_histset[ 11] = 1;		// "trkQual" tracks with 100 <= P < 110, and CalTrkFit e- with chi2 < 3 & d0 in [-20, 20] cm &  nhits >=10
  book_track_histset[ 12] = 1;		// "trkQual" tracks with 100 <= P < 110, and CalTrkFit e- with chi2 < 4 & d0 in [-20, 20] cm &  nhits >=10

  book_track_histset[ 13] = 1;		// "trkQual" tracks with 100 <= P < 110 found by CalPatRec , and CalTrkFit e- with chi2 < 3 & d0 in [-20, 20] cm &  nhits >=10
  book_track_histset[ 14] = 1;		// "trkQual" tracks with 100 <= P < 110 found by CalPatRec , and CalTrkFit e- with chi2 < 4 & d0 in [-20, 20] cm &  nhits >=10
  book_track_histset[ 15] = 1;		
  book_track_histset[ 16] = 1;		
  book_track_histset[ 17] = 1;		
  book_track_histset[ 18] = 1;		
  book_track_histset[ 19] = 1;		// track passedd trkQual

  book_track_histset[ 100] = 1;		// CalTrkFit e-
  book_track_histset[ 101] = 1;		// CalTrkFit e- with nhits > 15
  book_track_histset[ 102] = 1;		// CalTrkFit e- with chi2 < 3
  book_track_histset[ 103] = 1;		// CalTrkFit e- with chi2 < 4
  book_track_histset[ 104] = 1;		// CalTrkFit e- with chi2 < 3 & nhits >15
  book_track_histset[ 105] = 1;		// CalTrkFit e- with chi2 < 4 & nhits >15
  book_track_histset[ 106] = 1;		// CalTrkFit e- with chi2 < 3 & d0 in [-20, 20] cm
  book_track_histset[ 107] = 1;		// CalTrkFit e- with chi2 < 4 & d0 in [-20, 20] cm
  book_track_histset[ 108] = 1;		// CalTrkFit e- with chi2 < 3 & d0 in [-20, 20] cm &  nhits >15
  book_track_histset[ 109] = 1;		// CalTrkFit e- with chi2 < 4 & d0 in [-20, 20] cm &  nhits >15

  book_track_histset[ 110] = 1;		// CalTrkFit e- when track passedd trkQual

  book_track_histset[ 111] = 1;		// CalTrkFit e- when track passedd trkQual and TrackSeed with: nhits>=15 and chi2XY(ZPhi)<4
  book_track_histset[ 112] = 1;		// CalTrkFit e- when TrackSeed with: nhits>=15 and chi2XY(ZPhi)<4
  
  book_track_histset[ 113] = 1;		// CalTrkFit e- when track passedd trkQual and TrackSeed with: nhits>=15 and chi2XY(ZPhi)<4 and with CalTrkFit Track  with chi2 < 4 & d0 in [-20, 20] cm &  nhits >15
  book_track_histset[ 114] = 1;		// CalTrkFit e- when TrackSeed with: nhits>=15 and chi2XY(ZPhi)<4, and with CalTrkFit Track  with chi2 < 4 & d0 in [-20, 20] cm &  nhits >15

  book_track_histset[ 115] = 1;		// CalTrkFit e- when track has: N(active-hits)>=20 and chi2/ndof<5
  book_track_histset[ 116] = 1;		// CalTrkFit e- when track has: N(active-hits)>=20 and chi2/ndof<3
  book_track_histset[ 117] = 1;		// CalTrkFit e- when track has: N(active-hits)>=25 and chi2/ndof<5
  book_track_histset[ 118] = 1;		// CalTrkFit e- when track has: N(active-hits)>=25 and chi2/ndof<3
  
  book_track_histset[ 200] = 1;		// TrkPatRec e-
  book_track_histset[ 210] = 1;		// TrkPatRec passing quality cuts
  book_track_histset[ 215] = 1;		// TrkPAtRec e- when track has: N(active-hits)>=20 and chi2/ndof<5
  book_track_histset[ 216] = 1;		// TrkPAtRec e- when track has: N(active-hits)>=20 and chi2/ndof<3
  book_track_histset[ 217] = 1;		// TrkPAtRec e- when track has: N(active-hits)>=25 and chi2/ndof<5
  book_track_histset[ 218] = 1;		// TrkPAtRec e- when track has: N(active-hits)>=25 and chi2/ndof<3
 
  book_track_histset[ 300] = 1;		// Calpatrec e+
  book_track_histset[ 310] = 1;		// Calpatrec passing quality cuts
  book_track_histset[ 311] = 1;		// Calpatrec passing quality cuts
  book_track_histset[ 312] = 1;		// Calpatrec e+ when track has: N(active-hits)>=20 and chi2/ndof<5
  book_track_histset[ 313] = 1;		// Calpatrec passing quality cuts
  book_track_histset[ 314] = 1;		// Calpatrec e+ when track has: N(active-hits)>=20 and chi2/ndof<5
  book_track_histset[ 315] = 1;		// Calpatrec passing quality cuts
  book_track_histset[ 316] = 1;		// Calpatrec e+ when track has: N(active-hits)>=20 and chi2/ndof<3
  book_track_histset[ 317] = 1;		// Calpatrec e+ when track has: N(active-hits)>=25 and chi2/ndof<5
  book_track_histset[ 318] = 1;		// Calpatrec e+ when track has: N(active-hits)>=25 and chi2/ndof<3
 
  book_track_histset[ 400] = 1;		// TrkPatRec e+
  book_track_histset[ 410] = 1;		// TrkPatRec passing quality cuts
  book_track_histset[ 415] = 1;		// TrkPAtRec e+ when track has: N(active-hits)>=20 and chi2/ndof<5
  book_track_histset[ 416] = 1;		// TrkPAtRec e+ when track has: N(active-hits)>=20 and chi2/ndof<3
  book_track_histset[ 417] = 1;		// TrkPAtRec e+ when track has: N(active-hits)>=25 and chi2/ndof<5
  book_track_histset[ 418] = 1;		// TrkPAtRec e+ when track has: N(active-hits)>=25 and chi2/ndof<3
 

  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track_histset[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
      fBookedTracks[i] = 1;
    }else {
      fBookedTracks[i] = 0;
    }
    
  }
  //-----------------------------------------------------------------------------
  // book cluster histograms
  //-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[0] = 1;		// all clusters
  book_cluster_histset[1] = 1;		// clusters in events with the reconstructed e-
  book_cluster_histset[2] = 1;		// clusters in events with the track passing SetC cuts
  book_cluster_histset[3] = 1;		// clusters in events w/track passing SetC cuts and |dt|<2.5ns 
  book_cluster_histset[4] = 1;		// clusters > 10 MeV
  book_cluster_histset[5] = 1;		// clusters > 60 MeV
  book_cluster_histset[6] = 1;		// clusters disk#0
  book_cluster_histset[7] = 1;		// clusters disk#1

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
  //-----------------------------------------------------------------------------
  // book calorimeter histograms
  //-----------------------------------------------------------------------------
  int book_calo_histset[kNCaloHistSets];
  for (int i=0; i<kNCaloHistSets; i++) book_calo_histset[i] = 0;

  book_calo_histset[0] = 1;		// all crystals
  book_calo_histset[1] = 1;		// all crystals, e > 0
  book_calo_histset[2] = 1;		// all crystals, e > 0.1
  book_calo_histset[3] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNCaloHistSets; i++) {
    if (book_calo_histset[i] != 0) {
      sprintf(folder_name,"cal_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCalo[i] = new CaloHist_t;
      BookCaloHistograms(fHist.fCalo[i],Form("Hist/%s",folder_name));
    }
  }
  //-----------------------------------------------------------------------------
  // book Genp histograms
  //-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles
  book_genp_histset[1] = 1;		// electrons
  book_genp_histset[2] = 1;		// positrons
  book_genp_histset[3] = 1;		// mu minus
  book_genp_histset[4] = 1;		// mu plus
  book_genp_histset[5] = 1;		// pi minus
  book_genp_histset[6] = 1;		// pi plus
  book_genp_histset[7] = 1;		// kaon plus
  book_genp_histset[8] = 1;		// kaon minus
  book_genp_histset[9] = 1;		// photons
 

  //   book_genp_histset[1] = 1;		// all crystals, e > 0
  //   book_genp_histset[2] = 1;		// all crystals, e > 0.1
  //   book_genp_histset[3] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
      fBookedGenps[i] = 1;
    }else {
      fBookedGenps[i] = 0;
    }
  }

}

//--------------------------------------------------------------------------------
// function to fill TimeCluster block
//--------------------------------------------------------------------------------
void TTrigAna001Module::FillTimeClusterHistograms(TimeClusterHist_t*   Hist, TStnTimeCluster*    TimeCluster){
  
  int         nhits      = TimeCluster->NHits      ();
  int         ncombohits = TimeCluster->NComboHits ();
  double      time       = TimeCluster->T0();
  double      clusterE   = TimeCluster->ClusterEnergy();
  
  Hist->fNHits         ->Fill(nhits);	 
  Hist->fNComboHits    ->Fill(ncombohits);	 
  Hist->fT0            ->Fill(time);
  Hist->fClusterEnergy ->Fill(clusterE);
}

//--------------------------------------------------------------------------------
// function to fill Helix block
//--------------------------------------------------------------------------------
void TTrigAna001Module::FillHelixHistograms(HelixHist_t*   Hist, TStnHelix*    Helix){
  
  int         nhits    = Helix->NHits      ();
  double      clusterT = Helix->ClusterTime();
  double      clusterE = Helix->ClusterEnergy();
  
  double      radius   = Helix->Radius();

  double      lambda   = Helix->Lambda();  
  double      tanDip   = lambda/radius;
  double      mm2MeV   = 3/10.;
  double      pT       = radius*mm2MeV;
  double      p        = pT/std::cos( std::atan(tanDip));

  double      p1MC     = Helix->Mom1().Vect().Mag();
  double      pT1MC    = Helix->Mom1().Vect().Perp();
  double      x1Origin = Helix->Origin1().Vect().x()+3904;
  double      y1Origin = Helix->Origin1().Vect().y();
  double      z1Origin = Helix->Origin1().Vect().z();
  double      r1Origin = sqrt(x1Origin*x1Origin + y1Origin*y1Origin);
  double      deltaP1  = p - p1MC;
  double      pzMC1    = Helix->Mom1().Vect().Z();

  // 2nd most frequent particle origin data
  double      p2MC     = Helix->Mom2().Vect().Mag();
  double      pT2MC    = Helix->Mom2().Vect().Perp();
  double      x2Origin = Helix->Origin2().Vect().x()+3904;
  double      y2Origin = Helix->Origin2().Vect().y();
  double      z2Origin = Helix->Origin2().Vect().z();
  double      r2Origin = sqrt(x2Origin*x2Origin + y2Origin*y2Origin);
  double      deltaP2  = p - p2MC;

  // if (GetDebugBit(43) && (nhits > 55)) {
  //   GetHeaderBlock()->Print(Form(" bit:043 nHits = %i p = %2.3f E = %2.3f chi2-XY = %2.3f chi2-ZPhi = %2.3f", 
  // 				 nhits, p, clusterE, Helix->Chi2XY(), Helix->Chi2ZPhi()));
  // }
  Hist->fNHits         ->Fill(nhits);	 
  Hist->fClusterTime   ->Fill(clusterT);
  Hist->fClusterEnergy ->Fill(clusterE);
  
  Hist->fRadius        ->Fill(radius);    
  Hist->fMom           ->Fill(p);	 
  Hist->fPt            ->Fill(pT);	 
  Hist->fLambda        ->Fill(fabs(lambda));    
  		       
  Hist->fT0            ->Fill(Helix->T0());
  Hist->fT0Err         ->Fill(Helix->T0Err());
  if (fTrackT0 >  0.) Hist->fDT0->Fill(Helix->T0() - fTrackT0);
  if (fTrackCaloE  > 10. && clusterE > 10.) Hist->fDE->Fill(clusterE - fTrackCaloE);

  Hist->fAlg           ->Fill(Helix->AlgMask());
  Hist->fBestAlg       ->Fill(Helix->BestAlg());
  Hist->fD0            ->Fill(Helix->D0());
  Hist->fChi2XY        ->Fill(Helix->Chi2XY());
  Hist->fPhiZ          ->Fill(Helix->Chi2ZPhi());

  Hist->fNHitsVsChi2PhiZ ->Fill(nhits, Helix->Chi2ZPhi());

  //MC info
  Hist->fPDGMother1    ->Fill(Helix->PDGMother1());
  Hist->fPDG1          ->Fill(Helix->PDG1());
  Hist->fFrNhits1      ->Fill((double)Helix->ComboHitsFrom1()/Helix->NComboHits());
  Hist->fp1MC          ->Fill(p1MC);
  Hist->fz1Origin      ->Fill(z1Origin);
  Hist->fr1Origin      ->Fill(r1Origin);
  Hist->fdeltaP1       ->Fill(deltaP1);
  Hist->fpT1MC         ->Fill(pT1MC);
  Hist->fpzMC1         ->Fill(pzMC1);

  Hist->fPDGMother2    ->Fill(Helix->PDGMother2());
  Hist->fPDG2          ->Fill(Helix->PDG2());
  Hist->fFrNhits2      ->Fill((double)Helix->ComboHitsFrom2()/Helix->NComboHits());
  Hist->fz2Origin      ->Fill(z2Origin);
  Hist->fr2Origin      ->Fill(r2Origin);
  Hist->fdeltaP2       ->Fill(deltaP2);
  Hist->fpT2MC         ->Fill(pT2MC);

  
}

//--------------------------------------------------------------------------------
// function to fill TrasckSeedHit block
//--------------------------------------------------------------------------------
void TTrigAna001Module::FillTrackSeedHistograms(TrackSeedHist_t*   Hist, TStnTrackSeed*    TrkSeed){
  
  int         nhits    = TrkSeed->NHits      ();
  double      clusterT = TrkSeed->ClusterTime();
  double      clusterE = TrkSeed->ClusterEnergy();
  
  double      mm2MeV   = 3/10.;
  double      pT       = TrkSeed->Pt();
  double      p        = TrkSeed->P();
  double      radius   = pT/mm2MeV;

  double      tanDip   = TrkSeed->TanDip();  
  double      chi2d    = TrkSeed->Chi2()/(nhits - 5.);

  double      p1MC     = TrkSeed->Mom1().Vect().Mag();
  double      pT1MC    = TrkSeed->Mom1().Vect().Perp();
  double      x1Origin = TrkSeed->Origin1().Vect().x()+3904;
  double      y1Origin = TrkSeed->Origin1().Vect().y();
  double      z1Origin = TrkSeed->Origin1().Vect().z();
  double      r1Origin = sqrt(x1Origin*x1Origin + y1Origin*y1Origin);
  double      deltaP1  = p - p1MC;
  double      pzMC1    = TrkSeed->Mom1().Vect().Z();

  // 2nd most frequent particle origin data
  double      p2MC     = TrkSeed->Mom2().Vect().Mag();
  double      pT2MC    = TrkSeed->Mom2().Vect().Perp();
  double      x2Origin = TrkSeed->Origin2().Vect().x()+3904;
  double      y2Origin = TrkSeed->Origin2().Vect().y();
  double      z2Origin = TrkSeed->Origin2().Vect().z();
  double      r2Origin = sqrt(x2Origin*x2Origin + y2Origin*y2Origin);
  double      deltaP2  = p - p2MC;

  Hist->fNHits      ->Fill(nhits);	 
  Hist->fClusterTime->Fill(clusterT);
  Hist->fClusterEnergy->Fill(clusterE);
  
  Hist->fRadius     ->Fill(radius);    
  Hist->fMom        ->Fill(p);	 
  Hist->fPt         ->Fill(pT);	 
  Hist->fTanDip     ->Fill(fabs(tanDip));    
  
  Hist->fChi2       ->Fill(chi2d);
  Hist->fFitCons    ->Fill(TrkSeed->FitCons());
  Hist->fD0         ->Fill(TrkSeed->D0());

  //MC info
  Hist->fPDGMother1    ->Fill(TrkSeed->PDGMother1());
  Hist->fPDG1          ->Fill(TrkSeed->PDG1());
  Hist->fFrNhits1      ->Fill((double)TrkSeed->NHitsFrom1()/TrkSeed->NHits());
  Hist->fp1MC          ->Fill(p1MC);
  Hist->fz1Origin      ->Fill(z1Origin);
  Hist->fr1Origin      ->Fill(r1Origin);
  Hist->fdeltaP1       ->Fill(deltaP1);
  Hist->fpT1MC         ->Fill(pT1MC);
  Hist->fpzMC1         ->Fill(pzMC1);

  Hist->fPDGMother2    ->Fill(TrkSeed->PDGMother2());
  Hist->fPDG2          ->Fill(TrkSeed->PDG2());
  Hist->fFrNhits2      ->Fill((double)TrkSeed->NHitsFrom2()/TrkSeed->NHits());
  Hist->fz2Origin      ->Fill(z2Origin);
  Hist->fr2Origin      ->Fill(r2Origin);
  Hist->fdeltaP2       ->Fill(deltaP2);
  Hist->fpT2MC         ->Fill(pT2MC);

}


//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TTrigAna001Module::FillEventHistograms(EventHist_t* Hist) {

  //follow the loop pver the tracks from the calo-seeded algoriothm
  //we will evaluate the efficiency vs the E_calo threshold for electrons and positrons

  //Electrons
  TStnTrack*     trk;
  bool           demTracksFromTpr(false);
  
  for (int i=0; i<fNTracks[2]; i++ ) {
    trk = fTprDemTrackBlock->Track(i);
    double   chi2d = trk->Chi2Dof();
    double p = trk->P();
    if (chi2d > 8.)        continue;
    if (p<100.)            continue;
    demTracksFromTpr = true;
  }

  for (int i=0; i<fNTracks[1]; i++ ) {
    trk            = fCprDemTrackBlock->Track(i);
    double   chi2d = trk->Chi2Dof();
    double p = trk->P();
    if (chi2d > 8.)        continue;
    if (p<100.)            continue;

    double energy  = trk->ClusterE ();
    double binwidth = Hist->fThreshEffDem->GetBinWidth(1);
    double init = Hist->fThreshEffDem->GetBinLowEdge(1);
    int nbins = Hist->fThreshEffDem->GetNbinsX();
    for (int j=0; j<nbins; j++){
      if (energy>=(init+j*binwidth)){
	Hist->fThreshEffDem->Fill(init+j*binwidth);
	if (demTracksFromTpr==false){
	  Hist->fTPRnoneDem->Fill(init+j*binwidth);
	} 
      }
      if (energy>=(init+j*binwidth) || demTracksFromTpr==true){
	Hist->fcalThreshNormDem->Fill(init+j*binwidth);
      }
    }
  }

  bool           depTracksFromTpr(false);
  //Positrons
  for (int i=0; i<fNTracks[4]; i++ ) {
    trk = fTprDepTrackBlock->Track(i);
    double   chi2d = trk->Chi2Dof();
    double p = trk->P();
    if (chi2d > 8.)        continue;
    if (p<85.)            continue;
    depTracksFromTpr = true;
  }


  for (int i=0; i<fNTracks[3]; i++ ) {
    trk            = fCprDepTrackBlock->Track(i);
    double   chi2d = trk->Chi2Dof();
    double p = trk->P();
    if (chi2d > 8.)        continue;
    if (p<85.)            continue;

    double energy  = trk->ClusterE ();
    double binwidth = Hist->fThreshEffDep->GetBinWidth(1);
    double init = Hist->fThreshEffDep->GetBinLowEdge(1);
    int nbins = Hist->fThreshEffDep->GetNbinsX();
    for (int j=0; j<nbins; j++){
      if (energy>=(init+j*binwidth)){
	Hist->fThreshEffDep->Fill(init+j*binwidth);
	if (depTracksFromTpr==false){
	  Hist->fTPRnoneDep->Fill(init+j*binwidth);
	} 
      }
      if (energy>=(init+j*binwidth) || depTracksFromTpr==true){
	Hist->fcalThreshNormDep->Fill(init+j*binwidth);
      }
    }
  }


  Hist->fInstLum ->Fill(GetHeaderBlock()->InstLum());

  double            cos_th, xv, yv, rv, zv, p;
  //  double            e, m, r;
  TLorentzVector    mom;

  if (fParticle != NULL){
    fParticle->Momentum(mom);
  }
  
  p      = mom.P();

  cos_th = mom.Pz()/p;

  if (fParticle != NULL){
    xv = fParticle->Vx()+3904.;
    yv = fParticle->Vy();
    rv = sqrt(xv*xv+yv*yv);
    zv = fParticle->Vz();
  } else {
    xv = -9999.;
    yv = -9999.;
    rv = -9999.;
    zv = -9999.;
  }

  Hist->fEleMom->Fill(p);
  Hist->fNormMom->Fill(p);
  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  Hist->fNClusters[0]->Fill(fNClusters);
  
  double       edep_min(50);
  int          nCl_trig(0);
  TStnCluster* cluster(0);
  for (int i=0; i<fNClusters; ++i) {
    cluster = fClusterBlock->Cluster(i);
    double edep = cluster->Energy();
    if (edep >= edep_min) ++nCl_trig;
  }
  Hist->fNClusters[1]->Fill(nCl_trig);

  Hist->fNTracks->Fill  (fNTracks[0]);
  Hist->fNStrawHits[0]->Fill(fNStrawHits);
  Hist->fNStrawHits[1]->Fill(fNStrawHits);


  TStnTrack* track(0);
  int        hasQualityTrack(0),hasQualityTprTrack(0) ;
  int        hasQualityTrackDep(0),hasQualityTprTrackDep(0) ;

  //  double  fMinFitCons      = 2.e-3;
  double  fMinNActive      = 15;//25
  //  double  fMaxT0Err        = 0.9;  		// in ns
  // double  fMaxFitMomErr    = 0.25;  		// in MeV
  // double  fMinTanDip       = tan(M_PI/6.);	// 0.5773
  // double  fMaxTanDip       = 1.0;  
  // double  fMinD1           = -80.;		// in mm
  // double  fMaxD1           = 105.;
  // double  fMinT0           = 700.; // ns
  double  fChi2Max         = 5.;
  //  double  fD0Max           = 200.;
  double  fPMin            = 100.;
  double  fPMax            = 105.;
  double  fPMinDep         = 89.;
  double  fPMaxDep         = 94.;
  //  TStnTrackSeed*    trkSeed(0);

  int     tclFalgs[2]={0};

  for (int i=0; i<fNTracks[1]; ++i ) {
    track = fCprDemTrackBlock->Track(i);

    fTrackTrkSeedIndex = track->TrackSeedIndex();
    
    double chi2dof    = track->Chi2Dof();
    double nactive    = track->NActive();
    double p          = track->P();
    double dpf        = track->fP - track->fPFront;
    //assume fTrackBlock has only CalPatRec tracks
    if ( (chi2dof < fChi2Max   ) &&  
	 // (tan_dip < fMaxTanDip ) &&
	 // (tan_dip > fMinTanDip ) &&
	 (nactive >= fMinNActive) &&
	 // (fitmom_err< fMaxFitMomErr) &&
	 // (d0 < fMaxD1)           &&
	 // (d0 > fMinD1)           && 
	 // (t0 > fMinT0 )          &&
	 //	 ( fabs(d0) <= fD0Max)   &&
	 ((p >= fPMin) && (p <= fPMax)) && (fabs(dpf)<10.)
	 ){
      hasQualityTrack = 1;
    }

    if ( (nactive>=20) && (chi2dof<5)) tclFalgs[0] = 1;
    if ( (nactive>=25) && (chi2dof<5)) tclFalgs[0] = 1;
  }

  for (int i=0; i<fNTracks[3]; ++i ) {
    track = fCprDepTrackBlock->Track(i);

    double chi2dof    = track->Chi2Dof();
    double nactive    = track->NActive();
    double p          = track->P();
    double dpf        = track->fP - track->fPFront;

    //assume fTrackBlock has only CalPatRec tracks
    if ( (chi2dof < fChi2Max   ) &&  
	 // (tan_dip < fMaxTanDip ) &&
	 // (tan_dip > fMinTanDip ) &&
	 (nactive >= fMinNActive) &&
	 // (fitmom_err< fMaxFitMomErr) &&
	 // (d0 < fMaxD1)           &&
	 // (d0 > fMinD1)           && 
	 // (t0 > fMinT0 )          &&
	 //	 ( fabs(d0) <= fD0Max)   &&
	 ((p >= fPMinDep) && (p <= fPMaxDep)) && (fabs(dpf)<10.)
	 ){
      hasQualityTrackDep = 1;
    }
  }

  for (int i=0; i<fNTracks[2]; ++i ) {
    track = fTprDemTrackBlock->Track(i);

    double chi2dof    = track->Chi2Dof();
    double nactive    = track->NActive();
    double p          = track->P();
    double dpf        = track->fP - track->fPFront;

    //assume fTrackBlock has only CalPatRec tracks
    if ( (chi2dof < fChi2Max   ) &&  
	 // (tan_dip < fMaxTanDip ) &&
	 // (tan_dip > fMinTanDip ) &&
	 (nactive >= fMinNActive) &&
	 // (fitmom_err< fMaxFitMomErr) &&
	 // (d0 < fMaxD1)           &&
	 // (d0 > fMinD1)           && 
	 // (t0 > fMinT0 )          &&
	 //	 ( fabs(d0) <= fD0Max)   &&
	 ((p >= fPMin) && (p <= fPMax))  && (fabs(dpf)<10.)
	 ){
      hasQualityTprTrack = 1;
    }

    // if ( (nactive>=20) && (chi2dof<5)) tclFalgs[0] = 1;
    // if ( (nactive>=25) && (chi2dof<5)) tclFalgs[0] = 1;
  }

  for (int i=0; i<fNTracks[4]; ++i ) {
    track = fTprDepTrackBlock->Track(i);

    double chi2dof    = track->Chi2Dof();
    double nactive    = track->NActive();
    double p          = track->P();
    double dpf        = track->fP - track->fPFront;

    //assume fTrackBlock has only CalPatRec tracks
    if ( (chi2dof < fChi2Max   ) &&  
	 // (tan_dip < fMaxTanDip ) &&
	 // (tan_dip > fMinTanDip ) &&
	 (nactive >= fMinNActive) &&
	 // (fitmom_err< fMaxFitMomErr) &&
	 // (d0 < fMaxD1)           &&
	 // (d0 > fMinD1)           && 
	 // (t0 > fMinT0 )          &&
	 //	 ( fabs(d0) <= fD0Max)   &&
	 ((p >= fPMinDep) && (p <= fPMaxDep))  && (fabs(dpf)<10.)
	 ){
      hasQualityTprTrackDep = 1;
    }

    // if ( (nactive>=20) && (chi2dof<5)) tclFalgs[0] = 1;
    // if ( (nactive>=25) && (chi2dof<5)) tclFalgs[0] = 1;
  }


  //loop over the TimeClusters
  TStnTimeCluster* tcl(0);
  for (int i=0; i<fNTimeClusters [0]; ++i){
    tcl = fTimeClusterBlock->TimeCluster(i);
    for (int j=0; j<3; ++j){
      if (tclFalgs[j] == 1) Hist->fNTClHits[j]->Fill(tcl->NHits());
    }
  }


  FillTprDemTriggerHist(Hist, hasQualityTprTrack);
  FillTprDepTriggerHist(Hist, hasQualityTprTrackDep);


  FillCprDemTriggerHist(Hist, hasQualityTrack);
  FillCprDepTriggerHist(Hist, hasQualityTrackDep);

  //--------------------------------------------------------------------------------//

  //compare the the number of hits associated to the helices from TrkPatRec and CalPatRec
  int nhelices_tpr = fNHelices[0];
  int nhelices_cpr = fNHelices[1];
  TStnHelix* tpr_helix(0), *cpr_helix(0);
  
  if ( (nhelices_tpr > 0) && 
       (nhelices_cpr > 0)   ) {
    tpr_helix = fTprDemHelixBlock->Helix(0);
    cpr_helix = fCprDemHelixBlock->Helix(0);
    double     dNHits = cpr_helix->NHits() -  tpr_helix->NHits();
    Hist->fDNHitsHelix->Fill(dNHits);

    if (GetDebugBit(42) && (dNHits < 0) ) {
      double     mm2MeV   = 3/10.;
      double     cpr_chi2xy   = cpr_helix->Chi2XY();				       
      double     cpr_chi2zphi = cpr_helix->Chi2ZPhi();				       
      double     tanDip       = cpr_helix->Lambda()/cpr_helix->Radius();		       
      double     cpr_mom      = cpr_helix->Radius()*mm2MeV/std::cos( std::atan(tanDip));
      int        cpr_nhits    = cpr_helix->NHits();                            

      double     tpr_chi2xy   = tpr_helix->Chi2XY();				       
      double     tpr_chi2zphi = tpr_helix->Chi2ZPhi();				       
      tanDip       = tpr_helix->Lambda()/tpr_helix->Radius();		       
      double     tpr_mom      = tpr_helix->Radius()*mm2MeV/std::cos( std::atan(tanDip));
      int        tpr_nhits    = tpr_helix->NHits();   
      
      GetHeaderBlock()->Print(Form(" bit:042 Cpr: chi2XY = %10.3f chi2ZPhi = %10.3f p = %10.3f  nhits = %5i; Tpr: chi2XY = %10.3f chi2ZPhi = %10.3f p = %10.3f nhits = %5i",
				   cpr_chi2xy, cpr_chi2zphi, cpr_mom, cpr_nhits,
				   tpr_chi2xy, tpr_chi2zphi, tpr_mom, tpr_nhits));
    }

  }


  double emax   = -1;
  double t0_cls = -1;
  double dt     = 9999.;
  
  if (fNClusters > 0) cluster = fClusterBlock->Cluster(0);

  //  TStnTrack* track(0);
  if (fNTracks[1] > 0) track = fCprDemTrackBlock->Track(0);

  if (cluster) {
    emax   = cluster->Energy();
    t0_cls = cluster->Time();
  }

  double t0_trk = -1;
  if (track) {
    t0_trk = track->fT0;
  }

  if (track && cluster) {
    dt = t0_cls-t0_trk;
  }

  Hist->fDtClT->Fill(dt);
  Hist->fEMax->Fill(emax);

  TStrawHitData*  sh;
  int n_good_hits = 0;
  for (int i=0; i<fNStrawHits; i++ ) {
    sh  = fStrawDataBlock->Hit(i);
    dt  = t0_cls - sh->Time() + 15;
    Hist->fDtClS->Fill(dt);
    Hist->fSHTime->Fill(sh->Time());

    if (fabs(dt+15.)< 50) n_good_hits += 1;
  }

  Hist->fNGoodSH->Fill(n_good_hits);

  Hist->fNHyp->Fill(fNHyp);
  Hist->fBestHyp[0]->Fill(fBestHyp[0]);
  Hist->fBestHyp[1]->Fill(fBestHyp[1]);
  Hist->fNGenp->Fill(fNGenp);
  //-----------------------------------------------------------------------------
  // crystals - count crystals with E > 1MeV
  //-----------------------------------------------------------------------------
  //   TCalHitData* cch;

  //   int n_cch_1mev = 0;

  //   if (fCalorimeterType == 1) {
  // //-----------------------------------------------------------------------------
  // // vane calorimeter
  // //-----------------------------------------------------------------------------
  //     int  nhits_vane[2][4], nhits_row [2][20], nhits_col[2][50];
  //     int  crystal_id, vane_id, local_id, vane_row, vane_col;

  //     for (int i=0; i<4; i++) {
  //       nhits_vane[0][i] = 0;
  //       nhits_vane[1][i] = 0;
  //     }
      
  //     for (int i=0; i<20; i++) {
  //       nhits_row[0][i] = 0;
  //       nhits_row[1][i] = 0;
  //     }

  //     for (int i=0; i<50; i++) {
  //       nhits_col[0][i] = 0;
  //       nhits_col[1][i] = 0;
  //     }
      
  //     for (int ic=0; ic<fNCalHits; ic++) {
  //       cch        = fCalDataBlock->CalHitData(ic);
  //       crystal_id = cch->ID();

  //       if (cch->Energy() > 1.) {
  // 	n_cch_1mev += 1;
  //       }
  //       // for each crystal determine its row and column
  //       // the following is for vanes
  //       vane_id  = crystal_id/484.;
  //       local_id = crystal_id-vane_id*484;
  //       vane_row = local_id/44;
  //       vane_col = local_id-vane_row*44;
      
  //       nhits_vane[0][vane_id ] += 1;
  //       nhits_row [0][vane_row] += 1;
  //       nhits_col [0][vane_col] += 1;
      
  //       if (cch->Energy() > 1.) {
  // 	nhits_row [1][vane_row] += 1;
  // 	nhits_col [1][vane_col] += 1;
  // 	nhits_vane[1][vane_id ] += 1;
  //       }
  //     }

  //     Hist->fNCaloCrystalHits[0]->Fill(fNCalHits);
  //     Hist->fNCaloCrystalHits[1]->Fill(n_cch_1mev);

  //     for (int iv=0; iv<4; iv++) {
  //       Hist->fNCaloHitsVsVane[0]->Fill(iv,nhits_vane[0][iv]);
  //       Hist->fNCaloHitsVsVane[1]->Fill(iv,nhits_vane[1][iv]);
  //     }

  //     for (int ir=0; ir<20; ir++) {
  //       Hist->fNCaloHitsVsRow[0]->Fill(ir,nhits_row[0][ir]);
  //       Hist->fNCaloHitsVsRow[1]->Fill(ir,nhits_row[1][ir]);
  //     }

  //     for (int ic=0; ic<50; ic++) {
  //       Hist->fNCaloHitsVsCol[0]->Fill(ic,nhits_col[0][ic]);
  //       Hist->fNCaloHitsVsCol[1]->Fill(ic,nhits_col[1][ic]);
  //     }
  //   }
  //   else if (fCalorimeterType == 2) {
  // //-----------------------------------------------------------------------------
  // // disk calorimeter
  // //-----------------------------------------------------------------------------
  //     int      ndisks, n_hit_crystals[4], n_hit_crystals_tot;
  //     double   etot[4];

  //     TCalHitData* hit;

  //     //    TDisk*       disk;
  //     //    TEvdCrystal* cr;

  //     ndisks = fDiskCalorimeter->NDisks();

  //     int   bin, hit_id, idisk, nhits, nhits_r[4][100], n_hit_crystals_r[4][100];

  //     for (int id=0; id<ndisks; id++) {
  //       n_hit_crystals[id] = 0;
  //       etot[id]           = 0;

  //       for (int ib=0; ib<100; ib++) {
  // 	nhits_r         [id][ib] = 0;
  // 	n_hit_crystals_r[id][ib] = 0;
  //       }
  //     }

  //     nhits = fCalDataBlock->NHits();

  //     for (int i=0; i< nhits; i++) {
  //       hit    = fCalDataBlock->CalHitData(i);

  //       hit_id = hit->ID();
  //       idisk  = fDiskCalorimeter->DiskNumber(hit_id);
  //       r      = fDiskCalorimeter->CrystalRadius(hit_id);
  //       e      = hit->Energy(); 

  //       etot          [idisk] += e;
  //       n_hit_crystals[idisk] += 1;

  //       Hist->fECrVsR[idisk]->Fill(r,e);
  //       Hist->fNCrVsR[idisk]->Fill(r,1);

  //       bin  = (int) (r/10.);

  //       nhits_r         [idisk][bin] += 1;
  // //-----------------------------------------------------------------------------
  // // this is not correct, one needs to check whether this crystal has been hit,
  // // for the moment, to get going, ignore that
  // //-----------------------------------------------------------------------------
  //       n_hit_crystals_r[idisk][bin] += 1;
  //     }

  //     n_hit_crystals_tot = 0;

  //     double ecal = 0;
  //     for (int id=0; id<ndisks; id++) {
  //       n_hit_crystals_tot += n_hit_crystals[id];
  //       ecal += etot[id];
  // //-----------------------------------------------------------------------------
  // // fill 'per-disk' histograms
  // //-----------------------------------------------------------------------------
  //       Hist->fETot[id]->Fill(etot[id]);

  // //-----------------------------------------------------------------------------
  // // 100 is fixed by the number of bins in the radial distributions
  // //-----------------------------------------------------------------------------
  //       for (int ib=0; ib<100; ib++) {
  // 	r = (ib+0.5)*10.;
  // 	Hist->fNCrystalHitsVsR[id]->Fill(r,nhits_r         [id][ib]);
  // 	Hist->fNHitCrystalsVsR[id]->Fill(r,n_hit_crystals_r[id][ib]);
  //       }
  //     }

  //     Hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
  //     Hist->fECal->Fill(ecal);

  //     double ekin(-1.);
  // //-----------------------------------------------------------------------------
  // // there is an inconsistency in the SIMP block filling - in Mu2e offline 
  // // the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
  // // thus the energy is screwed up... kludge around
  // // assign muon mass
  // //-----------------------------------------------------------------------------
  //     if (fSimp) {
  //       p    = fSimp->fStartMom.P();
  //       m    = 105.658; // in MeV
  //       ekin = sqrt(p*p+m*m)-m;
  //     }
  //     Hist->fECalOverEKin->Fill(ecal/ekin);
  //   }
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::FillCaloHistograms(CaloHist_t* Hist, TStnCrystal* Cr) {

  int                    nhits;
  float                  t, e, r, e700, n700;
  TCalHitData*           hit;

  // determine crystal coordinates

  TDisk* disk = Cr->Disk();

  int idisk = disk->SectionID();
  // time needs to be defiend
  //    t  = Cr->Time();
  e     = Cr->Energy();
  r     = Cr->Radius();
  nhits = Cr->NHits();

  Hist->fVaneID->Fill(idisk);

  Hist->fEnergy  [idisk]->Fill(e);
  Hist->fNHits   [idisk]->Fill(nhits);
  //    Hist->fTime    [idisk]->Fill(t);
  Hist->fRadius  [idisk]->Fill(r);
  Hist->fRadiusWE[idisk]->Fill(r,e);
    
  e700 = 0;
  n700 = 0;
  for (int i=0; i<nhits; i++) {
    hit  = Cr->CalHitData(i);
    t   = hit->Time();
    Hist->fTime[idisk]->Fill(t);
    if (t > 700.) {
      n700 += 1;
      e700 += hit->Energy();
      Hist->fT700[idisk]->Fill(t);
    }
  }

  Hist->fE700   [idisk]->Fill(e700);
  Hist->fN700   [idisk]->Fill(n700);

  if (n700 > 0) {
    Hist->fR700  [idisk]->Fill(r);
    Hist->fRWE700[idisk]->Fill(r,e700);
  }
}


//-----------------------------------------------------------------------------
void TTrigAna001Module::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX+3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  Hist->fVaneID->Fill(Cluster->DiskID());
  Hist->fEnergy->Fill(Cluster->Energy());
  Hist->fT0->Fill(Cluster->Time());
  Hist->fRow->Fill(row);
  Hist->fCol->Fill(col);
  Hist->fX->Fill(x);
  Hist->fY->Fill(y);
  Hist->fZ->Fill(z);
  Hist->fR->Fill(r);

  Hist->fYMean->Fill(Cluster->fYMean);
  Hist->fZMean->Fill(Cluster->fZMean);
  Hist->fSigY->Fill(Cluster->fSigY);
  Hist->fSigZ->Fill(Cluster->fSigZ);
  Hist->fSigR->Fill(Cluster->fSigR);
  Hist->fNCr0->Fill(Cluster->fNCrystals);
  Hist->fNCr1->Fill(Cluster->fNCr1);
  Hist->fFrE1->Fill(Cluster->fFrE1);
  Hist->fFrE2->Fill(Cluster->fFrE2);
  Hist->fSigE1->Fill(Cluster->fSigE1);
  Hist->fSigE2->Fill(Cluster->fSigE2);
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
  int    gen_id;
  float  p, cos_th, z0, t0, r0, x0, y0;

  TLorentzVector mom, v;

  Genp->Momentum(mom);
  //  Genp->ProductionVertex(v);

  p      = mom.P();
  cos_th = mom.CosTheta();

  x0     = Genp->Vx()+3904.;
  y0     = Genp->Vy();

  z0     = Genp->Vz();
  t0     = Genp->T();
  r0     = sqrt(x0*x0+y0*y0);
  gen_id = Genp->GetStatusCode();

  Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
  Hist->fGenID->Fill(gen_id);
  Hist->fZ0->Fill(z0);
  Hist->fT0->Fill(t0);
  Hist->fR0->Fill(r0);
  Hist->fP->Fill(p);
  Hist->fCosTh->Fill(cos_th);
  Hist->fPNorm->Fill(p);
  Hist->fCosThNorm->Fill(cos_th);
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TTrigAna001Module::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

  TLorentzVector  mom;
  double          chi2c, r;
  //  int             itrk;
  //  TrackPar_t*     tp;
  // pointer to local track parameters
  // itrk = Track->Number();
  // tp   = fTrackPar+itrk;

  VirtualPar_t*   vp;
  int vDetId0, vDetId1;
  double dt0(-1e10);
  TStnTrack::InterData_t*    vt = Track->fVMaxEp;//fVMinS;
  bool vHitsFilled(false);
  double dx,dy,dz,dt,dr, tof;
  double dtOffset = 0.;//1.4;
  double deBuncherPeriod = 1695.;
  double tmp;
  double vdetT0(-9999.);
  
  double dpf = Track->fP - Track->fPFront;

  if(vt==NULL) goto NEXT_STEP;

  if(vt->fID == 0) {
    vDetId0 = 73; 
    vDetId1 = 77;
  }else {
    vDetId0 = 75; 
    vDetId1 = 79;
  }

  for (int i=0; i<fNVirtualHits; ++i) {
    vp  = fVirtualPar+i;
    
    if( (vp->fIndex == 11) ||
	(vp->fIndex == 12) ){         //middle of the tracker
      vdetT0 = vp->fTime;
      dt0    = vdetT0 - Track->fT0 + dtOffset;

      double t0Pull = dt0/Track->T0Err();
      Hist->fvDt0Norm->Fill (dt0);
      Hist->fvDt0    ->Fill (dt0);
      Hist->fT0Pull  ->Fill (t0Pull);
    }
    if ( ( (vp->fIndex == vDetId0) || (vp->fIndex == vDetId1) ) &&
	 ( !vHitsFilled ) ) {
      vHitsFilled = true;
      dx  = vp->fPosX - vt->fXTrk;
      dy  = vp->fPosY - vt->fYTrk;
      dz  = vp->fPosZ - vt->fZTrk;
      dt  = vp->fTime - vt->fTime + dtOffset;
      tmp = vp->fTime;
      while (tmp >= deBuncherPeriod) {
	dt  -= deBuncherPeriod;
	tmp -= deBuncherPeriod;
      }
      dr  = sqrt(dx*dx + dy*dy);
      tof = vp->fTime - Track->fT0 + dtOffset;
      Hist->fvDx->Fill  (dx);
      Hist->fvDy->Fill  (dy);
      Hist->fvDz->Fill  (dz);
      Hist->fvDt->Fill  (dt);
      Hist->fvDp->Fill  (Track->fP - Track->fPFront);
      Hist->fvtof->Fill(tof);
      Hist->fvDr->Fill  (dr);
      
      Hist->fvDtVsDt0->Fill  (dt, dt0);
      Hist->fvDtVsT0err->Fill(dt, Track->T0Err());
      Hist->fvDtVsDr->Fill   (dt, dr);
      Hist->fvDtVsDp->Fill   (dt, dpf); //tp->fDpF);
      Hist->fvtofVsDp->Fill   (tof, dpf); //tp->fDpF);

      //now fill with the hist to be normalized
      Hist->fvDxNorm->Fill  (dx);
      Hist->fvDyNorm->Fill  (dy);
      Hist->fvDzNorm->Fill  (dz);
      Hist->fvDtNorm->Fill  (dt);
      Hist->fvDpNorm->Fill  (dpf); //Track->fP - Track->fPFront);
      Hist->fvtofNorm->Fill(tof);
      Hist->fvDrNorm->Fill  (dr);
      
      Hist->fvDtVsDt0Norm->Fill  (dt, dt0);
      Hist->fvDtVsT0errNorm->Fill(dt, Track->T0Err());
      Hist->fvDtVsDrNorm->Fill   (dt, dr);
      Hist->fvDtVsDpNorm->Fill   (dt, dpf); //tp->fDpF);
      Hist->fvtofVsDpNorm->Fill   (tof, dpf); //tp->fDpF);
    }
  }
 NEXT_STEP:;
  Hist->fP[0]->Fill (Track->fP);
  Hist->fP[1]->Fill (Track->fP);
  Hist->fP[2]->Fill (Track->fP);
  Hist->fP0->  Fill (Track->fP0);
  Hist->fP2->  Fill (Track->fP2);

  Hist->fPNorm->Fill(Track->fP);

  Hist->fFitMomErr->Fill(Track->fFitMomErr);
  Hist->fFitMomErrNorm->Fill(Track->fFitMomErr);

  Hist->fPt    ->Fill(Track->fPt    );
  Hist->fPtNorm ->Fill(Track->fPt);
  Hist->fPFront->Fill(Track->fPFront);
  Hist->fPStOut->Fill(Track->fPStOut);
  // dp: Tracker-only resolution

  Hist->fDpFront ->Fill(dpf); //tp->fDpF);
  Hist->fDpFront0->Fill(Track->fP0    - Track->fPFront);//tp->fDp0);
  Hist->fDpFront2->Fill(Track->fP2    - Track->fPFront);//tp->fDp2);
  Hist->fDpFSt   ->Fill(Track->fPFront - Track->fPStOut);//tp->fDpFSt);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1,dpf); //tp->fDpF);

  double   nActive = Track->NActive();
  
  Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  Hist->fChi2->Fill (Track->fChi2);
  Hist->fChi2Norm->Fill (Track->fChi2);
  Hist->fNDof->Fill(nActive-5.);
  Hist->fChi2Dof->Fill(Track->fChi2/(nActive-5.));
  Hist->fNActive->Fill(nActive);
  Hist->fT0->Fill(Track->fT0);
  Hist->fT0Err->Fill(Track->fT0Err);
  Hist->fT0Norm->Fill(Track->fT0);
  Hist->fT0ErrNorm->Fill(Track->fT0Err);
  //  printf("TTrigAna001Module::FillTrackHistograms: track charge is not defined yet\n");
  Hist->fQ->Fill(-1);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);
  Hist->fFitConsNorm[0]->Fill(Track->fFitCons);
  Hist->fFitConsNorm[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);
  Hist->fTanDipNorm->Fill(Track->fTanDip);
  Hist->fAlgMask->Fill(Track->AlgMask());

  int          bestAlg = Track->BestAlg();
  Hist->fBestAlg->Fill(bestAlg);

  double       trkQual = Track->fTrkQual;
  if (bestAlg ==0 ){
    Hist->fTrkQual[0]->Fill(trkQual);
  }else if (bestAlg ==1 ){
    Hist->fTrkQual[1]->Fill(trkQual);
  }
  

  chi2c = Track->Chi2Dof();//fChi2C/(Track->fNActive-5.);
  Hist->fChi2DofC->Fill(chi2c);

  //  int nh, nst_with_nh[10];
  // 2014-04-29: currently not saved

  //  for (int i=0; i<10; i++) nst_with_nh[i] = 0;

  //   for (int i=0; i<40; i++) {
  //     Hist->fNHVsStation->Fill(i,Track->fNHPerStation[i]);
  //     nh = Track->fNHPerStation[i];
  //     if ((nh >= 0) && (nh < 10)) {
  //       nst_with_nh[nh] += 1;
  //     }
  //     else {
  //       printf(">>> ERROR : nh = %20i, IGNORE \n",nh);
  //     }
  //  }

  //   for (int i=0; i<10; i++) {
  //     Hist->fNHVsNSt->Fill(i,nst_with_nh[i]);
  //   }
  //-----------------------------------------------------------------------------
  // track-cluster matching part: 
  // - for residuals, determine intersection with the most energetic cluster
  // - for track -only parameters use intersection with lowest trajectory length
  //-----------------------------------------------------------------------------
  //  TStnTrack::InterData_t*    vt = Track->fVMinS;  // track-only
  //  TStnTrack::InterData_t*    vr = Track->fVMaxEp; // residuals

  if (vt) {
    Hist->fVaneID->Fill(vt->fID  );
    Hist->fXTrk->Fill  (vt->fXTrk);
    Hist->fYTrk->Fill  (vt->fYTrk);

    r = sqrt(vt->fXTrk*vt->fXTrk+vt->fYTrk*vt->fYTrk);
    Hist->fRTrk->Fill  (r);

    Hist->fZTrk->Fill  (vt->fZTrk);
  }
  else {
    //-----------------------------------------------------------------------------
    // fill histograms with numbers easy to recognize as dummy
    //-----------------------------------------------------------------------------
    Hist->fVaneID->Fill(-1.);
    Hist->fXTrk->Fill  (999.);
    Hist->fYTrk->Fill  (999.);
    Hist->fRTrk->Fill  (999.);
    Hist->fZTrk->Fill  (-1. );
  }

  //-----------------------------------------------------------------------------
  // there is an inconsistency in the SIMP block filling - in Mu2e offline 
  // the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
  // thus the energy is screwed up... kludge around
  // assign muon mass
  //-----------------------------------------------------------------------------
  double ekin(-1.);
  if (fSimp) {
    double p, m;
    //    p    = fSimp->fStartMom.P();
    p = Track->fP;
    m    = 105.658; // in MeV
    ekin = sqrt(p*p+m*m)-m;
  }
  
  double   ep  = vt != NULL ? vt->fEnergy/Track->fP : -1.;
  double   ecl = vt != NULL ? vt->fEnergy : -1.;
  double   vt_dx  = vt != NULL ? vt->fDx : -9999.;
  double   vt_dy  = vt != NULL ? vt->fDy : -9999.;
  double   vt_dz  = vt != NULL ? vt->fDz : -9999.;
  double   vt_dt  = vt != NULL ? vt->fDt : -9999.;

  double   chi2m  = vt != NULL ? vt->fChi2Match : -1.;
  double   path   = vt != NULL ? vt->fPath : -1.;

  Hist->fECl    ->Fill(ecl);
  Hist->fEClEKin->Fill(ecl/ekin);
  Hist->fEp[0]  ->Fill(ep);//tp->fEp);
  Hist->fEp[1]  ->Fill(ep);//tp->fEp);
  
  
  Hist->fEpVsPath->Fill(path, ep);

  Hist->fDx->Fill(vt_dx);
  Hist->fDy->Fill(vt_dy);
  Hist->fDz->Fill(vt_dz);

  Hist->fDt->Fill(vt_dt);
  Hist->fChi2Match->Fill(chi2m);

  double nx  = 0.;
  double ny  = 0.;
  double du  = 0.;
  double dv  = 0.;

  if (vt){
    nx  = vt->fNxTrk/sqrt(vt->fNxTrk*vt->fNxTrk+vt->fNyTrk*vt->fNyTrk);
    ny  = vt->fNyTrk/sqrt(vt->fNxTrk*vt->fNxTrk+vt->fNyTrk*vt->fNyTrk);
    du  = vt->fDx*nx+vt->fDy*ny;
    dv  = vt->fDx*ny-vt->fDy*nx;
  }
  
  Hist->fDu->Fill    (du);
  Hist->fDv->Fill    (dv);
  Hist->fDvVsDu->Fill(du, dv);

  Hist->fPath->Fill(path);
  Hist->fDuVsPath->Fill(path, du);//tp->fDu);

  //  double cu[4] = { -60.3049, -0.749111,    0.00522242,  -7.52018e-06};
  double cu[4] = { -59.5174, -0.541226, 0.00414309, -5.84989e-06 };
  double cv[4] = {  6.44161, 0.0722353, -0.000653084, 1.14054e-06};

  double x = path;//vt->fPath;
  double corr_u = cu[0]+cu[1]*x+cu[2]*x*x+cu[3]*x*x*x;
  double corr_v = cv[0]+cv[1]*x+cv[2]*x*x+cv[3]*x*x*x;

  //  double duc = tp->fDu-0.34*(tp->fPath-350.);
  double duc = du-corr_u;
  double dvc = dv-corr_v;
  
  Hist->fDucVsPath->Fill(path,duc);
  Hist->fDvcVsPath->Fill(path,dvc);

  Hist->fDvVsPath->Fill(path, dv);
  Hist->fDtVsPath->Fill(path, vt_dt);

  Hist->fDuVsTDip->Fill(Track->fTanDip, du);
  Hist->fDvVsTDip->Fill(Track->fTanDip, dv);

  Hist->fZ1->Fill(Track->fZ1);

  int ncl = Track->NClusters();
  Hist->fNClusters->Fill(ncl);

  Hist->fRSlope->Fill(Track->RSlope());
  Hist->fXSlope->Fill(Track->XSlope());

  double llhr_dedx, llhr_xs, llhr_cal, llhr_trk, llhr;

  Hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
  Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());

  llhr_cal = Track->LogLHRCal();
  Hist->fLogLHRCal->Fill(llhr_cal);

  llhr_dedx = Track->LogLHRDeDx();
  llhr_xs   = Track->LogLHRXs();
  llhr_trk  = Track->LogLHRTrk();
  llhr      = llhr_cal+llhr_trk;

  Hist->fEpVsDt->Fill(vt_dt, ep);
  Hist->fLogLHRDeDx->Fill(llhr_dedx);
  Hist->fLogLHRXs->Fill(llhr_xs);
  Hist->fLogLHRTrk->Fill(llhr_trk);
  Hist->fLogLHR->Fill(llhr);

  Hist->fPdgCode->Fill(Track->fPdgCode);
  Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->NActive()+1.e-5));

  int i1, i2, n1(0) ,n2(0), ndiff(0);
  int nbits = Track->fHitMask.GetNBits();
  for (int i=0; i<nbits; i++) {
    i1 = Track->HitMask()->GetBit(i);
    i2 = Track->ExpectedHitMask()->GetBit(i);
    n1 += i1;
    n2 += i2;
    if (i1 != i2) ndiff += 1;
  }
  

  Hist->fNEPlVsNHPl->Fill(n2   , n1);
  Hist->fNDPlVsNHPl->Fill(ndiff, n1);
  Hist->fChi2dVsNDPl->Fill(ndiff, Track->Chi2Dof());
  Hist->fDpFVsNDPl  ->Fill(ndiff, dpf); //tp->fDpF);

  float        fre1(-1), fre2(-1);
  int          icl;
  TStnCluster* cl;

  if (Track->fVMinS) {
    icl = Track->fVMinS->fClusterIndex;
    if (icl >= 0) {
      cl = fClusterBlock->Cluster(icl);
      fre1 = cl->fFrE1;
      fre2 = cl->fFrE2;
    }
  }

  Hist->fFrE1->Fill(fre1);
  Hist->fFrE2->Fill(fre2);
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrigAna001Module::BeginJob() {
  //-----------------------------------------------------------------------------
  // register data blocks
  //-----------------------------------------------------------------------------
  // // RegisterDataBlock("CprTimeClusterBlock","TStnTimeClusterBlock" ,&fTimeClusterBlock    );
  if (fIsBackground){
    RegisterDataBlock("TrackSeedBlockDemCpr","TStnTrackSeedBlock"     ,&fCprDemTrackSeedBlock);
    RegisterDataBlock("TrackSeedBlockDemTpr","TStnTrackSeedBlock"     ,&fTprDemTrackSeedBlock);
    RegisterDataBlock("TrackSeedBlockDepCpr","TStnTrackSeedBlock"     ,&fCprDepTrackSeedBlock);
    RegisterDataBlock("TrackSeedBlockDepTpr","TStnTrackSeedBlock"     ,&fTprDepTrackSeedBlock);
    RegisterDataBlock("HelixBlockDemTpr"    ,"TStnHelixBlock"         ,&fTprDemHelixBlock    );
    RegisterDataBlock("HelixBlockDepTpr"    ,"TStnHelixBlock"         ,&fTprDepHelixBlock    );
    RegisterDataBlock("HelixBlockDemCpr"    ,"TStnHelixBlock"         ,&fCprDemHelixBlock    );
    RegisterDataBlock("HelixBlockDepCpr"    ,"TStnHelixBlock"         ,&fCprDepHelixBlock    );
    RegisterDataBlock("TrackBlockDemTpr"    ,"TStnTrackBlock"         ,&fTprDemTrackBlock       );
    RegisterDataBlock("TrackBlockDemCpr"    ,"TStnTrackBlock"         ,&fCprDemTrackBlock       );
    RegisterDataBlock("TrackBlockDepTpr"    ,"TStnTrackBlock"         ,&fTprDepTrackBlock    );
    RegisterDataBlock("TrackBlockDepCpr"    ,"TStnTrackBlock"         ,&fCprDepTrackBlock    );
    // RegisterDataBlock("CprTimeClusterBlockDem","TStnTimeClusterBlockDem" ,&fTimeClusterBlock );
  } else {
    RegisterDataBlock("CprTrackSeedBlock","TStnTrackSeedBlock"     ,&fCprDemTrackSeedBlock);
    RegisterDataBlock("TprTrackSeedBlock","TStnTrackSeedBlock"     ,&fTprDemTrackSeedBlock);
    RegisterDataBlock("TprHelixBlock"    ,"TStnHelixBlock"         ,&fTprDemHelixBlock    );
    RegisterDataBlock("TprNegHelixBlock" ,"TStnHelixBlock"         ,&fTprDepHelixBlock );
    RegisterDataBlock("CprHelixBlock"    ,"TStnHelixBlock"         ,&fCprDemHelixBlock    );
    //RegisterDataBlock("HelixBlock"    ,"TStnHelixBlock"    ,&fHelixBlock       );
    RegisterDataBlock("TprTrackBlock"    ,"TStnTrackBlock"         ,&fTprDemTrackBlock    );
    RegisterDataBlock("CprTrackBlock"    ,"TStnTrackBlock"         ,&fCprDemTrackBlock    );
  
    RegisterDataBlock("TrackBlock"       ,"TStnTrackBlock"         ,&fTrackBlock       );
  }
  RegisterDataBlock("ClusterBlock"     ,"TStnClusterBlock"       ,&fClusterBlock     );
  RegisterDataBlock("CalDataBlock"     ,"TCalDataBlock"          ,&fCalDataBlock     );
  RegisterDataBlock("StrawDataBlock"   ,"TStrawDataBlock"        ,&fStrawDataBlock   );
  RegisterDataBlock("GenpBlock"        ,"TGenpBlock"             ,&fGenpBlock        );
  RegisterDataBlock("SimpBlock"        ,"TSimpBlock"             ,&fSimpBlock        );
  RegisterDataBlock("VdetBlock"        ,"TVdetDataBlock"         ,&fVdetDataBlock    );

  //-----------------------------------------------------------------------------
  // book histograms
  //-----------------------------------------------------------------------------
  BookHistograms();

  fLogLH->Init("v4_2_4");

  return 0;
}


//_____________________________________________________________________________
int TTrigAna001Module::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TTrigAna001Module::FillHistograms() {

  double       cos_th (-2.),  cl_e(-1.);
  int          disk_id(-1), alg_mask, nsh, nactive;
  float        pfront, ce_pitch, reco_pitch, fcons, t0, sigt, sigp, p; 
  TStnCluster  *cl0;

  //  cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

  if (fNClusters > 0) {
    cl0     = fClusterBlock->Cluster(0);
    cl_e    = cl0->Energy();
    disk_id = cl0->DiskID();
  }
  //-----------------------------------------------------------------------------
  // event histograms
  //-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fNTracks[1]> 0) FillEventHistograms(fHist.fEvent[1]);
  else                FillEventHistograms(fHist.fEvent[2]);

  if (fNClusters > 0) FillEventHistograms(fHist.fEvent[3]);
  else                FillEventHistograms(fHist.fEvent[4]);

  if ((fNTracks[1] == 0) && (fabs(cos_th) < 0.4)) {
    FillEventHistograms(fHist.fEvent[5]); 
  }

  if (fNGoodTracks > 0) {
    FillEventHistograms(fHist.fEvent[6]); 

    TLorentzVector    mom;
    
    if(fParticle != NULL)
      fParticle->Momentum(mom);

    double p, cos_th;

    p      = mom.P();
    cos_th = mom.Pz()/p;

    if (GetDebugBit(31) && (cos_th > 0.8)) {
      GetHeaderBlock()->Print(Form(" bit:031 cos_th = %10.3f p = %10.3f ntrk = %5i",
				   cos_th, p, fNTracks[0]));
    }
  }

  if (cl_e > 60.) {
    FillEventHistograms(fHist.fEvent[7]); 
    if (GetDebugBit(34)) {
      if (fNTracks[0] <= 0) {
	GetHeaderBlock()->Print(Form(" bit:034 cl_e = %10.3f",cl_e));
      }
    }
  }

  if      (disk_id == 0) FillEventHistograms(fHist.fEvent[8]);
  else if (disk_id == 1) FillEventHistograms(fHist.fEvent[9]);
  //-----------------------------------------------------------------------------
  // Dave's ladder for all tracks
  // 1. N(straw hits) > 20
  //-----------------------------------------------------------------------------
  if (fSimp) {
    nsh    = fSimp->NStrawHits();
    pfront = fSimp->fMomTrackerFront;
  }
  else {
    nsh    = -1;
    pfront = -1.e6;
  }
  
  if (nsh >= 20) {
    FillEventHistograms(fHist.fEvent[11]);
    if (pfront > 100.) {
      FillEventHistograms(fHist.fEvent[12]);
      
      ce_pitch = 0.7; // kludge
      if ((ce_pitch > 0.577) && (ce_pitch < 1.)) {
	FillEventHistograms(fHist.fEvent[13]);

	if (fNTracks[0] > 0) {
	  FillEventHistograms(fHist.fEvent[14]);

	  // here we have a track reconstructed

	  TStnTrack* trk = fTrackBlock->Track(0);

	  fcons = trk->fFitCons;
	  t0    = trk->T0();
	  reco_pitch = trk->fTanDip;
	  sigp       = trk->fFitMomErr;
	  sigt       = trk->fT0Err;
	  nactive    = trk->NActive();
	  p          = trk->fP;
	  // fit quality
	  if ((nactive > 25) && (fcons > 2.e-3) && (sigp < 0.25) && (sigt < 1.0))  {
	    FillEventHistograms(fHist.fEvent[15]);
	    if (t0 > 700) {
	      FillEventHistograms(fHist.fEvent[16]);
	      if ((reco_pitch > 0.577) && (reco_pitch < 1.)) {
		FillEventHistograms(fHist.fEvent[17]);
		if (p > 103.5) {
		  FillEventHistograms(fHist.fEvent[18]);
		}
	      }
	    }
	  }

	  alg_mask = trk->AlgMask();

	  if ((alg_mask == 1) || (alg_mask == 3)) {
	    //-----------------------------------------------------------------------------
	    // track reconstructed with TrkPatRec 
	    //-----------------------------------------------------------------------------
	    FillEventHistograms(fHist.fEvent[24]);
	    if ((nactive > 25) && (fcons > 2.e-3) && (sigp < 0.25) && (sigt < 1.0))  {
	      FillEventHistograms(fHist.fEvent[25]);
	      if (t0 > 700) {
		FillEventHistograms(fHist.fEvent[26]);
		if ((reco_pitch > 0.577) && (reco_pitch < 1.)) {
		  FillEventHistograms(fHist.fEvent[27]);
		  if (p > 103.5) {
		    FillEventHistograms(fHist.fEvent[28]);
		  }
		}
	      }
	    }
	  }
	  else if (alg_mask == 2) {
	    //-----------------------------------------------------------------------------
	    // track reconstructed with CalPatRec, but not with TrkPatRec
	    //-----------------------------------------------------------------------------
	    //int x=0;
	  }
	}
      }
    }
  }
  //-----------------------------------------------------------------------------
  // the same ladder for TrkPatRec tracks 
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // Simp histograms
  //-----------------------------------------------------------------------------
  if (fSimp) {
    FillSimpHistograms(fHist.fSimp[0],fSimp);
  }



  //--------------------------------------------------------------------------------
  // MergePatRec tracks
  //--------------------------------------------------------------------------------
  TStnTrack*     trk;
  int            trigFlag0(0), trigFlag1(0);
  int            trigFlag2(0), trigFlag3(0);
  int            trigSeedFlag(0);

  double  fMinNActive      = 25;
  double  fMinChi2d        = 5.;
  double  fMaxD0           = 200.;
  double  fTrigPMin        = 80.;
  //  double  fMaxT0Err        = 0.9;  		// in ns
  double  fMaxFitMomErr    = 0.25;  		// in MeV
  double  fMinTanDip       = tan(M_PI/6.);	// 0.5773
  double  fMaxTanDip       = 1.0;  
  double  fMinD1           = -80.;		// in mm
  double  fMaxD1           = 105.;
  double  fMinT0           = 700.; // ns
  double  fPMin            = 100.; //MeV/c
  double  fPMax            = 105.;
  double  fPMinDep         = 89.; //MeV/c
  double  fPMaxDep         = 94.;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock->Track(i);

    //    int   trackQual(0);
    fTrackCaloE = trk->Ep()*trk->fP;

    FillTrackHistograms(fHist.fTrack[0],trk);

    if (trk->fIDWord == 0) {	// track passes selection "C" 
      
      FillTrackHistograms(fHist.fTrack[1],trk);
      
      if (trk->fP >= 100. && trk->fP < 110.) {
	FillTrackHistograms(fHist.fTrack[3],trk);
      }

    }

    
    double     d0         = trk->fD0;
    double     chi2dof    = trk->Chi2Dof();
    double     t0         = trk->fT0;
    //    double t0_err     = trk->fT0Err;
    double     nactive    = trk->NActive();
    double     tan_dip    = trk->fTanDip;
    double     fitmom_err = trk->fFitMomErr;
    double     p          = trk->P();

    if ( (chi2dof < fMinChi2d  ) &&  
	 (tan_dip < fMaxTanDip ) &&
	 (tan_dip > fMinTanDip ) &&
	 (nactive > fMinNActive) &&
	 (fitmom_err< fMaxFitMomErr) &&
	 (d0 < fMaxD1)           &&
	 (d0 > fMinD1)           && 
	 (t0 > fMinT0 )          &&
	 ((p>=fPMin) && (p<=fPMax)) ){
      FillTrackHistograms(fHist.fTrack[19],trk);
    }
    //-----------------------------------------------------------------------------
    // split tracks by the algorithm mask: 1 , 2 , or 3
    //-----------------------------------------------------------------------------
    int     alg = trk->AlgMask();
    if ((alg == 2) || (alg == 3)){
      if (trk->fTrkQual > 0.2) {
	FillTrackHistograms(fHist.fTrack[7],trk);

	if (trk->fP >= 100. && trk->fP < 110.) {
	  FillTrackHistograms(fHist.fTrack[8],trk);
	  
	  if ( trigSeedFlag == 1){

	    if (trigFlag1 == 1){
	      FillTrackHistograms(fHist.fTrack[9],trk);
	    }
	    if (trigFlag0 == 1){
	      FillTrackHistograms(fHist.fTrack[10],trk);
	    }
	    
	    if (trigFlag3 == 1){
	      FillTrackHistograms(fHist.fTrack[13],trk);
	    }
	    if (trigFlag2 == 1){
	      FillTrackHistograms(fHist.fTrack[14],trk);
	    }
	  }
	}
      }
    }


    alg_mask = trk->BestAlg();//AlgMask();
    if      (alg_mask == 0) {
      //-----------------------------------------------------------------------------
      // TrkPatRec tracks
      //-----------------------------------------------------------------------------
      // double chi2dof    = trk->Chi2Dof();
      // double t0         = trk->fT0;
      // //    double t0_err     = trk->fT0Err;
      // double nactive    = trk->NActive();
      // double tan_dip    = trk->fTanDip;
      // double fitmom_err = trk->fFitMomErr;
      // double d0         = trk->fD0;
      // double p          = trk->P();
      
      // if ( (chi2dof < fMinChi2d  ) &&  
      // 	   (tan_dip < fMaxTanDip ) &&
      // 	   (tan_dip > fMinTanDip ) &&
      // 	   (nactive > fMinNActive) &&
      // 	   (fitmom_err< fMaxFitMomErr) &&
      // 	   (d0 < fMaxD1)           &&
      // 	   (d0 > fMinD1)           && 
      // 	   (t0 > fMinT0 )          &&
      // 	   ((p>=fPMin) && (p<=fPMax)) ){
      // 	trackQual = 1;
      // }
      
      if (trk->fTrkQual > 0.4) {
	FillTrackHistograms(fHist.fTrack[2],trk);
	
	if (trk->fP >= 100. && trk->fP < 110.) {
	  FillTrackHistograms(fHist.fTrack[4],trk);

	  //	  trackPassingTrkQualCut = 1;
	  
	  if ( trigSeedFlag == 1){
	    
	    if (trigFlag1 == 1){
	      FillTrackHistograms(fHist.fTrack[5],trk);
	    }
	    if (trigFlag0 == 1){
	      FillTrackHistograms(fHist.fTrack[6],trk);
	    }
	    
	    if (trigFlag3 == 1){
	      FillTrackHistograms(fHist.fTrack[11],trk);
	    }
	    if (trigFlag2 == 1){
	      FillTrackHistograms(fHist.fTrack[12],trk);
	    }
	  }
	}
      }
      
    }
    //    else if  (alg_mask == 1) {
    if ( (trk->fTrkQual > 0.2) && (alg_mask == 1)) {
      FillTrackHistograms(fHist.fTrack[2],trk);
    }

    if (trk->fP >= 100. && trk->fP < 110.) {
      FillTrackHistograms(fHist.fTrack[4],trk);
    }
    //	  trackPassingTrkQualCut = 1;
	
    if ((p>=fPMin) && (p<=fPMax)) {
      if ( (nactive >=20) && (chi2dof < 5) ) FillTrackHistograms(fHist.fTrack[15],trk);
      if ( (nactive >=20) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[16],trk);
      if ( (nactive >=25) && (chi2dof < 5) ) {
	FillTrackHistograms(fHist.fTrack[17],trk);
      }
      if ( (nactive >=25) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[18],trk);
    }
	

  }
    

  //--------------------------------------------------------------------------------  
  // TrkPatRec 
  //--------------------------------------------------------------------------------
  for (int i=0; i<fNTracks[2]; ++i){
    trk  = fTprDemTrackBlock->Track(i);

    FillTrackHistograms(fHist.fTrack[200],trk);
    int trkSeedIndex = trk->TrackSeedIndex();
    if ( trkSeedIndex >= 0 && ( trkSeedIndex < fNTrackSeeds[1])){
      TStnTrackSeed* trkSeed = fTprDemTrackSeedBlock->TrackSeed(trkSeedIndex);
      fTprTrackHelixIndex   = trkSeed->HelixIndex();
    }


    double chi2dof    = trk->Chi2Dof();
    double t0         = trk->fT0;
    //    double t0_err   = track->fT0Err;
    double nactive    = trk->NActive();
    double tan_dip    = trk->fTanDip;
    double fitmom_err = trk->fFitMomErr;
    double d0         = trk->fD0;
    double p          = trk->P();

    if (GetDebugBit(51))   GetHeaderBlock()->Print(Form(" bit:051 trkPatRec track"));// p = %2.3f MeV/c nActive = %2.0f chi2d = %2.3f ", p, nactive, chi2dof));

    if ( (chi2dof < fMinChi2d  ) &&  
	 (tan_dip < fMaxTanDip ) &&
	 (tan_dip > fMinTanDip ) &&
	 (nactive > fMinNActive) &&
	 (fitmom_err< fMaxFitMomErr) &&
	 (d0 < fMaxD1)           &&
	 (d0 > fMinD1)           && 
	 (t0 > fMinT0 )          &&
	 ((p>=fPMin) && (p<=fPMax)) ){
      FillTrackHistograms(fHist.fTrack[210],trk);
    }
    
    if ((p>=fPMin) && (p<=fPMax)) {
      if ( (nactive >=20) && (chi2dof < 5) ) FillTrackHistograms(fHist.fTrack[215],trk);
      if ( (nactive >=20) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[216],trk);
      if ( (nactive >=25) && (chi2dof < 5) ) {
	FillTrackHistograms(fHist.fTrack[217],trk);
	if (GetDebugBit(52))   GetHeaderBlock()->Print(Form(" bit:052 trkPatRec track"));
      }
      if ( (nactive >=25) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[218],trk);
    }

  }

  for (int i=0; i<fNTracks[4]; ++i){
    trk  = fTprDepTrackBlock->Track(i);

    FillTrackHistograms(fHist.fTrack[400],trk);
    int trkSeedIndex = trk->TrackSeedIndex();
    if ( trkSeedIndex >= 0 && ( trkSeedIndex < fNTrackSeeds[3])){
      TStnTrackSeed* trkSeed = fTprDepTrackSeedBlock->TrackSeed(trkSeedIndex);
      fTprDepTrackHelixIndex   = trkSeed->HelixIndex();
    }


    double chi2dof    = trk->Chi2Dof();
    double t0         = trk->fT0;
    //    double t0_err   = track->fT0Err;
    double nactive    = trk->NActive();
    double tan_dip    = trk->fTanDip;
    double fitmom_err = trk->fFitMomErr;
    double d0         = trk->fD0;
    double p          = trk->P();

    if ( (chi2dof < fMinChi2d  ) &&  
	 (tan_dip < fMaxTanDip ) &&
	 (tan_dip > fMinTanDip ) &&
	 (nactive > fMinNActive) &&
	 (fitmom_err< fMaxFitMomErr) &&
	 (d0 < fMaxD1)           &&
	 (d0 > fMinD1)           && 
	 (t0 > fMinT0 )          &&
	 ((p>=fPMinDep) && (p<=fPMaxDep)) ){
      FillTrackHistograms(fHist.fTrack[410],trk);
    }
    
    if ((p>=fPMinDep) && (p<=fPMaxDep)) {
      if ( (nactive >=20) && (chi2dof < 5) ) FillTrackHistograms(fHist.fTrack[415],trk);
      if ( (nactive >=20) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[416],trk);
      if ( (nactive >=25) && (chi2dof < 5) ) {
	FillTrackHistograms(fHist.fTrack[417],trk);
	if (GetDebugBit(52))   GetHeaderBlock()->Print(Form(" bit:052 trkPatRec track"));
      }
      if ( (nactive >=25) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[418],trk);
    }

  }

  //--------------------------------------------------------------------------------  
  // CalTrkFit 
  //--------------------------------------------------------------------------------
  int            trackPassingTrkQualCut(0);
  TStnTrackSeed* trkSeed;

  for (int i=0; i<fNTracks[1]; ++i ) {
    trk                = fCprDemTrackBlock->Track(i);

    fTrackT0           = trk->T0();
    fTrackTrkSeedIndex = trk->TrackSeedIndex();
    
    if (fTrackTrkSeedIndex>=0 && ( fTrackTrkSeedIndex < fNTrackSeeds[0] )) {
      trkSeed            = fCprDemTrackSeedBlock->TrackSeed(fTrackTrkSeedIndex);
      fTrackHelixIndex   = trkSeed->HelixIndex();
      if (GetDebugBit(49))   GetHeaderBlock()->Print(Form(" bit:049 track-seed index = %i helix-index = %i", fTrackTrkSeedIndex, fTrackHelixIndex));
    }

    FillTrackHistograms(fHist.fTrack[100],trk);

    if (GetDebugBit(40)) {
      double     p     = trk->fP;
      double     chi2c = trk->fChi2/(trk->NActive()-5.);
      double     d0    = trk->fD0;

      if( p < 50) {
	GetHeaderBlock()->Print(Form(" bit:040 p = %5.3f chi2c = %5.3f d0 = %5.3f", p, chi2c, d0));
      }
    }

    int        nhits = trk->NActive();
    double     chi2c = trk->fChi2/(nhits-5.);
    double      d0   = trk->fD0;

    //check if the trk passes the "set C" cuts
    double chi2dof    = trk->Chi2Dof();
    //    double t0         = trk->fT0;
    //    double t0_err     = trk->fT0Err;
    double nactive    = trk->NActive();
    double tan_dip    = trk->fTanDip;
    //    double fitmom_err = trk->fFitMomErr;
    d0                = trk->fD0;
    double p          = trk->P();

    if ( (nhits >=10) && (chi2c < 4) && (fabs(d0) <200.)){
      trigFlag2 = 1;

      if (chi2c < 3){
	trigFlag3 = 1;
      }
    }
    
    if ( (nhits >=15) && (chi2c < 4) && (fabs(d0) <200.)){
      trigFlag0 = 1;
      
      if (chi2c < 3){
	trigFlag1 = 1;
      }
    }
    
    if (nhits > 15) {
      FillTrackHistograms(fHist.fTrack[101],trk);
    }

    if (chi2c < 3){
      FillTrackHistograms(fHist.fTrack[102],trk);    
    }

    if (chi2c < 4){
      FillTrackHistograms(fHist.fTrack[103],trk);    
    }

    if ((chi2c < 3) && (nhits > 15)){
      FillTrackHistograms(fHist.fTrack[104],trk); 
    }

    if ( (chi2c < 4) && (nhits > 15)){
      FillTrackHistograms(fHist.fTrack[105],trk);    
      if (GetDebugBit(45)) {
	GetHeaderBlock()->Print(Form("bit:045 nHits = %i p = %2.3f tan_dip = %2.3f chi2 = %2.3f", 
				     nhits, p, tan_dip, chi2dof));
      }
    }
    
    if ((chi2c < 3) && (fabs(d0) < 200.)){
      FillTrackHistograms(fHist.fTrack[106],trk);    
    }
    
    if ( (chi2c < 4) && (fabs(d0) < 200.) ){
      FillTrackHistograms(fHist.fTrack[107],trk);    
    }
    
    if ((chi2c < 3) && (fabs(d0) < 200.) && (nhits > 15)){
      FillTrackHistograms(fHist.fTrack[108],trk);    
    }
    
    if ( (chi2c < 4) && (fabs(d0) < 200.) && (nhits > 15)){
      FillTrackHistograms(fHist.fTrack[109],trk);    
    }

  

    if ( (chi2dof < fMinChi2d  )  &&  
	 //	 (tan_dip < fMaxTanDip ) &&
	 //	 (tan_dip > fMinTanDip ) &&
	 (nactive >= fMinNActive) &&
	 //	 (fitmom_err< fMaxFitMomErr) &&
	 (fabs(d0) < fMaxD0)      &&
	 // (d0 < fMaxD1)           &&
	 // (d0 > fMinD1)           && 
	 // (t0 > fMinT0 )          &&
	 //	 ((p>=fPMin) && (p<=fPMax)) ){
	 (p>=fTrigPMin)  ){
      trackPassingTrkQualCut = 1;
    }

 
  }


  //CalPatRecDep
  for (int i=0; i<fNTracks[3]; ++i){
    trk  = fCprDepTrackBlock->Track(i);

    FillTrackHistograms(fHist.fTrack[300],trk);
    int trkSeedIndex = trk->TrackSeedIndex();
    if ( trkSeedIndex >= 0 && ( trkSeedIndex < fNTrackSeeds[3])){
      TStnTrackSeed* trkSeed = fCprDepTrackSeedBlock->TrackSeed(trkSeedIndex);
      fTrackDepHelixIndex   = trkSeed->HelixIndex();
    }


    double chi2dof    = trk->Chi2Dof();
    double t0         = trk->fT0;
    //    double t0_err   = track->fT0Err;
    double nactive    = trk->NActive();
    double tan_dip    = trk->fTanDip;
    double fitmom_err = trk->fFitMomErr;
    double d0         = trk->fD0;
    double p          = trk->P();

    if ( (chi2dof < fMinChi2d  ) &&  
	 (tan_dip < fMaxTanDip ) &&
	 (tan_dip > fMinTanDip ) &&
	 (nactive > fMinNActive) &&
	 (fitmom_err< fMaxFitMomErr) &&
	 (d0 < fMaxD1)           &&
	 (d0 > fMinD1)           && 
	 (t0 > fMinT0 )          &&
	 ((p>=fPMinDep) && (p<=fPMaxDep)) ){
      FillTrackHistograms(fHist.fTrack[310],trk);
    }
    
    if ((p>=fPMinDep) && (p<=fPMaxDep)) {
      if ( (nactive >=20) && (chi2dof < 5) ) FillTrackHistograms(fHist.fTrack[315],trk);
      if ( (nactive >=20) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[316],trk);
      if ( (nactive >=25) && (chi2dof < 5) ) {
	FillTrackHistograms(fHist.fTrack[317],trk);
	//	if (GetDebugBit(52))   GetHeaderBlock()->Print(Form(" bit:052 trkPatRec track"));
      }
      if ( (nactive >=25) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[318],trk);
    }

  }

  //--------------------------------------------------------------------------------
  bool cprDemHelixTrig(false),  cprDepHelixTrig(false), tprDemHelixTrig(false), tprDepHelixTrig(false);

  //--------------------------------------------------------------------------------
  // helix histograms from the mergere module
  //--------------------------------------------------------------------------------
  TStnHelix* helix;

  for (int i=0; i<fNHelices[0]; ++i){
    
    helix = fTprDemHelixBlock->Helix(i);
    
    FillHelixHistograms(fHist.fHelix[0], helix);
    
    int         nhits    = helix->NHits();
    double      radius   = helix->Radius();
    double      lambda   = helix->Lambda();
    double      tanDip   = lambda/radius;
    double      mm2MeV   = 3/10.;
    double      p        = radius*mm2MeV/std::cos( std::atan(tanDip));
    double      chi2xy   = helix->Chi2XY();
    double      chi2zphi = helix->Chi2ZPhi();
    

    if (GetDebugBit(43)) {
      // int    id1   = helix->fSimpId1;
      int    nh1   = helix->fSimpId1Hits;
      // double p1    = helix->fSimp1P;
      // double pt1   = helix->fSimp1Pt;
      
      //      TSimParticle*sim1(0);
      int    pdg1  = helix->fSimpPDG1;
      int    pdgM1 = helix->fSimpPDGM1;
      TLorentzVector p1 = helix->fMom1;
      TLorentzVector x1 = helix->fOrigin1;

      int    nh2   = helix->fSimpId2Hits;
      int    pdg2  = helix->fSimpPDG2;
      int    pdgM2 = helix->fSimpPDGM2;
      TLorentzVector p2 = helix->fMom2;
      TLorentzVector x2 = helix->fOrigin2;
      
      GetHeaderBlock()->Print(Form(" bit:043 TPR p = %5.3f  nhits = %3i chi2XY = %5.3f  chi2PhiZ = %5.3f nh1 = %i p1 = %2.3f pt1 = %2.3f pdg1 = %i pdgM1 = %i x1=  %2.3f y1 = %2.3f z1 = %2.3f pdg2 = %i pdgM2 = %i nh2 = %i p2 = %2.3f pt2 = %2.3f", 
				   p, nhits, chi2xy, chi2zphi, nh1, p1.Vect().Mag(), p1.Vect().Perp(),
				   pdg1,pdgM1, 
				   x1.Vect().x(), x1.Vect().y(), x1.Vect().z(),
				   pdg2, pdgM2, nh2, p2.Vect().Mag(), p2.Vect().Perp()));
    }

    if (i  == fTprTrackHelixIndex){
      FillHelixHistograms(fHist.fHelix[8], helix);
    }

    if (p > 80.) {
      FillHelixHistograms(fHist.fHelix[1], helix);
    }
    
    if (p > 90) {
      FillHelixHistograms(fHist.fHelix[2], helix);
    }

    if (p > 100) {
      FillHelixHistograms(fHist.fHelix[3], helix);
    }
    
    if ( (nhits>10) && (nhits<15)) {
      FillHelixHistograms(fHist.fHelix[4], helix);
    }

    if ( nhits>=15) {
      // tprDemHelixTrig = true;
      FillHelixHistograms(fHist.fHelix[5], helix);
    }
    
    if ( (chi2xy <= 5) && (chi2zphi <= 5) && (nhits>=15)){
      tprDemHelixTrig = true;
      FillHelixHistograms(fHist.fHelix[6], helix);
    }
  
    double      pMC     = helix->Mom1().Vect().Mag();
    int         pdgm    = helix->PDGMother1();
    int         pdg     = helix->PDG1();
    if ( (pMC>=100.) && (pMC<=105) && (pdgm == 11) && (pdg == 11)) {
      FillHelixHistograms(fHist.fHelix[10], helix);

      if (i  == fTprTrackHelixIndex){
	FillHelixHistograms(fHist.fHelix[9], helix);
      }

      if (GetDebugBit(53))   GetHeaderBlock()->Print(Form(" bit:053 trkPatRec helix"));
    }
    
    if ( GetDebugBit(41) ) {
      

      double    cpr_chi2xy(-1), cpr_chi2zphi(-1), cpr_mom (-1);
      int       cpr_nhits(-1);
      double    tpr_chi2xy(-1), tpr_chi2zphi(-1), tpr_mom (-1);
      int       tpr_nhits(-1);
      double    tanDip(-1);
      double    mm2MeV   = 3/10.;


      if (fNHelices[1] > 0) {
	TStnHelix* cpr_helix = fCprDemHelixBlock->Helix(0);
	cpr_chi2xy   = cpr_helix->Chi2XY();				       
	cpr_chi2zphi = cpr_helix->Chi2ZPhi();				       
	tanDip       = cpr_helix->Lambda()/cpr_helix->Radius();		       
	cpr_mom      = cpr_helix->Radius()*mm2MeV/std::cos( std::atan(tanDip));
	cpr_nhits    = cpr_helix->NHits();                                       
      }
      if (fNHelices[2] > 0) {
	TStnHelix* tpr_helix = fTprDemHelixBlock->Helix(0);
	tpr_chi2xy   = tpr_helix->Chi2XY();				       
	tpr_chi2zphi = tpr_helix->Chi2ZPhi();				       
	tanDip       = tpr_helix->Lambda()/tpr_helix->Radius();		       
	tpr_mom      = tpr_helix->Radius()*mm2MeV/std::cos( std::atan(tanDip));
      	tpr_nhits    = tpr_helix->NHits();                                     
      }
      
      GetHeaderBlock()->Print(Form(" bit:041 algMask = %i Cpr: chi2XY = %10.3f chi2ZPhi = %2.3e p = %10.3f  nhits = %5i; Tpr: chi2XY = %10.3f chi2ZPhi = %10.3f p = %10.3f nhits = %5i",
				   helix->AlgMask(), 
				   cpr_chi2xy, cpr_chi2zphi, cpr_mom, cpr_nhits,
				   tpr_chi2xy, tpr_chi2zphi, tpr_mom, tpr_nhits));
    }

  }


  for (int i=0; i<fNHelices[3]; ++i){
    
    helix = fTprDepHelixBlock->Helix(i);
    
    FillHelixHistograms(fHist.fHelix[300], helix);
    
    int         nhits    = helix->NHits();
    double      radius   = helix->Radius();
    double      lambda   = helix->Lambda();
    double      tanDip   = lambda/radius;
    double      mm2MeV   = 3/10.;
    double      p        = radius*mm2MeV/std::cos( std::atan(tanDip));
    double      chi2xy   = helix->Chi2XY();
    double      chi2zphi = helix->Chi2ZPhi();
    

    if (GetDebugBit(43)) {
      // int    id1   = helix->fSimpId1;
      int    nh1   = helix->fSimpId1Hits;
      // double p1    = helix->fSimp1P;
      // double pt1   = helix->fSimp1Pt;
      
      //      TSimParticle*sim1(0);
      int    pdg1  = helix->fSimpPDG1;
      int    pdgM1 = helix->fSimpPDGM1;
      TLorentzVector p1 = helix->fMom1;
      TLorentzVector x1 = helix->fOrigin1;

      int    nh2   = helix->fSimpId2Hits;
      int    pdg2  = helix->fSimpPDG2;
      int    pdgM2 = helix->fSimpPDGM2;
      TLorentzVector p2 = helix->fMom2;
      TLorentzVector x2 = helix->fOrigin2;
      
      GetHeaderBlock()->Print(Form(" bit:043 TPR p = %5.3f  nhits = %3i chi2XY = %5.3f  chi2PhiZ = %5.3f nh1 = %i p1 = %2.3f pt1 = %2.3f pdg1 = %i pdgM1 = %i x1=  %2.3f y1 = %2.3f z1 = %2.3f pdg2 = %i pdgM2 = %i nh2 = %i p2 = %2.3f pt2 = %2.3f", 
				   p, nhits, chi2xy, chi2zphi, nh1, p1.Vect().Mag(), p1.Vect().Perp(),
				   pdg1,pdgM1, 
				   x1.Vect().x(), x1.Vect().y(), x1.Vect().z(),
				   pdg2, pdgM2, nh2, p2.Vect().Mag(), p2.Vect().Perp()));
    }

    if (p > 80.) {
      FillHelixHistograms(fHist.fHelix[301], helix);
    }
    
    if (p > 90) {
      FillHelixHistograms(fHist.fHelix[302], helix);
    }

    if (p > 100) {
      FillHelixHistograms(fHist.fHelix[303], helix);
    }
    
    if ( (nhits>10) && (nhits<15)) {
      FillHelixHistograms(fHist.fHelix[304], helix);
    }

    if ( nhits>=15) {
      // tprDepHelixTrig = true;
      FillHelixHistograms(fHist.fHelix[305], helix);
    }
    
    if ( (chi2xy <= 5) && (chi2zphi <= 5) && (nhits>=15)){
      tprDepHelixTrig = true;
      FillHelixHistograms(fHist.fHelix[306], helix);
    }

    if (i  == fTprDepTrackHelixIndex){
      FillHelixHistograms(fHist.fHelix[308], helix);
    }

  
    double      pMC     = helix->Mom1().Vect().Mag();
    int         pdgm    = helix->PDGMother1();
    int         pdg     = helix->PDG1();
    if ( (pMC>=90.) && (pMC<=100) && (pdgm == -11) && (pdg == -11)) {
      FillHelixHistograms(fHist.fHelix[310], helix);
      if (i  == fTprDepTrackHelixIndex){
	FillHelixHistograms(fHist.fHelix[309], helix);
      }
    }
    
  }



  //--------------------------------------------------------------------------
  // calo-seeded algorithm
  //--------------------------------------------------------------------------
  for (int i=0; i<fNHelices[1]; ++i){
    
    helix = fCprDemHelixBlock->Helix(i);
    
    FillHelixHistograms(fHist.fHelix[100], helix);
    
    int         nhits    = helix->NHits();
    double      radius   = helix->Radius();
    double      lambda   = helix->Lambda();
    double      tanDip   = lambda/radius;
    double      mm2MeV   = 3/10.;
    double      p        = radius*mm2MeV/std::cos( std::atan(tanDip));
    double      chi2xy   = helix->Chi2XY();
    double      chi2zphi = helix->Chi2ZPhi();
    
    if (chi2zphi > 5) {
      FillHelixHistograms(fHist.fHelix[120], helix);
    }


    if (GetDebugBit(43)) {
      //      int    id1   = helix->fSimpId1;
      int    nh1   = helix->fSimpId1Hits;
      int    pdg1  = helix->fSimpPDG1;
      int    pdgM1 = helix->fSimpPDGM1;
      TLorentzVector p1 = helix->fMom1;
      TLorentzVector x1 = helix->fOrigin1;

      int    nh2   = helix->fSimpId2Hits;
      int    pdg2  = helix->fSimpPDG2;
      int    pdgM2 = helix->fSimpPDGM2;
      TLorentzVector p2 = helix->fMom2;
      TLorentzVector x2 = helix->fOrigin2;
      GetHeaderBlock()->Print(Form(" bit:043 CPR p = %5.3f  nhits = %3i chi2XY = %5.3f  chi2PhiZ = %5.3f nh1 = %i p1 = %2.3f pt1 = %2.3f pdg1 = %i pdgM1 = %i x1=  %2.3f y1 = %2.3f z1 = %2.3f  nh2 = %i p2 = %2.3f pt2 = %2.3f pdg2 = %i pdgM2 = %i", 
				   p, nhits, chi2xy, chi2zphi, nh1, p1.Vect().Mag(), p1.Vect().Perp(),
				   pdg1,pdgM1, 
				   x1.Vect().x(), x1.Vect().y(), x1.Vect().z(),
				   nh2, p2.Vect().Mag(), p2.Vect().Perp(), pdg2, pdgM2));
    }
    
    if (p > 80.) {
      FillHelixHistograms(fHist.fHelix[101], helix);
    }
    
    if (p > 90) {
      FillHelixHistograms(fHist.fHelix[102], helix);
    }

    if (p > 100) {
      FillHelixHistograms(fHist.fHelix[103], helix);
    }
    
    if ( (nhits>10) && (nhits<15)) {
      FillHelixHistograms(fHist.fHelix[104], helix);
    }

    if ( nhits>=15) {
      cprDemHelixTrig = true;
      FillHelixHistograms(fHist.fHelix[105], helix);
    }
    
    if ( (chi2xy <= 5.) && (chi2zphi <= 5.) && (nhits>=15)){
      FillHelixHistograms(fHist.fHelix[106], helix);
    }
  
    if ( trackPassingTrkQualCut == 1){
      FillHelixHistograms(fHist.fHelix[107], helix);
    }
  
    // double clusterE = helix->ClusterEnergy();
    // double dE       = (clusterE - fTrackCaloE);
    // if ( (fTrackCaloE  >= 50.) && (clusterE >= 50.) && 
    // 	 (dE > -3.)            && (dE < 0)             ) {
    if (i  == fTrackHelixIndex){
      FillHelixHistograms(fHist.fHelix[108], helix);

      // if (trackPassingTrkQualCut ==1){
      // 	FillHelixHistograms(fHist.fHelix[109], helix);
      // }
      
    }
    
    double      pMC     = helix->Mom1().Vect().Mag();
    int         pdgm    = helix->PDGMother1();
    int         pdg     = helix->PDG1();
    if ( (pMC>=100.) && (pMC<=105) && (pdgm == 11) && (pdg == 11)) {
      FillHelixHistograms(fHist.fHelix[110], helix);

      if (i  == fTrackHelixIndex){
	FillHelixHistograms(fHist.fHelix[109], helix);
      }
      if (GetDebugBit(54))   GetHeaderBlock()->Print(Form(" bit:054 calPatRec helix"));
    }
    
  }
  
  for (int i=0; i<fNHelices[2]; ++i){
    
    helix = fCprDepHelixBlock->Helix(i);
    
    FillHelixHistograms(fHist.fHelix[200], helix);
    
    int         nhits    = helix->NHits();
    double      radius   = helix->Radius();
    double      lambda   = helix->Lambda();
    double      tanDip   = lambda/radius;
    double      mm2MeV   = 3/10.;
    double      p        = radius*mm2MeV/std::cos( std::atan(tanDip));
    double      chi2xy   = helix->Chi2XY();
    double      chi2zphi = helix->Chi2ZPhi();
    
    if (p > 80.) {
      FillHelixHistograms(fHist.fHelix[201], helix);
    }
    
    if (p > 90) {
      FillHelixHistograms(fHist.fHelix[202], helix);
    }

    if (p > 100) {
      FillHelixHistograms(fHist.fHelix[203], helix);
    }
    
    if ( (nhits>10) && (nhits<15)) {
      FillHelixHistograms(fHist.fHelix[204], helix);
    }

    if ( nhits>=15) {
      cprDepHelixTrig = true;
      FillHelixHistograms(fHist.fHelix[205], helix);
    }
    
    if ( (chi2xy <= 5) && (chi2zphi <= 5) && (nhits>=15)){
      FillHelixHistograms(fHist.fHelix[206], helix);
    }
  
    if (i  == fTrackDepHelixIndex){
      FillHelixHistograms(fHist.fHelix[208], helix);
    }


    double      pMC     = helix->Mom1().Vect().Mag();
    int         pdgm    = helix->PDGMother1();
    int         pdg     = helix->PDG1();
    if ( (pMC>=90.) && (pMC<=100) && (pdgm == -11) && (pdg == -11)) {
      FillHelixHistograms(fHist.fHelix[210], helix);
      
      if (i  == fTrackDepHelixIndex){
	FillHelixHistograms(fHist.fHelix[209], helix);
      }

      if (GetDebugBit(54))   GetHeaderBlock()->Print(Form(" bit:054 calPatRec helix"));
    }
  }
  //-----------------------------------------------------------------------------
  // track histograms, fill them only for the downstream e- hypothesis
  //-----------------------------------------------------------------------------
 
  //  TrackPar_t*    tp;

  //--------------------------------------------------------------------------------
  // timecluster histograms
  //--------------------------------------------------------------------------------
  TStnTimeCluster*tCluster(0);
  for (int i=0; i<fNTimeClusters[0]; ++i){
    
    tCluster = fTimeClusterBlock->TimeCluster(i);
    
    FillTimeClusterHistograms(fHist.fTimeCluster[0], tCluster);
    
    int         nhits    = tCluster->NHits();
    if (nhits >= 10 ) FillTimeClusterHistograms(fHist.fTimeCluster[1], tCluster);

    if (nhits >= 15 ) FillTimeClusterHistograms(fHist.fTimeCluster[2], tCluster);
  }  
  
  //--------------------------------------------------------------------------------
  // trackseed histograms
  //--------------------------------------------------------------------------------
  //TprSeedFit Dem
  for (int i=0; i<fNTrackSeeds[1]; ++i){
    
    trkSeed = fTprDemTrackSeedBlock->TrackSeed(i);
    if ( tprDemHelixTrig)     FillTrackSeedHistograms(fHist.fTrackSeed[115], trkSeed);

    FillTrackSeedHistograms(fHist.fTrackSeed[100], trkSeed);

    int         pdg1     = trkSeed->fSimpPDG1;
    int         pdgM1    = trkSeed->fSimpPDGM1;
    TLorentzVector p1 = trkSeed->fMom1;
    TLorentzVector x1 = trkSeed->fOrigin1;
    double      pMC     = p1.Vect().Mag();
    if ( (pMC>=100.) && (pMC<=105) && (pdgM1 == 11) && (pdg1 == 11)) {
      FillTrackSeedHistograms(fHist.fTrackSeed[120], trkSeed);
    }
  }

  //TprSeedFit Dep
  for (int i=0; i<fNTrackSeeds[3]; ++i){
    
    trkSeed = fTprDepTrackSeedBlock->TrackSeed(i);
    if ( tprDepHelixTrig)     FillTrackSeedHistograms(fHist.fTrackSeed[315], trkSeed);

    FillTrackSeedHistograms(fHist.fTrackSeed[300], trkSeed);
    double      pMC     = trkSeed->Mom1().Vect().Mag();
    int         pdgm    = trkSeed->PDGMother1();
    int         pdg     = trkSeed->PDG1();
    if ( (pMC>=90.) && (pMC<=100) && (pdgm == -11) && (pdg == -11)) {
      FillTrackSeedHistograms(fHist.fTrackSeed[320], trkSeed);
    }
  }


  //CalSeedFit  Dem
  for (int i=0; i<fNTrackSeeds[0]; ++i){
    
    trkSeed = fCprDemTrackSeedBlock->TrackSeed(i);
    
    FillTrackSeedHistograms(fHist.fTrackSeed[0], trkSeed);

    if ( cprDemHelixTrig)     FillTrackSeedHistograms(fHist.fTrackSeed[15], trkSeed);

    int         nhits    = trkSeed->NHits();
    double      p        = trkSeed->P();
    double      chi2     = trkSeed->Chi2()/(nhits-5.);
    int         pdg1     = trkSeed->fSimpPDG1;
    int         pdgM1    = trkSeed->fSimpPDGM1;
    TLorentzVector p1 = trkSeed->fMom1;
    TLorentzVector x1 = trkSeed->fOrigin1;
    double      pMC     = p1.Vect().Mag();
    if ( (pMC>=100.) && (pMC<=105) && (pdgM1 == 11) && (pdg1 == 11)) {
      FillTrackSeedHistograms(fHist.fTrackSeed[20], trkSeed);
    }

    if (GetDebugBit(43)) {
      //      int    id1   = trkSeed->fSimpId1;
      int    nh1   = trkSeed->fSimpId1Hits;
   
      int    nh2   = trkSeed->fSimpId2Hits;
      int    pdg2  = trkSeed->fSimpPDG2;
      int    pdgM2 = trkSeed->fSimpPDGM2;
      TLorentzVector p2 = trkSeed->fMom2;
      TLorentzVector x2 = trkSeed->fOrigin2;
      GetHeaderBlock()->Print(Form(" bit:043 CPR-TRKSEED p = %5.3f  nhits = %3i chi2 = %5.3f nh1 = %i p1 = %2.3f pt1 = %2.3f pdg1 = %i pdgM1 = %i x1=  %2.3f y1 = %2.3f z1 = %2.3f  nh2 = %i p2 = %2.3f pt2 = %2.3f pdg2 = %i pdgM2 = %i", 
				   p, nhits, chi2, nh1, p1.Vect().Mag(), p1.Vect().Perp(),
				   pdg1,pdgM1, 
				   x1.Vect().x(), x1.Vect().y(), x1.Vect().z(),
				   nh2, p2.Vect().Mag(), p2.Vect().Perp(), pdg2, pdgM2));
    }

    if (GetDebugBit(48) && (p>58.5) && (p<61.5)) {
      GetHeaderBlock()->Print(Form(" bit:48 nHits = %i p = %2.3f tanDip = %2.3f d0 = %2.3f chi2/ndof = %5.3f", 
				   nhits , p, trkSeed->D0(),  trkSeed->TanDip(), chi2));
    }
    
    if (GetDebugBit(47)) {
      GetHeaderBlock()->Print(Form(" bit:47 nHits = %i p = %2.3f tanDip = %2.3f d0 = %2.3f chi2/ndof = %5.3f", 
				   nhits , p, trkSeed->D0(),  trkSeed->TanDip(), chi2));
    }

    if (GetDebugBit(44) && (nhits > 30)) {
      double      clusterE = trkSeed->ClusterEnergy();
      GetHeaderBlock()->Print(Form(" bit:044 nHits = %i p = %2.3f E = %2.3f chi2/ndof = %2.3f", 
				   nhits, p, clusterE, chi2));
    }
    // if ( (chi2xy < 4) && (chi2zphi < 4)){
    //       trigSeedFlag = 1;
    //     }
    
    if ( (chi2 < 4) /*&& (chi2zphi < 4)*/ && (nhits>=15)){
      trigSeedFlag = 1;
      FillTrackSeedHistograms(fHist.fTrackSeed[7], trkSeed);
    }
    
    if (p > 80.) {
      FillTrackSeedHistograms(fHist.fTrackSeed[1], trkSeed);
    }
    
    if (p > 90) {
      FillTrackSeedHistograms(fHist.fTrackSeed[2], trkSeed);
    }

    if (p > 100) {
      FillTrackSeedHistograms(fHist.fTrackSeed[3], trkSeed);
    }
    
    if ( (nhits>10) && (nhits<15)) {
      FillTrackSeedHistograms(fHist.fTrackSeed[4], trkSeed);
    }

    if ( nhits>=15) {
      FillTrackSeedHistograms(fHist.fTrackSeed[5], trkSeed);
    }

    //    double clusterE = trkSeed->ClusterEnergy();
    //    double dE       = (clusterE - fTrackCaloE);
    // if ( (fTrackCaloE  >= 50.) && (clusterE >= 50.) && 
    // 	 (dE > -3.)            && (dE < 0)             ) {
    if (i == fTrackTrkSeedIndex){
      FillTrackSeedHistograms(fHist.fTrackSeed[9], trkSeed);
      if (trackPassingTrkQualCut ==1){
	FillTrackSeedHistograms(fHist.fTrackSeed[6], trkSeed);
	double      d0       = trkSeed->fD0;

	if ( (d0 > -200) && (d0< 200)) {
	  FillTrackSeedHistograms(fHist.fTrackSeed[10], trkSeed);
	} 

	if ( (chi2 < 4) /*&& (chi2zphi < 4)*/ && (nhits>=15)){
	  FillTrackSeedHistograms(fHist.fTrackSeed[8], trkSeed);
	}
      }
    }

  }


  //CalSeedFit Dep
  for (int i=0; i<fNTrackSeeds[2]; ++i){
    
    trkSeed = fCprDepTrackSeedBlock->TrackSeed(i);
    
    FillTrackSeedHistograms(fHist.fTrackSeed[200], trkSeed);
    if ( cprDepHelixTrig)     FillTrackSeedHistograms(fHist.fTrackSeed[215], trkSeed);

    double      pMC     = trkSeed->Mom1().Vect().Mag();
    int         pdgm    = trkSeed->PDGMother1();
    int         pdg     = trkSeed->PDG1();
    if ( (pMC>=90.) && (pMC<=100) && (pdgm == -11) && (pdg == -11)) {
      FillTrackSeedHistograms(fHist.fTrackSeed[220], trkSeed);
    }
  }
  // if (trackPassingTrkQualCut ==1){
  //   for (int i=0; i<fNTrackSeeds[0]; ++i){
  //     trkSeed = fCprDemTrackSeedBlock->TrackSeed(i);
  //     FillTrackSeedHistograms(fHist.fTrackSeed[6], trkSeed);

  //     //      int         nhits    = trkSeed->NHits();
  //     //      double      chi2     = trkSeed->Chi2()/(nhits-5.);
  //     //      double      chi2zphi = trkSeed->Chi2ZPhi();
  //     double      d0       = trkSeed->fD0;

  //     if ( (d0 > -200) && (d0< 200)) {
  // 	FillTrackSeedHistograms(fHist.fTrackSeed[10], trkSeed);
  //     } 
  //   }
  // }
  //       FillTrackHistograms(fHist.fTrack[1],trk);
      
  //       if (trk->fP >= 100. && trk->fP < 110.) {
  // 	FillTrackHistograms(fHist.fTrack[3],trk);
  //       }

  //     }

  // //-----------------------------------------------------------------------------
  // // split tracks by the algorithm mask: 1 , 2 , or 3
  // //-----------------------------------------------------------------------------
  // int     alg = trk->AlgMask();
  // if ((alg == 2) || (alg == 3)){
  //       if (trk->fTrkQual > 0.2) {
  // 	FillTrackHistograms(fHist.fTrack[17],trk);

  // 	if (trk->fP >= 100. && trk->fP < 110.) {
  // 	  FillTrackHistograms(fHist.fTrack[18],trk);

  // 	  if ( trigSeedFlag == 1){

  // 	    if (trigFlag1 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[19],trk);
  // 	    }
  // 	    if (trigFlag0 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[20],trk);
  // 	    }
	    
  // 	    if (trigFlag3 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[24],trk);
  // 	    }
  // 	    if (trigFlag2 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[25],trk);
  // 	    }
  // 	  }
  // 	}
  //       }
  //     }


  //     alg_mask = trk->BestAlg();//AlgMask();
  //     if      (alg_mask == 0) {
  // //-----------------------------------------------------------------------------
  // // TrkPatRec tracks
  // //-----------------------------------------------------------------------------
  //       if (trk->fTrkQual > 0.4) {
  // 	FillTrackHistograms(fHist.fTrack[2],trk);
	
  // 	if (trk->fP >= 100. && trk->fP < 110.) {
  // 	  FillTrackHistograms(fHist.fTrack[4],trk);

  // 	  trackPassingTrkQualCut = 1;
	  
  // 	  if ( trigSeedFlag == 1){
	    
  // 	    if (trigFlag1 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[15],trk);
  // 	    }
  // 	    if (trigFlag0 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[16],trk);
  // 	    }
	    
  // 	    if (trigFlag3 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[22],trk);
  // 	    }
  // 	    if (trigFlag2 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[23],trk);
  // 	    }
  // 	  }
  // 	}
  //       }
      
  //     }
  //     else if  (alg_mask == 1) {
  // //-----------------------------------------------------------------------------
  // // CalPatRec  tracks
  // //-----------------------------------------------------------------------------
  //       if (trk->fTrkQual > 0.2) {
  // 	FillTrackHistograms(fHist.fTrack[2],trk);

  // 	if (trk->fP >= 100. && trk->fP < 110.) {
  // 	  FillTrackHistograms(fHist.fTrack[4],trk);
	  
  // 	  trackPassingTrkQualCut = 1;

  // 	  if ( trigSeedFlag == 1){

  // 	    if (trigFlag1 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[15],trk);
  // 	    }
  // 	    if (trigFlag0 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[16],trk);
  // 	    }
	    
  // 	    if (trigFlag3 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[22],trk);
  // 	    }
  // 	    if (trigFlag2 == 1){
  // 	      FillTrackHistograms(fHist.fTrack[23],trk);
  // 	    }
  // 	  }
  // 	}
  //       }
  //     }
  //   }
    

  // if (trackPassingTrkQualCut ==1){
  //   for (int i=0; i<fNTrackSeeds[0]; ++i){
  //     trkSeed = fTrackSeedBlock->TrackSeed(i);
  //     FillTrackSeedHistograms(fHist.fTrackSeed[6], trkSeed);

  //     int         nhits    = trkSeed->NHits();
  //     double      chi2xy   = trkSeed->Chi2();//FIX ME!
  //     double      chi2zphi = trkSeed->Chi2();//FIX ME!

  //     if ( (chi2xy < 4) && (chi2zphi < 4) && (nhits>=15)){
  // 	FillTrackSeedHistograms(fHist.fTrackSeed[8], trkSeed);
  //     }
  //   }
  // }



  //--------------------------------------------------------------------------------  
  // CalTrkFit (trigger "seddFit")
  //--------------------------------------------------------------------------------
  for (int i=0; i<fNTracks[1]; ++i ) {
    trk = fCprDemTrackBlock->Track(i);

    int        nhits = trk->NHits();
    double     chi2c = trk->fChi2/(nhits-5.);
    double      d0   = trk->fD0;

    if  ( trigSeedFlag == 1){
      FillTrackHistograms(fHist.fTrack[112],trk);
      if ( (nhits >=15) && (chi2c < 4) && (fabs(d0) <200.)){
	FillTrackHistograms(fHist.fTrack[114],trk);
      }
    }
    
    double chi2dof    = trk->Chi2Dof();
    double t0         = trk->fT0;
    //    double t0_err     = trk->fT0Err;
    double nactive    = trk->NActive();
    double tan_dip    = trk->fTanDip;
    double fitmom_err = trk->fFitMomErr;
    double p          = trk->P();
    
    if ( (chi2dof < fMinChi2d  ) &&  
	 (tan_dip < fMaxTanDip ) &&
	 (tan_dip > fMinTanDip ) &&
	 (nactive > fMinNActive) &&
	 (fitmom_err< fMaxFitMomErr) &&
	 (d0 < fMaxD1)           &&
	 (d0 > fMinD1)           && 
	 (t0 > fMinT0 )          &&
	 ((p>=fPMin) && (p<=fPMax)) ){
      FillTrackHistograms(fHist.fTrack[110],trk);
           
      if  ( trigSeedFlag == 1){
	FillTrackHistograms(fHist.fTrack[111],trk);
	if ( (nhits >=15) && (chi2c < 4) && (fabs(d0) <200.)){
	  FillTrackHistograms(fHist.fTrack[113],trk);
	}
      }
    }
    
    if ((p>=fPMin) && (p<=fPMax)) {
      if ( (nactive >=20) && (chi2dof < 5) ) FillTrackHistograms(fHist.fTrack[115],trk);
      if ( (nactive >=20) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[116],trk);
      if ( (nactive >=25) && (chi2dof < 5) ) {
	if (GetDebugBit(47)) {
	  GetHeaderBlock()->Print(Form(" bit:47 nactive = %3.0f p = %5.3f tanDip = %5.3f chi2/ndof = %5.3f", 
				       nactive , p, tan_dip, chi2dof));
	}
	FillTrackHistograms(fHist.fTrack[117],trk);
      }
      if ( (nactive >=25) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[118],trk);
    }
  }
  for (int i=0; i<fNTracks[3]; ++i ) {
    trk = fCprDepTrackBlock->Track(i);

    int        nhits = trk->NHits();
    double     chi2c = trk->fChi2/(nhits-5.);
    double      d0   = trk->fD0;

    if  ( trigSeedFlag == 1){
      FillTrackHistograms(fHist.fTrack[312],trk);
      if ( (nhits >=15) && (chi2c < 4) && (fabs(d0) <200.)){
	FillTrackHistograms(fHist.fTrack[314],trk);
      }
    }
    
    double chi2dof    = trk->Chi2Dof();
    double t0         = trk->fT0;
    //    double t0_err     = trk->fT0Err;
    double nactive    = trk->NActive();
    double tan_dip    = trk->fTanDip;
    double fitmom_err = trk->fFitMomErr;
    double p          = trk->P();
    
    if ( (chi2dof < fMinChi2d  ) &&  
	 (tan_dip < fMaxTanDip ) &&
	 (tan_dip > fMinTanDip ) &&
	 (nactive > fMinNActive) &&
	 (fitmom_err< fMaxFitMomErr) &&
	 (d0 < fMaxD1)           &&
	 (d0 > fMinD1)           && 
	 (t0 > fMinT0 )          &&
	 ((p>=fPMinDep) && (p<=fPMaxDep)) ){
      FillTrackHistograms(fHist.fTrack[310],trk);
           
      if  ( trigSeedFlag == 1){
	FillTrackHistograms(fHist.fTrack[311],trk);
	if ( (nhits >=15) && (chi2c < 4) && (fabs(d0) <200.)){
	  FillTrackHistograms(fHist.fTrack[313],trk);
	}
      }
    }
    
    if ((p>=fPMinDep) && (p<=fPMaxDep)) {
      if ( (nactive >=20) && (chi2dof < 5) ) FillTrackHistograms(fHist.fTrack[315],trk);
      if ( (nactive >=20) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[316],trk);
      if ( (nactive >=25) && (chi2dof < 5) ) {
	FillTrackHistograms(fHist.fTrack[317],trk);
      }
      if ( (nactive >=25) && (chi2dof < 3) ) FillTrackHistograms(fHist.fTrack[318],trk);
    }
  }
  //-----------------------------------------------------------------------------
  // cluster histograms
  //-----------------------------------------------------------------------------
  TStnCluster* cl;
  int id;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    id = cl->DiskID();
    FillClusterHistograms(fHist.fCluster[0],cl);

    if (fNTracks[0]     >  0 ) FillClusterHistograms(fHist.fCluster[1],cl);
    if (fNGoodTracks    >  0 ) FillClusterHistograms(fHist.fCluster[2],cl);
    if (fNMatchedTracks >  0 ) FillClusterHistograms(fHist.fCluster[3],cl);
    if (cl->Energy()    > 10.) FillClusterHistograms(fHist.fCluster[4],cl);
    if (cl->Energy()    > 60.) FillClusterHistograms(fHist.fCluster[5],cl);

    if      (id == 0         ) FillClusterHistograms(fHist.fCluster[6],cl);
    else if (id == 1         ) FillClusterHistograms(fHist.fCluster[7],cl);
  }
  //-----------------------------------------------------------------------------
  // calorimeter histograms
  //-----------------------------------------------------------------------------
  TDisk*         disk;
  TStnCrystal*   cr;

  if (fCalorimeterType == 2) {
    int nd = fDiskCalorimeter->NDisks();

    for (int i=0; i<nd; i++) {
      disk = fDiskCalorimeter->Disk(i);
      for (int ic=0; ic<disk->NCrystals(); ic++) {
	cr = disk->Crystal(ic);
	FillCaloHistograms(fHist.fCalo[0],cr);

	if (cr->Energy() > 0) {
	  FillCaloHistograms(fHist.fCalo[1],cr);
	}
	if (cr->Energy() > 0.1) {
	  FillCaloHistograms(fHist.fCalo[2],cr);
	}
	if (cr->Energy() > 1.0) {
	  FillCaloHistograms(fHist.fCalo[3],cr);
	}
      }
    }
  }
  //-----------------------------------------------------------------------------
  // radial distributions for crystals
  //-----------------------------------------------------------------------------
  static int first_entry(1);

  if (first_entry == 1) {
    if (fCalorimeterType == 2) {
      int nd = fDiskCalorimeter->NDisks();
	
      for (int i=0; i<nd; i++) {
	disk = fDiskCalorimeter->Disk(i);
	for (int ic=0; ic<disk->NCrystals(); ic++) {
	  cr = disk->Crystal(ic);

	  fHist.fCrystalR[i]->Fill(cr->Radius());
	}
      }
    }
  }

  //-----------------------------------------------------------------------------
  // fill GENP histograms
  // GEN_0: all particles
  //-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
    if (genp->GetPdgCode() == 11)  // track is made by an electron 
      FillGenpHistograms(fHist.fGenp[1],genp);
    if (genp->GetPdgCode() == -11)  // track is made by a positron
      FillGenpHistograms(fHist.fGenp[2],genp);
    if (genp->GetPdgCode() == 13)  // track is made by a mu minus 
      FillGenpHistograms(fHist.fGenp[3],genp);
    if (genp->GetPdgCode() == -13)  // track is made by a mu plus 
      FillGenpHistograms(fHist.fGenp[4],genp);
    if (genp->GetPdgCode() == -211)  // track is made by a pi minus 
      FillGenpHistograms(fHist.fGenp[5],genp);
    if (genp->GetPdgCode() == 211)  // track is made by a pi plus
      FillGenpHistograms(fHist.fGenp[6],genp);
    if (genp->GetPdgCode() == -321)  // track is made by a kaon minus	
      FillGenpHistograms(fHist.fGenp[7],genp);
    if (genp->GetPdgCode() == 321)  // track is made by a kaon plus
      FillGenpHistograms(fHist.fGenp[8],genp);
    if (genp->GetPdgCode() == 22)  // track is made by a photon
      FillGenpHistograms(fHist.fGenp[9],genp);
  }
  first_entry = 0;
  
}


//----------------------------------------------------------------------
// 2014-07-21: normalize the histograms for getting the
//             expected numbers after three years of run
//----------------------------------------------------------------------
void    TTrigAna001Module::Normalize() {
  double simEvents= fEventCounter;//fHist.fGenp[0]->fP->GetEntries();// = 9995.*1e3; //each data file used 1e3 events
 
  double norm = 1./fEventCounter;

  printf("[TTrigAna001Module::Normalize] simulated events = %3.3e  \n", simEvents);
  printf("[TTrigAna001Module::Normalize] normalization factor = %5.5e\n", norm);

  // Normalize the GenBlock histograms
  for (int i=0; i< kNEventHistSets; ++i) {
    if (fBookedEvent[i] == 0) goto NEXT_EVN;
    fHist.fEvent[i]->fCprDemTrigger[0]->Scale(norm);
    fHist.fEvent[i]->fCprDemTrigger[1]->Scale(norm);
    fHist.fEvent[i]->fCprDemTrigger[2]->Scale(norm);
    fHist.fEvent[i]->fCprDemTrigger[3]->Scale(norm);
    fHist.fEvent[i]->fCprDemTrigger[4]->Scale(norm);
    fHist.fEvent[i]->fCprDemTrigger[5]->Scale(norm);
    fHist.fEvent[i]->fCprDemTrigger[6]->Scale(norm);
    fHist.fEvent[i]->fCprDemTrigger[7]->Scale(norm);

  NEXT_EVN:;
  }

  // Normalize the GenBlock histograms
  for (int i=0; i< kNGenpHistSets; ++i) {
    if (fBookedGenps[i] == 0) goto NEXT_GEN;
    fHist.fGenp[i]->fPNorm->Scale(norm);
    fHist.fGenp[i]->fCosThNorm->Scale(norm);
  NEXT_GEN:;
  }


  //Normalize the TrackBlock
  for (int i=0; i< kNTrackHistSets ; ++i) {
    if (fBookedTracks[i] == 0) goto NEXT_HISTO;
    fHist.fTrack[i]->fPtNorm->Scale(norm);
    fHist.fTrack[i]->fPNorm->Scale(norm);
    fHist.fTrack[i]->fFitMomErrNorm->Scale(norm);
    fHist.fTrack[i]->fChi2Norm->Scale(norm);
    fHist.fTrack[i]->fT0Norm->Scale(norm);
    fHist.fTrack[i]->fT0ErrNorm->Scale(norm);
    fHist.fTrack[i]->fFitConsNorm[0]->Scale(norm);
    fHist.fTrack[i]->fFitConsNorm[1]->Scale(norm);
    fHist.fTrack[i]->fTanDipNorm->Scale(norm);
    fHist.fTrack[i]->fvDxNorm->Scale(norm);
    fHist.fTrack[i]->fvDyNorm->Scale(norm);
    fHist.fTrack[i]->fvDzNorm->Scale(norm);
    fHist.fTrack[i]->fvDpNorm->Scale(norm);
    fHist.fTrack[i]->fvDtNorm->Scale(norm);
    fHist.fTrack[i]->fvDt0Norm->Scale(norm);
    fHist.fTrack[i]->fvtofNorm->Scale(norm);
    fHist.fTrack[i]->fvDrNorm->Scale(norm);
    fHist.fTrack[i]->fvDtVsDt0Norm->Scale(norm);
    fHist.fTrack[i]->fvDtVsT0errNorm->Scale(norm);
    fHist.fTrack[i]->fvDtVsDrNorm->Scale(norm);
    fHist.fTrack[i]->fvDtVsDpNorm->Scale(norm);
    fHist.fTrack[i]->fvtofVsDpNorm->Scale(norm);
  NEXT_HISTO:;
  }
  printf("[TTrigAna001Module::Normalize] normalization done...\n");

}


// funtion that evaluates the background rejection and trigger efficiency
// from the Calo-seeded algorithm
void TTrigAna001Module::FillTprDemTriggerHist(EventHist_t*   Hist, int &hasQualityTprTrack){
  TStnHelix*      helix(0);
  int  nCuts(11);
  int  nhits_cut       [nCuts] = {0};
  int  triggerCounter  [nCuts] = {0};
  int  triggerCounter1 [nCuts] = {0};
  int  triggerCounter2 [nCuts] = {0};
  int  triggerCounter3 [nCuts] = {0};
  int  triggerCounter4 [nCuts] = {0};
  int  triggerCounter5 [nCuts] = {0};
  int  triggerCounter6 [nCuts] = {0};
  int  triggerCounter8 [nCuts] = {0};
  int  triggerCounter9 [nCuts] = {0};
  int  triggerCounter10[nCuts] = {0};
  int  triggerCounter11[nCuts] = {0};
  int  triggerCounter12[nCuts] = {0};

  for (int i=0; i<nCuts; ++i){ nhits_cut[i] = 10 + i; }

  int  nCuts2(8);
  int  triggerCounter7 [nCuts2] = {0};
  int  triggerCounter13[nCuts] = {0};
  int  triggerCounter14[nCuts] = {0};
  int  triggerCounter15[nCuts] = {0};
		       
  double clE_cut[nCuts2]       = {0};
  
  bool  hasTrigHelix(false);
  
  for (int i=0; i<nCuts2; ++i){ clE_cut[i] = 40 + i*5.; }

  TStnTrackSeed*    trkSeed(0);

  //effciecy of the tracker only pattern-recognition
  for (int i=0; i<fNHelices[0]; ++i){
    
    helix   = fTprDemHelixBlock->Helix(i);
    int         nhits    = helix->NHits();
    double      chi2XY   = helix->Chi2XY();
    double      chi2ZPhi = helix->Chi2ZPhi();
    
    if ( (chi2XY > 5.) || (chi2ZPhi>5.))        continue;
    
    if (nhits >= fMinTrigHelixNSh)   hasTrigHelix = true;

    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	triggerCounter13[j] = 1;
	
	if (hasQualityTprTrack == 1){
	  triggerCounter14[j] = 1;
	}
	if ( (chi2XY < 4.) && ( chi2ZPhi < 4)){
	  //	  chi2Condition      = 1;
	  
	  if (hasQualityTprTrack == 1){
	    triggerCounter15[j] = 1;
	  }
	}
      }
    }//end loop on cuts
  }//end loop over the helices

  //fill counter used for normalization
  for (int i=0; i<fNTrackSeeds[1]; ++i){
    
    trkSeed   = fTprDemTrackSeedBlock->TrackSeed(i);
    int         nhits    = trkSeed->NHits();

    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	if (hasQualityTprTrack == 1){
	  triggerCounter12[j] = 1;
	}	
      }
    }    
  }
  
  if (hasTrigHelix){
    //TrlPatRec Dem
    for (int i=0; i<fNTrackSeeds[1]; ++i){
    
      trkSeed   = fTprDemTrackSeedBlock->TrackSeed(i);
      int         nhits    = trkSeed->NHits();
      double      chi2     = trkSeed->Chi2()/(nhits-5.);
    
      int         helixIndex  = trkSeed->HelixIndex();

      if (helixIndex>=0 && helixIndex<fNHelices[0] ){
	helix                   = fTprDemHelixBlock->Helix(helixIndex);
	int         helix_nhits = helix->NHits();
	int         delta_nhits = nhits - helix_nhits;
	Hist->fDNHits->Fill(delta_nhits);
      }

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter[j] = 1;
	  // if (hasQualityTprTrack == 1){
	  //   triggerCounter12[j] = 1;
	  // }
	
	  if ( (chi2 < 5.) /*&& ( chi2ZPhi < 4)*/){
	    triggerCounter1[j] = 1;
	  }
	}
      }    
    }


    //TrkPatRec Dem
    int ntrk = fNTrackSeeds[1];
    for (int i=0; i<ntrk; ++i){
      trkSeed = fTprDemTrackSeedBlock->TrackSeed(i);
    
      int         nhits = trkSeed->NHits();
      double      d0    = trkSeed->fD0;
      double      chi2  = trkSeed->fChi2/(nhits-5.);
      double      p     = trkSeed->P();
      double      p_min(80.);

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter2[j] = 1;
	  
	  // if (hasQualityTrack == 1){
	  //   triggerCounter11[j] = 1;
	  // }
	  
	  if ( (d0>-200) && (d0<200)){
	    triggerCounter3[j] = 1;
	  }

	  if ( (d0>-200) && (d0<200) &&
	       (hasQualityTprTrack == 1)){
	    triggerCounter9[j] = 1;
	  }	

	  if (p > p_min){
	    triggerCounter4[j] = 1;
	  
	    if (hasQualityTprTrack == 1){
	      triggerCounter8[j] = 1;
	    }
	    if ( (d0>=-200) && (d0<=200)){
	      triggerCounter5[j] = 1;

	      if (hasQualityTprTrack == 1){
		triggerCounter11[j] = 1;
	      }

	      if (chi2 < 3.){
		triggerCounter6[j] = 1;
		if (hasQualityTprTrack == 1){
		  triggerCounter10[j] = 1;
		}
	      }
	    }
	  }
	}
      }
    }

  }

  for (int i=0; i<nCuts; ++i){
    double   content = Hist->fTprDemTrigger[0]->GetBinContent(i+1);
    content += triggerCounter[i];
    Hist->fTprDemTrigger[0]->SetBinContent(i+1, content);
    
    content = Hist->fTprDemTrigger[1]->GetBinContent(i+1);
    content += triggerCounter1[i];
    Hist->fTprDemTrigger[1]->SetBinContent(i+1, content);
    
    content = Hist->fTprDemTrigger[2]->GetBinContent(i+1);
    content += triggerCounter2[i];
    Hist->fTprDemTrigger[2]->SetBinContent(i+1, content);
    
    content = Hist->fTprDemTrigger[3]->GetBinContent(i+1);
    content += triggerCounter3[i];
    Hist->fTprDemTrigger[3]->SetBinContent(i+1, content);

    content = Hist->fTprDemTrigger[4]->GetBinContent(i+1);
    content += triggerCounter4[i];
    Hist->fTprDemTrigger[4]->SetBinContent(i+1, content);
    
    content = Hist->fTprDemTrigger[5]->GetBinContent(i+1);
    content += triggerCounter5[i];
    Hist->fTprDemTrigger[5]->SetBinContent(i+1, content);
    
    content = Hist->fTprDemTrigger[6]->GetBinContent(i+1);
    content += triggerCounter6[i];
    Hist->fTprDemTrigger[6]->SetBinContent(i+1, content);
    
    content = Hist->fTprDemTrigger[8]->GetBinContent(i+1);
    content += triggerCounter8[i];
    Hist->fTprDemTrigger[8]->SetBinContent(i+1, content);

    content = Hist->fTprDemTrigger[9]->GetBinContent(i+1);
    content += triggerCounter9[i];
    Hist->fTprDemTrigger[9]->SetBinContent(i+1, content);

    content = Hist->fTprDemTrigger[10]->GetBinContent(i+1);
    content += triggerCounter10[i];
    Hist->fTprDemTrigger[10]->SetBinContent(i+1, content);

    content = Hist->fTprDemTrigger[11]->GetBinContent(i+1);
    content += triggerCounter11[i];
    Hist->fTprDemTrigger[11]->SetBinContent(i+1, content);

    content = Hist->fTprDemTrigger[12]->GetBinContent(i+1);
    content += triggerCounter12[i];
    Hist->fTprDemTrigger[12]->SetBinContent(i+1, content);


    content = Hist->fTprDemTrigger[13]->GetBinContent(i+1);
    content += triggerCounter13[i];
    Hist->fTprDemTrigger[13]->SetBinContent(i+1, content);


    content = Hist->fTprDemTrigger[14]->GetBinContent(i+1);
    content += triggerCounter14[i];
    Hist->fTprDemTrigger[14]->SetBinContent(i+1, content);


    content = Hist->fTprDemTrigger[15]->GetBinContent(i+1);
    content += triggerCounter15[i];
    Hist->fTprDemTrigger[15]->SetBinContent(i+1, content);

  }

  
  for (int i=0; i<nCuts2; ++i){
    double   content = Hist->fTprDemTrigger[7]->GetBinContent(i+1);
    content += triggerCounter7[i];
    Hist->fTprDemTrigger[7]->SetBinContent(i+1, content);
    
  }

}



// funtion that evaluates the background rejection and trigger efficiency
// from the Calo-seeded algorithm
void TTrigAna001Module::FillCprDemTriggerHist(EventHist_t*   Hist, int &hasQualityTrack){
  TStnHelix*      helix(0);
  TStnTrackSeed*    trkSeed(0);

  int  nCuts(11);
  int  nhits_cut       [nCuts] = {0};
  int  triggerCounter  [nCuts] = {0};
  int  triggerCounter1 [nCuts] = {0};
  int  triggerCounter2 [nCuts] = {0};
  int  triggerCounter3 [nCuts] = {0};
  int  triggerCounter4 [nCuts] = {0};
  int  triggerCounter5 [nCuts] = {0};
  int  triggerCounter6 [nCuts] = {0};
  int  triggerCounter8 [nCuts] = {0};
  int  triggerCounter9 [nCuts] = {0};
  int  triggerCounter10[nCuts] = {0};
  int  triggerCounter11[nCuts] = {0};
  int  triggerCounter12[nCuts] = {0};

  for (int i=0; i<nCuts; ++i){ nhits_cut[i] = 10 + i; }

  int  nCuts2(8);
  int  triggerCounter7 [nCuts2] = {0};
  int  triggerCounter13[nCuts] = {0};
  int  triggerCounter14[nCuts] = {0};
  int  triggerCounter15[nCuts] = {0};
		       
  double clE_cut[nCuts2]       = {0};
  bool   hasTrigHelix(false);


  for (int i=0; i<nCuts2; ++i){ clE_cut[i] = 40 + i*5.; }



  //  int  chi2Condition(0);

  //effciecy of the pattern-recognition
  for (int i=0; i<fNHelices[1]; ++i){
    
    helix   = fCprDemHelixBlock->Helix(i);
    int         nhits    = helix->NHits();
    double      chi2XY   = helix->Chi2XY();
    double      chi2ZPhi = helix->Chi2ZPhi();

    if (nhits >= fMinTrigHelixNSh)   hasTrigHelix = true;

    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	triggerCounter13[j] = 1;
	
	if (hasQualityTrack == 1){
	  triggerCounter14[j] = 1;
	}
	if ( (chi2XY < 4.) && ( chi2ZPhi < 4)){
	  //	  chi2Condition      = 1;
	  
	  if (hasQualityTrack == 1){
	    triggerCounter15[j] = 1;
	  }
	}
      }
    }//end loop
  }

  //eval normalzation counter
  for (int i=0; i<fNTrackSeeds[0]; ++i){
    
    trkSeed   = fCprDemTrackSeedBlock->TrackSeed(i);
    int         nhits    = trkSeed->NHits();    

    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	if (hasQualityTrack == 1){
	  triggerCounter12[j] = 1;
	}
	
      }

    }
  }
  

  //--------------------------------------------------------------------------------
  if (hasTrigHelix) {
    
    //calPatREc
    for (int i=0; i<fNTrackSeeds[0]; ++i){
    
      trkSeed   = fCprDemTrackSeedBlock->TrackSeed(i);
      int         nhits    = trkSeed->NHits();
      double      chi2     = trkSeed->Chi2()/(nhits-5.);
    
      int         helixIndex  = trkSeed->HelixIndex();

      if (helixIndex>=0 && helixIndex<fNHelices[1] ){
	helix                   = fCprDemHelixBlock->Helix(helixIndex);
	int         helix_nhits = helix->NHits();
	int         delta_nhits = nhits - helix_nhits;
	Hist->fDNHits->Fill(delta_nhits);
      }

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter[j] = 1;
	  // if (hasQualityTrack == 1){
	  //   triggerCounter12[j] = 1;
	  // }
	
	  if ( (chi2 < 5.) /*&& ( chi2ZPhi < 4)*/){
	    triggerCounter1[j] = 1;
	    //	  chi2Condition      = 1;
	  
	    // if (hasQualityTrack == 1){
	    //   triggerCounter8[j] = 1;
	    // }
	  }
	}

      }



      double clE =  trkSeed->ClusterEnergy();
      for (int j=0; j<nCuts2; ++j){
	if (clE>= clE_cut[j]){
	  triggerCounter7[j] = 1;
	}
      }
    }


    //CalPatRec
    int ntrk = fNTrackSeeds[0];
    for (int i=0; i<ntrk; ++i){
      trkSeed = fCprDemTrackSeedBlock->TrackSeed(i);
    
      int         nhits = trkSeed->NHits();
      double      d0    = trkSeed->fD0;
      double      chi2  = trkSeed->fChi2/(nhits-5.);
      double      p     = trkSeed->P();
      double      p_min(80.);

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter2[j] = 1;
	  
	  // if (hasQualityTrack == 1){
	  //   triggerCounter11[j] = 1;
	  // }
	  
	  if ( (d0>-200) && (d0<200)){
	    triggerCounter3[j] = 1;
	  }

	  if ( (d0>-200) && (d0<200) &&
	       (hasQualityTrack == 1)){
	    triggerCounter9[j] = 1;
	  }	

	  if (p > p_min){
	    triggerCounter4[j] = 1;
	  
	    if (hasQualityTrack == 1){
	      triggerCounter8[j] = 1;
	    }
	    if ( (d0>=-200) && (d0<=200)){
	      triggerCounter5[j] = 1;

	      if (hasQualityTrack == 1){
		triggerCounter11[j] = 1;
	      }

	      if (chi2 < 3.){
		triggerCounter6[j] = 1;
		if (hasQualityTrack == 1){
		  triggerCounter10[j] = 1;
		}
	      }
	    }
	  }
	}
      }
    }
  
  }

  for (int i=0; i<nCuts; ++i){
    double   content = Hist->fCprDemTrigger[0]->GetBinContent(i+1);
    content += triggerCounter[i];
    Hist->fCprDemTrigger[0]->SetBinContent(i+1, content);
    
    content = Hist->fCprDemTrigger[1]->GetBinContent(i+1);
    content += triggerCounter1[i];
    Hist->fCprDemTrigger[1]->SetBinContent(i+1, content);
    
    content = Hist->fCprDemTrigger[2]->GetBinContent(i+1);
    content += triggerCounter2[i];
    Hist->fCprDemTrigger[2]->SetBinContent(i+1, content);
    
    content = Hist->fCprDemTrigger[3]->GetBinContent(i+1);
    content += triggerCounter3[i];
    Hist->fCprDemTrigger[3]->SetBinContent(i+1, content);

    content = Hist->fCprDemTrigger[4]->GetBinContent(i+1);
    content += triggerCounter4[i];
    Hist->fCprDemTrigger[4]->SetBinContent(i+1, content);
    
    content = Hist->fCprDemTrigger[5]->GetBinContent(i+1);
    content += triggerCounter5[i];
    Hist->fCprDemTrigger[5]->SetBinContent(i+1, content);
    
    content = Hist->fCprDemTrigger[6]->GetBinContent(i+1);
    content += triggerCounter6[i];
    Hist->fCprDemTrigger[6]->SetBinContent(i+1, content);
    
    content = Hist->fCprDemTrigger[8]->GetBinContent(i+1);
    content += triggerCounter8[i];
    Hist->fCprDemTrigger[8]->SetBinContent(i+1, content);

    content = Hist->fCprDemTrigger[9]->GetBinContent(i+1);
    content += triggerCounter9[i];
    Hist->fCprDemTrigger[9]->SetBinContent(i+1, content);

    content = Hist->fCprDemTrigger[10]->GetBinContent(i+1);
    content += triggerCounter10[i];
    Hist->fCprDemTrigger[10]->SetBinContent(i+1, content);

    content = Hist->fCprDemTrigger[11]->GetBinContent(i+1);
    content += triggerCounter11[i];
    Hist->fCprDemTrigger[11]->SetBinContent(i+1, content);

    content = Hist->fCprDemTrigger[12]->GetBinContent(i+1);
    content += triggerCounter12[i];
    Hist->fCprDemTrigger[12]->SetBinContent(i+1, content);


    content = Hist->fCprDemTrigger[13]->GetBinContent(i+1);
    content += triggerCounter13[i];
    Hist->fCprDemTrigger[13]->SetBinContent(i+1, content);


    content = Hist->fCprDemTrigger[14]->GetBinContent(i+1);
    content += triggerCounter14[i];
    Hist->fCprDemTrigger[14]->SetBinContent(i+1, content);


    content = Hist->fCprDemTrigger[15]->GetBinContent(i+1);
    content += triggerCounter15[i];
    Hist->fCprDemTrigger[15]->SetBinContent(i+1, content);

  }

  
  for (int i=0; i<nCuts2; ++i){
    double   content = Hist->fCprDemTrigger[7]->GetBinContent(i+1);
    content += triggerCounter7[i];
    Hist->fCprDemTrigger[7]->SetBinContent(i+1, content);
    
  }
}

void TTrigAna001Module::FillTprDepTriggerHist(EventHist_t*   Hist, int &hasQualityTprTrack){
  TStnHelix*      helix(0);
  int  nCuts(11);
  int  nhits_cut       [nCuts] = {0};
  int  triggerCounter  [nCuts] = {0};
  int  triggerCounter1 [nCuts] = {0};
  int  triggerCounter2 [nCuts] = {0};
  int  triggerCounter3 [nCuts] = {0};
  int  triggerCounter4 [nCuts] = {0};
  int  triggerCounter5 [nCuts] = {0};
  int  triggerCounter6 [nCuts] = {0};
  int  triggerCounter8 [nCuts] = {0};
  int  triggerCounter9 [nCuts] = {0};
  int  triggerCounter10[nCuts] = {0};
  int  triggerCounter11[nCuts] = {0};
  int  triggerCounter12[nCuts] = {0};

  for (int i=0; i<nCuts; ++i){ nhits_cut[i] = 10 + i; }

  int  nCuts2(8);
  int  triggerCounter7 [nCuts2] = {0};
  int  triggerCounter13[nCuts] = {0};
  int  triggerCounter14[nCuts] = {0};
  int  triggerCounter15[nCuts] = {0};
		       
  double clE_cut[nCuts2]       = {0};
  
  bool   hasTrigHelix(false);

  for (int i=0; i<nCuts2; ++i){ clE_cut[i] = 40 + i*5.; }

  TStnTrackSeed*    trkSeed(0);

  //effciecy of the tracker only pattern-recognition
  for (int i=0; i<fNHelices[3]; ++i){
    
    helix   = fTprDepHelixBlock->Helix(i);
    int         nhits    = helix->NHits();
    double      chi2XY   = helix->Chi2XY();
    double      chi2ZPhi = helix->Chi2ZPhi();
    if ( (chi2XY > 5.) || (chi2ZPhi>5.))        continue;

    if (nhits >= fMinTrigHelixNSh)   hasTrigHelix = true;
  
    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	triggerCounter13[j] = 1;
	
	if (hasQualityTprTrack == 1){
	  triggerCounter14[j] = 1;
	}
	if ( (chi2XY < 4.) && ( chi2ZPhi < 4)){
	  //	  chi2Condition      = 1;
	  
	  if (hasQualityTprTrack == 1){
	    triggerCounter15[j] = 1;
	  }
	}
      }
    }//end loop
  }

  for (int i=0; i<fNTrackSeeds[3]; ++i){
    
    trkSeed   = fTprDepTrackSeedBlock->TrackSeed(i);
    int         nhits    = trkSeed->NHits();

    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	if (hasQualityTprTrack == 1){
	  triggerCounter12[j] = 1;
	}
	
      }
    }    
  }

  if (hasTrigHelix){

    //TrlPatRec Dep
    for (int i=0; i<fNTrackSeeds[3]; ++i){
    
      trkSeed   = fTprDepTrackSeedBlock->TrackSeed(i);
      int         nhits    = trkSeed->NHits();
      double      chi2     = trkSeed->Chi2()/(nhits-5.);
    
      int         helixIndex  = trkSeed->HelixIndex();

      if (helixIndex>=0 && helixIndex<fNHelices[3] ){
	helix                   = fTprDepHelixBlock->Helix(helixIndex);
	int         helix_nhits = helix->NHits();
	int         delta_nhits = nhits - helix_nhits;
	Hist->fDNHits->Fill(delta_nhits);
      }

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter[j] = 1;
	  // if (hasQualityTprTrack == 1){
	  //   triggerCounter12[j] = 1;
	  // }
	
	  if ( (chi2 < 5.) /*&& ( chi2ZPhi < 4)*/){
	    triggerCounter1[j] = 1;
	  }
	}
      }    
    }


    //TrkPatRec Dep
    int ntrk = fNTrackSeeds[3];
    for (int i=0; i<ntrk; ++i){
      trkSeed = fTprDepTrackSeedBlock->TrackSeed(i);
    
      int         nhits = trkSeed->NHits();
      double      d0    = trkSeed->fD0;
      double      chi2  = trkSeed->fChi2/(nhits-5.);
      double      p     = trkSeed->P();
      double      p_min(80.);

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter2[j] = 1;
	  
	  // if (hasQualityTrack == 1){
	  //   triggerCounter11[j] = 1;
	  // }
	  
	  if ( (d0>-200) && (d0<200)){
	    triggerCounter3[j] = 1;
	  }

	  if ( (d0>-200) && (d0<200) &&
	       (hasQualityTprTrack == 1)){
	    triggerCounter9[j] = 1;
	  }	

	  if (p > p_min){
	    triggerCounter4[j] = 1;
	  
	    if (hasQualityTprTrack == 1){
	      triggerCounter8[j] = 1;
	    }
	    if ( (d0>=-200) && (d0<=200)){
	      triggerCounter5[j] = 1;

	      if (hasQualityTprTrack == 1){
		triggerCounter11[j] = 1;
	      }

	      if (chi2 < 3.){
		triggerCounter6[j] = 1;
		if (hasQualityTprTrack == 1){
		  triggerCounter10[j] = 1;
		}
	      }
	    }
	  }
	}
      }
    }

  }

  for (int i=0; i<nCuts; ++i){
    double   content = Hist->fTprDepTrigger[0]->GetBinContent(i+1);
    content += triggerCounter[i];
    Hist->fTprDepTrigger[0]->SetBinContent(i+1, content);
    
    content = Hist->fTprDepTrigger[1]->GetBinContent(i+1);
    content += triggerCounter1[i];
    Hist->fTprDepTrigger[1]->SetBinContent(i+1, content);
    
    content = Hist->fTprDepTrigger[2]->GetBinContent(i+1);
    content += triggerCounter2[i];
    Hist->fTprDepTrigger[2]->SetBinContent(i+1, content);
    
    content = Hist->fTprDepTrigger[3]->GetBinContent(i+1);
    content += triggerCounter3[i];
    Hist->fTprDepTrigger[3]->SetBinContent(i+1, content);

    content = Hist->fTprDepTrigger[4]->GetBinContent(i+1);
    content += triggerCounter4[i];
    Hist->fTprDepTrigger[4]->SetBinContent(i+1, content);
    
    content = Hist->fTprDepTrigger[5]->GetBinContent(i+1);
    content += triggerCounter5[i];
    Hist->fTprDepTrigger[5]->SetBinContent(i+1, content);
    
    content = Hist->fTprDepTrigger[6]->GetBinContent(i+1);
    content += triggerCounter6[i];
    Hist->fTprDepTrigger[6]->SetBinContent(i+1, content);
    
    content = Hist->fTprDepTrigger[8]->GetBinContent(i+1);
    content += triggerCounter8[i];
    Hist->fTprDepTrigger[8]->SetBinContent(i+1, content);

    content = Hist->fTprDepTrigger[9]->GetBinContent(i+1);
    content += triggerCounter9[i];
    Hist->fTprDepTrigger[9]->SetBinContent(i+1, content);

    content = Hist->fTprDepTrigger[10]->GetBinContent(i+1);
    content += triggerCounter10[i];
    Hist->fTprDepTrigger[10]->SetBinContent(i+1, content);

    content = Hist->fTprDepTrigger[11]->GetBinContent(i+1);
    content += triggerCounter11[i];
    Hist->fTprDepTrigger[11]->SetBinContent(i+1, content);

    content = Hist->fTprDepTrigger[12]->GetBinContent(i+1);
    content += triggerCounter12[i];
    Hist->fTprDepTrigger[12]->SetBinContent(i+1, content);


    content = Hist->fTprDepTrigger[13]->GetBinContent(i+1);
    content += triggerCounter13[i];
    Hist->fTprDepTrigger[13]->SetBinContent(i+1, content);


    content = Hist->fTprDepTrigger[14]->GetBinContent(i+1);
    content += triggerCounter14[i];
    Hist->fTprDepTrigger[14]->SetBinContent(i+1, content);


    content = Hist->fTprDepTrigger[15]->GetBinContent(i+1);
    content += triggerCounter15[i];
    Hist->fTprDepTrigger[15]->SetBinContent(i+1, content);

  }

  
  for (int i=0; i<nCuts2; ++i){
    double   content = Hist->fTprDepTrigger[7]->GetBinContent(i+1);
    content += triggerCounter7[i];
    Hist->fTprDepTrigger[7]->SetBinContent(i+1, content);
    
  }

}



// funtion that evaluates the background rejection and trigger efficiency
// from the Calo-seeded algorithm
void TTrigAna001Module::FillCprDepTriggerHist(EventHist_t*   Hist, int &hasQualityTrack){
  TStnHelix*      helix(0);
  TStnTrackSeed*    trkSeed(0);

  int  nCuts(11);
  int  nhits_cut       [nCuts] = {0};
  int  triggerCounter  [nCuts] = {0};
  int  triggerCounter1 [nCuts] = {0};
  int  triggerCounter2 [nCuts] = {0};
  int  triggerCounter3 [nCuts] = {0};
  int  triggerCounter4 [nCuts] = {0};
  int  triggerCounter5 [nCuts] = {0};
  int  triggerCounter6 [nCuts] = {0};
  int  triggerCounter8 [nCuts] = {0};
  int  triggerCounter9 [nCuts] = {0};
  int  triggerCounter10[nCuts] = {0};
  int  triggerCounter11[nCuts] = {0};
  int  triggerCounter12[nCuts] = {0};

  for (int i=0; i<nCuts; ++i){ nhits_cut[i] = 10 + i; }

  int  nCuts2(8);
  int  triggerCounter7 [nCuts2] = {0};
  int  triggerCounter13[nCuts] = {0};
  int  triggerCounter14[nCuts] = {0};
  int  triggerCounter15[nCuts] = {0};
		       
  double clE_cut[nCuts2]       = {0};
  
  bool   hasTrigHelix(false);

  for (int i=0; i<nCuts2; ++i){ clE_cut[i] = 40 + i*5.; }



  //  int  chi2Condition(0);

  //effciecy of the pattern-recognition
  for (int i=0; i<fNHelices[2]; ++i){
    
    helix   = fCprDepHelixBlock->Helix(i);
    int         nhits    = helix->NHits();
    double      chi2XY   = helix->Chi2XY();
    double      chi2ZPhi = helix->Chi2ZPhi();
    
    if (nhits >= fMinTrigHelixNSh)   hasTrigHelix = true;
   
    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	triggerCounter13[j] = 1;
	
	if (hasQualityTrack == 1){
	  triggerCounter14[j] = 1;
	}
	if ( (chi2XY < 4.) && ( chi2ZPhi < 4)){
	  //	  chi2Condition      = 1;
	  
	  if (hasQualityTrack == 1){
	    triggerCounter15[j] = 1;
	  }
	}
      }
    }//end loop
  }

  for (int i=0; i<fNTrackSeeds[2]; ++i){
    
    trkSeed   = fCprDepTrackSeedBlock->TrackSeed(i);
    int         nhits    = trkSeed->NHits();

    for (int j=0; j<nCuts; ++j){
      if (nhits >= nhits_cut[j])  {
	if (hasQualityTrack == 1){
	  triggerCounter12[j] = 1;
	}
      }
    }
  }
  //--------------------------------------------------------------------------------
  if (hasTrigHelix) { 
    
    //calPatREc
    for (int i=0; i<fNTrackSeeds[2]; ++i){
    
      trkSeed   = fCprDepTrackSeedBlock->TrackSeed(i);
      int         nhits    = trkSeed->NHits();
      double      chi2     = trkSeed->Chi2()/(nhits-5.);
    
      int         helixIndex  = trkSeed->HelixIndex();

      if (helixIndex>=0 && helixIndex<fNHelices[1] ){
	helix                   = fCprDepHelixBlock->Helix(helixIndex);
	int         helix_nhits = helix->NHits();
	int         delta_nhits = nhits - helix_nhits;
	Hist->fDNHits->Fill(delta_nhits);
      }

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter[j] = 1;
	  // if (hasQualityTrack == 1){
	  //   triggerCounter12[j] = 1;
	  // }
	
	  if ( (chi2 < 5.) /*&& ( chi2ZPhi < 4)*/){
	    triggerCounter1[j] = 1;
	    //	  chi2Condition      = 1;
	  
	    // if (hasQualityTrack == 1){
	    //   triggerCounter8[j] = 1;
	    // }
	  }
	}

      }



      double clE =  trkSeed->ClusterEnergy();
      for (int j=0; j<nCuts2; ++j){
	if (clE>= clE_cut[j]){
	  triggerCounter7[j] = 1;
	}
      }
    }


    //CalPatRec
    int ntrk = fNTrackSeeds[2];
    for (int i=0; i<ntrk; ++i){
      trkSeed = fCprDepTrackSeedBlock->TrackSeed(i);
    
      int         nhits = trkSeed->NHits();
      double      d0    = trkSeed->fD0;
      double      chi2  = trkSeed->fChi2/(nhits-5.);
      double      p     = trkSeed->P();
      double      p_min(80.);

      for (int j=0; j<nCuts; ++j){
	if (nhits >= nhits_cut[j])  {
	  triggerCounter2[j] = 1;
	  
	  // if (hasQualityTrack == 1){
	  //   triggerCounter11[j] = 1;
	  // }
	  
	  if ( (d0>-200) && (d0<200)){
	    triggerCounter3[j] = 1;
	  }

	  if ( (d0>-200) && (d0<200) &&
	       (hasQualityTrack == 1)){
	    triggerCounter9[j] = 1;
	  }	

	  if (p > p_min){
	    triggerCounter4[j] = 1;
	  
	    if (hasQualityTrack == 1){
	      triggerCounter8[j] = 1;
	    }
	    if ( (d0>=-200) && (d0<=200)){
	      triggerCounter5[j] = 1;

	      if (hasQualityTrack == 1){
		triggerCounter11[j] = 1;
	      }

	      if (chi2 < 3.){
		triggerCounter6[j] = 1;
		if (hasQualityTrack == 1){
		  triggerCounter10[j] = 1;
		}
	      }
	    }
	  }
	}
      }
    }
  
  }

  for (int i=0; i<nCuts; ++i){
    double   content = Hist->fCprDepTrigger[0]->GetBinContent(i+1);
    content += triggerCounter[i];
    Hist->fCprDepTrigger[0]->SetBinContent(i+1, content);
    
    content = Hist->fCprDepTrigger[1]->GetBinContent(i+1);
    content += triggerCounter1[i];
    Hist->fCprDepTrigger[1]->SetBinContent(i+1, content);
    
    content = Hist->fCprDepTrigger[2]->GetBinContent(i+1);
    content += triggerCounter2[i];
    Hist->fCprDepTrigger[2]->SetBinContent(i+1, content);
    
    content = Hist->fCprDepTrigger[3]->GetBinContent(i+1);
    content += triggerCounter3[i];
    Hist->fCprDepTrigger[3]->SetBinContent(i+1, content);

    content = Hist->fCprDepTrigger[4]->GetBinContent(i+1);
    content += triggerCounter4[i];
    Hist->fCprDepTrigger[4]->SetBinContent(i+1, content);
    
    content = Hist->fCprDepTrigger[5]->GetBinContent(i+1);
    content += triggerCounter5[i];
    Hist->fCprDepTrigger[5]->SetBinContent(i+1, content);
    
    content = Hist->fCprDepTrigger[6]->GetBinContent(i+1);
    content += triggerCounter6[i];
    Hist->fCprDepTrigger[6]->SetBinContent(i+1, content);
    
    content = Hist->fCprDepTrigger[8]->GetBinContent(i+1);
    content += triggerCounter8[i];
    Hist->fCprDepTrigger[8]->SetBinContent(i+1, content);

    content = Hist->fCprDepTrigger[9]->GetBinContent(i+1);
    content += triggerCounter9[i];
    Hist->fCprDepTrigger[9]->SetBinContent(i+1, content);

    content = Hist->fCprDepTrigger[10]->GetBinContent(i+1);
    content += triggerCounter10[i];
    Hist->fCprDepTrigger[10]->SetBinContent(i+1, content);

    content = Hist->fCprDepTrigger[11]->GetBinContent(i+1);
    content += triggerCounter11[i];
    Hist->fCprDepTrigger[11]->SetBinContent(i+1, content);

    content = Hist->fCprDepTrigger[12]->GetBinContent(i+1);
    content += triggerCounter12[i];
    Hist->fCprDepTrigger[12]->SetBinContent(i+1, content);


    content = Hist->fCprDepTrigger[13]->GetBinContent(i+1);
    content += triggerCounter13[i];
    Hist->fCprDepTrigger[13]->SetBinContent(i+1, content);


    content = Hist->fCprDepTrigger[14]->GetBinContent(i+1);
    content += triggerCounter14[i];
    Hist->fCprDepTrigger[14]->SetBinContent(i+1, content);


    content = Hist->fCprDepTrigger[15]->GetBinContent(i+1);
    content += triggerCounter15[i];
    Hist->fCprDepTrigger[15]->SetBinContent(i+1, content);

  }

  
  for (int i=0; i<nCuts2; ++i){
    double   content = Hist->fCprDepTrigger[7]->GetBinContent(i+1);
    content += triggerCounter7[i];
    Hist->fCprDepTrigger[7]->SetBinContent(i+1, content);
    
  }
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TTrigAna001Module::Event(int ientry) {

  double                xs, p;
  TEmuLogLH::PidData_t  dat;
  TStnTrack*            track;
  TVDetHitData*         vhit;
  int                   id_word;
  TLorentzVector        mom;

  TDiskCalorimeter::GeomData_t disk_geom;

  fTprDemTrackBlock    ->GetEntry(ientry);
  fCprDemTrackBlock    ->GetEntry(ientry);
  if (fIsBackground ==0) fTrackBlock    ->GetEntry(ientry);
  // fTimeClusterBlock ->GetEntry(ientry);
  fTprDemTrackSeedBlock->GetEntry(ientry);
  fCprDemTrackSeedBlock->GetEntry(ientry);

  if (fIsBackground){
    fTprDepTrackBlock->GetEntry(ientry);
    fCprDepTrackBlock->GetEntry(ientry);
    
    fTprDepTrackSeedBlock->GetEntry(ientry);
    fCprDepTrackSeedBlock->GetEntry(ientry);
    
    fTprDepHelixBlock ->GetEntry(ientry);
    fCprDepHelixBlock ->GetEntry(ientry);
  }
  fTprDemHelixBlock ->GetEntry(ientry);

  fCprDemHelixBlock ->GetEntry(ientry);
  //
  fClusterBlock  ->GetEntry(ientry);
  //  fStrawDataBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVdetDataBlock->GetEntry(ientry);
  //-----------------------------------------------------------------------------
  // assume electron in the first particle, otherwise the logic will need to 
  // be changed
  //-----------------------------------------------------------------------------
  fNGenp    = fGenpBlock->NParticles();
  
  //  GetHeaderBlock()->Print("DEBUG");

  TGenParticle* genp;
  int           /*pdg_code,*/ generator_code;

  fParticle = NULL;
  for (int i=fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    //  pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if (/*(abs(pdg_code) == fPdgCode) &&*/ (generator_code == fGeneratorCode)) {
      fParticle = genp;
      break;
    }
  }
	
  if (fSimp != NULL){ // may want to revisit the definition of fSimp
    fSimp     = fSimpBlock->Particle(0);
  }
  
  if (fParticle != NULL){
    fParticle->Momentum(mom);
  }				// this is a kludge, to be removed at the next 
  // ntupling 
  //  fEleE     = fParticle->Energy();
  p         = mom.P();
  fEMom     = p;
  fEleE     = sqrt(p*p+0.511*0.511);


  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
      //      disk_geom.fNCrystals[i] = fCalDataBlock->fNCrystals[i];
      disk_geom.fRMin[i]      = fCalDataBlock->fRMin[i];
      disk_geom.fRMax[i]      = fCalDataBlock->fRMax[i];
      disk_geom.fZ0  [i]      = fCalDataBlock->fZ0  [i];
    }

    disk_geom.fHexSize          = fCalDataBlock->CrystalSize()*2;
    // kludge , so far
    disk_geom.fMinFraction      = 1.; // fCalDataBlock->MinFraction();
    disk_geom.fWrapperThickness = fCalDataBlock->WrapperThickness();
    disk_geom.fShellThickness   = fCalDataBlock->ShellThickness();

    fDiskCalorimeter->Init(&disk_geom);
  }

  if (fIsBackground) {
    fNTracks[0]     = 0;
    fNTracks[3]     = fCprDepTrackBlock    ->NTracks();
    fNTracks[4]     = fTprDepTrackBlock    ->NTracks();
  } else {
    fNTracks[0]     = fTrackBlock    ->NTracks();
    fNTracks[3]     = 0;
    fNTracks[4]     = 0;
  }
  fNTracks[1]     = fCprDemTrackBlock    ->NTracks();
  fNTracks[2]     = fTprDemTrackBlock    ->NTracks();


  if ( (fNTracks[2]>0) && (fNTracks[1] == 0) ){
    TStnTrack* trk        = fTprDemTrackBlock->Track(0);
    double     ecl        = trk->ClusterE();
    if (ecl>50.){
      double     fcons      = trk->fFitCons;
      double     tanDip     = trk->fTanDip;
      int        nactive    = trk->NActive();
      double     p          = trk->fP;   
      double     chi2dof    = trk->Chi2Dof();
      
      if (GetDebugBit(46)) {
	GetHeaderBlock()->Print(Form(" bit:46 ecl = %5.3f nactive = %i p = %5.3f tanDip = %5.3f chi2/ndof = %5.3f fitcons = %5.3f", 
				     ecl, nactive , p, tanDip, chi2dof, fcons ));
      }
    }
  }


  fNTimeClusters [0] = 0;//fTimeClusterBlock->NTimeClusters();

  fNTrackSeeds[0] = fCprDemTrackSeedBlock->NTrackSeeds();
  fNTrackSeeds[1] = fTprDemTrackSeedBlock->NTrackSeeds();

  fNTrackSeeds[2] = 0;

  //fNTrackSeeds[3] = fTprDepTrackSeedBlock->NTrackSeeds();
  fNTrackSeeds[3] = 0;

  fNHelices[0]    = fTprDemHelixBlock   ->NHelices();
  fNHelices[1]    = fCprDemHelixBlock->NHelices();
  fNHelices[2]    = 0;//fCprDepHelixBlock->NHelices();
  fNHelices[3]    = 0;//fTprDepHelixBlock->NHelices();

  if(fIsBackground){
    fNTrackSeeds[2] = fCprDepTrackSeedBlock->NTrackSeeds();
    fNTrackSeeds[3] = fTprDepTrackSeedBlock->NTrackSeeds();

    fNHelices[2]    = fCprDepHelixBlock->NHelices();
    fNHelices[3]    = fTprDepHelixBlock->NHelices();
  }
  

  if ( (fNHelices[1] == 1) && (fNHelices[0] == 0) && GetDebugBit(50)) {
    TStnHelix*  helix    = fCprDemHelixBlock->Helix(0);
    int         nhits    = helix->NHits();
    double      radius   = helix->Radius();
    double      lambda   = helix->Lambda();
    double      tanDip   = lambda/radius;
    double      mm2MeV   = 3/10.;
    double      p        = radius*mm2MeV/std::cos( std::atan(tanDip));
    double      chi2xy   = helix->Chi2XY();
    double      chi2zphi = helix->Chi2ZPhi();
    GetHeaderBlock()->Print(Form(" bit:050 CPR p = %5.3f  nhits = %3i chi2XY = %5.3f  chi2PhiZ = %5.3f", 
				 p, nhits, chi2xy, chi2zphi));
  }
  //reset track t0 reference
  fTrackT0        = -99999;
  fTrackCaloE     = -99999; 

  fTrackTrkSeedIndex    = -1;
  fTrackHelixIndex      = -1;
  fTprTrackHelixIndex   = -1;
  fTrackDepHelixIndex      = -1;
  fTprDepTrackHelixIndex   = -1;

  if (fNTrackSeeds[0] > 0){
    TStnTrackSeed*    TrkSeed;
    for (int i=0; i<fNTrackSeeds[0]; ++i){
      TrkSeed = fCprDemTrackSeedBlock->TrackSeed(i);

      int         nhits    = TrkSeed->NHits      ();
      double      pT       = TrkSeed->Pt();
      
      double      tanDip   = TrkSeed->TanDip();  
      double      p        = TrkSeed->P();//radius*mm2MeV/std::cos( std::atan(tanDip));
      double      chi2     = TrkSeed->Chi2()/(nhits-5.);
      double      fitcons  = TrkSeed->FitCons();
      if (GetDebugBit(38)) {
	GetHeaderBlock()->Print(Form(" bit:38 nhits = %i p = %5.3f pt = %5.3f tanDip = %5.3f chi2 = %5.3f fitcons = %5.3f", 
				     nhits, p, pT, tanDip, chi2, fitcons ));
      }

      if (GetDebugBit(39)) {
	if (/*chi2zphi >5 ||*/ chi2 >5 ){
	  GetHeaderBlock()->Print(Form(" bit:39 nhits = %i p = %5.3f pt = %5.3f tanDip = %5.3f chi2 = %5.3f fitcons = %5.3f", 
				       nhits, p, pT, tanDip, chi2, fitcons));
	}
      }
    }
  }


  fNClusters    = fClusterBlock->NClusters();
  fNCalHits     = fCalDataBlock->NHits();
  fNStrawHits   = fStrawDataBlock->NHits();
  fNVirtualHits = fVdetDataBlock->NHits();
    
  fDiskCalorimeter->InitEvent(fCalDataBlock);

  fNHyp       = -1;
  fBestHyp[0] = -1;
  fBestHyp[1] = -1;

  fNGoodTracks    = 0;
  fNMatchedTracks = 0;

  if (fIsBackground) {
    fNTracks[0] = 0;
    
  }else {
    fNTracks[0] = fTrackBlock    ->NTracks();
  }

  fNTracks[1] = fCprDemTrackBlock->NTracks();
  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fCprDemTrackBlock->Track(0);

  int ntrk = fNTracks[0];

  TrackPar_t*   tp;

  for (int itrk=0; itrk<ntrk; itrk++) {
    // assume less 20 tracks
    tp             = fTrackPar+itrk;

    track          = fTrackBlock->Track(itrk);
    id_word        = fTrackID->IDWord(track);
    track->fIDWord = id_word;
    if (id_word == 0) {
      fNGoodTracks += 1;
      if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	fNMatchedTracks += 1;
      }
    }
    //-----------------------------------------------------------------------------
    // process hit masks
    //-----------------------------------------------------------------------------
    int i1, i2, n1(0) ,n2(0), ndiff(0);
    int nbits = track->fHitMask.GetNBits();
    for (int i=0; i<nbits; i++) {
      i1 = track->HitMask()->GetBit(i);
      i2 = track->ExpectedHitMask()->GetBit(i);
      n1 += i1;
      n2 += i2;
      if (i1 != i2) ndiff += 1;
    }
    //-----------------------------------------------------------------------------
    // define additional parameters
    //-----------------------------------------------------------------------------
    tp->fNHPl = n1;
    tp->fNEPl = n2;
    tp->fNDPl = ndiff;

    tp->fDpF   = track->fP     -track->fPFront;
    tp->fDp0   = track->fP0    -track->fPFront;
    tp->fDp2   = track->fP2    -track->fPFront;
    tp->fDpFSt = track->fPFront-track->fPStOut;

    tp->fT0    = track->fT0;
    tp->fT0err = track->fT0Err;
    tp->fXTrk  = track->fX1;
    tp->fYTrk  = track->fY1;
    tp->fZTrk  = track->fZ1;


    //-----------------------------------------------------------------------------
    // track residuals
    //-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMaxEp; 
    double    nx, ny;

    tp->fEcl       = -1.e6;
    tp->fEp        = -1.e6;

    tp->fDu        = -1.e6;
    tp->fDv        = -1.e6;
    tp->fDx        = -1.e6;
    tp->fDy        = -1.e6;
    tp->fDz        = -1.e6;
    tp->fDt        = -1.e6;

    tp->fChi2Match = -1.e6;
    tp->fPath      = -1.e6;

    if (vr) {
      tp->fEcl = vr->fEnergy;
      tp->fEp  = tp->fEcl/track->fP;

      tp->fDx  = vr->fDx;
      tp->fDy  = vr->fDy;
      tp->fDz  = vr->fDz;
      //-----------------------------------------------------------------------------
      // v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
      //-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt;

      nx  = vr->fNxTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);
      ny  = vr->fNyTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);

      tp->fDu        = vr->fDx*nx+vr->fDy*ny;
      tp->fDv        = vr->fDx*ny-vr->fDy*nx;
      tp->fChi2Match = vr->fChi2Match;
      tp->fPath      = vr->fPath;
    }

    // if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
    //   GetHeaderBlock()->Print(Form(" TTrigAna001Module ERROR: tp->fEp = %10.5f  track->fEp = %10.5f\n ",tp->fEp,track->fEp));
    // }
    //-----------------------------------------------------------------------------
    // PID likelihoods
    //-----------------------------------------------------------------------------
    dat.fDt   = tp->fDt;
    dat.fEp   = tp->fEp;
    dat.fPath = tp->fPath;
      
    xs = track->XSlope();

    track->fEleLogLHCal = fLogLH->LogLHCal(&dat,11);
    track->fMuoLogLHCal = fLogLH->LogLHCal(&dat,13);

    double llhr_cal = track->fEleLogLHCal-track->fMuoLogLHCal;

    if (GetDebugBit(7)) {
      if ((id_word == 0) && (llhr_cal > 20)) {
	GetHeaderBlock()->Print(Form("bit:007: dt = %10.3f ep = %10.3f",track->Dt(),tp->fEp));
      }
    }

    if (GetDebugBit(8)) {
      if ((id_word == 0) && (llhr_cal < -20)) {
	GetHeaderBlock()->Print(Form("bit:008: p = %10.3f dt = %10.3f ep = %10.3f",
				     track->P(),track->Dt(),tp->fEp));
      }
    }

    track->fLogLHRXs    = fLogLH->LogLHRXs(xs);
  }

  //////////////////////////////////////////////////////////////////////
  // Now fill the virtual detectors information 
  //////////////////////////////////////////////////////////////////////
  VirtualPar_t*   vp;

  for (int ivhit=0; ivhit<fNVirtualHits; ivhit++) {
    
           
    vhit           = fVdetDataBlock->Hit(ivhit);
    if (vhit == NULL){
      printf("[filling Virtual det hits] step %i of %i returns NULL pointer. Max limit = %i\n",
	     ivhit, fNVirtualHits, fVdetDataBlock->NHits());
      break;
    }
    vp             = fVirtualPar+ivhit;

    //-----------------------------------------------------------------------------
    // define additional parameters
    //-----------------------------------------------------------------------------
    vp->fIndex     = vhit->Index();
    vp->fPdgCode   = vhit->PdgCode();
    vp->fGenCode   = vhit->GeneratorCode();
    vp->fTime      = vhit->Time();
    vp->fMass      = vhit->Mass();
    vp->fEnergyKin = vhit->EnergyKin();
    vp->fEnergy    = vhit->Energy();
    vp->fMom       = vhit->McMomentum();			// MC particle momentum
    vp->fMomX      = vhit->McMomentumX();			// MC particle momentum X - component
    vp->fMomY      = vhit->McMomentumY();			// MC particle momentum Y - component
    vp->fMomZ      = vhit->McMomentumZ();			// MC particle momentum Z - component
    vp->fPosX      = vhit->McPositionX();			// MC particle position X - component
    vp->fPosY      = vhit->McPositionY();			// MC particle position Y - component
    vp->fPosZ      = vhit->McPositionZ();
  }
  ////////////////////////////////////////////////////////////////////////////////



  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

  //  fDiskCalorimeter->InitEvent(fCalDataBlock);

  FillHistograms();

  //  ++fEventCounter;

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrigAna001Module::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp;
  int ntrk = fCprDemTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fCprDemTrackBlock->Track(itrk);
    tp  = &fTrackPar[itrk];
    //-----------------------------------------------------------------------------
    // bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
    //-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }
    //-----------------------------------------------------------------------------
    // bit 4: tracks with DpF > 1MeV - positive tail...
    //-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if (tp->fDpF > 1.) {
	GetHeaderBlock()->Print(Form("pF pRec, fDpf = %10.3f  %10.3f  %10.3f",
				     trk->fPFront, trk->Momentum()->P(),tp->fDpF));
      }
    }
    //-----------------------------------------------------------------------------
    // bit 9: Set C tracks with DpF > 1MeV - positive tail...
    //-----------------------------------------------------------------------------
    if (GetDebugBit(9) == 1) {
      double ep = trk->Ep();
      if (trk->fIDWord == 0) { 
	if (((ep > 0.42) && (ep < 0.46)) || ((ep > 0.35) && (ep < 0.39))) {
	  GetHeaderBlock()->Print(Form("bit:009 ep = %10.3f e = %10.3f p = %10.3f",
				       trk->fEp,trk->fEp*trk->fP,trk->fP));
	}
      }
    }
    //-----------------------------------------------------------------------------
    // bit 10: Set C tracks with Ecl > 80
    //-----------------------------------------------------------------------------
    if (GetDebugBit(10) == 1) {
      double ecl = trk->ClusterE();
      if (trk->fIDWord == 0) { 
	if (ecl > 60) {
	  GetHeaderBlock()->Print(Form("bit:010 e = %10.3f p = %10.3f",
				       ecl,trk->fP));
	}
      }
    }
  }

  //-----------------------------------------------------------------------------
  // bit 5: events with N(tracks) > 1
  //-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    int ntrk = fCprDemTrackBlock->NTracks();
    if (ntrk > 1) {
      GetHeaderBlock()->Print(Form("NTracks = %i5",ntrk));
    }
  }
}

//_____________________________________________________________________________
int TTrigAna001Module::EndJob() {
  if (fEventCounter>0)  Normalize();
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrigAna001Module::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

