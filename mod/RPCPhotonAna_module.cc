//
// An EDAnalyzer module that is used to study the RPC photon
//
// $Id:  $
// $Author:  $
// $Date: $
//
// Original author Heather Harrington
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"


#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"


#include "MCDataProducts/inc/EventWeight.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "CaloMC/inc/ClusterContentMC.hh"
#include "CaloMC/inc/CrystalContentMC.hh"

#include "CaloCluster/inc/ClusterMoments.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

// #include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>




using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;



namespace mu2e {


  class RPCPhotonAna : public art::EDFilter {
    
    struct  particleHist_ {
      TH1F* energy[3];
      TH1F* time  [3];
      TH1F* pdgId [3];
      TH1F* cosTh [3];
    };

    struct  clusterHist_ {
      TH1F* energy[3];   // reco energy in the cluster 
      TH1F* time  [3];     // reco cluster time
      TH1F* dt    [3];
      TH1F* nCr   [3];      // cluster size in crystal units
      TH1F* rDist [3];    // cluster radial distance from the DS axis

      TH1F* e1eRatio [3];  // most  energetic crystal energy over total cluster energy
      TH1F* e2eRatio [3];  // second most  energetic crystal energy over total cluster energy
      TH1F* e2e1Ratio[3];  // secondd most  energetic crystal energy over most energetic crystal energy
      TH1F* e3eRatio [3];  // third most  energetic crystal energy over total cluster energy
      TH1F* e4eRatio [3];  // fourth most  energetic crystal energy over total cluster energy

       //correlation w.r.t. cluster radial distance
      TH2F* eRDist    [3];  // cluster energy
      TH2F* e1eRDist  [3];  // most  energetic crystal energy over total cluster energy
      TH2F* e2eRDist  [3];  // second most  energetic crystal energy over total cluster energy
      TH2F* e2e1RDist [3];  // secondd most  energetic crystal energy over most energetic crystal energy
      TH2F* e3eRDist  [3];  // third most  energetic crystal energy over total cluster energy
      TH2F* e4eRDist  [3];  // fourth most  energetic crystal energy over total cluster energy
      TH2F* nCrRDist  [3];  // cluster size vs cluster cog radial distance

    };

    struct  eventHist_ {
      TH1F* nCl[3];
    };

    struct  lhHist_ {
      TH1F* hEnergy    ;  
      TH1F* hE1E       ;
      TH1F* hNCr       ;
      TH1F* hRDist     ;
      TH1F* hE2E       ;
      TH1F* hE2E1      ;
      TH1F* hE3E       ;
      TH1F* hE4E       ;
      
      TH1F* hCorEnergy ;
      TH1F* hCorE1E    ;
      TH1F* hCorE2E    ;
      TH1F* hCorE2E1   ;
      TH1F* hCorE3E    ;
      TH1F* hCorE4E    ;
      TH1F* hCorNCr    ;

      TH1F* hTot       ;
      TH1F* hCorTot    ;
    };

     struct triggeredHist_ {
      TH1F* energy;    // reco energy in the cluster 
      TH1F* time  ;    // reco cluster time
      TH1F* nCr   ;    // cluster size in crystal units
      TH1F* rDist ;    // cluster radial distance from the DS axis

      TH1F* e1eRatio ;  // most  energetic crystal energy over total cluster energy
      TH1F* e2eRatio ;  // second most  energetic crystal energy over total cluster energy
      TH1F* e2e1Ratio;  // secondd most  energetic crystal energy over most energetic crystal energy
      TH1F* e3eRatio ;  // third most  energetic crystal energy over total cluster energy
      TH1F* e4eRatio ;  // fourth most  energetic crystal energy over total cluster energy

      TH2F* eRDist   ;  // cluster energy
      TH2F* e1eRDist ;  // most  energetic crystal energy over total cluster energy
      TH2F* e2eRDist ;  // second most  energetic crystal energy over total cluster energy
      TH2F* e2e1RDist;  // secondd most  energetic crystal energy over most energetic crystal energy
      TH2F* e3eRDist ;  // third most  energetic crystal energy over total cluster energy
      TH2F* e4eRDist ;  // fourth most  energetic crystal energy over total cluster energy
      TH2F* nCrRDist ;  // cluster size vs cluster cog radial distance
     };
      
      
  public:
    enum {
      kN1DVar    = 10,
      kN2DVar    = 10,
      kNCorHist  = 10
    };

    typedef art::Ptr<StepPointMC> StepPtr;
    typedef std::vector<StepPtr>  StepPtrs;
    typedef std::map<int,StepPtrs > HitMap;

    virtual ~RPCPhotonAna() { }

    virtual void beginJob();
    virtual bool beginRun(art::Run&   run   );
    virtual void endJob();



    // This is called for each event.
    //     virtual void analyze(const art::Event& e);

    virtual bool filter(art::Event& event) override;
    explicit RPCPhotonAna(const fhicl::ParameterSet& PSet);

    void     bookHistograms();

    void     fillGenHistograms    (const GenParticle*           Gen    );
    void     fillEventHistograms  (const CaloClusterCollection* ClCol  );
    void     fillClusterHistograms(const CaloCluster*           Cluster);
    void     fillVDetHistograms   (const StepPointMC*           Step   );
    void     fillLHHistograms     (const CaloClusterCollection* ClCol  );

  private:
       
    int                        _diagLevel;
    int                        _nProcess;
    
    double                     _evtWeight;
    double                     _vdHitsSize;
    double                     _genSize;
    double                     _vdetE;
    double                     _vdetX, _vdetY;
    double                     _vdetT;
    double                     _vdetE_e;
    double                     _vdetX_e, _vdetY_e;
    double                     _vdetT_e;

    bool                       _evtWithPhoton;
    bool                       _evtWithElectron;

    ClusterMoments::cogtype    _cogType;
    bool                       _dropSecondDisk;
    std::string                _generatorModuleLabel;
    std::string                _weightModuleLabel;

    std::string                _caloClusterModuleLabel;

    std::string                _virtualDetectorLabel;

    SimParticleTimeOffset      _toff;     // time offset smearing
    double                     _mbtime;
    double                     _mbbuffer;
    double                     _blindTime;
    
    std::string                _signalTemplateFile;
    std::string                _bkgTemplateFile;

    double                     _minRDist   , _rDistStep;
   
    //histograms collectors
    particleHist_   _hGen;
    particleHist_   _hVDet;
    clusterHist_    _hCl;
    eventHist_      _hEvt;
    lhHist_         _hLh;
    triggeredHist_  _hTrig;

    GlobalConstantsHandle<ParticleDataTable> _pdt;
    std::vector<TString>                     _clNames;

    //template histograms
    TH1F* _bkgEnergy;   
    TH1F* _bkgTime   ;  
    TH1F* _bkgnCr     ; 
    TH1F* _bkgrDist    ;
    TH1F* _bkge1eRatio ;
    TH1F* _bkge2eRatio ;
    TH1F* _bkge2e1Ratio;
    TH1F* _bkge3eRatio ;
    TH1F* _bkge4eRatio ;

    TH2F* _bkgeRDist    ; 
    TH2F* _bkge1eRDist   ;
    TH2F* _bkge2eRDist   ;
    TH2F* _bkge2e1RDist   ;
    TH2F* _bkge3eRDist   ;
    TH2F* _bkge4eRDist   ;
    TH2F* _bkgnCrRDist  ;
                           
    TH1F* _signalEnergy ;  
    TH1F* _signalTime   ;  
    TH1F* _signalnCr    ;  
    TH1F* _signalrDist  ;  
    TH1F* _signale1eRatio; 
    TH1F* _signale2eRatio ;
    TH1F* _signale2e1Ratio;
    TH1F* _signale3eRatio ;
    TH1F* _signale4eRatio ;
                           
    TH2F* _signaleRDist   ;
    TH2F* _signale1eRDist ;
    TH2F* _signale2eRDist ;
    TH2F* _signale2e1RDist ;
    TH2F* _signale3eRDist ;
    TH2F* _signale4eRDist ;
    TH2F* _signalnCrRDist ;

    TH1F* _signalCorrHistEnergy[10];
    TH1F* _bkgCorrHistEnergy[10];
    TH1F* _signalCorrHiste1e[10];
    TH1F* _bkgCorrHiste1e[10];
    TH1F* _signalCorrHiste2e[10];
    TH1F* _bkgCorrHiste2e[10];
    TH1F* _signalCorrHiste2e1[10];
    TH1F* _bkgCorrHiste2e1[10];
    TH1F* _signalCorrHiste3e[10];
    TH1F* _bkgCorrHiste3e[10];
    TH1F* _signalCorrHiste4e[10];
    TH1F* _bkgCorrHiste4e[10];
    TH1F* _signalCorrHistnCr[10];
    TH1F* _bkgCorrHistnCr[10];

    void        buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize);
    
    //initialize calclate prob in private
    double      calculateProb  (double &variable, TH1F* templates);
    double      calculate2DProb(double& Ref, double&Variable, TH1F** Template, double MinX, double Step);

  };


  RPCPhotonAna::RPCPhotonAna(fhicl::ParameterSet const & pset) :
    art::EDFilter{pset},
    _diagLevel                   (pset.get<int>("diagLevel",0)),
    _nProcess                    (0),
    _cogType                     (ClusterMoments::Linear),                
    _generatorModuleLabel        (pset.get<string>("generatorModuleLabel")),
    _weightModuleLabel           (pset.get<string>("weightModuleLabel")),
    _caloClusterModuleLabel      (pset.get<string>("caloClusterModuleLabel")),
    _virtualDetectorLabel        (pset.get<string>("virtualDetectorName")),
    _toff                        (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _mbbuffer                    (pset.get<double>("TimeFoldingBuffer")),  // ns
    _blindTime                   (pset.get<double>("blindTime" )) ,         // ns
    _signalTemplateFile          (pset.get<string>("signalTemplates")),
    _bkgTemplateFile             (pset.get<string>("backgroundTemplates")),
    _minRDist                    (pset.get<double>("minClusterRadialDist" ,  350.)),   // mm
    _rDistStep                   (pset.get<double>("clusterRadialDistStep",   50.)){

      ConfigFileLookupPolicy configFile;
      _signalTemplateFile = configFile(_signalTemplateFile);
      _bkgTemplateFile    = configFile(_bkgTemplateFile);

      TFile* signal = TFile::Open(_signalTemplateFile.c_str());
      TFile* bkg    = TFile::Open(_bkgTemplateFile.c_str());

      //get the templates histograms
      _bkgEnergy    = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_E0");
      _bkgTime      = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_time0");
      _bkgnCr       = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_nCr0");
      _bkgrDist     = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_rDist0");
      _bkge1eRatio  = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e1eRatio0");
      _bkge2eRatio  = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e2eRatio0");
      _bkge2e1Ratio = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e2e1Ratio0");
      _bkge3eRatio  = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e3eRatio0");
      _bkge4eRatio  = (TH1F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e4eRatio0");
    
      _bkgeRDist    =  (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/cl_eRDist0");
      _bkge1eRDist  =  (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e1eRDist0");
      _bkge2eRDist  =  (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e2eRDist0");
      _bkge2e1RDist =  (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e2e1RDist0");
      _bkge3eRDist  =  (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e3eRDist0");
      _bkge4eRDist  =  (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/cl_e4eRDist0");
      _bkgnCrRDist  =  (TH2F*)bkg->Get("RPCPhotonAna/cluster_all/cl_nCrRDist0");
                
      _signalEnergy     = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_E1");
      _signalTime       = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_time1");
      _signalnCr        = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_nCr1");
      _signalrDist      = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_rDist1");
      _signale1eRatio   = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e1eRatio1");
      _signale2eRatio   = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e2eRatio1");
      _signale2e1Ratio  = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e2e1Ratio1");
      _signale3eRatio   = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e3eRatio1");
      _signale4eRatio   = (TH1F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e4eRatio1");
                
      _signaleRDist     = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/cl_eRDist1");
      _signale1eRDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e1eRDist1");
      _signale2eRDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e2eRDist1");
      _signale2e1RDist  = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e2e1RDist1");
      _signale3eRDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e3eRDist1");
      _signale4eRDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/cl_e4eRDist1");
      _signalnCrRDist   = (TH2F*)signal->Get("RPCPhotonAna/cluster_photon/cl_nCrRDist1");

      //make correlation templates  
      //cluster energy vs cluster radial distance
      buildTemplateHist(_signaleRDist, _signalCorrHistEnergy, Form("SignaleRDist"), _minRDist, _rDistStep);
      buildTemplateHist(_bkgeRDist, _bkgCorrHistEnergy, Form("BkgeRDist"), _minRDist, _rDistStep);
      
      //E1 (seedHitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signale1eRDist, _signalCorrHiste1e, Form("Signale1eRDist"), _minRDist, _rDistStep);
      buildTemplateHist(_bkge1eRDist, _bkgCorrHiste1e, Form("Bkge1eRDist"), _minRDist, _rDistStep);
      
      //E2 (2nd_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signale2eRDist, _signalCorrHiste2e, Form("Signale2eRDist"), _minRDist, _rDistStep);
      buildTemplateHist(_bkge2eRDist, _bkgCorrHiste2e, Form("Bkge2eRDist"), _minRDist, _rDistStep);

      //E2/E1 (2nd_HitEnergy/seed_energy) vs cluster radial distance
      buildTemplateHist(_signale2e1RDist, _signalCorrHiste2e1, Form("Signale2e1RDist"), _minRDist, _rDistStep);
      buildTemplateHist(_bkge2e1RDist, _bkgCorrHiste2e1, Form("Bkge2e1RDist"), _minRDist, _rDistStep);

      //E3 (3rd_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signale3eRDist, _signalCorrHiste3e, Form("Signale3eRDist"), _minRDist, _rDistStep);
      buildTemplateHist(_bkge3eRDist, _bkgCorrHiste3e, Form("Bkge3eRDist"), _minRDist, _rDistStep);

      //E3 (4th_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signale4eRDist, _signalCorrHiste4e, Form("Signale4eRDist"), _minRDist, _rDistStep);
      buildTemplateHist(_bkge4eRDist, _bkgCorrHiste4e, Form("Bkge4eRDist"), _minRDist, _rDistStep);

      //nCrystals vs cluster radial distance
      buildTemplateHist(_signalnCrRDist, _signalCorrHistnCr, Form("SignalnCrRDist"), _minRDist, _rDistStep);
      buildTemplateHist(_bkgnCrRDist, _bkgCorrHistnCr, Form("BkgnCrRDist"), _minRDist, _rDistStep);
    }
  //--------------------------------------------------------------------------------
  // routine function used to produce the template histograms 
  // from the 2D distributions
  //--------------------------------------------------------------------------------
  void   RPCPhotonAna::buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize){
    double    binsize = Input->GetXaxis()->GetBinWidth(1);
    double    binlow  = Input->GetXaxis()->GetBinLowEdge(1);
    
    for (int i=0; i<kNCorHist; ++i){
      double  start      = (MinX + i*StepSize);
      int     binstart   = (start - binlow)/binsize;
      int     binend     = binstart + (StepSize/binsize);
      FinalHist[i]       = (TH1F*)Input->ProjectionY(Form("%s_%i", Label.Data(), i), binstart, binend);
      double  area       = FinalHist[i]->Integral();
      FinalHist[i]->Scale(1/area);
    }
    
  }
  //--------------------------------------------------------------------------------//
  void RPCPhotonAna::bookHistograms(){
    
    art::ServiceHandle<art::TFileService> tfs;

    art::TFileDirectory genDir = tfs->mkdir("gen");
    _hGen.energy[0]  = genDir.make<TH1F>("gen_E0"    ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hGen.time  [0]  = genDir.make<TH1F>("gen_time0" ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hGen.pdgId [0]  = genDir.make<TH1F>("gen_pdgId0","PdgID; pdgId"               , 6000, -3000, 3000);
    _hGen.cosTh [0]  = genDir.make<TH1F>("gen_costh0","cos(#theta); cos(#theta)"   , 2000, -1, 1);
    _hGen.energy[1]  = genDir.make<TH1F>("gen_E1"    ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hGen.time  [1]  = genDir.make<TH1F>("gen_time1" ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hGen.pdgId [1]  = genDir.make<TH1F>("gen_pdgId1","PdgID; pdgId"               , 6000, -3000, 3000);
    _hGen.cosTh [1]  = genDir.make<TH1F>("gen_costh1","cos(#theta); cos(#theta)"   , 2000, -1, 1);
    _hGen.energy[2]  = genDir.make<TH1F>("gen_E2"    ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hGen.time  [2]  = genDir.make<TH1F>("gen_time2" ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hGen.pdgId [2]  = genDir.make<TH1F>("gen_pdgId2","PdgID; pdgId"               , 6000, -3000, 3000);
    _hGen.cosTh [2]  = genDir.make<TH1F>("gen_costh2","cos(#theta); cos(#theta)"   , 2000, -1, 1);

    art::TFileDirectory vdetDir = tfs->mkdir("vdet");
    _hVDet.energy[0] = vdetDir.make<TH1F>("vdet_E0"    ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hVDet.time  [0] = vdetDir.make<TH1F>("vdet_time0" ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hVDet.pdgId [0] = vdetDir.make<TH1F>("vdet_pdgId0","PdgID; pdgId"               , 6000, -3000, 3000);
    _hVDet.cosTh [0] = vdetDir.make<TH1F>("vdet_costh0","cos(#theta); cos(#theta)"   , 2000, -1, 1);
    _hVDet.energy[1] = vdetDir.make<TH1F>("vdet_E1"    ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hVDet.time  [1] = vdetDir.make<TH1F>("vdet_time1" ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hVDet.pdgId [1] = vdetDir.make<TH1F>("vdet_pdgId1","PdgID; pdgId"               , 6000, -3000, 3000);
    _hVDet.cosTh [1] = vdetDir.make<TH1F>("vdet_costh1","cos(#theta); cos(#theta)"   , 2000, -1, 1);
    _hVDet.energy[2] = vdetDir.make<TH1F>("vdet_E2"    ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hVDet.time  [2] = vdetDir.make<TH1F>("vdet_time2" ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hVDet.pdgId [2] = vdetDir.make<TH1F>("vdet_pdgId2","PdgID; pdgId"               , 6000, -3000, 3000);
    _hVDet.cosTh [2] = vdetDir.make<TH1F>("vdet_costh2","cos(#theta); cos(#theta)"   , 2000, -1, 1);

    art::TFileDirectory clDir = tfs->mkdir("cluster_all");
    _hCl.energy   [0]     = clDir.make<TH1F>("cl_E0"       ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hCl.time     [0]     = clDir.make<TH1F>("cl_time0"    ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hCl.nCr      [0]     = clDir.make<TH1F>("cl_nCr0"     ,"Cluster size distribution; nCrystalHits"  , 21, -0.5, 20.5);
    _hCl.rDist    [0]     = clDir.make<TH1F>("cl_rDist0"   ,"Time distribution; cluster radial distance [mm]"  , 500, 300, 800);
    _hCl.e1eRatio [0]     = clDir.make<TH1F>("cl_e1eRatio0 ","most crystalHit energy/E_{cluster}; E^{max}_{hit}/E_{cluster} "  , 100, 0, 1.);
    _hCl.e2eRatio [0]     = clDir.make<TH1F>("cl_e2eRatio0 ","2^{nd} most crystalHit energy/E_{cluster}; E^{2nd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hCl.e2e1Ratio[0]     = clDir.make<TH1F>("cl_e2e1Ratio0","2^{nd} most crystalHit energy/E^{mx}_{hit}; E^{2nd-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
    _hCl.e3eRatio [0]     = clDir.make<TH1F>("cl_e3eRatio0 ","3^{rd} most crystalHit energy/E_{cluster}; E^{3rd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hCl.e4eRatio [0]     = clDir.make<TH1F>("cl_e4eRatio0","4^{th} most crystalHit energy/E^{mx}_{hit}; E^{4th-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
  
    _hCl.eRDist   [0]     = clDir.make<TH2F>("cl_eRDist0"   ,"r_{cog} vs E_{cluster}; r_{cog} [mm]; E_{cluster} [MeV]"        , 500, 300, 800, 400,   0,   200);
    _hCl.e1eRDist [0]     = clDir.make<TH2F>("cl_e1eRDist0"  ,"r_{cog} vs E1/E_{cluster}; r_{cog} [mm]; E^{max}_{hit}/E_{cluster} "        , 500, 300, 800, 100,   0,   1.2);
    _hCl.e2eRDist [0]     = clDir.make<TH2F>("cl_e2eRDist0"  ,"r_{cog} vs E2/E_{cluster}; r_{cog} [mm]; E^{2nd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hCl.e2e1RDist[0]     = clDir.make<TH2F>("cl_e2e1RDist0","r_{cog} vs E2/E1; r_{cog} [mm]; E^{2nd-max}_{hit}/E^{max}_{hit}", 500, 300, 800, 100,   0,   1.2);
    _hCl.e3eRDist [0]     = clDir.make<TH2F>("cl_e3eRDist0"  ,"r_{cog} vs E3/E_{cluster}; r_{cog} [mm]; E^{3rd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hCl.e4eRDist [0]     = clDir.make<TH2F>("cl_e4eRDist0","r_{cog} vs E4/E_{cluster}; r_{cog} [mm]; E^{4rd-max}_{hit}/E_{cluster}"   , 500, 300, 800, 100,   0,   1.2);
    _hCl.nCrRDist [0]     = clDir.make<TH2F>("cl_nCrRDist0" ,"r_{cog} vs nCrystals; r_{cog} [mm]; nCrystals"                  , 500, 300, 800,  20,   0,    20);


    art::TFileDirectory cl1Dir = tfs->mkdir("cluster_photon");
    _hCl.energy   [1]     = cl1Dir.make<TH1F>("cl_E1"       ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hCl.time     [1]     = cl1Dir.make<TH1F>("cl_time1"    ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hCl.nCr      [1]     = cl1Dir.make<TH1F>("cl_nCr1"     ,"Cluster size distribution; nCrystalHits"  , 21, -0.5, 20.5);
    _hCl.rDist    [1]     = cl1Dir.make<TH1F>("cl_rDist1"   ,"Time distribution; cluster radial distance [mm]"  , 500, 300, 800);
    _hCl.e1eRatio [1]     = cl1Dir.make<TH1F>("cl_e1eRatio1","most crystalHit energy/E_{cluster}; E^{max}_{hit}/E_{cluster} "  , 100, 0, 1.);
    _hCl.e2eRatio [1]     = cl1Dir.make<TH1F>("cl_e2eRatio1","2^{nd} most crystalHit energy/E_{cluster}; E^{2nd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hCl.e2e1Ratio[1]     = cl1Dir.make<TH1F>("cl_e2e1Ratio1","2^{nd} most crystalHit energy/E^{mx}_{hit}; E^{2nd-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
    _hCl.e3eRatio [1]     = cl1Dir.make<TH1F>("cl_e3eRatio1","3^{rd} most crystalHit energy/E_{cluster}; E^{3rd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hCl.e4eRatio [1]     = cl1Dir.make<TH1F>("cl_e4eRatio1","4^{th} most crystalHit energy/E^{mx}_{hit}; E^{4th-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
  
    _hCl.eRDist   [1]     = cl1Dir.make<TH2F>("cl_eRDist1"   ,"r_{cog} vs E_{cluster}; r_{cog} [mm]; E_{cluster} [MeV]"        , 500, 300, 800, 400,   0,   200);
    _hCl.e1eRDist [1]     = cl1Dir.make<TH2F>("cl_e1eRDist1"  ,"r_{cog} vs E1/E_{cluster}; r_{cog} [mm]; E^{max}_{hit}/E_{cluster} "        , 500, 300, 800, 100,   0,   1.2);
    _hCl.e2eRDist [1]     = cl1Dir.make<TH2F>("cl_e2eRDist1"  ,"r_{cog} vs E2/E_{cluster}; r_{cog} [mm]; E^{2nd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hCl.e2e1RDist[1]     = cl1Dir.make<TH2F>("cl_e2e1RDist1","r_{cog} vs E2/E1; r_{cog} [mm]; E^{2nd-max}_{hit}/E^{max}_{hit}", 500, 300, 800, 100,   0,   1.2);
    _hCl.e3eRDist [1]     = cl1Dir.make<TH2F>("cl_e3eRDist1"  ,"r_{cog} vs E3/E_{cluster}; r_{cog} [mm]; E^{3rd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hCl.e4eRDist [1]     = cl1Dir.make<TH2F>("cl_e4eRDist1","r_{cog} vs E4/E_{cluster}; r_{cog} [mm]; E^{4rd-max}_{hit}/E_{cluster}"   , 500, 300, 800, 100,   0,   1.2);
    _hCl.nCrRDist [1]     = cl1Dir.make<TH2F>("cl_nCrRDist1" ,"r_{cog} vs nCrystals; r_{cog} [mm]; nCrystals"                  , 500, 300, 800,  20,   0,    20);


    art::TFileDirectory cl2Dir = tfs->mkdir("cluster_electron");
    _hCl.energy   [2]     = cl2Dir.make<TH1F>("cl_E2"       ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hCl.time     [2]     = cl2Dir.make<TH1F>("cl_time2"    ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hCl.nCr      [2]     = cl2Dir.make<TH1F>("cl_nCr2"     ,"Cluster size distribution; nCrystalHits"  , 21, -0.5, 20.5);
    _hCl.rDist    [2]     = cl2Dir.make<TH1F>("cl_rDist2"   ,"Time distribution; cluster radial distance [mm]"  , 500, 300, 800);
    _hCl.e1eRatio [2]     = cl2Dir.make<TH1F>("cl_e1eRatio2","most crystalHit energy/E_{cluster}; E^{max}_{hit}/E_{cluster} "  , 100, 0, 1.);
    _hCl.e2eRatio [2]     = cl2Dir.make<TH1F>("cl_e2eRatio2","2^{nd} most crystalHit energy/E_{cluster}; E^{2nd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hCl.e2e1Ratio[2]     = cl2Dir.make<TH1F>("cl_e2e1Ratio2","2^{nd} most crystalHit energy/E^{mx}_{hit}; E^{2nd-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
    _hCl.e3eRatio [2]     = cl2Dir.make<TH1F>("cl_e3eRatio2","3^{rd} most crystalHit energy/E_{cluster}; E^{3rd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hCl.e4eRatio [2]     = cl2Dir.make<TH1F>("cl_e4eRatio2","4^{th} most crystalHit energy/E^{mx}_{hit}; E^{4th-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
  
    _hCl.eRDist   [2]     = cl2Dir.make<TH2F>("cl_eRDist2"   ,"r_{cog} vs E_{cluster}; r_{cog} [mm]; E_{cluster} [MeV]"        , 500, 300, 800, 400,   0,   200);
    _hCl.e1eRDist [2]     = cl2Dir.make<TH2F>("cl_e1eRDist2"  ,"r_{cog} vs E1/E_{cluster}; r_{cog} [mm]; E^{max}_{hit}/E_{cluster} "        , 500, 300, 800, 100,   0,   1.2);
    _hCl.e2eRDist [2]     = cl2Dir.make<TH2F>("cl_e2eRDist2"  ,"r_{cog} vs E2/E_{cluster}; r_{cog} [mm]; E^{2nd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hCl.e2e1RDist[2]     = cl2Dir.make<TH2F>("cl_e2e1RDist2","r_{cog} vs E2/E1; r_{cog} [mm]; E^{2nd-max}_{hit}/E^{max}_{hit}", 500, 300, 800, 100,   0,   1.2);
    _hCl.e3eRDist [2]     = cl2Dir.make<TH2F>("cl_e3eRDist2"  ,"r_{cog} vs E3/E_{cluster}; r_{cog} [mm]; E^{3rd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hCl.e4eRDist [2]     = cl2Dir.make<TH2F>("cl_e4eRDist2","r_{cog} vs E4/E_{cluster}; r_{cog} [mm]; E^{4rd-max}_{hit}/E_{cluster}"   , 500, 300, 800, 100,   0,   1.2);
    _hCl.nCrRDist [2]     = cl2Dir.make<TH2F>("cl_nCrRDist2" ,"r_{cog} vs nCrystals; r_{cog} [mm]; nCrystals"                  , 500, 300, 800,  20,   0,    20);


    art::TFileDirectory evtDir = tfs->mkdir("evt");
    _hEvt.nCl[0]       = evtDir.make<TH1F>("evt_nCl0"    ,"n-cluster distrbution; nClusters",  201, -0.5, 200.5);
    _hEvt.nCl[1]       = evtDir.make<TH1F>("evt_nCl1"    ,"n-cluster distrbution photons with E>70 MeV; nClusters",  201, -0.5, 200.5);
    _hEvt.nCl[2]       = evtDir.make<TH1F>("evt_nCl2"    ,"n-cluster distrbution electrons with E>70 MeV; nClusters",  201, -0.5, 200.5);
  
    art::TFileDirectory lhDir = tfs->mkdir("likelihood");
    _hLh.hEnergy     = lhDir.make<TH1F>("lh_clEnergy"    ,"lh from cluster energy; log(P_{signal}(E)/P_{bkg}(E))"           ,  400, -100, 100);
    _hLh.hE1E        = lhDir.make<TH1F>("lh_clE1E"        ,"lh from E1=E_{seed}/E_{cluster}; log(P_{signal}(E1)/P_{bkg}(E1))",  400, -100, 100);
    _hLh.hNCr        = lhDir.make<TH1F>("lh_clSize"      ,"lh from cluster size; log(P_{signal}(Size)/P_{bkg}(Size))"       ,  400, -100, 100);
    _hLh.hRDist      = lhDir.make<TH1F>("lh_rDist"       ,"lh from E1=E_{seed}/E_{cluster}; log(P_{signal}(nCr)/P_{bkg}(nCr))",  400, -100, 100);
    _hLh.hE2E        = lhDir.make<TH1F>("lh_clE2E"        ,"lh from E2=E_{2nd/E_{cluster}; log(P_{signal}(E2)/P_{bkg}(E2))"  ,  400, -100, 100);
    _hLh.hE2E1       = lhDir.make<TH1F>("lh_clE2E1"      ,"lh from E2/E1; log(P_{signal}(E2/E1)/P_{bkg}(E2/E1))"            ,  400, -100, 100);
    _hLh.hE3E        = lhDir.make<TH1F>("lh_clE3E"        ,"lh from E_{3rd/E_{cluster}; log(P_{signal}(E3/E)/P_{bkg}(E3/E))"  ,  400, -100, 100);
    _hLh.hE4E        = lhDir.make<TH1F>("lh_clE4E"      ,"lh from  E_{3rd/E_{cluster}; log(P_{signal}(E4/E)/P_{bkg}(E4/E))"            ,  400, -100, 100);
	            
    _hLh.hCorEnergy  = lhDir.make<TH1F>("lhCor_clEnergy" ,"lh from cluster energy vs rdist; log(P_{signal}(E)/P_{bkg}(E))"           ,  400, -100, 100);
    _hLh.hCorE1E     = lhDir.make<TH1F>("lhCor_clE1E"     ,"lh from rdist vs E1; log(P_{signal}(E1)/P_{bkg}(E1))",  4000, -1000, 1000);
    _hLh.hCorE2E     = lhDir.make<TH1F>("lhCor_clE2E"     ,"lh from rdist vs E2; log(P_{signal}(E2)/P_{bkg}(E2))"  ,  4000, -1000, 1000);
    _hLh.hCorE2E1    = lhDir.make<TH1F>("lhCor_clE2E1"   ,"lh from rdist vsE2/E1; log(P_{signal}(E2/E1)/P_{bkg}(E2/E1))"            ,  4000, -1000, 1000);
    _hLh.hCorE3E     = lhDir.make<TH1F>("lhCor_clE3E"     ,"lh from rdist vs E3; log(P_{signal}(E3)/P_{bkg}(E3))"            ,  4000, -1000, 1000);
    _hLh.hCorE4E     = lhDir.make<TH1F>("lhCor_clE4E"     ,"lh from rdist vs E4; log(P_{signal}(E4)/P_{bkg}(E4))"            ,  4000, -1000, 1000);
    _hLh.hCorNCr     = lhDir.make<TH1F>("lhCor_clNCr"    ,"lh from rdist vs nCrystal; log(P_{signal}(nCr)/P_{bkg}(nCr))"            ,  4000, -1000, 1000);
  
    _hLh.hTot     = lhDir.make<TH1F>("lh_tot" ,"Global likelihood; #sum log(P_{signal}(x_{i})/P_{bkg}(x_{i}))"  , 400, -200, 200);
    _hLh.hCorTot  = lhDir.make<TH1F>("lhCor_tot" ,"Global likelihood; #sum log(P_{signal}(x_{i})/P_{bkg}(x_{i}))"  , 4000, -1000, 1000);

    art::TFileDirectory trigDir = tfs->mkdir("triggered");
    _hTrig.energy        = trigDir.make<TH1F>("trigcl_E"       ,"Energy distrbution; E [MeV]",  400, 0, 200);
    _hTrig.time          = trigDir.make<TH1F>("trigcl_time"    ,"Time distribution; t [ns]"  , 1700, 0, 1700);
    _hTrig.nCr           = trigDir.make<TH1F>("trigcl_nCr"     ,"Cluster size distribution; nCrystalHits"  , 21, -0.5, 20.5);
    _hTrig.rDist         = trigDir.make<TH1F>("trigcl_rDist"   ,"Time distribution; cluster radial distance [mm]"  , 500, 300, 800);
    _hTrig.e1eRatio      = trigDir.make<TH1F>("trigcl_e1eRatio","most crystalHit energy/E_{cluster}; E^{max}_{hit}/E_{cluster} "  , 100, 0, 1.);
    _hTrig.e2eRatio      = trigDir.make<TH1F>("trigcl_e2eRatio","2^{nd} most crystalHit energy/E_{cluster}; E^{2nd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hTrig.e2e1Ratio     = trigDir.make<TH1F>("trigcl_e2e1Ratio","2^{nd} most crystalHit energy/E^{mx}_{hit}; E^{2nd-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
    _hTrig.e3eRatio      = trigDir.make<TH1F>("trigcl_e3eRatio","3^{rd} most crystalHit energy/E_{cluster}; E^{3rd-max}_{hit}/E_{cluster}"   , 100, 0, 1.);
    _hTrig.e4eRatio      = trigDir.make<TH1F>("trigcl_e4eRatio","4^{th} most crystalHit energy/E^{mx}_{hit}; E^{4th-max}_{hit}/E^{max}_{hit}", 100, 0, 1.);
  
    _hTrig.eRDist        = trigDir.make<TH2F>("trigcl_eRDist"   ,"r_{cog} vs E_{cluster}; r_{cog} [mm]; E_{cluster} [MeV]"        , 500, 300, 800, 400,   0,   200);
    _hTrig.e1eRDist      = trigDir.make<TH2F>("trigcl_e1eRDist"  ,"r_{cog} vs E1/E_{cluster}; r_{cog} [mm]; E^{max}_{hit}/E_{cluster} "        , 500, 300, 800, 100,   0,   1.2);
    _hTrig.e2eRDist      = trigDir.make<TH2F>("trigcl_e2eRDist"  ,"r_{cog} vs E2/E_{cluster}; r_{cog} [mm]; E^{2nd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hTrig.e2e1RDist     = trigDir.make<TH2F>("trigcl_e2e1RDist","r_{cog} vs E2/E1; r_{cog} [mm]; E^{2nd-max}_{hit}/E^{max}_{hit}", 500, 300, 800, 100,   0,   1.2);
    _hTrig.e3eRDist      = trigDir.make<TH2F>("trigcl_e3eRDist"  ,"r_{cog} vs E3/E_{cluster}; r_{cog} [mm]; E^{3rd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    _hTrig.e4eRDist      = trigDir.make<TH2F>("trigcl_e4eRDist","r_{cog} vs E4/E_{cluster}; r_{cog} [mm]; E^{4rd-max}_{hit}/E_{cluster}"   , 500, 300, 800, 100,   0,   1.2);
    _hTrig.nCrRDist      = trigDir.make<TH2F>("trigcl_nCrRDist" ,"r_{cog} vs nCrystals; r_{cog} [mm]; nCrystals"                  , 500, 300, 800,  20,   0,    20);
  }


  //--------------------------------------------------------------------------------//
  void RPCPhotonAna::beginJob(){

    bookHistograms();

  }

  bool     RPCPhotonAna::beginRun(art::Run& ) {
    // mu2e::GeomHandle<mu2e::Calorimeter> ch;
    // _calorimeter = ch.get();
    return true;
  }	

  void     RPCPhotonAna::endJob(){}
  
  //--------------------------------------------------------------------------------
  double   RPCPhotonAna::calculate2DProb(double& Ref, double&Variable, TH1F** Template, double MinX, double Step){
    double     thisprob = 1e-17;

    for (int j=0; j<kNCorHist; ++j){
      double   ref_min = (MinX + (double)j*Step);
      double   ref_max = (MinX + (double)(j+1)*Step);
      if ( ( Ref >= ref_min) && (Ref < ref_max) ){
	thisprob = calculateProb(Variable,  Template[j]);
	break;
      }
    }
    
    return thisprob;
  }


  //--------------------------------------------------------------------------------
  double   RPCPhotonAna::calculateProb(double&Variable, TH1F* Template){

    double     thisprob = 1;
    //not sure how to loop over variable value   i should be event number
    double     binSize   = Template->GetBinWidth(1);
    int        binIndex = (Variable - Template->GetBinLowEdge(1))/binSize;
    //what is the probability for this variable
    int           templateNBins = Template->GetNbinsX();
    if (binIndex > templateNBins) {
      thisprob = 1e-17;
      printf("[CaloAna::calculateProb] %s value = %2.3f index %i out of range %i\n", Template->GetName(), Variable,  binIndex, templateNBins);
    } else{
      thisprob = Template->GetBinContent(binIndex);
    }
      
    if (thisprob < 1e-10) thisprob = 1e-17;

    return thisprob;

  }


  void     RPCPhotonAna::fillGenHistograms    (const GenParticle*           Gen    ){
    double energy  = Gen->momentum().e();
    double timegen = Gen->time();
    double costh   = Gen->momentum().cosTheta();

    _hGen.energy[0]->Fill(energy , _evtWeight);
    _hGen.time  [0]->Fill(timegen, _evtWeight);
    _hGen.pdgId [0]->Fill(Gen->pdgId(), _evtWeight);
    _hGen.cosTh [0]->Fill(costh       , _evtWeight);

    if (_evtWithPhoton){
      _hGen.energy[1]->Fill(energy , _evtWeight);
      _hGen.time  [1]->Fill(timegen, _evtWeight);
      _hGen.pdgId [1]->Fill(Gen->pdgId(), _evtWeight);
      _hGen.cosTh [1]->Fill(costh       , _evtWeight);
    }

    if (_evtWithElectron){
      _hGen.energy[2]->Fill(energy , _evtWeight);
      _hGen.time  [2]->Fill(timegen, _evtWeight);
      _hGen.pdgId [2]->Fill(Gen->pdgId(), _evtWeight);
      _hGen.cosTh [2]->Fill(costh       , _evtWeight);
    }
  }

  //--------------------------------------------------------------------------------//
  void     RPCPhotonAna::fillEventHistograms  (const CaloClusterCollection* ClCol  ){
    int  nCl = (int)ClCol->size();
    
    //count the number of cluster above the threshold
    double  e_min(70.);//FIXME!
    int     nCl1(0);

    const CaloCluster* cl(0);
    
    for (int i=0; i<nCl; ++i){
      cl = &ClCol->at(i);
      double  edep = cl->energyDep();
      if (edep >= e_min) ++nCl1;
    }
    
    _hEvt.nCl[0]   ->Fill(nCl , _evtWeight);
    _hEvt.nCl[1]   ->Fill(nCl1, _evtWeight);

  }

  void     RPCPhotonAna::fillLHHistograms    (const CaloClusterCollection* ClCol ){
    int                nClFast = ClCol->size();

    const CaloCluster*    cluster(0);
    const CaloCrystalHit* crystalHit(0);
    const CaloCluster::CaloCrystalHitPtrVector* caloClusterHits(0);
  
    double    maxLhValue(-1e10), maxLh2DValue(-1e10);
    int       sectionMaxValue(-1), sectionMax2DValue(-1);

    //for loop over the clusters in the calorimeter
    for( int i=0; i<nClFast; ++i){
      cluster=&(ClCol->at(i));
      int clSection = cluster->diskId();
      
      //get the variables needed
      //  ClusterMoments cogCalculator(*_calorimeter, *cluster,clSection);
      // cogCalculator.calculate(_cogType);

      double                   clEnergy = cluster->energyDep();
      double                   clTime   = cluster->time();
      const CLHEP::Hep3Vector& cog      = cluster->cog3Vector();
      double                   xpos     = cog(0);
      double                   ypos     = cog(1);
      double                   rDist    = sqrt(xpos*xpos+ypos*ypos);

      //if (rDist < 550|| rDist > 600)                    continue; 

      caloClusterHits = &cluster->caloCrystalHitsPtrVector();

      int       nCrystalHits = caloClusterHits->size();
      double    maxECrystal(0);
      double    secondmaxECrystal(-1), thirdmaxECrystal(-1), fourthmaxECrystal(-1);
      double    indexMaxECrystal(0), index2ndMaxECrystal(-1), index3rdMaxECrystal(-1);

      //first loop to find the most energetic crystalHit
      for (int j=0; j<nCrystalHits; ++j){
	crystalHit = &(*caloClusterHits->at(j));
	double   crystalEnergy = crystalHit->energyDep();
	if (crystalEnergy > maxECrystal) {
	  maxECrystal        = crystalEnergy;
	  indexMaxECrystal   = j;
	}
      }

      //second loop to find the second most energetic crystalHit
      for (int j=0; j<nCrystalHits; ++j){
	if (j == indexMaxECrystal)              continue;
	crystalHit = &(*caloClusterHits->at(j));

	double crystalEnergy = crystalHit->energyDep();

	if (crystalEnergy > secondmaxECrystal){
	  secondmaxECrystal   = crystalEnergy;
	  index2ndMaxECrystal = j;
	}
      }
      
      for (int j=0; j<nCrystalHits; ++j){
	if (j == indexMaxECrystal || 
	    j == index2ndMaxECrystal  )         continue;
	crystalHit = &(*caloClusterHits->at(j));

	double crystalEnergy = crystalHit->energyDep();

	if (crystalEnergy > thirdmaxECrystal){
	  thirdmaxECrystal    = crystalEnergy;
	  index3rdMaxECrystal = j;
	}
      }
    
      for (int j=0; j<nCrystalHits; ++j){
	if (j == indexMaxECrystal  || 
	    j == index2ndMaxECrystal ||
	    j == index3rdMaxECrystal )          continue;
	crystalHit = &(*caloClusterHits->at(j));

	double crystalEnergy = crystalHit->energyDep();

	if (crystalEnergy > fourthmaxECrystal){
	  fourthmaxECrystal    = crystalEnergy;
	}
      }       
      
      double   e1Ratio   = maxECrystal/clEnergy;
      double   e2Ratio   = secondmaxECrystal/clEnergy;
      double   e3Ratio   = thirdmaxECrystal/clEnergy;
      double   e4Ratio   = fourthmaxECrystal/clEnergy;
      double   e2e1Ratio = secondmaxECrystal/maxECrystal;

      //evaluate the likelihood values using the 1D tempaltes
      double  clusterSize = (double) nCrystalHits;
      
      double  signalEnergyProb      = calculateProb(clEnergy, _signalEnergy);
      double  bkgEnergyProb         = calculateProb(clEnergy, _bkgEnergy);
      double  logEnergyRatio        = log(signalEnergyProb/bkgEnergyProb);
  
      double  signalE1Prob          = calculateProb(e1Ratio, _signale1eRatio);
      double  bkgE1Prob             = calculateProb(e1Ratio, _bkge1eRatio   );
      double  logE1Ratio            = log(signalE1Prob/bkgE1Prob);
  
      double  signalClusterSizeProb = calculateProb(clusterSize, _signalnCr);
      double  bkgClusterSizeProb    = calculateProb(clusterSize , _bkgnCr);
      double  logClusterSizeRatio   = log(signalClusterSizeProb/bkgClusterSizeProb);
      
      double  signalRDistProb       = calculateProb(rDist, _signalrDist);
      double  bkgRDistProb          = calculateProb(rDist, _bkgrDist);
      double  logRDistRatio         = log(signalRDistProb/bkgRDistProb);

      double  signalE2Prob          = calculateProb(e2Ratio, _signale2eRatio);
      double  bkgE2Prob             = calculateProb(e2Ratio, _bkge2eRatio);
      double  logE2Ratio            = log(signalE2Prob/bkgE2Prob);
  
      double  signalE2E1Prob        = calculateProb(e2e1Ratio, _signale2e1Ratio);
      double  bkgE2E1Prob           = calculateProb(e2e1Ratio, _bkge2e1Ratio);
      double  logE2E1Ratio          = log(signalE2E1Prob/bkgE2E1Prob);
 
      double  signalE3Prob          = calculateProb(e3Ratio, _signale3eRatio);
      double  bkgE3Prob             = calculateProb(e3Ratio, _bkge3eRatio);
      double  logE3Ratio            = log(signalE3Prob/bkgE3Prob);

      double  signalE4Prob          = calculateProb(e4Ratio, _signale4eRatio);
      double  bkgE4Prob             = calculateProb(e4Ratio, _bkge4eRatio);
      double  logE4Ratio            = log(signalE4Prob/bkgE4Prob);

      //fill the histograms
      _hLh.hEnergy->Fill(logEnergyRatio     , _evtWeight);
      _hLh.hE1E   ->Fill(logE1Ratio         , _evtWeight);
      _hLh.hNCr   ->Fill(logClusterSizeRatio, _evtWeight);
      _hLh.hRDist ->Fill(logRDistRatio      , _evtWeight);
      _hLh.hE2E   ->Fill(logE2Ratio         , _evtWeight);
      _hLh.hE2E1  ->Fill(logE2E1Ratio       , _evtWeight);
      _hLh.hE3E   ->Fill(logE3Ratio         , _evtWeight);
      _hLh.hE4E   ->Fill(logE4Ratio         , _evtWeight);

      double  lhValue = logEnergyRatio + logE1Ratio + logClusterSizeRatio + logRDistRatio;
      if (nCrystalHits>=2) lhValue += logE2Ratio + logE2E1Ratio;
      if (nCrystalHits>=3) lhValue += logE3Ratio;
      if (nCrystalHits>=4) lhValue += logE4Ratio;
      if (lhValue > maxLhValue) {
	maxLhValue      = lhValue;
	sectionMaxValue = clSection;
      }
 
      //now use the 2D histograms
      // double  signalRDistProb         = calculateProb(rDist, _signalHist1D[clSection][3]);
      // double  bkgRDistProb            = calculateProb(rDist, _bkgHist1D   [clSection][3]);
      // double  logRDistRatio           = log(signalRDistProb/bkgRDistProb);
      // double  signalRDist2DProb       = calculate2DProb(clEnergy, rDist, _signalCorrHist1D[clSection][0], _minClEnergy, _clEStep);
      // double  bkgRDist2DProb          = calculate2DProb(clEnergy, rDist, _bkgCorrHist1D   [clSection][0], _minClEnergy, _clEStep);
      // double  logRDist2DRatio         = log(signalRDist2DProb/bkgRDist2DProb);
      double  signalEnergy2DProb      = calculate2DProb(rDist, clEnergy, _signalCorrHistEnergy, _minRDist, _rDistStep);
      double  bkgEnergy2DProb         = calculate2DProb(rDist, clEnergy, _bkgCorrHistEnergy, _minRDist, _rDistStep);
      double  logEnergy2DRatio        = log(signalEnergy2DProb/bkgEnergy2DProb);

      double  signalE12DProb          = calculate2DProb(rDist, e1Ratio, _signalCorrHiste1e, _minRDist, _rDistStep);
      double  bkgE12DProb             = calculate2DProb(rDist, e1Ratio, _bkgCorrHiste1e, _minRDist, _rDistStep);
      double  logE12DRatio            = log(signalE12DProb/bkgE12DProb);
  
      double  signalE22DProb          = calculate2DProb(rDist, e2Ratio, _signalCorrHiste2e, _minRDist, _rDistStep);
      double  bkgE22DProb             = calculate2DProb(rDist, e2Ratio, _bkgCorrHiste2e, _minRDist, _rDistStep);
      double  logE22DRatio            = log(signalE22DProb/bkgE22DProb);
  
      double  signalE2E12DProb        = calculate2DProb(rDist, e2e1Ratio, _signalCorrHiste2e1, _minRDist, _rDistStep);
      double  bkgE2E12DProb           = calculate2DProb(rDist, e2e1Ratio, _bkgCorrHiste2e1, _minRDist, _rDistStep);
      double  logE2E12DRatio          = log(signalE2E12DProb/bkgE2E12DProb);
 
      double  signalE32DProb          = calculate2DProb(rDist, e3Ratio, _signalCorrHiste3e, _minRDist, _rDistStep);
      double  bkgE32DProb             = calculate2DProb(rDist, e3Ratio, _bkgCorrHiste3e, _minRDist, _rDistStep);
      double  logE32DRatio            = log(signalE32DProb/bkgE32DProb);
  
      double  signalE42DProb          = calculate2DProb(rDist, e4Ratio, _signalCorrHiste4e, _minRDist, _rDistStep);
      double  bkgE42DProb             = calculate2DProb(rDist, e4Ratio, _bkgCorrHiste4e, _minRDist, _rDistStep);
      double  logE42DRatio            = log(signalE42DProb/bkgE42DProb);
  
      double  signalClSize2DProb      = calculate2DProb(rDist, clusterSize, _signalCorrHistnCr, _minRDist, _rDistStep);
      double  bkgClSize2DProb         = calculate2DProb(rDist, clusterSize , _bkgCorrHistnCr, _minRDist, _rDistStep);
      double  logClSize2DRatio        = log(signalClSize2DProb/bkgClSize2DProb);

      //fill the histograms
      _hLh.hCorEnergy ->Fill(logEnergy2DRatio    , _evtWeight);
      _hLh.hCorE1E    ->Fill(logE12DRatio        , _evtWeight);
      _hLh.hCorE2E    ->Fill(logE22DRatio        , _evtWeight);
      _hLh.hCorE2E1   ->Fill(logE2E12DRatio      , _evtWeight);
      _hLh.hCorE3E    ->Fill(logE32DRatio        , _evtWeight);
      _hLh.hCorE4E    ->Fill(logE42DRatio        , _evtWeight);
      _hLh.hCorNCr    ->Fill(logClSize2DRatio    , _evtWeight);
      

      double  lh2DValue = logRDistRatio + logEnergy2DRatio + logE12DRatio + logClSize2DRatio;
      if (nCrystalHits>=2) lh2DValue += logE22DRatio + logE2E12DRatio;
      if (nCrystalHits>=3) lh2DValue += logE32DRatio;
      if (nCrystalHits>=4) lh2DValue += logE42DRatio;
      // printf("logE22DRatio = %2.3f logE2E12DRatio = %2.3f logE32DRatio = %2.3f logE42DRatio = %2.3f\n",
      // 	     logE22DRatio, logE2E12DRatio, logE32DRatio, logE42DRatio);
      if (lh2DValue > maxLh2DValue) {
	maxLh2DValue      = lh2DValue;
	sectionMax2DValue = clSection;
      }

      //Skip events that don't pass trigger
      if(maxLh2DValue < 40)   continue;
      //Fill historgams of events that pass trigger
      
      _hTrig.energy     ->Fill(clEnergy, _evtWeight);
      _hTrig.time       ->Fill(clTime, _evtWeight);
      _hTrig.nCr        ->Fill(nCrystalHits, _evtWeight);
      _hTrig.rDist      ->Fill(rDist, _evtWeight);
      _hTrig.nCrRDist   ->Fill(rDist, nCrystalHits, _evtWeight);
      //_hCl.timeRDist  ->Fill(clRDist, clTime, _evtWeight);

      if (clEnergy    > 1e-3){
	_hTrig.e1eRatio ->Fill(maxECrystal/clEnergy, _evtWeight);
	_hTrig.eRDist   ->Fill( rDist,                   clEnergy, _evtWeight);
	_hTrig.e1eRDist  ->Fill( rDist,       maxECrystal/clEnergy, _evtWeight);
	if (nCrystalHits>=2) {
	  _hTrig.e2eRatio ->Fill(secondmaxECrystal/clEnergy, _evtWeight);
	  _hTrig.e2eRDist ->Fill( rDist, secondmaxECrystal/maxECrystal, _evtWeight);
	  if (maxECrystal > 1e-3){
	    _hTrig.e2e1Ratio->Fill(secondmaxECrystal/maxECrystal, _evtWeight);
	    _hTrig.e2e1RDist->Fill( rDist, secondmaxECrystal/maxECrystal, _evtWeight);
	  }
	  if (nCrystalHits>=3) {
	    _hTrig.e3eRatio ->Fill(thirdmaxECrystal/clEnergy, _evtWeight);
	    _hTrig.e3eRDist  ->Fill( rDist,  thirdmaxECrystal/clEnergy, _evtWeight);
	  }
	  if (nCrystalHits>=4) {
	    _hTrig.e4eRatio ->Fill(fourthmaxECrystal/clEnergy, _evtWeight);
	    _hTrig.e4eRDist  ->Fill( rDist,  fourthmaxECrystal/clEnergy, _evtWeight);
	  }
	}
      }


    }
    
    if (sectionMaxValue>=0)    _hLh.hTot ->Fill(maxLhValue);

    if (sectionMax2DValue>=0)    _hLh.hCorTot ->Fill(maxLh2DValue);

  }

  //--------------------------------------------------------------------------------//
  void     RPCPhotonAna::fillClusterHistograms(const CaloCluster*           Cluster){
    
    const CLHEP::Hep3Vector &cog        = Cluster->cog3Vector();
    double                  xpos        = cog.x();
    double                  ypos        = cog.y();

    double                  clEnergy    = Cluster->energyDep();
    double                  clTime      = Cluster->time();
    double                  clRDist     = sqrt(xpos*xpos+ypos*ypos);

    const CaloCluster::CaloCrystalHitPtrVector caloClusterHits = Cluster->caloCrystalHitsPtrVector();
    int                     nCrystalHits = caloClusterHits.size();
    const CaloCrystalHit*   crystalHit(0);

    double                  maxECrystal(0);
    double    secondmaxECrystal(-1), thirdmaxECrystal(-1), fourthmaxECrystal(-1);
    double    indexMaxECrystal(0), index2ndMaxECrystal(-1), index3rdMaxECrystal(-1);  

    for (int j=0; j<nCrystalHits; ++j){
      crystalHit = &(*caloClusterHits.at(j));

      double crystalEnergy = crystalHit->energyDep();
      if (crystalEnergy > maxECrystal)
	{
	  maxECrystal      = crystalEnergy;
	  indexMaxECrystal = j;
	}
    }

    for (int j=0; j<nCrystalHits; ++j){
      if (j == indexMaxECrystal)             continue;
      crystalHit = &(*caloClusterHits.at(j));

      double crystalEnergy = crystalHit->energyDep();

      if (crystalEnergy > secondmaxECrystal){
	secondmaxECrystal = crystalEnergy;
	index2ndMaxECrystal = j;
      }
    }

    for (int j=0; j<nCrystalHits; ++j){
      if (j == indexMaxECrystal ||
          j == index2ndMaxECrystal)             continue;
      crystalHit = &(*caloClusterHits.at(j));

      double crystalEnergy = crystalHit->energyDep();

      if (crystalEnergy > thirdmaxECrystal){
	thirdmaxECrystal = crystalEnergy;
	index3rdMaxECrystal = j;
      }
    }

    for (int j=0; j<nCrystalHits; ++j){
      if (j == indexMaxECrystal  || 
	  j == index2ndMaxECrystal ||
	  j == index3rdMaxECrystal )          continue;
      crystalHit = &(*caloClusterHits.at(j));

      double crystalEnergy = crystalHit->energyDep();

      if (crystalEnergy > fourthmaxECrystal){
	fourthmaxECrystal    = crystalEnergy;
      }
    }
  
    //fill the cluster info
    _hCl.energy     [0]->Fill(clEnergy, _evtWeight);
    _hCl.time       [0]->Fill(clTime, _evtWeight);
    _hCl.nCr        [0]->Fill(nCrystalHits, _evtWeight);
    _hCl.rDist      [0]->Fill(clRDist, _evtWeight);
    _hCl.nCrRDist   [0]->Fill(clRDist, nCrystalHits, _evtWeight);
    //_hCl.timeRDist  [0]->Fill(clRDist, clTime, _evtWeight);

    if (clEnergy    > 1e-3){
      _hCl.e1eRatio [0]->Fill(maxECrystal/clEnergy, _evtWeight);
      _hCl.eRDist   [0]->Fill( clRDist,                   clEnergy, _evtWeight);
      _hCl.e1eRDist  [0]->Fill( clRDist,       maxECrystal/clEnergy, _evtWeight);
      if (nCrystalHits>=2) {
	_hCl.e2eRatio [0]->Fill(secondmaxECrystal/clEnergy, _evtWeight);
	_hCl.e2eRDist [0]->Fill( clRDist, secondmaxECrystal/clEnergy, _evtWeight);
       	if (maxECrystal > 1e-3){
	  _hCl.e2e1Ratio[0]->Fill(secondmaxECrystal/maxECrystal, _evtWeight);
	  _hCl.e2e1RDist[0]->Fill( clRDist, secondmaxECrystal/maxECrystal, _evtWeight);
	}
	if (nCrystalHits>=3) {
	  _hCl.e3eRatio [0]->Fill(thirdmaxECrystal/clEnergy, _evtWeight);
	  _hCl.e3eRDist  [0]->Fill( clRDist,  thirdmaxECrystal/clEnergy, _evtWeight);
	}
	if (nCrystalHits>=4) {
	  _hCl.e4eRatio [0]->Fill(fourthmaxECrystal/clEnergy, _evtWeight);
	  _hCl.e4eRDist  [0]->Fill( clRDist,  fourthmaxECrystal/clEnergy, _evtWeight);
	}
      }
    }

    //fill the cluster info only if this cluster is associated with the generated RPC photon
    if ( ( _evtWithPhoton )           && 
	 (fabs(_vdetT - clTime) < 6)){
      _hCl.energy     [1]->Fill(clEnergy, _evtWeight);
      _hCl.time       [1]->Fill(clTime, _evtWeight);
      _hCl.nCr        [1]->Fill(nCrystalHits, _evtWeight);
      _hCl.rDist      [1]->Fill(clRDist, _evtWeight);
      _hCl.nCrRDist   [1]->Fill(clRDist, nCrystalHits, _evtWeight);
      //_hCl.timeRDist  [1]->Fill(clRDist, clTime, _evtWeight);

      if (clEnergy    > 1e-3){
	_hCl.e1eRatio   [1]->Fill(maxECrystal/clEnergy, _evtWeight);
	_hCl.eRDist     [1]->Fill( clRDist,                   clEnergy, _evtWeight);
	_hCl.e1eRDist    [1]->Fill( clRDist,       maxECrystal/clEnergy, _evtWeight);
	if (nCrystalHits>=2) {
	  _hCl.e2eRatio   [1]->Fill(secondmaxECrystal/clEnergy, _evtWeight);
	  _hCl.e2eRDist  [1]->Fill( clRDist, secondmaxECrystal/clEnergy, _evtWeight);
	  if (maxECrystal > 1e-3){
	    _hCl.e2e1Ratio  [1]->Fill(secondmaxECrystal/maxECrystal, _evtWeight);
	    _hCl.e2e1RDist  [1]->Fill( clRDist, secondmaxECrystal/maxECrystal, _evtWeight);
	  }
	  if (nCrystalHits>=3) {
	    _hCl.e3eRatio   [1]->Fill(thirdmaxECrystal/clEnergy, _evtWeight);
	    _hCl.e3eRDist    [1]->Fill( clRDist,  thirdmaxECrystal/clEnergy, _evtWeight);
	  }
	  if (nCrystalHits>=4) {
	    _hCl.e4eRatio   [1]->Fill(fourthmaxECrystal/clEnergy, _evtWeight);
	    _hCl.e4eRDist    [1]->Fill( clRDist,  fourthmaxECrystal/clEnergy, _evtWeight);
	  }
	}
      }
    }
    //fill the cluster info only if this cluster is associated with an electron
    if ( ( _evtWithElectron )           && 
	 ( fabs(_vdetT_e - clTime) < 6)){
      _hCl.energy     [2]->Fill(clEnergy, _evtWeight);
      _hCl.time       [2]->Fill(clTime, _evtWeight);
      _hCl.nCr        [2]->Fill(nCrystalHits, _evtWeight);
      _hCl.rDist      [2]->Fill(clRDist, _evtWeight);
      _hCl.nCrRDist   [2]->Fill(clRDist, nCrystalHits, _evtWeight);
      // _hCl.timeRDist  [2]->Fill(clRDist, clTime, _evtWeight);

      if (clEnergy    > 1e-3){
	_hCl.e1eRatio   [2]->Fill(maxECrystal/clEnergy, _evtWeight);
	_hCl.eRDist     [2]->Fill( clRDist,                   clEnergy, _evtWeight);
	_hCl.e1eRDist    [2]->Fill( clRDist,       maxECrystal/clEnergy, _evtWeight);
	if (nCrystalHits>=2) {
	  _hCl.e2eRatio   [2]->Fill(secondmaxECrystal/clEnergy, _evtWeight);
	  _hCl.e2eRDist  [2]->Fill( clRDist, secondmaxECrystal/clEnergy, _evtWeight);
	  if (maxECrystal > 1e-3){
	    _hCl.e2e1Ratio  [2]->Fill(secondmaxECrystal/maxECrystal, _evtWeight);
	    _hCl.e2e1RDist  [2]->Fill( clRDist, secondmaxECrystal/maxECrystal, _evtWeight);
	  }
	  if (nCrystalHits>=3) {
	    _hCl.e3eRatio   [2]->Fill(thirdmaxECrystal/clEnergy, _evtWeight);
	    _hCl.e3eRDist    [2]->Fill( clRDist,  thirdmaxECrystal/clEnergy, _evtWeight);
	  }
	  if (nCrystalHits>=4) {
	    _hCl.e4eRatio   [2]->Fill(fourthmaxECrystal/clEnergy, _evtWeight);
	    _hCl.e4eRDist    [2]->Fill( clRDist,  fourthmaxECrystal/clEnergy, _evtWeight);
	  }
	}
      }
    }
  }


  //--------------------------------------------------------------------------------//
  void     RPCPhotonAna::fillVDetHistograms   (const StepPointMC*           Step   ){
    art::Ptr<SimParticle> const& simptr = Step->simParticle();

    //    SimParticle const& sim  = *simptr;
    
    //    if ( !sim.fromGenerator() )   return;
    
    double hitTime = fmod(Step->time() + _toff.totalTimeOffset(simptr), _mbtime);
    // if (hitTime < _mbbuffer) {
    //   if (hitTime+_mbtime > _blindTime) {
    // 	hitTime = hitTime + _mbtime;
    //   }
    // }
    // else {
    //   if (hitTime > (_mbtime - _mbbuffer)) {
    // 	if (hitTime - _mbtime > _blindTime) {
    // 	  hitTime =   hitTime - _mbtime;
    // 	}
    //   }
    // }

          	    
    double  hitP  = Step->momentum().mag();
    double  pdgId = Step->simParticle()->pdgId();
    double  hitM  = _pdt->particle(pdgId).ref().mass();
    double  hitE  = sqrt(hitP*hitP + hitM*hitM);
    double costh   = Step->momentum().cosTheta();
    
    //    if (hitE < 20.)                      return;

    //check if the particle is a photon
    //double radius = sqrt(xpos*xpos + ypos*ypos);
    // double vdcosthetahere = costheta;

    if (pdgId == 22 && (hitE > _vdetE)) {
      _evtWithPhoton    = true;
      _vdetE            = hitE;
      _vdetX            = Step->position().x()+3904.;
      _vdetY            = Step->position().y();
      _vdetT            = hitTime;
      _hVDet.energy[1]->Fill( hitE   , _evtWeight);
      _hVDet.time  [1]->Fill( hitTime, _evtWeight);
      _hVDet.cosTh [1]->Fill( costh , _evtWeight);
    }

    if (pdgId == 11 && (hitE > _vdetE_e)) {
      _evtWithElectron    = true;
      _vdetE_e            = hitE;
      _vdetX_e            = Step->position().x()+3904.;
      _vdetY_e            = Step->position().y();
      _vdetT_e            = hitTime;
      _hVDet.energy[2]->Fill( hitE   , _evtWeight);
      _hVDet.time  [2]->Fill( hitTime, _evtWeight);
      _hVDet.cosTh [2]->Fill( costh , _evtWeight);
    }
    
    _hVDet.energy[0]->Fill( hitE   , _evtWeight);
    _hVDet.time  [0]->Fill( hitTime, _evtWeight);
    _hVDet.pdgId [0]->Fill( pdgId  , _evtWeight);
    _hVDet.cosTh [0]->Fill( costh , _evtWeight);
    
    
    
  }


  //changed this from   void CaloAnalysis::analyze(const art::Event& event) {
  bool     RPCPhotonAna::filter(art::Event& event) {

    ++_nProcess;
    int  _iev            = event.id().event(); 
    //skip odd events to build templates. Change to ==1 to make liklyhood histograms
    if (_iev %2 == 1) 
      return false;

    //initialize
    _evtWithPhoton    = false;
    _evtWithElectron    = false;
    _vdetX            = -9999.;
    _vdetY            = -9999.;
    _vdetT            =-9999.;
    _vdetE            = 40.;//-9999.;
    _vdetX_e            = -9999.;
    _vdetY_e            = -9999.;
    _vdetT              = -9999.;
    _vdetE_e            = 40.;//-9999.;

    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from RPCPhotonAna =  "<<_nProcess <<std::endl;
   
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime=accPar->deBuncherPeriod;
    _toff.updateMap(event);

    //Handle to VD steps
    _vdHitsSize = 0;
    art::Handle<StepPointMCCollection> vdStepsHandle;
    event.getByLabel(_virtualDetectorLabel, vdStepsHandle);
    const StepPointMCCollection *vdHits(0);
    if (vdStepsHandle.isValid()) {
      vdHits = vdStepsHandle.product();
      _vdHitsSize = vdHits->size();}

    //get info from the generator
    _genSize = 0;
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel, gensHandle);
    const GenParticleCollection *genParticles(0);
    const GenParticle *gen(0);
    if (gensHandle.isValid()){
      genParticles = gensHandle.product();
      _genSize =  genParticles->size();
      if (_genSize > 0) gen = &genParticles->at(0);
    }

    //get the weight of the event    
    // _evtWeight  = 0;
    // double x, w{20.};
    // double K = gen->momentum().e();
    // double KMax = 90.1;
    // x = K/KMax;
    // if (x < 1)  _evtWeight = w*(1-2*x+2*x*x)*(1-x)*(1-x)*x;
    //Handle to EventWeight
    _evtWeight  = 1;
    art::Handle<EventWeight> weightHand;
    event.getByLabel(_weightModuleLabel, weightHand);
    if (weightHand.isValid()) _evtWeight = weightHand->weight()
				;

    //loop over the hits in the virtual detector in front of the first disk
    const StepPointMC* hit(0);
    for (size_t i=0; i<_vdHitsSize; ++i) {
      hit = &vdHits->at(i);

      int id = hit->volumeId();

      //check the location of the virtual detector: we want to study only the particles in the first disk (closer to the tracker)
      if (id == VirtualDetectorId::EMC_Disk_0_SurfIn  ||
	  //	  id == VirtualDetectorId::EMC_Disk_1_SurfIn  ||
	  id == VirtualDetectorId::EMC_Disk_0_EdgeIn//   ||
	  // id == VirtualDetectorId::EMC_Disk_1_EdgeIn    
	  ) {

	fillVDetHistograms(hit);
      }
    }

    for (unsigned i=0; i<_genSize; ++i){
      gen = &genParticles->at(i);
      fillGenHistograms(gen);
    }

    //Get calo cluster 2
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    CaloClusterCollection const& caloClusters(*caloClustersHandle);

    int                nClFast = caloClusters.size();
    const CaloCluster* cluster(0);
  
    //for loop over the clusters in the calorimeter
    for( int i=0; i<nClFast; ++i){
      cluster=&(caloClusters.at(i));
      int iSection = cluster->diskId();

      if (iSection == 1)                continue;

      fillClusterHistograms(cluster);
    }
     //now fill the likelihood distributions
    fillLHHistograms  (&caloClusters);
        
    //fill the event histograms
    fillEventHistograms(&caloClusters);
    
    return 1;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCPhotonAna);


