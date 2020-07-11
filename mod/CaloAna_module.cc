//
// An EDAnalyzer module that is used to study the RPC photon
//
// $Id: CaloAna_module.cc,v 1.1.2.1 2018/12/05 21:42:24 gianipez Exp $
// $Author: gianipez $
// $Date: 2018/12/05 21:42:24 $
//
// Original author Gianantonio Pezzullo 
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
#include "TString.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>




using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;



namespace mu2e {


  class CaloAna : public art::EDFilter {
    
    struct  vDetInfo_ {
      double                     _energy;
      double                     _x, _y, _time;
      double                     _energy_e;
      double                     _x_e, _y_e, _time_e;
      
      void   init () {
	_x            = -9999.;
	_y            = -9999.;
	_energy       = 50.;//-9999.;
	_x_e          = -9999.;
	_y_e          = -9999.;
	_energy_e     = 10.;//-9999.;
	_time         = -1e3;
	_time_e       = -1e3;
      }
      
      vDetInfo_(){
	init();
      };
    };

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
      TH2F* eRDist   [3];  // cluster energy
      TH2F* e1RDist  [3];  // most  energetic crystal energy over total cluster energy
      TH2F* e2RDist  [3];  // second most  energetic crystal energy over total cluster energy
      TH2F* e2e1RDist[3];  // secondd most  energetic crystal energy over most energetic crystal energy
      TH2F* e3RDist  [3];  // third most  energetic crystal energy over total cluster energy
      TH2F* e4RDist  [3];  // fourth most  energetic crystal energy over total cluster energy
      TH2F* nCrRDist [3];  // cluster size vs cluster cog radial distance

    };

    struct  eventHist_ {
      TH1F* nCl[3];
    };

    struct  lhHist_ {
      TH1F* hVar   [2][10];
      TH1F* hCorVar[2][10];
      TH1F* hTot   [2];
      TH1F* hCorTot[2];
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

    virtual ~CaloAna() { }

    virtual void beginJob();
    virtual bool beginRun(art::Run&   run   );
    virtual void endJob();



    // This is called for each event.
    //     virtual void analyze(const art::Event& e);

    virtual bool filter(art::Event& event);
    explicit CaloAna(const fhicl::ParameterSet& PSet);

    void     bookHistograms ();
    void     bookGenHist    (art::ServiceHandle<art::TFileService> & Tfs);
    void     bookEvtHist    (art::ServiceHandle<art::TFileService> & Tfs);
    void     bookVDetHist   (art::ServiceHandle<art::TFileService> & Tfs, particleHist_   &Hist, int Id);
    void     bookClusterHist(art::ServiceHandle<art::TFileService> & Tfs, clusterHist_    &Hist, int Id, int NameId );
    void     bookLHHist     (art::ServiceHandle<art::TFileService> & Tfs);

    void     fillGenHistograms    (const GenParticle*           Gen    );
    void     fillEventHistograms  (const CaloClusterCollection* ClCol  );
    void     fillClusterHistograms(const CaloCluster*           Cluster, clusterHist_    &Hist);
    void     fillVDetHistograms   (const StepPointMC*           Step   , particleHist_   &Hist);
    void     fillLHHistogrgams    (const CaloClusterCollection* ClCol  );
  private:
       
    int                        _diagLevel;
    int                        _nProcess;
    
    double                     _evtWeight;
    vDetInfo_                  _vdet;
    // double                     _vdetE;
    // double                     _vdetX, _vdetY, _vdetT;
    // double                     _vdetE_e;
    // double                     _vdetX_e, _vdetY_e, _vdetT_e;

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
    //    double                     _minClEnergy, _clEStep;
    double                     _minRDist   , _rDistStep;
 
    //copied these from private section of PID
    // std::string _TemplateFile;
   
   
    //histograms collectors
    particleHist_   _hGen;
    particleHist_   _hVDetDisk0, _hVDetDisk1;  
    clusterHist_    _hClDisk0  , _hClDisk1;
    eventHist_      _hEvt;
    lhHist_         _hLh;
    GlobalConstantsHandle<ParticleDataTable> _pdt;
    std::vector<TString>                     _clNames;

    TH1F*       _signalHist1D[2][kN1DVar];
    TH2F*       _signalHist2D[2][kN2DVar];

    TH1F*       _bkgHist1D   [2][kN1DVar];
    TH2F*       _bkgHist2D   [2][kN2DVar];

    TH1F*       _signalCorrHist1D[2][kN2DVar][kNCorHist];
    TH1F*       _bkgCorrHist1D   [2][kN2DVar][kNCorHist];
   
    void        buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize);
    
    //initialize calclate prob in private
    double      calculateProb  (double &variable, TH1F* templates);
    double      calculate2DProb(double& Ref, double&Variable, TH1F** Template, double MinX, double Step);

 
    const Calorimeter*                    _calorimeter ;
  };


  CaloAna::CaloAna(const fhicl::ParameterSet & pset) :
    art::EDFilter{pset},
    _diagLevel                   (pset.get<int>("diagLevel",0)),
    _nProcess                    (0),
    _cogType                     (ClusterMoments::Linear),                
    _dropSecondDisk              (pset.get<bool>  ("dropSecondDisk", false)),
    _generatorModuleLabel        (pset.get<string>("generatorModuleLabel")),
    _weightModuleLabel           (pset.get<string>("weightModuleLabel")),
    _caloClusterModuleLabel      (pset.get<string>("caloClusterModuleLabel")),
    _virtualDetectorLabel        (pset.get<string>("virtualDetectorName")),
    _toff                        (pset.get<fhicl::ParameterSet>("timeOffsets", fhicl::ParameterSet())),
    _mbbuffer                    (pset.get<double>("timeFoldingBuffer")),  // ns
    _blindTime                   (pset.get<double>("blindTime" )),         //ns
    _signalTemplateFile          (pset.get<string>("signalTemplates")),
    _bkgTemplateFile             (pset.get<string>("backgroundTemplates")),
    // _minClEnergy                 (pset.get<double>("minClusterEnergy"     ,   50.)),   // MeV
    // _clEStep                     (pset.get<double>("clusterEnergyStep"    ,   10.)),   // MeV
    _minRDist                    (pset.get<double>("minClusterRadialDist" ,  350.)),   // mm
    _rDistStep                   (pset.get<double>("clusterRadialDistStep",   50.)){   // mm

    ConfigFileLookupPolicy configFile;
    _signalTemplateFile = configFile(_signalTemplateFile);
    _bkgTemplateFile    = configFile(_bkgTemplateFile);
   
    _clNames.push_back("all"     );
    _clNames.push_back("photon"  );
    _clNames.push_back("electron");


    TFile* signalFile = TFile::Open(_signalTemplateFile.c_str());
    TFile* bkgFile    = TFile::Open(_bkgTemplateFile.c_str());

        //get the templates histograms
    //signal on calo-disk 0
    _signalHist1D[0][0]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_E2");	
    _signalHist1D[0][1]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e1eRatio2");   
    _signalHist1D[0][2]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_nCr2");   	
    _signalHist1D[0][3]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_rDist2");	
    _signalHist1D[0][4]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2eRatio2");	
    _signalHist1D[0][5]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2e1Ratio2");
    //signal on calo-disk 1
    _signalHist1D[1][0]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_E2");
    _signalHist1D[1][1]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e1eRatio2");   
    _signalHist1D[1][2]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_nCr2");   
    _signalHist1D[1][3]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_rDist2");
    _signalHist1D[1][4]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2eRatio2");
    _signalHist1D[1][5]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2e1Ratio2");

    //background on calo-disk 0
    _bkgHist1D   [0][0]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_E0");	
    _bkgHist1D   [0][1]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e1eRatio0");
    _bkgHist1D   [0][2]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_nCr0");   
    _bkgHist1D   [0][3]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_rDist0");
    _bkgHist1D   [0][4]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2eRatio0");
    _bkgHist1D   [0][5]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2e1Ratio0");
    //background on calo-disk 1
    _bkgHist1D   [1][0]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_E0");	
    _bkgHist1D   [1][1]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e1eRatio0");
    _bkgHist1D   [1][2]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_nCr0");   
    _bkgHist1D   [1][3]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_rDist0");
    _bkgHist1D   [1][4]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2eRatio0");
    _bkgHist1D   [1][5]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2e1Ratio0");
    

    //get the 2D histograms
    // signal
    //disk 0
    _signalHist2D[0][0]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_eRDist2"); 
    _signalHist2D[0][1]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e1RDist2");
    _signalHist2D[0][2]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2RDist2");
    _signalHist2D[0][3]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2e1RDist2");
    _signalHist2D[0][4]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e3RDist2");
    _signalHist2D[0][5]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e4RDist2");
    _signalHist2D[0][6]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_nCrRDist2");

    //disk 1
    _signalHist2D[1][0]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_eRDist2");
    _signalHist2D[1][1]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e1RDist2");
    _signalHist2D[1][2]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2RDist2");
    _signalHist2D[1][3]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2e1RDist2");
    _signalHist2D[1][4]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e3RDist2");
    _signalHist2D[1][5]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e4RDist2");
    _signalHist2D[1][6]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_nCrRDist2");

    //background
    //disk 0
    _bkgHist2D   [0][0]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_eRDist0"); 
    _bkgHist2D   [0][1]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e1RDist0");
    _bkgHist2D   [0][2]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2RDist0");
    _bkgHist2D   [0][3]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2e1RDist0");
    _bkgHist2D   [0][4]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e3RDist0");
    _bkgHist2D   [0][5]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e4RDist0");
    _bkgHist2D   [0][6]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_nCrRDist0");

    //disk 1
    _bkgHist2D   [1][0]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_eRDist0"); 
    _bkgHist2D   [1][1]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e1RDist0");
    _bkgHist2D   [1][2]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2RDist0");
    _bkgHist2D   [1][3]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2e1RDist0");
    _bkgHist2D   [1][4]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e3RDist0");
    _bkgHist2D   [1][5]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e4RDist0");
    _bkgHist2D   [1][6]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_nCrRDist0");

    //make correlation templates
    int     nCaloDisks(2);
    for (int i=0; i<nCaloDisks; ++i){    
      //cluster energy vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][0], _signalCorrHist1D[i][0], Form("SignalDisk%i_eRDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][0], _bkgCorrHist1D   [i][0], Form("BkgDisk%i_eRDist"   , i), _minRDist, _rDistStep);
      
      //E1 (seedHitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][1], _signalCorrHist1D[i][1], Form("SignalDisk%i_e1RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][1], _bkgCorrHist1D   [i][1], Form("BkgDisk%i_e1RDist"   , i), _minRDist, _rDistStep);
      
      //E2 (2nd_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][2], _signalCorrHist1D[i][2], Form("SignalDisk%i_e2RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][2], _bkgCorrHist1D   [i][2], Form("BkgDisk%i_e2RDist"   , i), _minRDist, _rDistStep);

      //E2/E1 vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][3], _signalCorrHist1D[i][3], Form("SignalDisk%i_e2e1RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][3], _bkgCorrHist1D   [i][3], Form("BkgDisk%i_e2e1RDist"   , i), _minRDist, _rDistStep);

      //E3 (3rd_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][4], _signalCorrHist1D[i][4], Form("SignalDisk%i_e3RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][4], _bkgCorrHist1D   [i][4], Form("BkgDisk%i_e3RDist"   , i), _minRDist, _rDistStep);

      //E4 (4th_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][5], _signalCorrHist1D[i][5], Form("SignalDisk%i_e4RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][5], _bkgCorrHist1D   [i][5], Form("BkgDisk%i_e4RDist"   , i), _minRDist, _rDistStep);

      //nCrystals vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][6], _signalCorrHist1D[i][6], Form("SignalDisk%i_nCrRDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][6], _bkgCorrHist1D   [i][6], Form("BkgDisk%i_nCrRDist"   , i), _minRDist, _rDistStep);

    }
  }
  
  //--------------------------------------------------------------------------------
  // routine function used to produce the template histograms 
  // from the 2D distributions
  //--------------------------------------------------------------------------------
  void   CaloAna::buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize){
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
  void CaloAna::bookHistograms(){
    
    art::ServiceHandle<art::TFileService> tfs;

    bookGenHist(tfs);

    bookVDetHist(tfs, _hVDetDisk0, 0);

    bookVDetHist(tfs, _hVDetDisk1, 1);
    
    //all clusters
    bookClusterHist(tfs, _hClDisk0, 0, 0);
    bookClusterHist(tfs, _hClDisk1, 1, 0);

    //clusters from  photons
    bookClusterHist(tfs, _hClDisk0, 0, 1);
    bookClusterHist(tfs, _hClDisk1, 1, 1);

    //clusters from electrons
    bookClusterHist(tfs, _hClDisk0, 0, 2);
    bookClusterHist(tfs, _hClDisk1, 1, 2);

    bookLHHist(tfs);

    bookEvtHist(tfs);
  }

  void     CaloAna::bookLHHist(art::ServiceHandle<art::TFileService> & Tfs){
    art::TFileDirectory lhDir = Tfs->mkdir("likelihood");
    _hLh.hVar    [0][0]  = lhDir.make<TH1F>("lh_disk0clEnergy"    ,"lh disk0 from cluster energy; log(P_{signal}(E)/P_{bkg}(E))"           ,  400, -100, 100);
    _hLh.hVar    [0][1]  = lhDir.make<TH1F>("lh_disk0clE1"        ,"lh disk0 from E1=E_{seed}/E_{cluster}; log(P_{signal}(E1)/P_{bkg}(E1))",  400, -100, 100);
    _hLh.hVar    [0][2]  = lhDir.make<TH1F>("lh_disk0clSize"      ,"lh disk0 from cluster size; log(P_{signal}(Size)/P_{bkg}(Size))"       ,  400, -100, 100);
    _hLh.hVar    [0][3]  = lhDir.make<TH1F>("lh_disk0rDist"       ,"lh disk0 from E1=E_{seed}/E_{cluster}; log(P_{signal}(E1)/P_{bkg}(E1))",  400, -100, 100);
    _hLh.hVar    [0][4]  = lhDir.make<TH1F>("lh_disk0clE2"        ,"lh disk0 from E2=E_{2nd/E_{cluster}; log(P_{signal}(E2)/P_{bkg}(E2))"  ,  400, -100, 100);
    _hLh.hVar    [0][5]  = lhDir.make<TH1F>("lh_disk0clE2E1"      ,"lh disk0 from E2/E1; log(P_{signal}(E2/E1)/P_{bkg}(E2/E1))"            ,  400, -100, 100);

    _hLh.hVar    [1][0]  = lhDir.make<TH1F>("lh_disk1clEnergy"    ,"lh disk1 from cluster energy; log(P_{signal}(E)/P_{bkg}(E))"           ,  400, -100, 100);
    _hLh.hVar    [1][1]  = lhDir.make<TH1F>("lh_disk1clE1"        ,"lh disk1 from E1=E_{seed}/E_{cluster}; log(P_{signal}(E1)/P_{bkg}(E1))",  400, -100, 100);
    _hLh.hVar    [1][2]  = lhDir.make<TH1F>("lh_disk1clSize"      ,"lh disk1 from cluster size; log(P_{signal}(Size)/P_{bkg}(Size))"       ,  400, -100, 100);
    _hLh.hVar    [1][3]  = lhDir.make<TH1F>("lh_disk1rDist"       ,"lh disk1 from E1=E_{seed}/E_{cluster}; log(P_{signal}(E1)/P_{bkg}(E1))",  400, -100, 100);
    _hLh.hVar    [1][4]  = lhDir.make<TH1F>("lh_disk1clE2"        ,"lh disk1 from E2=E_{2nd/E_{cluster}; log(P_{signal}(E2)/P_{bkg}(E2))"  ,  400, -100, 100);
    _hLh.hVar    [1][5]  = lhDir.make<TH1F>("lh_disk1clE2E1"      ,"lh disk1 from E2/E1; log(P_{signal}(E2/E1)/P_{bkg}(E2/E1))"            ,  400, -100, 100);

    _hLh.hCorVar [0][0]  = lhDir.make<TH1F>("lhCor_disk0clEnergy" ,"lh disk0 from cluster energy vs rdist; log(P_{signal}(E)/P_{bkg}(E))"           ,  400, -100, 100);
    _hLh.hCorVar [0][1]  = lhDir.make<TH1F>("lhCor_disk0clE1"     ,"lh disk0 from rdist vs E1; log(P_{signal}(E1)/P_{bkg}(E1))",  4000, -1000, 1000);
    //    _hLh.hVar  [2]  = lhDir.make<TH1F>("lh_clSize"   ,"lh disk0 from cluster size; log(P_{signal}(Size)/P_{bkg}(Size))"       ,  4000, -1000, 1000);
    //    _hLh.hVar  [3]  = lhDir.make<TH1F>("lh_rDist"    ,"lh disk0 from E1=E_{seed}/E_{cluster}; log(P_{signal}(E1)/P_{bkg}(E1))",  4000, -1000, 1000);
    _hLh.hCorVar [0][2]  = lhDir.make<TH1F>("lhCor_disk0clE2"     ,"lh disk0 from rdist vs E2; log(P_{signal}(E2)/P_{bkg}(E2))"  ,  4000, -1000, 1000);
    _hLh.hCorVar [0][3]  = lhDir.make<TH1F>("lhCor_disk0clE2E1"   ,"lh disk0 from rdist vsE2/E1; log(P_{signal}(E2/E1)/P_{bkg}(E2/E1))"            ,  4000, -1000, 1000);
    _hLh.hCorVar [0][4]  = lhDir.make<TH1F>("lhCor_disk0clE3"     ,"lh disk0 from rdist vs E3; log(P_{signal}(E3)/P_{bkg}(E3))"            ,  4000, -1000, 1000);
    _hLh.hCorVar [0][5]  = lhDir.make<TH1F>("lhCor_disk0clE4"     ,"lh disk0 from rdist vs E4; log(P_{signal}(E4)/P_{bkg}(E4))"            ,  4000, -1000, 1000);
    _hLh.hCorVar [0][6]  = lhDir.make<TH1F>("lhCor_disk0clNCr"    ,"lh disk0 from rdist vs nCrystal; log(P_{signal}(nCr)/P_{bkg}(nCr))"            ,  4000, -1000, 1000);
  
    _hLh.hCorVar [1][0]  = lhDir.make<TH1F>("lhCor_disk1clEnergy" ,"lh disk1 from cluster energy vs rdist; log(P_{signal}(E)/P_{bkg}(E))"           ,  400, -100, 100);
    _hLh.hCorVar [1][1]  = lhDir.make<TH1F>("lhCor_disk1clE1"     ,"lh disk1 from rdist vs E1; log(P_{signal}(E1)/P_{bkg}(E1))",  4000, -1000, 1000);
    //    _hLh.hVar  [2]  = lhDir.make<TH1F>("lh_clSize"   ,"lh disk1 from cluster size; log(P_{signal}(Size)/P_{bkg}(Size))"       ,  4000, -1000, 1000);
    //    _hLh.hVar  [3]  = lhDir.make<TH1F>("lh_rDist"    ,"lh disk1 from E1=E_{seed}/E_{cluster}; log(P_{signal}(E1)/P_{bkg}(E1))",  4000, -1000, 1000);
    _hLh.hCorVar [1][2]  = lhDir.make<TH1F>("lhCor_disk1clE2"     ,"lh disk1 from rdist vs E2; log(P_{signal}(E2)/P_{bkg}(E2))"  ,  4000, -1000, 1000);
    _hLh.hCorVar [1][3]  = lhDir.make<TH1F>("lhCor_disk1clE2E1"   ,"lh disk1 from rdist vs E2/E1; log(P_{signal}(E2/E1)/P_{bkg}(E2/E1))"            ,  4000, -1000, 1000);
    _hLh.hCorVar [1][4]  = lhDir.make<TH1F>("lhCor_disk1clE3"     ,"lh disk1 from rdist vs E3; log(P_{signal}(E3)/P_{bkg}(E3))"            ,  4000, -1000, 1000);
    _hLh.hCorVar [1][5]  = lhDir.make<TH1F>("lhCor_disk1clE4"     ,"lh disk1 from rdist vs E4; log(P_{signal}(E4)/P_{bkg}(E4))"            ,  4000, -1000, 1000);
    _hLh.hCorVar [1][6]  = lhDir.make<TH1F>("lhCor_disk1clNCr"    ,"lh disk1 from rdist vs nCrystal; log(P_{signal}(nCr)/P_{bkg}(nCr))"            ,  4000, -1000, 1000);
    //    _hLh.hVar  [5]  = lhDir.make<TH1F>("lh_clE2E1"      ,"likelihood from E2/E1; log(P_{signal}(E2/E1)/P_{bkg}(E2/E1))"            ,  400, -100, 100);


    _hLh.hTot    [0]     = lhDir.make<TH1F>("lh_disk0tot" ,"Global likelihood disk 0; #sum log(P_{signal}(x_{i})/P_{bkg}(x_{i}))"  , 400, -200, 200);
    _hLh.hTot    [1]     = lhDir.make<TH1F>("lh_disk1tot" ,"Global likelihood disk 1; #sum log(P_{signal}(x_{i})/P_{bkg}(x_{i}))"  , 400, -200, 200);

    _hLh.hCorTot    [0]     = lhDir.make<TH1F>("lhCor_disk0tot" ,"Global likelihood disk 0; #sum log(P_{signal}(x_{i})/P_{bkg}(x_{i}))"  , 4000, -1000, 1000);
    _hLh.hCorTot    [1]     = lhDir.make<TH1F>("lhCor_disk1tot" ,"Global likelihood disk 1; #sum log(P_{signal}(x_{i})/P_{bkg}(x_{i}))"  , 4000, -1000, 1000);
  }

  //--------------------------------------------------------------------------------
  // Book  histograms associated with the generation stage
  //--------------------------------------------------------------------------------
  void     CaloAna::bookGenHist    (art::ServiceHandle<art::TFileService> & Tfs){
    art::TFileDirectory genDir = Tfs->mkdir("gen");
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

    
  }
  
  
 //--------------------------------------------------------------------------------
  // Book  histograms associated with the generation stage
  //--------------------------------------------------------------------------------
  void     CaloAna::bookEvtHist    (art::ServiceHandle<art::TFileService> & Tfs){
    art::TFileDirectory evtDir = Tfs->mkdir("evt");
    _hEvt.nCl[0]       = evtDir.make<TH1F>("evt_nCl0"    ,"n-cluster distrbution; nClusters",  201, -0.5, 200.5);
    _hEvt.nCl[1]       = evtDir.make<TH1F>("evt_nCl1"    ,"n-cluster distrbution photons with E>70 MeV; nClusters",  201, -0.5, 200.5);
    _hEvt.nCl[2]       = evtDir.make<TH1F>("evt_nCl2"    ,"n-cluster distrbution electrons with E>70 MeV; nClusters",  201, -0.5, 200.5);

  }
  //--------------------------------------------------------------------------------
  // Book  histograms associated with the generation stage
  //--------------------------------------------------------------------------------
  void     CaloAna::bookVDetHist   (art::ServiceHandle<art::TFileService> & Tfs, particleHist_   &Hist, int Id){
    art::TFileDirectory vdetDir = Tfs->mkdir("vdet");
    Hist.energy[0] = vdetDir.make<TH1F>(Form("vdetDisk%i_E0"    , Id),"Energy distrbution; E [MeV]",  400, 0, 200);
    Hist.time  [0] = vdetDir.make<TH1F>(Form("vdetDisk%i_time0" , Id),"Time distribution; t [ns]"  , 1700, 0, 1700);
    Hist.pdgId [0] = vdetDir.make<TH1F>(Form("vdetDisk%i_pdgId0", Id),"PdgID; pdgId"               , 6000, -3000, 3000);
    Hist.cosTh [0] = vdetDir.make<TH1F>(Form("vdetDisk%i_costh0", Id),"cos(#theta); cos(#theta)"   , 2000, -1, 1);
    Hist.energy[1] = vdetDir.make<TH1F>(Form("vdetDisk%i_E1"    , Id),"Energy distrbution; E [MeV]",  400, 0, 200);
    Hist.time  [1] = vdetDir.make<TH1F>(Form("vdetDisk%i_time1" , Id),"Time distribution; t [ns]"  , 1700, 0, 1700);
    Hist.pdgId [1] = vdetDir.make<TH1F>(Form("vdetDisk%i_pdgId1", Id),"PdgID; pdgId"               , 6000, -3000, 3000);
    Hist.cosTh [1] = vdetDir.make<TH1F>(Form("vdetDisk%i_costh1", Id),"cos(#theta); cos(#theta)"   , 2000, -1, 1);
    Hist.energy[2] = vdetDir.make<TH1F>(Form("vdetDisk%i_E2"    , Id),"Energy distrbution; E [MeV]",  400, 0, 200);
    Hist.time  [2] = vdetDir.make<TH1F>(Form("vdetDisk%i_time2" , Id),"Time distribution; t [ns]"  , 1700, 0, 1700);
    Hist.pdgId [2] = vdetDir.make<TH1F>(Form("vdetDisk%i_pdgId2", Id),"PdgID; pdgId"               , 6000, -3000, 3000);
    Hist.cosTh [2] = vdetDir.make<TH1F>(Form("vdetDisk%i_costh2", Id),"cos(#theta); cos(#theta)"   , 2000, -1, 1);

  }
  //--------------------------------------------------------------------------------
  // Book  histograms associated with the generation stage
  //--------------------------------------------------------------------------------
  void     CaloAna::bookClusterHist(art::ServiceHandle<art::TFileService> & Tfs, clusterHist_    &Hist, int Id, int NameId){
    
    TString             name  = _clNames.at(NameId);
    art::TFileDirectory clDir = Tfs->mkdir(Form("cluster_%s", name.Data()));

    Hist.energy   [NameId]     = clDir.make<TH1F>(Form("clDisk%i_E%i"        , Id, NameId),"Energy distrbution; E [MeV]",  400, 0, 200);
    Hist.time     [NameId]     = clDir.make<TH1F>(Form("clDisk%i_time%i"     , Id, NameId),"Time distribution; t [ns]"  , 1700, 0, 1700);
    Hist.dt       [NameId]     = clDir.make<TH1F>(Form("clDisk%i_dt%i"       , Id, NameId),"dT=t_{calo}-t_{vdet} [ns] distribution; #Delta t = t_{calo}-t_{vdet} [ns]"  , 400, -100, 100);
    Hist.nCr      [NameId]     = clDir.make<TH1F>(Form("clDisk%i_nCr%i"      , Id, NameId),"Cluster size distribution; nCrystalHits"  , 21, -0.5, 20.5);
    Hist.rDist    [NameId]     = clDir.make<TH1F>(Form("clDisk%i_rDist%i"    , Id, NameId),"Time distribution; cluster radial distance [mm]"  , 500, 300, 800);
    Hist.e1eRatio [NameId]     = clDir.make<TH1F>(Form("clDisk%i_e1eRatio%i" , Id, NameId),"most crystalHit energy/E_{cluster}; E^{max}_{hit}/E_{cluster} "  , 100, 0, 1.);
    Hist.e2eRatio [NameId]     = clDir.make<TH1F>(Form("clDisk%i_e2eRatio%i" , Id, NameId),"2^{nd} most crystalHit energy/E_{cluster}; E^{2nd-max}_{hit}/E_{cluster}"   , 100, 0, 1.2);
    Hist.e2e1Ratio[NameId]     = clDir.make<TH1F>(Form("clDisk%i_e2e1Ratio%i", Id, NameId),"2^{nd} most crystalHit energy/E^{mx}_{hit}; E^{2nd-max}_{hit}/E^{max}_{hit}", 100, 0, 1.2);
    Hist.e3eRatio [NameId]     = clDir.make<TH1F>(Form("clDisk%i_e3eRatio%i" , Id, NameId),"3rd crystalHit energy/E_{cluster}; E^{3rd}_{hit}/E_{cluster} "  , 100, 0, 1.2);
    Hist.e4eRatio [NameId]     = clDir.make<TH1F>(Form("clDisk%i_e4eRatio%i" , Id, NameId),"4th crystalHit energy/E_{cluster}; E^{4th}_{hit}/E_{cluster} "  , 100, 0, 1.2);

    Hist.eRDist   [NameId]     = clDir.make<TH2F>(Form("clDisk%i_eRDist%i"   , Id, NameId),"r_{cog} vs E_{cluster}; r_{cog} [mm]; E_{cluster} [MeV]"        , 500, 300, 800, 400,   0,   200);
    Hist.e1RDist  [NameId]     = clDir.make<TH2F>(Form("clDisk%i_e1RDist%i"  , Id, NameId),"r_{cog} vs E1; r_{cog} [mm]; E^{max}_{hit}/E_{cluster} "        , 500, 300, 800, 100,   0,   1.2);
    Hist.e2RDist  [NameId]     = clDir.make<TH2F>(Form("clDisk%i_e2RDist%i"  , Id, NameId),"r_{cog} vs E2; r_{cog} [mm]; E^{2nd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    Hist.e2e1RDist[NameId]     = clDir.make<TH2F>(Form("clDisk%i_e2e1RDist%i", Id, NameId),"r_{cog} vs E2/E1; r_{cog} [mm]; E^{2nd-max}_{hit}/E^{max}_{hit}", 500, 300, 800, 100,   0,   1.2);
    Hist.e3RDist  [NameId]     = clDir.make<TH2F>(Form("clDisk%i_e3RDist%i"  , Id, NameId),"r_{cog} vs E3; r_{cog} [mm]; E^{3rd-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    Hist.e4RDist  [NameId]     = clDir.make<TH2F>(Form("clDisk%i_e4RDist%i"  , Id, NameId),"r_{cog} vs E4; r_{cog} [mm]; E^{4th-max}_{hit}/E_{cluster}"     , 500, 300, 800, 100,   0,   1.2);
    Hist.nCrRDist [NameId]     = clDir.make<TH2F>(Form("clDisk%i_nCrRDist%i" , Id, NameId),"r_{cog} vs nCrystals; r_{cog} [mm]; nCrystals"                  , 500, 300, 800,  20,   0,    20);

    
  }




  //--------------------------------------------------------------------------------//
  void     CaloAna::beginJob(){

    bookHistograms();

  }

  bool     CaloAna::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    return true;
  }

  void     CaloAna::endJob(){}

  //--------------------------------------------------------------------------------
  double   CaloAna::calculate2DProb(double& Ref, double&Variable, TH1F** Template, double MinX, double Step){
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
  double   CaloAna::calculateProb(double&Variable, TH1F* Template){

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

  //------------------------------------------------------------------------------------------
  void     CaloAna::fillGenHistograms    (const GenParticle*           Gen    ){
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
  void     CaloAna::fillEventHistograms  (const CaloClusterCollection* ClCol  ){
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

  //--------------------------------------------------------------------------------
  void     CaloAna::fillLHHistogrgams    (const CaloClusterCollection* ClCol  ){
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
      ClusterMoments cogCalculator(*_calorimeter, *cluster,clSection);
      cogCalculator.calculate(_cogType);

      double                   clEnergy = cluster->energyDep();
      const CLHEP::Hep3Vector& cog      = cluster->cog3Vector();
      double                   xpos     = cog(0);
      double                   ypos     = cog(1);
      double                   rDist    = sqrt(xpos*xpos+ypos*ypos);

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
      
      double  signalEnergyProb      = calculateProb(clEnergy, _signalHist1D[clSection][0]);
      double  bkgEnergyProb         = calculateProb(clEnergy, _bkgHist1D   [clSection][0]);
      double  logEnergyRatio        = log(signalEnergyProb/bkgEnergyProb);
  
      double  signalE1Prob          = calculateProb(e1Ratio, _signalHist1D[clSection][1]);
      double  bkgE1Prob             = calculateProb(e1Ratio, _bkgHist1D   [clSection][1]);
      double  logE1Ratio            = log(signalE1Prob/bkgE1Prob);
  
      double  signalClusterSizeProb = calculateProb(clusterSize, _signalHist1D[clSection][2]);
      double  bkgClusterSizeProb    = calculateProb(clusterSize , _bkgHist1D  [clSection][2]);
      double  logClusterSizeRatio   = log(signalClusterSizeProb/bkgClusterSizeProb);
      
      double  signalRDistProb       = calculateProb(rDist, _signalHist1D[clSection][3]);
      double  bkgRDistProb          = calculateProb(rDist, _bkgHist1D   [clSection][3]);
      double  logRDistRatio         = log(signalRDistProb/bkgRDistProb);

      double  signalE2Prob          = calculateProb(e2Ratio, _signalHist1D[clSection][1]);
      double  bkgE2Prob             = calculateProb(e2Ratio, _bkgHist1D   [clSection][1]);
      double  logE2Ratio            = log(signalE2Prob/bkgE2Prob);
  
      double  signalE2E1Prob        = calculateProb(e2e1Ratio, _signalHist1D[clSection][5]);
      double  bkgE2E1Prob           = calculateProb(e2e1Ratio, _bkgHist1D   [clSection][5]);
      double  logE2E1Ratio          = log(signalE2E1Prob/bkgE2E1Prob);
 
      //fill the histograms
      _hLh.hVar[clSection][0]->Fill(logEnergyRatio     );
      _hLh.hVar[clSection][1]->Fill(logE1Ratio         );
      _hLh.hVar[clSection][2]->Fill(logClusterSizeRatio);
      _hLh.hVar[clSection][3]->Fill(logRDistRatio      );
      _hLh.hVar[clSection][4]->Fill(logE2Ratio         );
      _hLh.hVar[clSection][5]->Fill(logE2E1Ratio       );

      double  lhValue = logEnergyRatio + logE1Ratio + logClusterSizeRatio + logRDistRatio;
      if (nCrystalHits>=2) lhValue += logE2Ratio + logE2E1Ratio;
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
      double  signalEnergy2DProb      = calculate2DProb(rDist, clEnergy, _signalCorrHist1D[clSection][0], _minRDist, _rDistStep);
      double  bkgEnergy2DProb         = calculate2DProb(rDist, clEnergy, _bkgCorrHist1D   [clSection][0], _minRDist, _rDistStep);
      double  logEnergy2DRatio        = log(signalEnergy2DProb/bkgEnergy2DProb);

      double  signalE12DProb          = calculate2DProb(rDist, e1Ratio, _signalCorrHist1D[clSection][1], _minRDist, _rDistStep);
      double  bkgE12DProb             = calculate2DProb(rDist, e1Ratio, _bkgCorrHist1D   [clSection][1], _minRDist, _rDistStep);
      double  logE12DRatio            = log(signalE12DProb/bkgE12DProb);
  
      double  signalE22DProb          = calculate2DProb(rDist, e2Ratio, _signalCorrHist1D[clSection][2], _minRDist, _rDistStep);
      double  bkgE22DProb             = calculate2DProb(rDist, e2Ratio, _bkgCorrHist1D   [clSection][2], _minRDist, _rDistStep);
      double  logE22DRatio            = log(signalE22DProb/bkgE22DProb);
  
      double  signalE2E12DProb        = calculate2DProb(rDist, e2e1Ratio, _signalCorrHist1D[clSection][3], _minRDist, _rDistStep);
      double  bkgE2E12DProb           = calculate2DProb(rDist, e2e1Ratio, _bkgCorrHist1D   [clSection][3], _minRDist, _rDistStep);
      double  logE2E12DRatio          = log(signalE2E12DProb/bkgE2E12DProb);
 
      double  signalE32DProb          = calculate2DProb(rDist, e3Ratio, _signalCorrHist1D[clSection][4], _minRDist, _rDistStep);
      double  bkgE32DProb             = calculate2DProb(rDist, e3Ratio, _bkgCorrHist1D   [clSection][4], _minRDist, _rDistStep);
      double  logE32DRatio            = log(signalE32DProb/bkgE32DProb);
  
      double  signalE42DProb          = calculate2DProb(rDist, e4Ratio, _signalCorrHist1D[clSection][5], _minRDist, _rDistStep);
      double  bkgE42DProb             = calculate2DProb(rDist, e4Ratio, _bkgCorrHist1D   [clSection][5], _minRDist, _rDistStep);
      double  logE42DRatio            = log(signalE42DProb/bkgE42DProb);
  
      double  signalClSize2DProb      = calculate2DProb(rDist, clusterSize, _signalCorrHist1D[clSection][6], _minRDist, _rDistStep);
      double  bkgClSize2DProb         = calculate2DProb(rDist, clusterSize , _bkgCorrHist1D  [clSection][6], _minRDist, _rDistStep);
      double  logClSize2DRatio        = log(signalClSize2DProb/bkgClSize2DProb);

      //fill the histograms
      _hLh.hCorVar[clSection][0]->Fill(logEnergy2DRatio     );
      _hLh.hCorVar[clSection][1]->Fill(logE12DRatio         );
      // _hLh.hCorVar[clSection][2]->Fill(logClusterSize2DRatio);
      // _hLh.hCorVar[clSection][3]->Fill(logRDist2DRatio      );
      _hLh.hCorVar[clSection][2]->Fill(logE22DRatio         );
      _hLh.hCorVar[clSection][3]->Fill(logE2E12DRatio       );
      _hLh.hCorVar[clSection][4]->Fill(logE32DRatio         );
      _hLh.hCorVar[clSection][5]->Fill(logE42DRatio         );
      _hLh.hCorVar[clSection][6]->Fill(logClSize2DRatio     );
      

      double  lh2DValue = logRDistRatio + logEnergy2DRatio + logE12DRatio + logRDistRatio + logClSize2DRatio;
      if (nCrystalHits>=2) lh2DValue += logE22DRatio + logE2E12DRatio;
      if (nCrystalHits>=3) lh2DValue += logE32DRatio;
      if (nCrystalHits>=4) lh2DValue += logE42DRatio;
      
      if (lh2DValue > maxLh2DValue) {
	maxLh2DValue      = lh2DValue;
	sectionMax2DValue = clSection;
      }
      
    }
    
    if (sectionMaxValue>=0)    _hLh.hTot[sectionMaxValue] ->Fill(maxLhValue);

    if (sectionMax2DValue>=0)    _hLh.hCorTot[sectionMax2DValue] ->Fill(maxLh2DValue);

  }

  //--------------------------------------------------------------------------------//
  void     CaloAna::fillClusterHistograms(const CaloCluster*           Cluster, clusterHist_   &Hist ){
    
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
    double                  secondmaxECrystal(0), thirdmaxECrystal(0), fourthmaxECrystal(0);
    double                  indexMaxECrystal(0), index2ndMaxECrystal(-1), index3rdMaxECrystal(-1);

    double                  toll(6.);//ns

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
      if (j == indexMaxECrystal)              continue;
      crystalHit = &(*caloClusterHits.at(j));

      double crystalEnergy = crystalHit->energyDep();

      if (crystalEnergy > secondmaxECrystal){
	secondmaxECrystal   = crystalEnergy;
	index2ndMaxECrystal = j;
      }
    }
    
    for (int j=0; j<nCrystalHits; ++j){
      if (j == indexMaxECrystal || 
	  j == index2ndMaxECrystal  )         continue;
      crystalHit = &(*caloClusterHits.at(j));

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
      crystalHit = &(*caloClusterHits.at(j));

      double crystalEnergy = crystalHit->energyDep();

      if (crystalEnergy > fourthmaxECrystal){
	fourthmaxECrystal    = crystalEnergy;
      }
    }
    
    double     dt = clTime - _vdet._time_e;
 
    //fill the cluster info
    Hist.energy   [0]->Fill(clEnergy, _evtWeight);
    Hist.time     [0]->Fill(clTime, _evtWeight);
    Hist.dt       [0]->Fill(dt, _evtWeight);      
    Hist.nCr      [0]->Fill(nCrystalHits, _evtWeight);
    Hist.rDist    [0]->Fill(clRDist, _evtWeight);
    Hist.nCrRDist [0]->Fill(clRDist, nCrystalHits, _evtWeight);

    if (clEnergy    > 1e-3) {
      Hist.e1eRatio [0]->Fill(maxECrystal/clEnergy, _evtWeight);
      if (nCrystalHits>=2) {
	Hist.e2eRatio [0]->Fill(secondmaxECrystal/clEnergy, _evtWeight);
	Hist.e2RDist  [0]->Fill( clRDist, secondmaxECrystal/clEnergy, _evtWeight);

	if (maxECrystal > 1e-3){
	  Hist.e2e1Ratio[0]->Fill(secondmaxECrystal/maxECrystal, _evtWeight);
	  Hist.e2e1RDist[0]->Fill( clRDist, secondmaxECrystal/maxECrystal, _evtWeight);
	}
	
	if (nCrystalHits>=3){
	  Hist.e3eRatio [0]->Fill( thirdmaxECrystal/clEnergy, _evtWeight);
	  Hist.e3RDist  [0]->Fill( clRDist,  thirdmaxECrystal/clEnergy, _evtWeight);
	  if (nCrystalHits>=4) {		    
	    Hist.e4eRatio [0]->Fill(fourthmaxECrystal/clEnergy, _evtWeight);
	    Hist.e4RDist  [0]->Fill( clRDist, fourthmaxECrystal/clEnergy, _evtWeight);	
	  }
	}
      }
      Hist.eRDist   [0]->Fill( clRDist,                   clEnergy, _evtWeight);
      Hist.e1RDist  [0]->Fill( clRDist,       maxECrystal/clEnergy, _evtWeight);
    }
    
    //fill the cluster info only if this cluster is associated with the generated RPC photon
    double dt_photon = clTime - _vdet._time;
    if ( ( _evtWithPhoton )           && 
	 ( fabs(dt_photon) < toll)){
      Hist.energy   [1]->Fill(clEnergy    , _evtWeight);
      Hist.time     [1]->Fill(clTime      , _evtWeight);
      Hist.dt       [1]->Fill(dt_photon   , _evtWeight);      
      Hist.nCr      [1]->Fill(nCrystalHits, _evtWeight);
      Hist.rDist    [1]->Fill(clRDist     , _evtWeight);
      Hist.nCrRDist [1]->Fill(clRDist, nCrystalHits, _evtWeight);

      if (clEnergy    > 1e-3){
	Hist.e1eRatio [1]->Fill(maxECrystal/clEnergy, _evtWeight);

	if (nCrystalHits>=2) {	
	  Hist.e2eRatio [1]->Fill(secondmaxECrystal/clEnergy, _evtWeight);
	  Hist.e2RDist  [1]->Fill( clRDist, secondmaxECrystal/clEnergy, _evtWeight);

	  if (maxECrystal > 1e-3) {
	    Hist.e2e1Ratio[1]->Fill(secondmaxECrystal/maxECrystal, _evtWeight);
	    Hist.e2e1RDist[1]->Fill( clRDist, secondmaxECrystal/maxECrystal, _evtWeight);
	  }
	  
	  if (nCrystalHits>=3) {	
	    Hist.e3eRatio [1]->Fill( thirdmaxECrystal/clEnergy, _evtWeight);
	    Hist.e3RDist  [1]->Fill( clRDist,  thirdmaxECrystal/clEnergy, _evtWeight);
	    if (nCrystalHits>=4) {		    
	      Hist.e4eRatio [1]->Fill(fourthmaxECrystal/clEnergy, _evtWeight);
	      Hist.e4RDist  [1]->Fill( clRDist, fourthmaxECrystal/clEnergy, _evtWeight);
	    }
	  }
	}
	
	Hist.eRDist   [1]->Fill( clRDist,                   clEnergy, _evtWeight);
	Hist.e1RDist  [1]->Fill( clRDist,       maxECrystal/clEnergy, _evtWeight);

      }
      
    }
    //fill the cluster info only if this cluster is associated with an electron
    if ( ( _evtWithElectron )           && 
	 ( fabs(dt) < toll)  ){
      Hist.energy   [2]->Fill(clEnergy, _evtWeight);
      Hist.time     [2]->Fill(clTime, _evtWeight);
      Hist.dt       [2]->Fill(dt, _evtWeight);      
      Hist.nCr      [2]->Fill(nCrystalHits, _evtWeight);
      Hist.rDist    [2]->Fill(clRDist, _evtWeight);
      Hist.nCrRDist [2]->Fill(clRDist, nCrystalHits, _evtWeight);

      if (clEnergy    > 1e-3) {
	Hist.e1eRatio [2]->Fill(maxECrystal/clEnergy, _evtWeight);
	
	if (nCrystalHits>=2) {	
	  Hist.e2eRatio [2]->Fill(secondmaxECrystal/clEnergy, _evtWeight);
	  Hist.e2RDist  [2]->Fill( clRDist, secondmaxECrystal/clEnergy, _evtWeight);
	  
	  if (maxECrystal > 1e-3){
	    Hist.e2e1Ratio[2]->Fill(secondmaxECrystal/maxECrystal, _evtWeight);
	    Hist.e2e1RDist[2]->Fill( clRDist, secondmaxECrystal/maxECrystal, _evtWeight);
	  }	

	  if (nCrystalHits>=3) {	
	    Hist.e3eRatio [2]->Fill( thirdmaxECrystal/clEnergy, _evtWeight);
	    Hist.e3RDist  [2]->Fill( clRDist,  thirdmaxECrystal/clEnergy, _evtWeight);
	    if (nCrystalHits>=4) {		    
	      Hist.e4RDist  [2]->Fill( clRDist, fourthmaxECrystal/clEnergy, _evtWeight);
	      Hist.e4eRatio [2]->Fill(fourthmaxECrystal/clEnergy, _evtWeight);
	    }
	  }
	}
	Hist.eRDist   [2]->Fill( clRDist,                   clEnergy, _evtWeight);
	Hist.e1RDist  [2]->Fill( clRDist,       maxECrystal/clEnergy, _evtWeight);
       }
      
    }
  }


  //--------------------------------------------------------------------------------//
  void     CaloAna::fillVDetHistograms   (const StepPointMC*           Step, particleHist_   &Hist   ){
    art::Ptr<SimParticle> const& simptr = Step->simParticle();

    //    SimParticle const& sim  = *simptr;
    
    //    if ( !sim.fromGenerator() )   return;
    
    double  hitTime = fmod(Step->time() + _toff.totalTimeOffset(simptr), _mbtime);
    double  hitP    = Step->momentum().mag();
    double  pdgId   = Step->simParticle()->pdgId();
    double  hitM    = _pdt->particle(pdgId).ref().mass();
    double  hitE    = sqrt(hitP*hitP + hitM*hitM);
    double  costh   = Step->momentum().cosTheta();
    
    //check if the particle is a photon
    //double radius = sqrt(xpos*xpos + ypos*ypos);
    // double vdcosthetahere = costheta;

    if (pdgId == 22 && (hitE > _vdet._energy)) {
      _evtWithPhoton    = true;
      _vdet._energy     = hitE;
      _vdet._x          = Step->position().x()+3904.;
      _vdet._y          = Step->position().y();
      _vdet._time       = hitTime;
      Hist.energy[1]->Fill( hitE   , _evtWeight);
      Hist.time  [1]->Fill( hitTime, _evtWeight);
      Hist.cosTh [1]->Fill( costh , _evtWeight);
    }

    if (pdgId == 11 && (hitE > _vdet._energy_e)) {
      _evtWithElectron    = true;
      _vdet._energy_e     = hitE;
      _vdet._x_e          = Step->position().x()+3904.;
      _vdet._y_e          = Step->position().y();
      _vdet._time_e       = hitTime;
      Hist.energy[2]->Fill( hitE   , _evtWeight);
      Hist.time  [2]->Fill( hitTime, _evtWeight);
      Hist.cosTh [2]->Fill( costh  , _evtWeight);
    }
    
    Hist.energy[0]->Fill( hitE   , _evtWeight);
    Hist.time  [0]->Fill( hitTime, _evtWeight);
    Hist.pdgId [0]->Fill( pdgId  , _evtWeight);
    Hist.cosTh [0]->Fill( costh  , _evtWeight);
    
    
    
  }


  //changed this from   void CaloAnalysis::analyze(const art::Event& event) {
  bool     CaloAna::filter(art::Event& event) {

    ++_nProcess;

    //initialize
    _evtWithPhoton    = false;
    _evtWithElectron  = false;
    _vdet.init();
    // _vdetX            = -9999.;
    // _vdetY            = -9999.;
    // _vdetE            = 50.;//-9999.;
    // _vdetX_e          = -9999.;
    // _vdetY_e          = -9999.;
    // _vdetE_e          = 10.;//-9999.;
    // _vdetT            = -1e3;
    // _vdetT_e          = -1e3;

    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloAna =  "<<_nProcess <<std::endl;

    //get the weight of the event    
    _evtWeight  = 1.;
    art::Handle<EventWeight> evtWeightH;
    event.getByLabel(_weightModuleLabel, evtWeightH);
    if (evtWeightH.isValid()){
      _evtWeight  = evtWeightH->weight();
    }
   
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime=accPar->deBuncherPeriod;
    _toff.updateMap(event);

    //Handle to VD steps
    art::Handle<StepPointMCCollection> vdStepsHandle;
    event.getByLabel(_virtualDetectorLabel, vdStepsHandle);
    int   nVdHits(0);
    const StepPointMCCollection *vdHits(0);
    if (vdStepsHandle.isValid()){
      vdHits   = vdStepsHandle.product();
      nVdHits  = vdHits->size(); 
    }

    //get info from the generator
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel, gensHandle);
    GenParticleCollection const& genParticles(*gensHandle);
    const GenParticle*gen(0);

    //loop over the hits in the virtual detector in front of the first disk
    const StepPointMC* hit(0);
    for (int i=0; i<nVdHits; ++i) {
      hit = &vdHits->at(i);

      int id = hit->volumeId();
      bool conditionDisk0 = (id == VirtualDetectorId::EMC_Disk_0_SurfIn  ||
			     id == VirtualDetectorId::EMC_Disk_0_EdgeIn);
      bool conditionDisk1 = (id == VirtualDetectorId::EMC_Disk_1_SurfIn || 
			     id == VirtualDetectorId::EMC_Disk_1_EdgeIn);
      //check the location of the virtual detector: we want to study only the particles in the first disk (closer to the tracker)
      if (conditionDisk0) {
	fillVDetHistograms(hit, _hVDetDisk0);
      }else if (conditionDisk1) {
	fillVDetHistograms(hit, _hVDetDisk1);
      }
    }

    for (unsigned i=0; i<genParticles.size(); ++i){
      gen = &genParticles[i];
      fillGenHistograms(gen);
    }

    //Get calo cluster 2
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    const CaloClusterCollection* caloClusters = caloClustersHandle.product();

    int                nClFast = caloClusters->size();
    const CaloCluster* cluster(0);
  
    //for loop over the clusters in the calorimeter
    for( int i=0; i<nClFast; ++i){
      cluster=&(caloClusters->at(i));
      int iSection = cluster->diskId();

      if (iSection == 0) {
	fillClusterHistograms(cluster, _hClDisk0);
      }else if (iSection == 1) {
	fillClusterHistograms(cluster, _hClDisk1);
      } 
    }
    
    //now fill the likelihood distributions
    fillLHHistogrgams  (caloClusters);
        
    //fill the event histograms
    fillEventHistograms(caloClusters);

    return 1;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloAna);


