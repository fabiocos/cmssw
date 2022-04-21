#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TMath.h"
#include "TLorentzVector.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepMC/GenRanges.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Validation/MtdValidation/interface/Utils.h"
#include "Validation/MtdValidation/interface/B0KstMuMuTreeContent.h"
#include "Validation/MtdValidation/interface/B0Isolation.h"
#include "Validation/MtdValidation/interface/B0ImpactPars.h"

class MtdSecondaryPvValidation : public DQMEDAnalyzer {
public:
  explicit MtdSecondaryPvValidation(const edm::ParameterSet&);
  ~MtdSecondaryPvValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const bool mvaGenSel(const HepMC::GenParticle&, const float&);
  const bool mvaRecSel(const reco::TrackBase&, const reco::Vertex&, const double&, const double&);
  const bool mvaGenRecMatch(const HepMC::GenParticle&, const double&, const reco::TrackBase&);

  const edm::Ref<std::vector<TrackingParticle>>* getMatchedTP(const reco::TrackBaseRef&);
  const bool genSelBPHkstar(const HepMC::GenParticle&);
  const bool genSelBPHmu(const HepMC::GenParticle&);
  std::string getMuCat(reco::Muon const&);
  void MonteCarloStudies(const edm::Event&);
  bool isAncestor(const reco::Candidate*, const reco::Candidate*);
  bool skipOscillations(const reco::GenParticle&, edm::Handle<reco::GenParticleCollection>);
  bool genSelB0kstarmumu(const reco::GenParticle&);
  const void pidFromTrackRef(double, double, double, unsigned int&, bool&, bool&, bool&, bool&);

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;
  const float minProbHeavy_;

  bool printMsg_;

  static constexpr double etacutGEN_ = 4.;     // |eta| < 4;
  static constexpr double etacutREC_ = 3.;     // |eta| < 3;
  static constexpr double pTcut_ = 0.7;        // PT > 0.7 GeV
  static constexpr double deltaZcut_ = 0.1;    // dz separation 1 mm
  static constexpr double deltaPTcut_ = 0.05;  // dPT < 5%
  static constexpr double deltaDRcut_ = 0.03;  // DeltaR separation
  static constexpr double tol_ = 1.e-4;        // tolerance on reconstructed track time, [ns]
  static constexpr double mvaSel_ = 0.8;       // minimum MVA value for PID analysis

  static constexpr double c_cm_ns = geant_units::operators::convertMmToCm(CLHEP::c_light);  // [mm/ns] -> [cm/ns]
  static constexpr double c_inv = 1.0 / c_cm_ns;
  static constexpr double m_pi = 0.13957018;
  static constexpr double m_pi_inv2 = 1.0 / m_pi / m_pi;
  static constexpr double m_k = 0.493677;
  static constexpr double m_k_inv2 = 1.0 / m_k / m_k;
  static constexpr double m_p = 0.9382720813;
  static constexpr double m_p_inv2 = 1.0 / m_p / m_p;

  // ####################
  // # HLT-trigger cuts #
  // ####################
  static constexpr double CLMUMUVTX = 0.1;       // mu-mu Vtx CL [0.1]
  static constexpr double LSMUMUBS = 3;          //mu-mu L/sigma w/respect to BS [3.0]
  static constexpr double DCAMUMU = 0.5;         //mu-mu DCA w/respect to each other [0.5 cm]
  static constexpr double DCAMUBS = 2.0;         //mu DCA w/respect to BS [2.0 cm]
  static constexpr double COSALPHAMUMUBS = 0.9;  // mu-mu cos(alpha) w/respect to BS [0.9]
  static constexpr double MUMINPT = 4.0;         // mu min pT [4.0 GeV/c]
  static constexpr double MUMAXETA = 2.4;        // mu max eta [2.4]
  static constexpr double MINMUMUPT = 6.9;       // mu-mu min pT [6.9 GeV/c]
  static constexpr double MINMUMUINVMASS = 1.0;  // mu-mu min inv. mass [1.0 GeV/c2]
  static constexpr double MAXMUMUINVMASS = 4.8;  // mu-mu max inv. mass [4.8 GeV/c2]

  // ######################
  // # Pre-selection cuts #
  // ######################
  static constexpr double B0MASSLOWLIMIT = 4.5;  // B0 mass lower limit [4.5 GeV/c2]
  static constexpr double B0MASSUPLIMIT = 6.5;   // B0 mass upper limit [6.5 GeV/c2]
  static constexpr double CLB0VTX = 0.01;        // B0 Vtx CL [0.01]
  static constexpr double KSTMASSWINDOW = 3.0;   // K*0 (OR K*0bar) mass window sigma [3.0]
  static constexpr double HADDCASBS = 0.8;       // hadron DCA/sigma w/respect to BS [0.8] (also in HLT, now is 2)
  static constexpr double MINHADPT = .8;         // hadron min pT [0.8 GeV/c] (also in HLT)
  static constexpr double MAXB0PREMASS =
      25.;  // B0 mass upper limit  before performing the fit#     electrons = cms.InputTag("slimmedElectrons"),
  static constexpr double TRKMAXR = 110.0;  // [cm]
  static constexpr double TRKMAXZ = 280.0;  // [cm]

  static constexpr double MUVARTOLE = 0.01;   // [From 0 to 1]
  static constexpr double HADVARTOLE = 0.10;  // [From 0 to 1]

  static constexpr double PRIVTXNDOF = 4.0;
  static constexpr double PRIVTXMAXZ = 50.0;  // [cm]
  static constexpr double PRIVTXMAXR = 2.0;   // [cm]

  // #######################
  // # Truth matching cuts #
  // #######################
  static constexpr double RCUTMU = 0.004;  // [eta-phi]
  static constexpr double RCUTTRK = 0.1;   // [eta-phi] // was 0.3

  static constexpr double mumasserr = 3.5e-9;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;

  edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPToken_;

  edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
  edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken_;

  const reco::RecoToSimCollection* r2s_;
  const reco::SimToRecoCollection* s2r_;

  MonitorElement* meKstar_;
  MonitorElement* meKstarBar_;
  MonitorElement* meKstarM_;
  MonitorElement* meKstarBarM_;
  MonitorElement* meKstarT_;
  MonitorElement* meKstarBarT_;

  MonitorElement* meKaonp_;
  MonitorElement* meKaoneta_;
  MonitorElement* meKaonIP_;
  MonitorElement* meKaonpM_;
  MonitorElement* meKaonetaM_;
  MonitorElement* meKaonIPM_;
  MonitorElement* meKaonpT_;
  MonitorElement* meKaonetaT_;
  MonitorElement* meKaonIPT_;

  MonitorElement* mePionp_;
  MonitorElement* mePioneta_;
  MonitorElement* mePionIP_;
  MonitorElement* mePionpM_;
  MonitorElement* mePionetaM_;
  MonitorElement* mePionIPM_;
  MonitorElement* mePionpT_;
  MonitorElement* mePionetaT_;
  MonitorElement* mePionIPT_;

  MonitorElement* mePIDKstar_;
  MonitorElement* mePIDKstarBar_;
  MonitorElement* mePIDKstarM_;
  MonitorElement* mePIDKstarBarM_;
  MonitorElement* mePIDKstarT_;
  MonitorElement* mePIDKstarBarT_;

  MonitorElement* meCombKstar_;
  MonitorElement* meCombKstarBar_;
  MonitorElement* meCombKstarM_;
  MonitorElement* meCombKstarBarM_;
  MonitorElement* meCombKstarT_;
  MonitorElement* meCombKstarBarT_;

  MonitorElement* mePIDKaonp_;
  MonitorElement* mePIDKaoneta_;
  MonitorElement* mePIDKaonIP_;
  MonitorElement* mePIDKaonpM_;
  MonitorElement* mePIDKaonetaM_;
  MonitorElement* mePIDKaonIPM_;
  MonitorElement* mePIDKaonpT_;
  MonitorElement* mePIDKaonetaT_;
  MonitorElement* mePIDKaonIPT_;

  MonitorElement* mePIDPionp_;
  MonitorElement* mePIDPioneta_;
  MonitorElement* mePIDPionIP_;
  MonitorElement* mePIDPionpM_;
  MonitorElement* mePIDPionetaM_;
  MonitorElement* mePIDPionIPM_;
  MonitorElement* mePIDPionpT_;
  MonitorElement* mePIDPionetaT_;
  MonitorElement* mePIDPionIPT_;

  MonitorElement* meKaonpR1_;
  MonitorElement* meKaonetaR1_;
  MonitorElement* meKaonpMR1_;
  MonitorElement* meKaonetaMR1_;
  MonitorElement* meKaonpTR1_;
  MonitorElement* meKaonetaTR1_;

  MonitorElement* mePionpR1_;
  MonitorElement* mePionetaR1_;
  MonitorElement* mePionpMR1_;
  MonitorElement* mePionetaMR1_;
  MonitorElement* mePionpTR1_;
  MonitorElement* mePionetaTR1_;

  MonitorElement* meKaonpR2_;
  MonitorElement* meKaonetaR2_;
  MonitorElement* meKaonpMR2_;
  MonitorElement* meKaonetaMR2_;
  MonitorElement* meKaonpTR2_;
  MonitorElement* meKaonetaTR2_;

  MonitorElement* mePionpR2_;
  MonitorElement* mePionetaR2_;
  MonitorElement* mePionpMR2_;
  MonitorElement* mePionetaMR2_;
  MonitorElement* mePionpTR2_;
  MonitorElement* mePionetaTR2_;

  MonitorElement* meSVip_;
  MonitorElement* meBarrelPCASVdiff_;
  MonitorElement* meEndcapPCASVdiff_;

  MonitorElement* meBarrelNoPIDtype_;
  MonitorElement* meEndcapNoPIDtype_;

  MonitorElement* meBarrelRecoPiNoPID_;
  MonitorElement* meBarrelRecoKNoPID_;
  MonitorElement* meEndcapRecoPiNoPID_;
  MonitorElement* meEndcapRecoKNoPID_;

  MonitorElement* meBarrelRecoPiAsPi_;
  MonitorElement* meBarrelRecoPiAsOth_;
  MonitorElement* meEndcapRecoPiAsPi_;
  MonitorElement* meEndcapRecoPiAsOth_;

  MonitorElement* meBarrelRecoKAsOth_;
  MonitorElement* meBarrelRecoKAsK_;
  MonitorElement* meEndcapRecoKAsOth_;
  MonitorElement* meEndcapRecoKAsK_;

  B0KstMuMuTreeContent* NTuple;
  Utils* Utility;
};

// ------------ constructor and destructor --------------
MtdSecondaryPvValidation::MtdSecondaryPvValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMaxBtlEta_(iConfig.getParameter<double>("trackMaximumBtlEta")),
      trackMinEtlEta_(iConfig.getParameter<double>("trackMinimumEtlEta")),
      trackMaxEtlEta_(iConfig.getParameter<double>("trackMaximumEtlEta")),
      minProbHeavy_(iConfig.getParameter<double>("minProbHeavy")),
      printMsg_(iConfig.getUntrackedParameter<bool>("printMsg")),
      magFieldToken_(esConsumes()) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));
  muonToken_ = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  beamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  HepMCProductToken_ = consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("inputTagH"));
  prunedGenToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned"));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  trackingVertexCollectionToken_ = consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  tofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofPi"));
  tofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofK"));
  tofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofP"));
  probPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPi"));
  probKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probK"));
  probPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probP"));
  mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();
  particleTableToken_ = esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>();

  NTuple = new B0KstMuMuTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();

  Utility = new Utils();
}

MtdSecondaryPvValidation::~MtdSecondaryPvValidation() {}

// ------------ method called for each event  ------------
void MtdSecondaryPvValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));
  auto RecVertexHandle = makeValid(iEvent.getHandle(RecVertexToken_));

  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);
  const auto& tofPi = iEvent.get(tofPiToken_);
  const auto& tofK = iEvent.get(tofKToken_);
  const auto& tofP = iEvent.get(tofPToken_);
  const auto& probPi = iEvent.get(probPiToken_);
  const auto& probK = iEvent.get(probKToken_);
  const auto& probP = iEvent.get(probPToken_);

  // generator level information (HepMC format)
  auto GenEventHandle = makeValid(iEvent.getHandle(HepMCProductToken_));
  const HepMC::GenEvent* mc = GenEventHandle->GetEvent();
  //double zsim = convertMmToCm((*(mc->vertices_begin()))->position().z());
  //double tsim = (*(mc->vertices_begin()))->position().t() * CLHEP::mm / CLHEP::c_light;

  // TrackingParticle collections and association maps
  //auto tpCollectionH = makeValid(iEvent.getHandle(trackingParticleCollectionToken_));
  //const auto& tpColl = tpCollectionH.product();
  //auto tvCollectionH = makeValid(iEvent.getHandle(trackingVertexCollectionToken_));

  //auto simToRecoH = makeValid(iEvent.getHandle(simToRecoAssociationToken_));
  //s2r_ = simToRecoH.product();

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  r2s_ = recoToSimH.product();

  auto pdt = iSetup.getHandle(particleTableToken_);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  //std::vector<reco::Track> candTrkMu;
  //std::vector<reco::TrackRef> candRefMu;
  //std::vector<const HepMC::GenParticle*> candGenMu;
  //std::vector<reco::Track> candTrk;
  //std::vector<reco::TrackRef> candRef;
  //std::vector<const HepMC::GenParticle*> candGen;

  // Get magnetic field
  auto bFieldHandle = iSetup.getHandle(magFieldToken_);
  // Get BeamSpot
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken_, beamSpotH);
  reco::BeamSpot beamSpot = *beamSpotH;

  reco::TrackRef muTrackm;
  reco::TrackRef muTrackp;

  double pT;
  double eta;
  double chi;
  double ndf;

  double LSBS;
  double LSBSErr;
  double cosAlphaBS;
  double cosAlphaBSErr;
  double bMcosAlphaBS;
  double bMcosAlphaBSErr;
  double bPcosAlphaBS;
  double bPcosAlphaBSErr;
  double bMinusVtxCL;
  double bPlusVtxCL;

  TrajectoryStateClosestToPoint theDCAXBS;
  ClosestApproachInRPhi ClosestApp;
  GlobalPoint XingPoint;

  float muonMassErr = Utility->muonMassErr;
  float pionMassErr = Utility->pionMassErr;
  float kaonMassErr = Utility->kaonMassErr;
  const ParticleMass muonMass = Utility->muonMass;
  const ParticleMass pionMass = Utility->pionMass;
  const ParticleMass kaonMass = Utility->kaonMass;

  std::pair<bool, Measurement1D> theDCAXVtx;
  std::vector<reco::CandidatePtr> footprint;

  std::string MuMCat, MuPCat;
  std::string tmpString1, tmpString2, tmpString3, tmpString4;
  std::stringstream myString;

  TLorentzVector mum_lv;
  TLorentzVector mup_lv;
  TLorentzVector jpsi_lv;

  TLorentzVector tkm_lv;
  TLorentzVector tkp_lv;
  TLorentzVector kst_lv;
  TLorentzVector kstbar_lv;

  KinematicParticleFactoryFromTransientTrack partFactory;
  AdaptiveVertexFitter theVtxFitter;            // Vertex fitter in nominal reconstruction
  KinematicParticleVertexFitter PartVtxFitter;  // Vertex fit with vtx constraint

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(RecVertexToken_, vertices);
  if (vertices->empty())
    return;  // skip the event if no PV found
  reco::Vertex bestVtx;
  for (std::vector<reco::Vertex>::const_iterator iVertex = vertices->begin(); iVertex != vertices->end(); iVertex++) {
    bestVtx = *(iVertex);
    if (bestVtx.isValid() == true)
      break;
  }

  //edm::Handle<pat::MuonCollection> muons;
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  // Get PAT Tracks
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(GenRecTrackToken_, tracks);

  //for (const pat::Muon &mum : *muons) {
  for (const reco::Muon& mum : *muons) {
    muTrackm = mum.innerTrack();
    if ((muTrackm.isNull() == true) || (muTrackm->charge() != -1))
      continue;

    pT = muTrackm->pt();
    eta = muTrackm->eta();
    //         if (muTrackm->hitPattern().trackerLayersWithMeasurement() < 6) continue;
    //         if (muTrackm->hitPattern().pixelLayersWithMeasurement()   < 1) continue;

    if ((pT < (MUMINPT * (1.0 - MUVARTOLE))) || (fabs(eta) > (MUMAXETA * (1.0 + MUVARTOLE)))) {
      if (printMsg_)
        std::cout << __LINE__ << " : break --> too low pT of mu- : " << pT << " or " << mum.pt()
                  << " or too high eta : " << eta << std::endl;
      break;
    }

    const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle));
    if (!muTrackmTT.isValid())
      continue;
    theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(
        GlobalPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()));
    if (theDCAXBS.isValid() == false) {
      if (printMsg_)
        std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu-" << std::endl;
      continue;
    }
    double DCAmumBS = theDCAXBS.perigeeParameters().transverseImpactParameter();
    double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    if (fabs(DCAmumBS) > DCAMUBS) {
      if (printMsg_)
        std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu- : " << DCAmumBS << std::endl;
      continue;
    }

    //for (const pat::Muon &mup : *muons) {
    for (const reco::Muon& mup : *muons) {
      muTrackp = mup.innerTrack();
      if ((muTrackp.isNull() == true) || (muTrackp->charge() != 1))
        continue;

      pT = muTrackp->pt();
      eta = muTrackp->eta();
      if ((pT < (MUMINPT * (1.0 - MUVARTOLE))) || (fabs(eta) > (MUMAXETA * (1.0 + MUVARTOLE)))) {
        if (printMsg_)
          std::cout << __LINE__ << " : break --> too low pT of mu+ : " << pT << " or too high eta : " << eta
                    << std::endl;
        break;
      }

      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));
      if (!muTrackpTT.isValid())
        continue;

      // ###############################
      // # Compute mu+ DCA to BeamSpot #
      // ###############################
      theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(
          GlobalPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()));
      if (theDCAXBS.isValid() == false) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu+" << std::endl;
        continue;
      }
      double DCAmupBS = theDCAXBS.perigeeParameters().transverseImpactParameter();
      double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
      if (fabs(DCAmupBS) > DCAMUBS) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu+: " << DCAmupBS
                    << std::endl;
        continue;
      }

      // ############################################
      // # Check goodness of muons closest approach #
      // ############################################
      ClosestApp.calculate(muTrackpTT.initialFreeState(), muTrackmTT.initialFreeState());
      if (ClosestApp.status() == false) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
        continue;
      }
      XingPoint = ClosestApp.crossingPoint();
      if ((sqrt(XingPoint.x() * XingPoint.x() + XingPoint.y() * XingPoint.y()) > TRKMAXR) ||
          (fabs(XingPoint.z()) > TRKMAXZ)) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume"
                    << std::endl;
        continue;
      }

      double mumuDCA = ClosestApp.distance();
      if (mumuDCA > DCAMUMU) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> bad 3D-DCA of mu+(-) with respect to mu-(+): " << mumuDCA
                    << std::endl;
        continue;
      }

      // ############################################
      // # Cut on the dimuon inviariant mass and pT #
      // ############################################
      jpsi_lv.SetPxPyPzE(muTrackmTT.track().px() + muTrackpTT.track().px(),
                         muTrackmTT.track().py() + muTrackpTT.track().py(),
                         muTrackmTT.track().pz() + muTrackpTT.track().pz(),
                         sqrt(pow(muTrackmTT.track().p(), 2) + pow(Utility->muonMass, 2)) +
                             sqrt(pow(muTrackpTT.track().p(), 2) + pow(Utility->muonMass, 2)));
      if ((jpsi_lv.Pt() < (MINMUMUPT * (1.0 - MUVARTOLE))) || (jpsi_lv.M() < (MINMUMUINVMASS * (1.0 - MUVARTOLE))) ||
          (jpsi_lv.M() > (MAXMUMUINVMASS * (1.0 + MUVARTOLE)))) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << jpsi_lv.Pt()
                    << "\tinv. mass: " << jpsi_lv.M() << std::endl;
        continue;
      }

      chi = 0.;
      ndf = 0.;
      // ####################################################
      // # Try to vertex the two muons to get dimuon vertex #
      // ####################################################
      std::vector<RefCountedKinematicParticle> muonParticles;
      muonParticles.push_back(partFactory.particle(muTrackmTT, muonMass, chi, ndf, muonMassErr));
      muonParticles.push_back(partFactory.particle(muTrackpTT, muonMass, chi, ndf, muonMassErr));

      RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
      if (mumuVertexFitTree->isValid() == false) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
        continue;
      }

      mumuVertexFitTree->movePointerToTheTop();
      RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();
      if (mumu_KV->vertexIsValid() == false) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
        continue;
      }
      if (TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) <
          CLMUMUVTX) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> bad vtx CL from mu+ mu- fit: "
                    << TMath::Prob(static_cast<double>(mumu_KV->chiSquared()),
                                   static_cast<int>(rint(mumu_KV->degreesOfFreedom())))
                    << std::endl;
        continue;
      }

      RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();

      // ######################################################
      // # Compute the distance between mumu vtx and BeamSpot #
      // ######################################################
      double MuMuLSBS;
      double MuMuLSBSErr;
      Utility->computeLS(mumu_KV->position().x(),
                         mumu_KV->position().y(),
                         0.0,
                         beamSpot.position().x(),
                         beamSpot.position().y(),
                         0.0,
                         mumu_KV->error().cxx(),
                         mumu_KV->error().cyy(),
                         0.0,
                         mumu_KV->error().matrix()(0, 1),
                         0.0,
                         0.0,
                         beamSpot.covariance()(0, 0),
                         beamSpot.covariance()(1, 1),
                         0.0,
                         beamSpot.covariance()(0, 1),
                         0.0,
                         0.0,
                         &MuMuLSBS,
                         &MuMuLSBSErr);
      if (MuMuLSBS / MuMuLSBSErr < LSMUMUBS) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> bad mumu L/sigma with respect to BeamSpot: " << MuMuLSBS << "+/-"
                    << MuMuLSBSErr << std::endl;
        continue;
      }
      // ###################################################################
      // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
      // ###################################################################
      double MuMuCosAlphaBS;
      double MuMuCosAlphaBSErr;
      Utility->computeCosAlpha(mumu_KP->currentState().globalMomentum().x(),
                               mumu_KP->currentState().globalMomentum().y(),
                               0.0,
                               mumu_KV->position().x() - beamSpot.position().x(),
                               mumu_KV->position().y() - beamSpot.position().y(),
                               0.0,
                               mumu_KP->currentState().kinematicParametersError().matrix()(3, 3),
                               mumu_KP->currentState().kinematicParametersError().matrix()(4, 4),
                               0.0,
                               mumu_KP->currentState().kinematicParametersError().matrix()(3, 4),
                               0.0,
                               0.0,
                               mumu_KV->error().cxx() + beamSpot.covariance()(0, 0),
                               mumu_KV->error().cyy() + beamSpot.covariance()(1, 1),
                               0.0,
                               mumu_KV->error().matrix()(0, 1) + beamSpot.covariance()(0, 1),
                               0.0,
                               0.0,
                               &MuMuCosAlphaBS,
                               &MuMuCosAlphaBSErr);
      if (MuMuCosAlphaBS < COSALPHAMUMUBS) {
        if (printMsg_)
          std::cout << __LINE__ << " : continue --> bad mumu cos(alpha) with respect to BeamSpot: " << MuMuCosAlphaBS
                    << "+/-" << MuMuCosAlphaBSErr << std::endl;
        continue;
      }

      // ########################### convert KinFit vertex to reco Vertex
      //             reco::Vertex::Point mumu_GP  = reco::Vertex::Point(mumu_KV->position().x(), mumu_KV->position().y(), mumu_KV->position().z());
      //             const reco::Vertex::Error mumu_error = mumu_KV->vertexState().error().matrix();
      //             float mumu_chi2      = mumu_KV -> chiSquared();
      //             float mumu_ndof      = mumu_KV -> degreesOfFreedom();
      //             reco::Vertex mumu_rv =  reco::Vertex( mumu_GP, mumu_error, mumu_chi2, mumu_ndof, 2 );

      for (uint itrkm = 0; itrkm < tracks->size(); itrkm++) {
        reco::TrackRef tkm(tracks, itrkm);
        if (tkm.isNull() == true)
          continue;
        //                 if (! tkm->quality(reco::TrackBase::highPurity))     continue;

        if (tkm->charge() != -1)
          continue;

        if (tkm->pt() < (MINHADPT * (1.0 - HADVARTOLE)))
          continue;
        if (fabs(tkm->eta()) > MUMAXETA)
          continue;

        const reco::TransientTrack TrackmTT((*tkm), &(*bFieldHandle));
        if (!TrackmTT.isValid())
          continue;

        // ######################################
        // # Compute K*0 track- DCA to BeamSpot #
        // ######################################
        theDCAXBS = TrackmTT.trajectoryStateClosestToPoint(
            GlobalPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()));
        if (theDCAXBS.isValid() == false) {
          if (printMsg_)
            std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track-" << std::endl;
          continue;
        }
        double DCAKstTrkmBS = theDCAXBS.perigeeParameters().transverseImpactParameter();
        double DCAKstTrkmBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
        if (printMsg_)
          std::cout << __LINE__ << " ---> DCAKstTrkmBS " << DCAKstTrkmBS << "+/-" << DCAKstTrkmBSErr << std::endl;

        if (!(DCAKstTrkmBSErr > 0))
          continue;
        if (fabs(DCAKstTrkmBS / DCAKstTrkmBSErr) < HADDCASBS) {
          if (printMsg_)
            std::cout << __LINE__
                      << " : continue --> track- DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkmBS
                      << "+/-" << DCAKstTrkmBSErr << std::endl;
          continue;
        }

        for (uint itrkp = 0; itrkp < tracks->size(); itrkp++) {
          reco::TrackRef tkp(tracks, itrkp);
          if (tkp.isNull() == true)
            continue;
          //                     if (! tkp->quality(reco::TrackBase::highPurity))     continue;

          if (tkp->charge() != +1)
            continue;

          // same sign!!
          //                     if ( (tkp->charge()) * (tkm->charge()) != 1 )                 continue;

          if (tkp->pt() < (MINHADPT * (1.0 - HADVARTOLE)))
            continue;
          if (fabs(tkp->eta()) > MUMAXETA)
            continue;
          const reco::TransientTrack TrackpTT((*tkp), &(*bFieldHandle));
          if (!TrackpTT.isValid())
            continue;

          // ######################################
          // # Compute K*0 track+ DCA to BeamSpot #
          // ######################################
          theDCAXBS = TrackpTT.trajectoryStateClosestToPoint(
              GlobalPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()));
          if (theDCAXBS.isValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track+" << std::endl;
            continue;
          }
          double DCAKstTrkpBS = theDCAXBS.perigeeParameters().transverseImpactParameter();
          double DCAKstTrkpBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
          if (!(DCAKstTrkpBSErr > 0))
            continue;
          if (fabs(DCAKstTrkpBS / DCAKstTrkpBSErr) < HADDCASBS) {
            if (printMsg_)
              std::cout << __LINE__
                        << " : continue --> track+ DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkpBS
                        << "+/-" << DCAKstTrkpBSErr << std::endl;
            continue;
          }

          // ##############################################
          // # Check goodness of hadrons closest approach #
          // ##############################################
          ClosestApp.calculate(TrackpTT.initialFreeState(), TrackmTT.initialFreeState());
          if (ClosestApp.status() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
            continue;
          }
          XingPoint = ClosestApp.crossingPoint();
          if ((sqrt(XingPoint.x() * XingPoint.x() + XingPoint.y() * XingPoint.y()) > TRKMAXR) ||
              (fabs(XingPoint.z()) > TRKMAXZ)) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume"
                        << std::endl;
            continue;
          }

          // ######################################################
          // # Check if K*0 (OR K*0bar) mass is within acceptance #
          // ######################################################
          tkm_lv.SetPxPyPzE(TrackmTT.track().momentum().x(),
                            TrackmTT.track().momentum().y(),
                            TrackmTT.track().momentum().z(),
                            sqrt(pow(TrackmTT.track().p(), 2) + pow(Utility->pionMass, 2)));
          tkp_lv.SetPxPyPzE(TrackpTT.track().momentum().x(),
                            TrackpTT.track().momentum().y(),
                            TrackpTT.track().momentum().z(),
                            sqrt(pow(TrackpTT.track().p(), 2) + pow(Utility->kaonMass, 2)));
          double kstInvMass = (tkm_lv + tkp_lv).M();

          tkm_lv.SetE(sqrt(pow(TrackmTT.track().p(), 2) + pow(Utility->kaonMass, 2)));
          tkp_lv.SetE(sqrt(pow(TrackpTT.track().p(), 2) + pow(Utility->pionMass, 2)));
          double kstBarInvMass = (tkm_lv + tkp_lv).M();

          if ((fabs(kstInvMass - Utility->kstMass) > (KSTMASSWINDOW * Utility->kstSigma * (1.0 + HADVARTOLE))) &&
              (fabs(kstBarInvMass - Utility->kstMass) > (KSTMASSWINDOW * Utility->kstSigma * (1.0 + HADVARTOLE)))) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass
                        << " AND K*0bar mass: " << kstBarInvMass << std::endl;
            continue;
          }

          // ####################################################
          // # @@@ Make K* and implement pre-selection cuts @@@ #
          // ####################################################
          if (printMsg_)
            std::cout << "\n"
                      << __LINE__ << " : @@@ I have 2 good oppositely-charged tracks. I'm trying to vertex them @@@"
                      << std::endl;

          chi = 0.;
          ndf = 0.;
          // ##############################################################################
          // # Try to vertex the two Tracks to get K*0 vertex: pion = track- | k = track+ #
          // ##############################################################################
          std::vector<RefCountedKinematicParticle> kstParticles;
          kstParticles.push_back(partFactory.particle(TrackmTT, pionMass, chi, ndf, pionMassErr));
          kstParticles.push_back(partFactory.particle(TrackpTT, kaonMass, chi, ndf, kaonMassErr));

          RefCountedKinematicTree kstVertexFitTree = PartVtxFitter.fit(kstParticles);
          if (kstVertexFitTree->isValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
            continue;
          }

          kstVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle kst_KP = kstVertexFitTree->currentParticle();
          RefCountedKinematicVertex kst_KV = kstVertexFitTree->currentDecayVertex();
          if (kst_KV->vertexIsValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
            continue;
          }

          chi = 0.;
          ndf = 0.;
          // #################################################################################
          // # Try to vertex the two Tracks to get K*0bar vertex: pion = track+ | k = track- #
          // #################################################################################
          std::vector<RefCountedKinematicParticle> kstBarParticles;
          kstBarParticles.push_back(partFactory.particle(TrackmTT, kaonMass, chi, ndf, kaonMassErr));
          kstBarParticles.push_back(partFactory.particle(TrackpTT, pionMass, chi, ndf, pionMassErr));

          RefCountedKinematicTree kstBarVertexFitTree = PartVtxFitter.fit(kstBarParticles);
          if (kstBarVertexFitTree->isValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
            continue;
          }

          kstBarVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle kstBar_KP = kstBarVertexFitTree->currentParticle();
          RefCountedKinematicVertex kstBar_KV = kstBarVertexFitTree->currentDecayVertex();
          if (kstBar_KV->vertexIsValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
            continue;
          }

          // ######################################################
          // # Check if K*0 (OR K*0bar) mass is within acceptance #
          // ######################################################
          kstInvMass = kst_KP->currentState().mass();
          kstBarInvMass = kstBar_KP->currentState().mass();
          if ((fabs(kstInvMass - Utility->kstMass) > KSTMASSWINDOW * Utility->kstSigma) &&
              (fabs(kstBarInvMass - Utility->kstMass) > KSTMASSWINDOW * Utility->kstSigma)) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass
                        << " AND K*0bar mass: " << kstBarInvMass << std::endl;
            continue;
          }

          // ####################################################
          // # @@@ Make B0 and implement pre-selection cuts @@@ #
          // ####################################################
          if (printMsg_)
            std::cout << "\n"
                      << __LINE__ << " : @@@ I have 4 good charged tracks. I'm trying to vertex them @@@" << std::endl;

          TLorentzVector a_lv, b_lv, c_lv, d_lv, tot_lv;
          a_lv.SetPxPyPzE(muTrackmTT.track().momentum().x(),
                          muTrackmTT.track().momentum().y(),
                          muTrackmTT.track().momentum().z(),
                          sqrt(pow(muTrackmTT.track().p(), 2) + pow(Utility->muonMass, 2)));
          b_lv.SetPxPyPzE(muTrackpTT.track().momentum().x(),
                          muTrackpTT.track().momentum().y(),
                          muTrackpTT.track().momentum().z(),
                          sqrt(pow(muTrackpTT.track().p(), 2) + pow(Utility->muonMass, 2)));
          c_lv.SetPxPyPzE(TrackmTT.track().momentum().x(),
                          TrackmTT.track().momentum().y(),
                          TrackmTT.track().momentum().z(),
                          sqrt(pow(TrackmTT.track().p(), 2) + pow(Utility->pionMass, 2)));
          d_lv.SetPxPyPzE(TrackpTT.track().momentum().x(),
                          TrackpTT.track().momentum().y(),
                          TrackpTT.track().momentum().z(),
                          sqrt(pow(TrackpTT.track().p(), 2) + pow(Utility->kaonMass, 2)));

          tot_lv = a_lv + b_lv + c_lv + d_lv;
          if (tot_lv.M() > MAXB0PREMASS) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> b0 mass before fit is > max value" << std::endl;
            continue;
          }

          // #################################################
          // # Check if the hadron tracks are actually muons #
          // #################################################
          MuMCat.clear();
          MuPCat.clear();
          MuMCat = "NotMatched";
          MuPCat = "NotMatched";
          bool foundTkmMum = false;
          bool foundTkpMup = false;

          //for (const pat::Muon &imutmp : *muons) {
          for (const reco::Muon& imutmp : *muons) {
            for (unsigned int i = 0; i < imutmp.numberOfSourceCandidatePtrs(); ++i) {
              const edm::Ptr<reco::Candidate>& source = imutmp.sourceCandidatePtr(i);
              if (!(imutmp.sourceCandidatePtr(i)).isNonnull())
                continue;
              if (!(imutmp.sourceCandidatePtr(i)).isAvailable())
                continue;

              const reco::Candidate& cand = *(source);
              if (cand.charge() == 0)
                continue;
              if (cand.bestTrack() == nullptr)
                continue;
              try {
                cand.bestTrack()->eta();
              } catch (...) {
                std::cout << "should continue: " << std::endl;
                continue;
              }

              if (imutmp.charge() == -1 &&
                  deltaR(tkm->eta(), tkm->phi(), cand.bestTrack()->eta(), cand.bestTrack()->phi()) < 0.00001) {
                MuMCat.clear();
                MuMCat.append(getMuCat(imutmp));
              } else if (imutmp.charge() == 1 &&
                         deltaR(tkp->eta(), tkp->phi(), cand.bestTrack()->eta(), cand.bestTrack()->phi()) < 0.00001) {
                MuPCat.clear();
                MuPCat.append(getMuCat(imutmp));
              }
            }
          }

          if (printMsg_)
            std::cout << __LINE__ << " : now vertexing" << std::endl;

          chi = 0.;
          ndf = 0.;
          // #####################################
          // # B0 vertex fit with vtx constraint #
          // #####################################
          std::vector<RefCountedKinematicParticle> bParticles;
          bParticles.push_back(partFactory.particle(muTrackmTT, muonMass, chi, ndf, muonMassErr));
          bParticles.push_back(partFactory.particle(muTrackpTT, muonMass, chi, ndf, muonMassErr));
          bParticles.push_back(partFactory.particle(TrackmTT, pionMass, chi, ndf, pionMassErr));
          bParticles.push_back(partFactory.particle(TrackpTT, kaonMass, chi, ndf, kaonMassErr));

          RefCountedKinematicTree bVertexFitTree = PartVtxFitter.fit(bParticles);
          if (bVertexFitTree->isValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
            continue;
          }

          bVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle b_KP = bVertexFitTree->currentParticle();
          RefCountedKinematicVertex b_KV = bVertexFitTree->currentDecayVertex();
          if (b_KV->vertexIsValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
            continue;
          }

          chi = 0.;
          ndf = 0.;
          // ########################################
          // # B0bar vertex fit with vtx constraint #
          // ########################################
          std::vector<RefCountedKinematicParticle> bBarParticles;
          bBarParticles.push_back(partFactory.particle(muTrackmTT, muonMass, chi, ndf, muonMassErr));
          bBarParticles.push_back(partFactory.particle(muTrackpTT, muonMass, chi, ndf, muonMassErr));
          bBarParticles.push_back(partFactory.particle(TrackmTT, kaonMass, chi, ndf, kaonMassErr));
          bBarParticles.push_back(partFactory.particle(TrackpTT, pionMass, chi, ndf, pionMassErr));

          RefCountedKinematicTree bBarVertexFitTree = PartVtxFitter.fit(bBarParticles);
          if (bBarVertexFitTree->isValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
            continue;
          }

          bBarVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle bBar_KP = bBarVertexFitTree->currentParticle();
          RefCountedKinematicVertex bBar_KV = bBarVertexFitTree->currentDecayVertex();
          if (bBar_KV->vertexIsValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
            continue;
          }

          // ###########################################################
          // # Extract the re-fitted tracks after the B0 vtx fit #
          // ###########################################################
          bVertexFitTree->movePointerToTheTop();

          bVertexFitTree->movePointerToTheFirstChild();
          const reco::TransientTrack refitMumTT = bVertexFitTree->currentParticle()->refittedTransientTrack();
          bVertexFitTree->movePointerToTheNextChild();
          const reco::TransientTrack refitMupTT = bVertexFitTree->currentParticle()->refittedTransientTrack();
          bVertexFitTree->movePointerToTheNextChild();
          const reco::TransientTrack refitTrkmTT = bVertexFitTree->currentParticle()->refittedTransientTrack();
          bVertexFitTree->movePointerToTheNextChild();
          const reco::TransientTrack refitTrkpTT = bVertexFitTree->currentParticle()->refittedTransientTrack();

          // ########################
          // # Muon pT and eta cuts #
          // ########################
          pT = refitMupTT.track().pt();
          eta = refitMupTT.track().eta();
          if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA)) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> too low pT of mu+ : " << pT << " or too high eta : " << eta
                        << std::endl;
            continue;
          }

          pT = refitMumTT.track().pt();
          eta = refitMumTT.track().eta();
          if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA)) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> too low pT of mu- : " << pT << " or too high eta : " << eta
                        << std::endl;
            continue;
          }

          // ############################################
          // # Cut on the dimuon invariant mass and pT #
          // ############################################
          pT = sqrt((refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) *
                        (refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) +
                    (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()) *
                        (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()));
          double MuMuInvMass = mumu_KP->currentState().mass();
          if ((pT < MINMUMUPT) || (MuMuInvMass < MINMUMUINVMASS) || (MuMuInvMass > MAXMUMUINVMASS)) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << pT << "\tinv. mass: " << MuMuInvMass
                        << std::endl;
            continue;
          }

          // ##########################
          // # Hadron pT and eta cuts #
          // ##########################
          pT = refitTrkpTT.track().pt();
          if (pT < MINHADPT) {
            if (printMsg_)
              std::cout << __LINE__ << " : break --> too low pT of track+ : " << pT << std::endl;
            continue;
          }

          pT = refitTrkmTT.track().pt();
          if (pT < MINHADPT) {
            if (printMsg_)
              std::cout << __LINE__ << " : break --> too low pT of track- : " << pT << std::endl;
            continue;
          }

          // ########################
          // # Cuts on B0 AND B0bar #
          // ########################
          if (((b_KP->currentState().mass() < B0MASSLOWLIMIT) || (b_KP->currentState().mass() > B0MASSUPLIMIT)) &&
              ((bBar_KP->currentState().mass() < B0MASSLOWLIMIT) || (bBar_KP->currentState().mass() > B0MASSUPLIMIT))) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> bad B0 mass: " << b_KP->currentState().mass()
                        << " AND B0bar mass: " << bBar_KP->currentState().mass() << std::endl;
            continue;
          }
          if ((TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))) <
               CLB0VTX) &&
              (TMath::Prob(static_cast<double>(bBar_KV->chiSquared()),
                           static_cast<int>(rint(bBar_KV->degreesOfFreedom()))) < CLB0VTX)) {
            if (printMsg_) {
              std::cout << __LINE__ << " : continue --> bad vtx CL from B0 fit: "
                        << TMath::Prob(static_cast<double>(b_KV->chiSquared()),
                                       static_cast<int>(rint(b_KV->degreesOfFreedom())));
              std::cout << " AND bad vtx CL from B0bar fit: "
                        << TMath::Prob(static_cast<double>(bBar_KV->chiSquared()),
                                       static_cast<int>(rint(bBar_KV->degreesOfFreedom())))
                        << std::endl;
            }
            continue;
          }

          // ###############################
          // # Cuts on B0 L/sigma BeamSpot #
          // ###############################
          Utility->computeLS(b_KV->position().x(),
                             b_KV->position().y(),
                             0.0,
                             beamSpot.position().x(),
                             beamSpot.position().y(),
                             0.0,
                             b_KV->error().cxx(),
                             b_KV->error().cyy(),
                             0.0,
                             b_KV->error().matrix()(0, 1),
                             0.0,
                             0.0,
                             beamSpot.covariance()(0, 0),
                             beamSpot.covariance()(1, 1),
                             0.0,
                             beamSpot.covariance()(0, 1),
                             0.0,
                             0.0,
                             &LSBS,
                             &LSBSErr);

          // ##############################
          // # Compute B0 DCA to BeamSpot #
          // ##############################
          theDCAXBS = b_KP->refittedTransientTrack().trajectoryStateClosestToPoint(
              GlobalPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()));
          if (theDCAXBS.isValid() == false) {
            if (printMsg_)
              std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for B0" << std::endl;
            continue;
          }
          double DCAB0BS = theDCAXBS.perigeeParameters().transverseImpactParameter();
          double DCAB0BSErr = theDCAXBS.perigeeError().transverseImpactParameterError();

          // #####################################
          // # Compute B0 cos(alpha) to BeamSpot #
          // #####################################
          Utility->computeCosAlpha(b_KP->currentState().globalMomentum().x(),
                                   b_KP->currentState().globalMomentum().y(),
                                   0.0,
                                   b_KV->position().x() - beamSpot.position().x(),
                                   b_KV->position().y() - beamSpot.position().y(),
                                   0.0,
                                   b_KP->currentState().kinematicParametersError().matrix()(3, 3),
                                   b_KP->currentState().kinematicParametersError().matrix()(4, 4),
                                   0.0,
                                   b_KP->currentState().kinematicParametersError().matrix()(3, 4),
                                   0.0,
                                   0.0,
                                   b_KV->error().cxx() + beamSpot.covariance()(0, 0),
                                   b_KV->error().cyy() + beamSpot.covariance()(1, 1),
                                   0.0,
                                   b_KV->error().matrix()(0, 1) + beamSpot.covariance()(0, 1),
                                   0.0,
                                   0.0,
                                   &cosAlphaBS,
                                   &cosAlphaBSErr);

          // #################################################
          // # Try to fit 2muon + 1 trk vertex (K mass hyp) ##
          // #################################################
          chi = 0.;
          ndf = 0.;
          bMinusVtxCL = -99;
          bMcosAlphaBS = -99;
          std::vector<RefCountedKinematicParticle> bMinusParticles;
          bMinusParticles.push_back(partFactory.particle(muTrackmTT, muonMass, chi, ndf, muonMassErr));
          bMinusParticles.push_back(partFactory.particle(muTrackpTT, muonMass, chi, ndf, muonMassErr));
          bMinusParticles.push_back(partFactory.particle(TrackmTT, kaonMass, chi, ndf, kaonMassErr));

          RefCountedKinematicTree bMinusVertexFitTree = PartVtxFitter.fit(bMinusParticles);
          if (bMinusVertexFitTree->isValid() == true) {
            bMinusVertexFitTree->movePointerToTheTop();
            RefCountedKinematicParticle bM_KP = bMinusVertexFitTree->currentParticle();
            RefCountedKinematicVertex bM_KV = bMinusVertexFitTree->currentDecayVertex();

            if (bM_KV->vertexIsValid() == true) {
              bMinusVtxCL = TMath::Prob(static_cast<double>(bM_KV->chiSquared()),
                                        static_cast<int>(rint(bM_KV->degreesOfFreedom())));
              Utility->computeCosAlpha(bM_KP->currentState().globalMomentum().x(),
                                       bM_KP->currentState().globalMomentum().y(),
                                       0.0,
                                       bM_KV->position().x() - beamSpot.position().x(),
                                       bM_KV->position().y() - beamSpot.position().y(),
                                       0.0,
                                       bM_KP->currentState().kinematicParametersError().matrix()(3, 3),
                                       bM_KP->currentState().kinematicParametersError().matrix()(4, 4),
                                       0.0,
                                       bM_KP->currentState().kinematicParametersError().matrix()(3, 4),
                                       0.0,
                                       0.0,
                                       bM_KV->error().cxx() + beamSpot.covariance()(0, 0),
                                       bM_KV->error().cyy() + beamSpot.covariance()(1, 1),
                                       0.0,
                                       bM_KV->error().matrix()(0, 1) + beamSpot.covariance()(0, 1),
                                       0.0,
                                       0.0,
                                       &bMcosAlphaBS,
                                       &bMcosAlphaBSErr);
            }
          }

          chi = 0.;
          ndf = 0.;
          bPlusVtxCL = -99;
          bPcosAlphaBS = -99;
          std::vector<RefCountedKinematicParticle> bPlusParticles;
          bPlusParticles.push_back(partFactory.particle(muTrackmTT, muonMass, chi, ndf, muonMassErr));
          bPlusParticles.push_back(partFactory.particle(muTrackpTT, muonMass, chi, ndf, muonMassErr));
          bPlusParticles.push_back(partFactory.particle(TrackpTT, kaonMass, chi, ndf, kaonMassErr));

          RefCountedKinematicTree bPlusVertexFitTree = PartVtxFitter.fit(bPlusParticles);
          if (bPlusVertexFitTree->isValid() == true) {
            bPlusVertexFitTree->movePointerToTheTop();
            RefCountedKinematicParticle bP_KP = bPlusVertexFitTree->currentParticle();
            RefCountedKinematicVertex bP_KV = bPlusVertexFitTree->currentDecayVertex();

            if (bP_KV->vertexIsValid() == true) {
              bPlusVtxCL = TMath::Prob(static_cast<double>(bP_KV->chiSquared()),
                                       static_cast<int>(rint(bP_KV->degreesOfFreedom())));
              Utility->computeCosAlpha(bP_KP->currentState().globalMomentum().x(),
                                       bP_KP->currentState().globalMomentum().y(),
                                       0.0,
                                       bP_KV->position().x() - beamSpot.position().x(),
                                       bP_KV->position().y() - beamSpot.position().y(),
                                       0.0,
                                       bP_KP->currentState().kinematicParametersError().matrix()(3, 3),
                                       bP_KP->currentState().kinematicParametersError().matrix()(4, 4),
                                       0.0,
                                       bP_KP->currentState().kinematicParametersError().matrix()(3, 4),
                                       0.0,
                                       0.0,
                                       bP_KV->error().cxx() + beamSpot.covariance()(0, 0),
                                       bP_KV->error().cyy() + beamSpot.covariance()(1, 1),
                                       0.0,
                                       bP_KV->error().matrix()(0, 1) + beamSpot.covariance()(0, 1),
                                       0.0,
                                       0.0,
                                       &bPcosAlphaBS,
                                       &bPcosAlphaBSErr);
            }
          }

          // #######################################
          // # @@@ Fill B0-candidate variables @@@ #
          // #######################################
          if (printMsg_)
            std::cout << __LINE__ << " : @@@ Filling B0 candidate variables @@@\n\n" << std::endl;

          //// ############
          //// # Save: B0 #
          //// ############
          NTuple->nB++;

          NTuple->bMass->push_back(b_KP->currentState().mass());
          NTuple->bMassE->push_back(sqrt(b_KP->currentState().kinematicParametersError().matrix()(6, 6)));
          NTuple->bBarMass->push_back(bBar_KP->currentState().mass());
          NTuple->bBarMassE->push_back(sqrt(bBar_KP->currentState().kinematicParametersError().matrix()(6, 6)));

          NTuple->bPx->push_back(b_KP->currentState().globalMomentum().x());
          NTuple->bPy->push_back(b_KP->currentState().globalMomentum().y());
          NTuple->bPz->push_back(b_KP->currentState().globalMomentum().z());

          NTuple->bVtxCL->push_back(
              TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))));
          NTuple->bVtxX->push_back(b_KV->position().x());
          NTuple->bVtxY->push_back(b_KV->position().y());
          NTuple->bVtxZ->push_back(b_KV->position().z());

          //edm::LogPrint("MtdSecondaryPvValidation") << "B0    vtx = " << b_KV->position();
          //edm::LogPrint("MtdSecondaryPvValidation") << "B0bar vtx = " << bBar_KV->position();

          NTuple->bCosAlphaBS->push_back(cosAlphaBS);
          NTuple->bCosAlphaBSE->push_back(cosAlphaBSErr);

          NTuple->bLBS->push_back(LSBS);
          NTuple->bLBSE->push_back(LSBSErr);

          NTuple->bDCABS->push_back(DCAB0BS);
          NTuple->bDCABSE->push_back(DCAB0BSErr);

          //// #############
          //// # Save: K*0 #
          //// #############
          NTuple->kstMass->push_back(kstInvMass);
          NTuple->kstMassE->push_back(sqrt(kst_KP->currentState().kinematicParametersError().matrix()(6, 6)));
          NTuple->kstBarMass->push_back(kstBarInvMass);
          NTuple->kstBarMassE->push_back(sqrt(kstBar_KP->currentState().kinematicParametersError().matrix()(6, 6)));

          NTuple->kstPx->push_back(kst_KP->currentState().globalMomentum().x());
          NTuple->kstPy->push_back(kst_KP->currentState().globalMomentum().y());
          NTuple->kstPz->push_back(kst_KP->currentState().globalMomentum().z());

          NTuple->kstVtxCL->push_back(TMath::Prob(static_cast<double>(kst_KV->chiSquared()),
                                                  static_cast<int>(rint(kst_KV->degreesOfFreedom()))));
          NTuple->kstVtxX->push_back(kst_KV->position().x());
          NTuple->kstVtxY->push_back(kst_KV->position().y());
          NTuple->kstVtxZ->push_back(kst_KV->position().z());

          //// #################
          //// # Save: mu+ mu- #
          //// #################
          NTuple->mumuMass->push_back(MuMuInvMass);
          NTuple->mumuMassE->push_back(sqrt(mumu_KP->currentState().kinematicParametersError().matrix()(6, 6)));

          NTuple->mumuPx->push_back(mumu_KP->currentState().globalMomentum().x());
          NTuple->mumuPy->push_back(mumu_KP->currentState().globalMomentum().y());
          NTuple->mumuPz->push_back(mumu_KP->currentState().globalMomentum().z());

          NTuple->mumuVtxCL->push_back(TMath::Prob(static_cast<double>(mumu_KV->chiSquared()),
                                                   static_cast<int>(rint(mumu_KV->degreesOfFreedom()))));
          NTuple->mumuVtxX->push_back(mumu_KV->position().x());
          NTuple->mumuVtxY->push_back(mumu_KV->position().y());
          NTuple->mumuVtxZ->push_back(mumu_KV->position().z());

          NTuple->mumuCosAlphaBS->push_back(MuMuCosAlphaBS);
          NTuple->mumuCosAlphaBSE->push_back(MuMuCosAlphaBSErr);
          NTuple->mumuLBS->push_back(MuMuLSBS);
          NTuple->mumuLBSE->push_back(MuMuLSBSErr);
          NTuple->mumuDCA->push_back(mumuDCA);

          //// #############
          //// # Save: mu- #
          //// #############
          NTuple->mumHighPurity->push_back((int)muTrackm->quality(reco::Track::highPurity));
          NTuple->mumCL->push_back(TMath::Prob(muTrackmTT.chi2(), static_cast<int>(rint(muTrackmTT.ndof()))));
          NTuple->mumNormChi2->push_back(muTrackm->normalizedChi2());
          NTuple->mumPx->push_back(refitMumTT.track().momentum().x());
          NTuple->mumPy->push_back(refitMumTT.track().momentum().y());
          NTuple->mumPz->push_back(refitMumTT.track().momentum().z());

          //                     NTuple->mumDCAVtx->push_back(DCAmumVtx);
          //                     NTuple->mumDCAVtxE->push_back(DCAmumVtxErr);
          NTuple->mumDCABS->push_back(DCAmumBS);
          NTuple->mumDCABSE->push_back(DCAmumBSErr);

          //                     NTuple->mumKinkChi2->push_back(iMuonM->combinedQuality().trkKink);
          NTuple->mumFracHits->push_back(
              static_cast<double>(muTrackm->hitPattern().numberOfValidHits()) /
              static_cast<double>(muTrackm->hitPattern().numberOfValidHits() +
                                  muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                  muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
                                  muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
          //                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackmTT, bestVtxReFit);
          //                     NTuple->mumdxyVtx->push_back(theDCAXVtx.second.value());
          //                     NTuple->mumdzVtx->push_back(muTrackmTT.track().dz( ));
          NTuple->mumdxyBS->push_back(muTrackmTT.track().dxy((beamSpot.position())));
          NTuple->mumdzBS->push_back(muTrackmTT.track().dz((beamSpot.position())));

          NTuple->mumCat->push_back(getMuCat(mum));

          NTuple->mumNPixHits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
          NTuple->mumNPixLayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());
          NTuple->mumNTrkHits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
          NTuple->mumNTrkLayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
          if (mum.isGlobalMuon() == true)
            NTuple->mumNMuonHits->push_back(mum.globalTrack()->hitPattern().numberOfValidMuonHits());
          else
            NTuple->mumNMuonHits->push_back(0);
          NTuple->mumNMatchStation->push_back(mum.numberOfMatchedStations());

          NTuple->trkrefMum->push_back(muTrackm);
          reco::TrackBaseRef tbrTrk1(muTrackm);
          auto tp_info = getMatchedTP(tbrTrk1);
          bool isMatched = false;
          if (tp_info != nullptr) {
            // select for BPH
            if ((*tp_info)->g4Tracks()[0].genpartIndex() != -1) {
              const HepMC::GenParticle* genP = mc->barcode_to_particle((*tp_info)->g4Tracks()[0].genpartIndex());
              if (genSelBPHmu(*genP)) {
                isMatched = true;
                if (printMsg_) {
                  edm::LogPrint("MtdSecondaryPvValidation")
                      << "Mu- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                      << genP->momentum().eta() << " phi " << genP->momentum().phi();
                }
              }
            }
          }
          if (isMatched) {
            NTuple->matchMum->push_back(1);
          } else {
            NTuple->matchMum->push_back(0);
          }

          //// #############
          //// # Save: mu+ #
          //// #############
          NTuple->mupHighPurity->push_back((int)muTrackp->quality(reco::Track::highPurity));
          NTuple->mupCL->push_back(TMath::Prob(muTrackpTT.chi2(), static_cast<int>(rint(muTrackpTT.ndof()))));
          NTuple->mupNormChi2->push_back(muTrackp->normalizedChi2());
          NTuple->mupPx->push_back(refitMupTT.track().momentum().x());
          NTuple->mupPy->push_back(refitMupTT.track().momentum().y());
          NTuple->mupPz->push_back(refitMupTT.track().momentum().z());

          NTuple->mupDCABS->push_back(DCAmupBS);
          NTuple->mupDCABSE->push_back(DCAmupBSErr);
          NTuple->mupdxyBS->push_back(muTrackpTT.track().dxy((beamSpot.position())));
          NTuple->mupdzBS->push_back(muTrackpTT.track().dz((beamSpot.position())));

          //                     NTuple->mupKinkChi2->push_back(iMuonP->combinedQuality().trkKink);
          NTuple->mupFracHits->push_back(
              static_cast<double>(muTrackp->hitPattern().numberOfValidHits()) /
              static_cast<double>(muTrackp->hitPattern().numberOfValidHits() +
                                  muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                  muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
                                  muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
          //                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackpTT, bestVtxReFit);
          //                     NTuple->mupdxyVtx->push_back(theDCAXVtx.second.value());
          //                     NTuple->mupdzVtx->push_back(muTrackpTT.track().dz(bestVtxReFit.position()));

          NTuple->mupCat->push_back(getMuCat(mup));

          NTuple->mupNPixHits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
          NTuple->mupNPixLayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());
          NTuple->mupNTrkHits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
          NTuple->mupNTrkLayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
          if (mup.isGlobalMuon() == true)
            NTuple->mupNMuonHits->push_back(mup.globalTrack()->hitPattern().numberOfValidMuonHits());
          else
            NTuple->mupNMuonHits->push_back(0);
          NTuple->mupNMatchStation->push_back(mup.numberOfMatchedStations());

          NTuple->trkrefMup->push_back(muTrackp);
          reco::TrackBaseRef tbrTrk2(muTrackp);
          tp_info = getMatchedTP(tbrTrk2);
          isMatched = false;
          if (tp_info != nullptr) {
            // select for BPH
            if ((*tp_info)->g4Tracks()[0].genpartIndex() != -1) {
              const HepMC::GenParticle* genP = mc->barcode_to_particle((*tp_info)->g4Tracks()[0].genpartIndex());
              if (genSelBPHmu(*genP)) {
                isMatched = true;
                if (printMsg_) {
                  edm::LogPrint("MtdSecondaryPvValidation")
                      << "Mu- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                      << genP->momentum().eta() << " phi " << genP->momentum().phi();
                }
              }
            }
          }
          if (isMatched) {
            NTuple->matchMup->push_back(1);
          } else {
            NTuple->matchMup->push_back(0);
          }

          //// ################
          //// # Save: Track- #
          //// ################
          NTuple->kstTrkmHighPurity->push_back((int)tkm->quality(reco::Track::highPurity));
          NTuple->kstTrkmCL->push_back(TMath::Prob(TrackmTT.chi2(), static_cast<int>(rint(TrackmTT.ndof()))));
          NTuple->kstTrkmNormChi2->push_back(tkm->normalizedChi2());
          NTuple->kstTrkmPx->push_back(refitTrkmTT.track().momentum().x());
          NTuple->kstTrkmPy->push_back(refitTrkmTT.track().momentum().y());
          NTuple->kstTrkmPz->push_back(refitTrkmTT.track().momentum().z());

          //                     NTuple->kstTrkmDCAVtx->push_back(DCAKstTrkmVtx);
          //                     NTuple->kstTrkmDCAVtxE->push_back(DCAKstTrkmVtxErr);
          NTuple->kstTrkmDCABS->push_back(DCAKstTrkmBS);
          NTuple->kstTrkmDCABSE->push_back(DCAKstTrkmBSErr);

          NTuple->kstTrkmFracHits->push_back(
              static_cast<double>(tkm->hitPattern().numberOfValidHits()) /
              static_cast<double>(tkm->hitPattern().numberOfValidHits() +
                                  tkm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                  tkm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
          //                     NTuple->kstTrkmFracHits->push_back(static_cast<double>(Trackm->hitPattern().numberOfValidHits()) / static_cast<double>(Trackm->hitPattern().numberOfValidHits() +
          //                                                                            Trackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
          //                                                                            Trackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
          theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackmTT, bestVtx);
          NTuple->kstTrkmdxyVtx->push_back(theDCAXVtx.second.value());
          NTuple->kstTrkmdzVtx->push_back(TrackmTT.track().dz(bestVtx.position()));

          NTuple->kstTrkmMuMatch->push_back(MuMCat);

          // I do NOT include the number of missing outer hits because the hadron might interact
          NTuple->kstTrkmNPixHits->push_back(tkm->hitPattern().numberOfValidPixelHits());
          NTuple->kstTrkmNPixLayers->push_back(tkm->hitPattern().pixelLayersWithMeasurement());
          NTuple->kstTrkmNTrkHits->push_back(tkm->hitPattern().numberOfValidTrackerHits());
          NTuple->kstTrkmNTrkLayers->push_back(tkm->hitPattern().trackerLayersWithMeasurement());

          NTuple->trkrefTkm->push_back(tkm);
          reco::TrackBaseRef tbrTrk3(tkm);
          tp_info = getMatchedTP(tbrTrk3);
          isMatched = false;
          if (tp_info != nullptr) {
            // select for BPH
            if ((*tp_info)->g4Tracks()[0].genpartIndex() != -1) {
              const HepMC::GenParticle* genP = mc->barcode_to_particle((*tp_info)->g4Tracks()[0].genpartIndex());
              if (genSelBPHkstar(*genP)) {
                isMatched = true;
                if (printMsg_) {
                  edm::LogPrint("MtdSecondaryPvValidation")
                      << "Trk- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                      << genP->momentum().eta() << " phi " << genP->momentum().phi();
                }
              }
            }
          }
          if (isMatched) {
            NTuple->matchTkm->push_back(1);
          } else {
            NTuple->matchTkm->push_back(0);
          }

          // ################
          // # Save: Track+ #
          // ################
          NTuple->kstTrkpHighPurity->push_back((int)tkp->quality(reco::Track::highPurity));
          NTuple->kstTrkpCL->push_back(TMath::Prob(TrackpTT.chi2(), static_cast<int>(rint(TrackpTT.ndof()))));
          NTuple->kstTrkpNormChi2->push_back(tkp->normalizedChi2());
          NTuple->kstTrkpPx->push_back(refitTrkpTT.track().momentum().x());
          NTuple->kstTrkpPy->push_back(refitTrkpTT.track().momentum().y());
          NTuple->kstTrkpPz->push_back(refitTrkpTT.track().momentum().z());

          //                     NTuple->kstTrkpDCAVtx->push_back(DCAKstTrkpVtx);
          //                     NTuple->kstTrkpDCAVtxE->push_back(DCAKstTrkpVtxErr);
          NTuple->kstTrkpDCABS->push_back(DCAKstTrkpBS);
          NTuple->kstTrkpDCABSE->push_back(DCAKstTrkpBSErr);

          NTuple->kstTrkpFracHits->push_back(
              static_cast<double>(tkp->hitPattern().numberOfValidHits()) /
              static_cast<double>(tkp->hitPattern().numberOfValidHits() +
                                  tkp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                  tkp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
          theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackpTT, bestVtx);
          NTuple->kstTrkpdxyVtx->push_back(theDCAXVtx.second.value());
          NTuple->kstTrkpdzVtx->push_back(TrackpTT.track().dz(bestVtx.position()));

          NTuple->kstTrkpMuMatch->push_back(MuPCat);

          // I do NOT include the number of missing outer hits because the hadron might interact
          NTuple->kstTrkpNPixHits->push_back(tkp->hitPattern().numberOfValidPixelHits());
          NTuple->kstTrkpNPixLayers->push_back(tkp->hitPattern().pixelLayersWithMeasurement());
          NTuple->kstTrkpNTrkHits->push_back(tkp->hitPattern().numberOfValidTrackerHits());
          NTuple->kstTrkpNTrkLayers->push_back(tkp->hitPattern().trackerLayersWithMeasurement());

          NTuple->trkrefTkp->push_back(tkp);
          reco::TrackBaseRef tbrTrk4(tkp);
          tp_info = getMatchedTP(tbrTrk4);
          isMatched = false;
          if (tp_info != nullptr) {
            // select for BPH
            if ((*tp_info)->g4Tracks()[0].genpartIndex() != -1) {
              const HepMC::GenParticle* genP = mc->barcode_to_particle((*tp_info)->g4Tracks()[0].genpartIndex());
              if (genSelBPHkstar(*genP)) {
                isMatched = true;
                if (printMsg_) {
                  edm::LogPrint("MtdSecondaryPvValidation")
                      << "Trk- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                      << genP->momentum().eta() << " phi " << genP->momentum().phi();
                }
              }
            }
          }
          if (isMatched) {
            NTuple->matchTkp->push_back(1);
          } else {
            NTuple->matchTkp->push_back(0);
          }

          NTuple->bPlusCosAlphaBS->push_back(bPcosAlphaBS);
          NTuple->bPlusVtxCL->push_back(bPlusVtxCL);
          NTuple->bMinusCosAlphaBS->push_back(bMcosAlphaBS);
          NTuple->bMinusVtxCL->push_back(bMinusVtxCL);

          // impact parameters from PVs
          B0ImpactPars B0IPs = B0ImpactPars(vertices, refitMumTT, refitMupTT, refitTrkmTT, refitTrkpTT);

          NTuple->mumMinIP2D->push_back(B0IPs.mumMind0.first);
          NTuple->mumMinIP2DE->push_back(B0IPs.mumMind0.second);
          NTuple->mupMinIP2D->push_back(B0IPs.mupMind0.first);
          NTuple->mupMinIP2DE->push_back(B0IPs.mupMind0.second);
          NTuple->kstTrkmMinIP2D->push_back(B0IPs.tkmMind0.first);
          NTuple->kstTrkmMinIP2DE->push_back(B0IPs.tkmMind0.second);
          NTuple->kstTrkpMinIP2D->push_back(B0IPs.tkpMind0.first);
          NTuple->kstTrkpMinIP2DE->push_back(B0IPs.tkpMind0.second);

          NTuple->mumMinIP->push_back(B0IPs.mumMinIP.first);
          NTuple->mumMinIPS->push_back(B0IPs.mumMinIP.second);
          NTuple->mupMinIP->push_back(B0IPs.mupMinIP.first);
          NTuple->mupMinIPS->push_back(B0IPs.mupMinIP.second);
          NTuple->kstTrkmMinIP->push_back(B0IPs.tkmMinIP.first);
          NTuple->kstTrkmMinIPS->push_back(B0IPs.tkmMinIP.second);
          NTuple->kstTrkpMinIP->push_back(B0IPs.tkpMinIP.first);
          NTuple->kstTrkpMinIPS->push_back(B0IPs.tkpMinIP.second);

          // isolation
          B0Isolation B0Iso =
              B0Isolation(bestVtx,
                          tracks,
                          bFieldHandle,
                          ClosestApp,
                          beamSpot,
                          mum,
                          mup,
                          //                                                     muTrackm    , muTrackp  ,
                          itrkm,
                          itrkp,
                          refitMumTT,
                          refitMupTT,
                          refitTrkmTT,
                          refitTrkpTT);

          NTuple->mumIso->push_back(B0Iso.mum_isovec);
          NTuple->mupIso->push_back(B0Iso.mup_isovec);
          NTuple->kstTrkmIso->push_back(B0Iso.trkm_isovec);
          NTuple->kstTrkpIso->push_back(B0Iso.trkp_isovec);

          NTuple->mumIsoPt->push_back(B0Iso.mum_isopts);
          NTuple->mupIsoPt->push_back(B0Iso.mup_isopts);
          NTuple->kstTrkmIsoPt->push_back(B0Iso.trkm_isopts);
          NTuple->kstTrkpIsoPt->push_back(B0Iso.trkp_isopts);
          NTuple->mumIsodR->push_back(B0Iso.mum_isodr);
          NTuple->mupIsodR->push_back(B0Iso.mup_isodr);
          NTuple->kstTrkmIsodR->push_back(B0Iso.trkm_isodr);
          NTuple->kstTrkpIsodR->push_back(B0Iso.trkp_isodr);

          NTuple->rawmumPt->push_back(muTrackm->pt());
          NTuple->rawmumPhi->push_back(muTrackm->phi());
          NTuple->rawmumEta->push_back(muTrackm->eta());
          NTuple->rawmupPt->push_back(muTrackp->pt());
          NTuple->rawmupPhi->push_back(muTrackp->phi());
          NTuple->rawmupEta->push_back(muTrackp->eta());
          NTuple->rawkstTrkmPt->push_back(tkm->pt());
          NTuple->rawkstTrkmPhi->push_back(tkm->phi());
          NTuple->rawkstTrkmEta->push_back(tkm->eta());
          NTuple->rawkstTrkpPt->push_back(tkp->pt());
          NTuple->rawkstTrkpPhi->push_back(tkp->phi());
          NTuple->rawkstTrkpEta->push_back(tkp->eta());

          bParticles.clear();
          bBarParticles.clear();
          kstParticles.clear();
          kstBarParticles.clear();
        }  // end for trackp
      }
      muonParticles.clear();
    }
  }

  NTuple->runN = iEvent.id().run();
  NTuple->eventN = iEvent.id().event();
  NTuple->bsX = beamSpot.position().x();
  NTuple->bsY = beamSpot.position().y();

  for (const reco::Vertex& iVertex : *vertices) {
    if (iVertex.ndof() < PRIVTXNDOF)
      continue;
    if (fabs(iVertex.z()) > PRIVTXMAXZ)
      continue;
    if (fabs(iVertex.position().rho()) > PRIVTXMAXR)
      continue;
    NTuple->recoVtxN++;
  }

  MonteCarloStudies(iEvent);

  for (size_t iB = 0; iB < NTuple->nB; iB++) {
    double pimom(0.);
    double kmom(0.);
    double pieta(0.);
    double keta(0.);
    double pionip(0.);
    double kaonip(0.);
    double distPi(0.);
    double distK(0.);

    unsigned int noPIDtype_Pi(0);
    bool isNoPID_Pi(false);
    bool isPi_Pi(false);
    bool isK_Pi(false);
    bool isP_Pi(false);
    unsigned int noPIDtype_K(0);
    bool isNoPID_K(false);
    bool isPi_K(false);
    bool isK_K(false);
    bool isP_K(false);

    bool isTrue(false);

    reco::Vertex::Point secVtx(NTuple->bVtxX->at(iB), NTuple->bVtxY->at(iB), NTuple->bPz->at(iB));

    double kstar_mass(0.);
    double kstarbar_mass(0.);

    bool pionFail1(false);
    bool kaonFail1(false);
    bool pionFail2(false);
    bool kaonFail2(false);
    if (std::abs(NTuple->kstMass->at(iB) - Utility->kstMass) <
        std::abs(NTuple->kstBarMass->at(iB) - Utility->kstMass)) {
      kstar_mass = NTuple->kstMass->at(iB);
      meKstar_->Fill(kstar_mass);
      pimom = std::sqrt(NTuple->kstTrkmPx->at(iB) * NTuple->kstTrkmPx->at(iB) +
                        NTuple->kstTrkmPy->at(iB) * NTuple->kstTrkmPy->at(iB) +
                        NTuple->kstTrkmPz->at(iB) * NTuple->kstTrkmPz->at(iB));
      kmom = std::sqrt(NTuple->kstTrkpPx->at(iB) * NTuple->kstTrkpPx->at(iB) +
                       NTuple->kstTrkpPy->at(iB) * NTuple->kstTrkpPy->at(iB) +
                       NTuple->kstTrkpPz->at(iB) * NTuple->kstTrkpPz->at(iB));
      pieta = std::abs(
          Utility->computeEta(NTuple->kstTrkmPx->at(iB), NTuple->kstTrkmPy->at(iB), NTuple->kstTrkmPz->at(iB)));
      keta = std::abs(
          Utility->computeEta(NTuple->kstTrkpPx->at(iB), NTuple->kstTrkpPy->at(iB), NTuple->kstTrkpPz->at(iB)));
      pionip = NTuple->kstTrkmMinIPS->at(iB);
      kaonip = NTuple->kstTrkpMinIPS->at(iB);
      auto trackref_pi = NTuple->trkrefTkm->at(iB);
      auto trackref_k = NTuple->trkrefTkp->at(iB);
      distPi = std::sqrt((secVtx - (*trackref_pi).referencePoint()).mag2());
      distK = std::sqrt((secVtx - (*trackref_k).referencePoint()).mag2());
      pidFromTrackRef((double)probPi[trackref_pi],
                      (double)probK[trackref_pi],
                      (double)probP[trackref_pi],
                      noPIDtype_Pi,
                      isNoPID_Pi,
                      isPi_Pi,
                      isK_Pi,
                      isP_Pi);
      pidFromTrackRef((double)probPi[trackref_k],
                      (double)probK[trackref_k],
                      (double)probP[trackref_k],
                      noPIDtype_K,
                      isNoPID_K,
                      isPi_K,
                      isK_K,
                      isP_K);
      pionFail1 = Sigmat0Safe[trackref_pi] == -1 || mtdQualMVA[trackref_pi] < 0.5;
      kaonFail1 = Sigmat0Safe[trackref_k] == -1 || mtdQualMVA[trackref_k] < 0.5;
      if (NTuple->truthMatchSignal->at(iB)) {
        meKstarM_->Fill(kstar_mass);
        if (std::abs(NTuple->genKstTrkmID) == 211 && std::abs(NTuple->genKstTrkpID) == 321) {
          isTrue = true;
          meKstarT_->Fill(kstar_mass);
        }
      }
    } else {
      kstarbar_mass = NTuple->kstBarMass->at(iB);
      meKstarBar_->Fill(kstarbar_mass);
      kmom = std::sqrt(NTuple->kstTrkmPx->at(iB) * NTuple->kstTrkmPx->at(iB) +
                       NTuple->kstTrkmPy->at(iB) * NTuple->kstTrkmPy->at(iB) +
                       NTuple->kstTrkmPz->at(iB) * NTuple->kstTrkmPz->at(iB));
      pimom = std::sqrt(NTuple->kstTrkpPx->at(iB) * NTuple->kstTrkpPx->at(iB) +
                        NTuple->kstTrkpPy->at(iB) * NTuple->kstTrkpPy->at(iB) +
                        NTuple->kstTrkpPz->at(iB) * NTuple->kstTrkpPz->at(iB));
      keta = Utility->computeEta(NTuple->kstTrkmPx->at(iB), NTuple->kstTrkmPy->at(iB), NTuple->kstTrkmPz->at(iB));
      pieta = Utility->computeEta(NTuple->kstTrkpPx->at(iB), NTuple->kstTrkpPy->at(iB), NTuple->kstTrkpPz->at(iB));
      kaonip = NTuple->kstTrkmMinIPS->at(iB);
      pionip = NTuple->kstTrkpMinIPS->at(iB);
      auto trackref_k = NTuple->trkrefTkm->at(iB);
      auto trackref_pi = NTuple->trkrefTkp->at(iB);
      distPi = std::sqrt((secVtx - (*trackref_pi).referencePoint()).mag2());
      distK = std::sqrt((secVtx - (*trackref_k).referencePoint()).mag2());
      pidFromTrackRef((double)probPi[trackref_pi],
                      (double)probK[trackref_pi],
                      (double)probP[trackref_pi],
                      noPIDtype_Pi,
                      isNoPID_Pi,
                      isPi_Pi,
                      isK_Pi,
                      isP_Pi);
      pidFromTrackRef((double)probPi[trackref_k],
                      (double)probK[trackref_k],
                      (double)probP[trackref_k],
                      noPIDtype_K,
                      isNoPID_K,
                      isPi_K,
                      isK_K,
                      isP_K);
      pionFail1 = Sigmat0Safe[trackref_pi] == -1. || mtdQualMVA[trackref_pi] < 0.5;
      kaonFail1 = Sigmat0Safe[trackref_k] == -1. || mtdQualMVA[trackref_k] < 0.5;
      if (NTuple->truthMatchSignal->at(iB)) {
        meKstarBarM_->Fill(kstarbar_mass);
        if (std::abs(NTuple->genKstTrkmID) == 321 && std::abs(NTuple->genKstTrkpID) == 211) {
          isTrue = true;
          meKstarBarT_->Fill(kstarbar_mass);
        }
      }
    }
    mePionp_->Fill(pimom);
    meKaonp_->Fill(kmom);
    mePioneta_->Fill(pieta);
    meKaoneta_->Fill(keta);
    if (pionFail1) {
      mePionpR1_->Fill(pimom);
      mePionetaR1_->Fill(pieta);
    }
    if (kaonFail1) {
      meKaonpR1_->Fill(kmom);
      meKaonetaR1_->Fill(keta);
    }
    mePionIP_->Fill(pionip);
    meKaonIP_->Fill(kaonip);
    if (NTuple->matchTkm->at(iB) && NTuple->matchTkp->at(iB)) {
      mePionpM_->Fill(pimom);
      meKaonpM_->Fill(kmom);
      mePionetaM_->Fill(pieta);
      meKaonetaM_->Fill(keta);
      mePionIPM_->Fill(pionip);
      meKaonIPM_->Fill(kaonip);
      if (pionFail1) {
        mePionpMR1_->Fill(pimom);
        mePionetaMR1_->Fill(pieta);
      }
      if (kaonFail1) {
        meKaonpMR1_->Fill(kmom);
        meKaonetaMR1_->Fill(keta);
      }
      if (isTrue) {
        mePionpT_->Fill(pimom);
        meKaonpT_->Fill(kmom);
        mePionetaT_->Fill(pieta);
        meKaonetaT_->Fill(keta);
        mePionIPT_->Fill(pionip);
        meKaonIPT_->Fill(kaonip);
        if (pionFail1) {
          mePionpTR1_->Fill(pimom);
          mePionetaTR1_->Fill(pieta);
        }
        if (kaonFail1) {
          meKaonpTR1_->Fill(kmom);
          meKaonetaTR1_->Fill(keta);
        }
      }
    }
    double dx = NTuple->bVtxX->at(iB) - bestVtx.x();
    double dy = NTuple->bVtxY->at(iB) - bestVtx.y();
    double dz = NTuple->bVtxZ->at(iB) - bestVtx.z();
    meSVip_->Fill(std::sqrt(dx * dx + dy * dy + dz * dz));

    if (pieta < 1.5) {
      meBarrelPCASVdiff_->Fill(distPi);
      meBarrelNoPIDtype_->Fill(noPIDtype_Pi);
      if (isNoPID_Pi) {
        meBarrelRecoPiNoPID_->Fill(pimom);
      } else if (isPi_Pi) {
        meBarrelRecoPiAsPi_->Fill(pimom);
      } else {
        meBarrelRecoPiAsOth_->Fill(pimom);
      }
    } else if (pieta > 1.6) {
      meEndcapPCASVdiff_->Fill(distPi);
      meEndcapNoPIDtype_->Fill(noPIDtype_Pi);
      if (isNoPID_Pi) {
        meEndcapRecoPiNoPID_->Fill(pimom);
      } else if (isPi_Pi) {
        meEndcapRecoPiAsPi_->Fill(pimom);
      } else {
        meEndcapRecoPiAsOth_->Fill(pimom);
      }
    }
    if (keta < 1.5) {
      meBarrelPCASVdiff_->Fill(distK);
      meBarrelNoPIDtype_->Fill(noPIDtype_K);
      if (isNoPID_K) {
        meBarrelRecoKNoPID_->Fill(kmom);
      } else if (isK_K) {
        meBarrelRecoKAsK_->Fill(kmom);
      } else {
        meBarrelRecoKAsOth_->Fill(kmom);
      }
    } else if (keta > 1.6) {
      meEndcapPCASVdiff_->Fill(distK);
      meEndcapNoPIDtype_->Fill(noPIDtype_K);
      if (isNoPID_K) {
        meEndcapRecoKNoPID_->Fill(kmom);
      } else if (isK_K) {
        meEndcapRecoKAsK_->Fill(kmom);
      } else {
        meEndcapRecoKAsOth_->Fill(kmom);
      }
    }

    // Tof for secondary vertex
    // compute time at secondary vertex for muon candidates

    auto trackref_mum = NTuple->trkrefMum->at(iB);
    //if (Sigmat0Safe[trackref_mum] == -1.) {
    //continue;
    //}
    //if (mtdQualMVA[trackref_mum] < 0.5) {
    //continue;
    //}
    auto trackref_mup = NTuple->trkrefMup->at(iB);
    //if (Sigmat0Safe[trackref_mup] == -1.) {
    //continue;
    //}
    //if (mtdQualMVA[trackref_mup] < 0.5) {
    //continue;
    //}
    auto trackref_tkm = NTuple->trkrefTkm->at(iB);
    auto trackref_tkp = NTuple->trkrefTkp->at(iB);
    if (Sigmat0Safe[trackref_tkm] == -1. || mtdQualMVA[trackref_tkm] < 0.5 || Sigmat0Safe[trackref_tkp] == -1. ||
        mtdQualMVA[trackref_tkp] < 0.5) {
      if (kstar_mass > 0.) {
        meCombKstar_->Fill(kstar_mass);
        if (NTuple->truthMatchSignal->at(iB)) {
          meCombKstarM_->Fill(kstar_mass);
          if (isTrue) {
            meCombKstarT_->Fill(kstar_mass);
          }
        }
      } else if (kstarbar_mass > 0.) {
        meCombKstarBar_->Fill(kstarbar_mass);
        if (NTuple->truthMatchSignal->at(iB)) {
          meCombKstarBarM_->Fill(kstarbar_mass);
          if (isTrue) {
            meCombKstarBarT_->Fill(kstarbar_mass);
          }
        }
      } else {
        edm::LogError("MtdSecondaryPvValidation") << "Not reconstructed K0*/bar mass when needed";
      }
      continue;
    }

    double mummom =
        std::sqrt(NTuple->mumPx->at(iB) * NTuple->mumPx->at(iB) + NTuple->mumPy->at(iB) * NTuple->mumPy->at(iB) +
                  NTuple->mumPz->at(iB) * NTuple->mumPz->at(iB));

    double gammasq_mum = 1. + mummom * mummom * m_pi_inv2;
    double beta_mum = std::sqrt(1. - 1. / gammasq_mum);
    double distMum = std::sqrt((secVtx - (*trackref_mum).referencePoint()).mag2());
    double corrtofMum = distMum / beta_mum * c_inv;
    double crossMum = (secVtx.x() - (*trackref_mum).referencePoint().x()) * NTuple->mumPx->at(iB) +
                      (secVtx.y() - (*trackref_mum).referencePoint().y()) * NTuple->mumPy->at(iB) +
                      (secVtx.z() - (*trackref_mum).referencePoint().z()) * NTuple->mumPz->at(iB);
    int signMum = crossMum > 0. ? 1. : -1;
    double t1 = tMtd[trackref_mum] - tofPi[trackref_mum] + signMum * corrtofMum;
    double w1 = 1. / SigmatMtd[trackref_mum] / SigmatMtd[trackref_mum];
    if (Sigmat0Safe[trackref_mum] == -1. || mtdQualMVA[trackref_mum] < 0.5) {
      t1 = 0.;
      w1 = 0.;
    }

    double mupmom =
        std::sqrt(NTuple->mupPx->at(iB) * NTuple->mupPx->at(iB) + NTuple->mupPy->at(iB) * NTuple->mupPy->at(iB) +
                  NTuple->mupPz->at(iB) * NTuple->mupPz->at(iB));
    double gammasq_mup = 1. + mupmom * mupmom * m_pi_inv2;
    double beta_mup = std::sqrt(1. - 1. / gammasq_mup);
    double distMup = std::sqrt((secVtx - (*trackref_mup).referencePoint()).mag2());
    double corrtofMup = distMup / beta_mup * c_inv;
    double crossMup = (secVtx.x() - (*trackref_mup).referencePoint().x()) * NTuple->mupPx->at(iB) +
                      (secVtx.y() - (*trackref_mup).referencePoint().y()) * NTuple->mupPy->at(iB) +
                      (secVtx.z() - (*trackref_mup).referencePoint().z()) * NTuple->mupPz->at(iB);
    int signMup = crossMup > 0. ? 1. : -1;
    double t2 = tMtd[trackref_mup] - tofPi[trackref_mup] + signMup * corrtofMup;
    double w2 = 1. / SigmatMtd[trackref_mup] / SigmatMtd[trackref_mup];
    if (Sigmat0Safe[trackref_mup] == -1. || mtdQualMVA[trackref_mup] < 0.5) {
      t2 = 0.;
      w2 = 0.;
    }

    // compute time at secondary vertex for negative and positive tracks in different hypothesis

    std::array<std::array<double, 3>, 2> hypTof;
    std::array<double, 3> tof;

    double tkmmom = std::sqrt(NTuple->kstTrkmPx->at(iB) * NTuple->kstTrkmPx->at(iB) +
                              NTuple->kstTrkmPy->at(iB) * NTuple->kstTrkmPy->at(iB) +
                              NTuple->kstTrkmPz->at(iB) * NTuple->kstTrkmPz->at(iB));
    double tkmeta =
        std::abs(Utility->computeEta(NTuple->kstTrkmPx->at(iB), NTuple->kstTrkmPy->at(iB), NTuple->kstTrkmPz->at(iB)));
    double tkmip = NTuple->kstTrkmMinIPS->at(iB);
    double distTkm = std::sqrt((secVtx - (*trackref_tkm).referencePoint()).mag2());
    double crossTkm = (secVtx.x() - (*trackref_tkm).referencePoint().x()) * NTuple->kstTrkmPx->at(iB) +
                      (secVtx.y() - (*trackref_tkm).referencePoint().y()) * NTuple->kstTrkmPy->at(iB) +
                      (secVtx.z() - (*trackref_tkm).referencePoint().z()) * NTuple->kstTrkmPz->at(iB);
    int signTkm = crossTkm > 0. ? 1. : -1;

    double gammasq_pi_tkm = 1. + tkmmom * tkmmom * m_pi_inv2;
    double beta_pi_tkm = std::sqrt(1. - 1. / gammasq_pi_tkm);
    double corrtof_pi_tkm = distTkm / beta_pi_tkm * c_inv;
    double gammasq_k_tkm = 1. + tkmmom * tkmmom * m_k_inv2;
    double beta_k_tkm = std::sqrt(1. - 1. / gammasq_k_tkm);
    double corrtof_k_tkm = distTkm / beta_k_tkm * c_inv;
    double gammasq_p_tkm = 1. + tkmmom * tkmmom * m_p_inv2;
    double beta_p_tkm = std::sqrt(1. - 1. / gammasq_p_tkm);
    double corrtof_p_tkm = distTkm / beta_p_tkm * c_inv;
    tof[0] = tMtd[trackref_tkm] - tofPi[trackref_tkm] + signTkm * corrtof_pi_tkm;
    tof[1] = tMtd[trackref_tkm] - tofK[trackref_tkm] + signTkm * corrtof_k_tkm;
    tof[2] = tMtd[trackref_tkm] - tofP[trackref_tkm] + signTkm * corrtof_p_tkm;
    hypTof[0] = tof;

    if (printMsg_) {
      edm::LogPrint("MtdSecondaryPvValidation")
          << "Track- tof pi/k/p " << tof[0] << " " << tof[1] << " " << tof[2] << " corr " << corrtof_pi_tkm << " "
          << corrtof_k_tkm << " " << corrtof_p_tkm << " sigmatMtd " << SigmatMtd[trackref_tkm];
    }

    double tkpmom = std::sqrt(NTuple->kstTrkpPx->at(iB) * NTuple->kstTrkpPx->at(iB) +
                              NTuple->kstTrkpPy->at(iB) * NTuple->kstTrkpPy->at(iB) +
                              NTuple->kstTrkpPz->at(iB) * NTuple->kstTrkpPz->at(iB));
    double tkpeta =
        std::abs(Utility->computeEta(NTuple->kstTrkpPx->at(iB), NTuple->kstTrkpPy->at(iB), NTuple->kstTrkpPz->at(iB)));
    double tkpip = NTuple->kstTrkpMinIPS->at(iB);
    double distTkp = std::sqrt((secVtx - (*trackref_tkp).referencePoint()).mag2());
    double crossTkp = (secVtx.x() - (*trackref_tkp).referencePoint().x()) * NTuple->kstTrkpPx->at(iB) +
                      (secVtx.y() - (*trackref_tkp).referencePoint().y()) * NTuple->kstTrkpPy->at(iB) +
                      (secVtx.z() - (*trackref_tkp).referencePoint().z()) * NTuple->kstTrkpPz->at(iB);
    int signTkp = crossTkp > 0. ? 1. : -1;

    double gammasq_pi_tkp = 1. + tkpmom * tkpmom * m_pi_inv2;
    double beta_pi_tkp = std::sqrt(1. - 1. / gammasq_pi_tkp);
    double corrtof_pi_tkp = distTkp / beta_pi_tkp * c_inv;
    double gammasq_k_tkp = 1. + tkpmom * tkpmom * m_k_inv2;
    double beta_k_tkp = std::sqrt(1. - 1. / gammasq_k_tkp);
    double corrtof_k_tkp = distTkp / beta_k_tkp * c_inv;
    double gammasq_p_tkp = 1. + tkpmom * tkpmom * m_p_inv2;
    double beta_p_tkp = std::sqrt(1. - 1. / gammasq_p_tkp);
    double corrtof_p_tkp = distTkp / beta_p_tkp * c_inv;
    tof[0] = tMtd[trackref_tkp] - tofPi[trackref_tkp] + signTkp * corrtof_pi_tkp;
    tof[1] = tMtd[trackref_tkp] - tofK[trackref_tkp] + signTkp * corrtof_k_tkp;
    tof[2] = tMtd[trackref_tkp] - tofP[trackref_tkp] + signTkp * corrtof_p_tkp;
    hypTof[1] = tof;

    if (printMsg_) {
      edm::LogPrint("MtdSecondaryPvValidation")
          << "Track+ tof pi/k/p " << tof[0] << " " << tof[1] << " " << tof[2] << " corr " << corrtof_pi_tkp << " "
          << corrtof_k_tkp << " " << corrtof_p_tkp << " sigmatMtd " << SigmatMtd[trackref_tkp];
    }

    // Use PID only if pion/kaon hypothesis at least 1 sigma away

    if (std::abs((hypTof[0])[0] - (hypTof[0])[1]) > SigmatMtd[trackref_tkm]) {
      if (std::abs(NTuple->kstMass->at(iB) - Utility->kstMass) <
          std::abs(NTuple->kstBarMass->at(iB) - Utility->kstMass)) {
        pionFail2 = true;
      } else {
        kaonFail2 = true;
      }
    }
    if (std::abs((hypTof[1])[0] - (hypTof[1])[1]) > SigmatMtd[trackref_tkp]) {
      if (std::abs(NTuple->kstMass->at(iB) - Utility->kstMass) <
          std::abs(NTuple->kstBarMass->at(iB) - Utility->kstMass)) {
        kaonFail2 = true;
      } else {
        pionFail2 = true;
      }
    }

    if (pionFail2) {
      mePionpR2_->Fill(pimom);
      mePionetaR2_->Fill(pieta);
      if (NTuple->truthMatchSignal->at(iB)) {
        mePionpMR2_->Fill(pimom);
        mePionetaMR2_->Fill(pieta);
        if (isTrue) {
          mePionpTR2_->Fill(pimom);
          mePionetaTR2_->Fill(pieta);
        }
      }
    }
    if (kaonFail2) {
      meKaonpR2_->Fill(kmom);
      meKaonetaR2_->Fill(keta);
      if (NTuple->truthMatchSignal->at(iB)) {
        meKaonpMR2_->Fill(kmom);
        meKaonetaMR2_->Fill(keta);
        if (isTrue) {
          meKaonpTR2_->Fill(kmom);
          meKaonetaTR2_->Fill(keta);
        }
      }
    }

    if (std::abs((hypTof[0])[0] - (hypTof[0])[1]) > SigmatMtd[trackref_tkm] ||
        std::abs((hypTof[1])[0] - (hypTof[1])[1]) > SigmatMtd[trackref_tkp]) {
      if (std::abs(NTuple->kstMass->at(iB) - Utility->kstMass) <
          std::abs(NTuple->kstBarMass->at(iB) - Utility->kstMass)) {
        meCombKstar_->Fill(kstar_mass);
        if (NTuple->truthMatchSignal->at(iB)) {
          meCombKstarM_->Fill(kstar_mass);
          if (isTrue) {
            meCombKstarT_->Fill(kstar_mass);
          }
        }
      } else if (std::abs(NTuple->kstMass->at(iB) - Utility->kstMass) >
                 std::abs(NTuple->kstBarMass->at(iB) - Utility->kstMass)) {
        meCombKstarBar_->Fill(kstarbar_mass);
        if (NTuple->truthMatchSignal->at(iB)) {
          meCombKstarBarM_->Fill(kstarbar_mass);
          if (isTrue) {
            meCombKstarBarT_->Fill(kstarbar_mass);
          }
        }
      } else {
        edm::LogError("MtdSecondaryPvValidation") << "Not reconstructed K0*/bar mass when needed";
      }
      continue;
    }

    double sigdtAtSV = 9999.;
    unsigned int id1 = 3;
    unsigned int id2 = 3;
    for (unsigned int ihyp1 = 0; ihyp1 < 3; ihyp1++) {
      for (unsigned int ihyp2 = 0; ihyp2 < 3; ihyp2++) {
        double w3 = 1. / SigmatMtd[trackref_tkm] / SigmatMtd[trackref_tkm];
        double t3 = (hypTof[0])[ihyp1];
        double w4 = 1. / SigmatMtd[trackref_tkp] / SigmatMtd[trackref_tkp];
        double t4 = (hypTof[1])[ihyp2];
        double wave = (w1 * t1 + w2 * t2 + w3 * t3 + w4 * t4) / (w1 + w2 + w3 + w4);
        double wstd = std::sqrt((w1 * (t1 - wave) * (t1 - wave) + w2 * (t2 - wave) * (t2 - wave) +
                                 w3 * (t3 - wave) * (t3 - wave) + w4 * (t4 - wave) * (t4 - wave)) /
                                (3. / 4. * (w1 + w2 + w3 + w4)));
        double sigdtTmp = wstd;
        if (sigdtTmp < sigdtAtSV) {
          sigdtAtSV = sigdtTmp;
          id1 = ihyp1;
          id2 = ihyp2;
        }
      }
    }

    // only pion - kaon pairs are considered

    if (id1 > 1 || id2 > 1 || id1 == id2) {
      continue;
    }

    // tkm == pion && tkp == kaon, i.e. K0*

    if (id1 == 0 && id2 == 1) {
      mePIDKstar_->Fill(NTuple->kstMass->at(iB));
      meCombKstar_->Fill(NTuple->kstMass->at(iB));
      mePIDPionp_->Fill(tkmmom);
      mePIDKaonp_->Fill(tkpmom);
      mePIDPioneta_->Fill(tkmeta);
      mePIDKaoneta_->Fill(tkpeta);
      mePIDPionIP_->Fill(tkmip);
      mePIDKaonIP_->Fill(tkpip);
      if (NTuple->truthMatchSignal->at(iB)) {
        mePIDKstarM_->Fill(NTuple->kstMass->at(iB));
        meCombKstarM_->Fill(NTuple->kstMass->at(iB));
        if (std::abs(NTuple->genKstTrkmID) == 211 && std::abs(NTuple->genKstTrkpID) == 321) {
          mePIDKstarT_->Fill(NTuple->kstMass->at(iB));
          meCombKstarT_->Fill(NTuple->kstMass->at(iB));
        }
      }
      if (NTuple->matchTkm->at(iB) && NTuple->matchTkp->at(iB)) {
        mePIDPionpM_->Fill(tkmmom);
        mePIDPionpM_->Fill(tkmmom);
        mePIDKaonpM_->Fill(tkpmom);
        mePIDPionetaM_->Fill(tkmeta);
        mePIDKaonetaM_->Fill(tkpeta);
        mePIDPionIPM_->Fill(tkmip);
        mePIDKaonIPM_->Fill(tkpip);
        if (std::abs(NTuple->genKstTrkmID) == 211 && std::abs(NTuple->genKstTrkpID) == 321) {
          mePIDPionpT_->Fill(tkmmom);
          mePIDKaonpT_->Fill(tkpmom);
          mePIDPionetaT_->Fill(tkmeta);
          mePIDKaonetaT_->Fill(tkpeta);
          mePIDPionIPT_->Fill(tkmip);
          mePIDKaonIPT_->Fill(tkpip);
        }
      }

      // tkm == kaon && tkp == pion, i.e. K0*bar

    } else if (id1 == 1 && id2 == 0) {
      mePIDKstarBar_->Fill(NTuple->kstMass->at(iB));
      meCombKstarBar_->Fill(NTuple->kstMass->at(iB));
      mePIDPionp_->Fill(tkpmom);
      mePIDKaonp_->Fill(tkmmom);
      mePIDPioneta_->Fill(tkpeta);
      mePIDKaoneta_->Fill(tkmeta);
      mePIDPionIP_->Fill(tkpip);
      mePIDKaonIP_->Fill(tkmip);
      if (NTuple->truthMatchSignal->at(iB)) {
        mePIDKstarBarM_->Fill(NTuple->kstBarMass->at(iB));
        meCombKstarBarM_->Fill(NTuple->kstBarMass->at(iB));
        if (std::abs(NTuple->genKstTrkmID) == 321 && std::abs(NTuple->genKstTrkpID) == 211) {
          mePIDKstarBarT_->Fill(NTuple->kstBarMass->at(iB));
          meCombKstarBarT_->Fill(NTuple->kstBarMass->at(iB));
        }
      }
      if (NTuple->matchTkm->at(iB) && NTuple->matchTkp->at(iB)) {
        mePIDPionpM_->Fill(tkpmom);
        mePIDKaonpM_->Fill(tkmmom);
        mePIDPionetaM_->Fill(tkpeta);
        mePIDKaonetaM_->Fill(tkmeta);
        mePIDPionIPM_->Fill(tkpip);
        mePIDKaonIPM_->Fill(tkmip);
        if (std::abs(NTuple->genKstTrkmID) == 321 && std::abs(NTuple->genKstTrkpID) == 211) {
          mePIDPionpT_->Fill(tkpmom);
          mePIDKaonpT_->Fill(tkmmom);
          mePIDPionetaT_->Fill(tkpeta);
          mePIDKaonetaT_->Fill(tkmeta);
          mePIDPionIPT_->Fill(tkpip);
          mePIDKaonIPT_->Fill(tkmip);
        }
      }

    } else {
      edm::LogError("MtdSecondaryPvValidation") << "Incorrect ID combination found!!! " << id1 << " " << id2;
      continue;
    }
  }

  NTuple->ClearNTuple();
}

// ------------ method for histogram booking ------------
void MtdSecondaryPvValidation::bookHistograms(DQMStore::IBooker& ibook,
                                              edm::Run const& run,
                                              edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  meKstar_ = ibook.book1D("Kstar", "Kstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meKstarBar_ = ibook.book1D("KstarBar", "KstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meKstarM_ = ibook.book1D("KstarM", "Matched Kstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meKstarBarM_ = ibook.book1D("KstarBarM", "Matched KstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meKstarT_ = ibook.book1D("KstarT", "True Kstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meKstarBarT_ = ibook.book1D("KstarBarT", "True KstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);

  mePionp_ = ibook.book1D("Pionp", "Pion candidate p; p [GeV]", 25, 0., 10.);
  mePionpM_ = ibook.book1D("PionpM", "Matched Pion candidate p; p [GeV]", 25, 0., 10.);
  mePionpT_ = ibook.book1D("PionpT", "True Pion candidate p; p [GeV]", 25, 0., 10.);
  meKaonp_ = ibook.book1D("Kaonp", "Kaon candidate p; p [GeV]", 25, 0., 10.);
  meKaonpM_ = ibook.book1D("KaonpM", "Matched Kaon candidate p; p [GeV]", 25, 0., 10.);
  meKaonpT_ = ibook.book1D("KaonpT", "True Kaon candidate p; p [GeV]", 25, 0., 10.);

  mePioneta_ = ibook.book1D("Pioneta", "Pion candidate eta; eta", 66, 0., 3.3);
  mePionetaM_ = ibook.book1D("PionetaM", "Matched Pion candidate eta; eta", 66, 0., 3.3);
  mePionetaT_ = ibook.book1D("PionetaT", "True Pion candidate eta; eta", 66, 0., 3.3);
  meKaoneta_ = ibook.book1D("Kaoneta", "Kaon candidate eta; eta", 66, 0., 3.3);
  meKaonetaM_ = ibook.book1D("KaonetaM", "Matched Kaon candidate eta; eta", 66, 0., 3.3);
  meKaonetaT_ = ibook.book1D("KaonetaT", "True Kaon candidate eta; eta", 66, 0., 3.3);

  mePionIP_ = ibook.book1D("PionIP", "Pion candidate min IP; IP [cm]", 100, 0., 20.);
  mePionIPM_ = ibook.book1D("PionIPM", "Matched Pion candidate min IP; IP [cm]", 100, 0., 20.);
  mePionIPT_ = ibook.book1D("PionIPT", "True Pion candidate min IP; IP [cm]", 100, 0., 20.);
  meKaonIP_ = ibook.book1D("KaonIP", "Kaon candidate min IP; IP [cm]", 100, 0., 20.);
  meKaonIPM_ = ibook.book1D("KaonIPM", "Matched Kaon candidate min IP; IP [cm]", 100, 0., 20.);
  meKaonIPT_ = ibook.book1D("KaonIPT", "True Kaon candidate min IP; IP [cm]", 100, 0., 20.);

  mePIDKstar_ = ibook.book1D("PIDKstar", "PIDKstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  mePIDKstarBar_ = ibook.book1D("PIDKstarBar", "PIDKstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  mePIDKstarM_ = ibook.book1D("PIDKstarM", "Matched PIDKstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  mePIDKstarBarM_ = ibook.book1D("PIDKstarBarM", "Matched PIDKstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  mePIDKstarT_ = ibook.book1D("PIDKstarT", "True PIDKstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  mePIDKstarBarT_ = ibook.book1D("PIDKstarBarT", "True PIDKstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);

  meCombKstar_ = ibook.book1D("CombKstar", "CombKstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meCombKstarBar_ = ibook.book1D("CombKstarBar", "CombKstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meCombKstarM_ = ibook.book1D("CombKstarM", "Matched CombKstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meCombKstarBarM_ = ibook.book1D("CombKstarBarM", "Matched CombKstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meCombKstarT_ = ibook.book1D("CombKstarT", "True CombKstar candidate mass; mass [GeV]", 50, 0.65, 1.15);
  meCombKstarBarT_ = ibook.book1D("CombKstarBarT", "True CombKstarBar candidate mass; mass [GeV]", 50, 0.65, 1.15);

  mePIDPionp_ = ibook.book1D("PIDPionp", "PIDPion candidate p; p [GeV]", 25, 0., 10.);
  mePIDPionpM_ = ibook.book1D("PIDPionpM", "Matched PIDPion candidate p; p [GeV]", 25, 0., 10.);
  mePIDPionpT_ = ibook.book1D("PIDPionpT", "True PIDPion candidate p; p [GeV]", 25, 0., 10.);
  mePIDKaonp_ = ibook.book1D("PIDKaonp", "PIDKaon candidate p; p [GeV]", 25, 0., 10.);
  mePIDKaonpM_ = ibook.book1D("PIDKaonpM", "Matched PIDKaon candidate p; p [GeV]", 25, 0., 10.);
  mePIDKaonpT_ = ibook.book1D("PIDKaonpT", "True PIDKaon candidate p; p [GeV]", 25, 0., 10.);

  mePIDPioneta_ = ibook.book1D("PIDPioneta", "PIDPion candidate eta; eta", 66, 0., 3.3);
  mePIDPionetaM_ = ibook.book1D("PIDPionetaM", "Matched PIDPion candidate eta; eta", 66, 0., 3.3);
  mePIDPionetaT_ = ibook.book1D("PIDPionetaT", "True PIDPion candidate eta; eta", 66, 0., 3.3);
  mePIDKaoneta_ = ibook.book1D("PIDKaoneta", "PIDKaon candidate eta; eta", 66, 0., 3.3);
  mePIDKaonetaM_ = ibook.book1D("PIDKaonetaM", "Matched PIDKaon candidate eta; eta", 66, 0., 3.3);
  mePIDKaonetaT_ = ibook.book1D("PIDKaonetaT", "True PIDKaon candidate eta; eta", 66, 0., 3.3);

  mePIDPionIP_ = ibook.book1D("PIDPionIP", "PIDPion candidate min IP; IP [cm]", 100, 0., 20.);
  mePIDPionIPM_ = ibook.book1D("PIDPionIPM", "Matched PIDPion candidate min IP; IP [cm]", 100, 0., 20.);
  mePIDPionIPT_ = ibook.book1D("PIDPionIPT", "True PIDPion candidate min IP; IP [cm]", 100, 0., 20.);
  mePIDKaonIP_ = ibook.book1D("PIDKaonIP", "PIDKaon candidate min IP; IP [cm]", 100, 0., 20.);
  mePIDKaonIPM_ = ibook.book1D("PIDKaonIPM", "Matched PIDKaon candidate min IP; IP [cm]", 100, 0., 20.);
  mePIDKaonIPT_ = ibook.book1D("PIDKaonIPT", "True PIDKaon candidate min IP; IP [cm]", 100, 0., 20.);

  mePionpR1_ = ibook.book1D("PionpR1", "Pion candidate p rejected 1; p [GeV]", 25, 0., 10.);
  mePionpMR1_ = ibook.book1D("PionpMR1", "Matched Pion candidate p rejected 1; p [GeV]", 25, 0., 10.);
  mePionpTR1_ = ibook.book1D("PionpTR1", "True Pion candidate p rejected 1; p [GeV]", 25, 0., 10.);
  meKaonpR1_ = ibook.book1D("KaonpR1", "Kaon candidate p rejected 1; p [GeV]", 25, 0., 10.);
  meKaonpMR1_ = ibook.book1D("KaonpMR1", "Matched Kaon candidate p rejected 1; p [GeV]", 25, 0., 10.);
  meKaonpTR1_ = ibook.book1D("KaonpTR1", "True Kaon candidate p rejected 1; p [GeV]", 25, 0., 10.);

  mePionetaR1_ = ibook.book1D("PionetaR1", "Pion candidate eta rejected 1; eta", 66, 0., 3.3);
  mePionetaMR1_ = ibook.book1D("PionetaMR1", "Matched Pion candidate eta rejected 1; eta", 66, 0., 3.3);
  mePionetaTR1_ = ibook.book1D("PionetaTR1", "True Pion candidate eta rejected 1; eta", 66, 0., 3.3);
  meKaonetaR1_ = ibook.book1D("KaonetaR1", "Kaon candidate eta rejected 1; eta", 66, 0., 3.3);
  meKaonetaMR1_ = ibook.book1D("KaonetaMR1", "Matched Kaon candidate eta rejected 1; eta", 66, 0., 3.3);
  meKaonetaTR1_ = ibook.book1D("KaonetaTR1", "True Kaon candidate eta rejected 1; eta", 66, 0., 3.3);

  mePionpR2_ = ibook.book1D("PionpR2", "Pion candidate p rejected 2; p [GeV]", 25, 0., 10.);
  mePionpMR2_ = ibook.book1D("PionpMR2", "Matched Pion candidate p rejected 2; p [GeV]", 25, 0., 10.);
  mePionpTR2_ = ibook.book1D("PionpTR2", "True Pion candidate p rejected 2; p [GeV]", 25, 0., 10.);
  meKaonpR2_ = ibook.book1D("KaonpR2", "Kaon candidate p rejected 2; p [GeV]", 25, 0., 10.);
  meKaonpMR2_ = ibook.book1D("KaonpMR2", "Matched Kaon candidate p rejected 2; p [GeV]", 25, 0., 10.);
  meKaonpTR2_ = ibook.book1D("KaonpTR2", "True Kaon candidate p rejected 2; p [GeV]", 25, 0., 10.);

  mePionetaR2_ = ibook.book1D("PionetaR2", "Pion candidate eta rejected 2; eta", 66, 0., 3.3);
  mePionetaMR2_ = ibook.book1D("PionetaMR2", "Matched Pion candidate eta rejected 2; eta", 66, 0., 3.3);
  mePionetaTR2_ = ibook.book1D("PionetaTR2", "True Pion candidate eta rejected 2; eta", 66, 0., 3.3);
  meKaonetaR2_ = ibook.book1D("KaonetaR2", "Kaon candidate eta rejected 2; eta", 66, 0., 3.3);
  meKaonetaMR2_ = ibook.book1D("KaonetaMR2", "Matched Kaon candidate eta rejected 2; eta", 66, 0., 3.3);
  meKaonetaTR2_ = ibook.book1D("KaonetaTR2", "True Kaon candidate eta rejected 2; eta", 66, 0., 3.3);

  meSVip_ = ibook.book1D("SVip", "Secondary vtx distance from PV; ip [cm]", 100, 0., 2.);
  meBarrelPCASVdiff_ = ibook.book1D("BarrelPCASVdiff", "PCA - SV distance, |eta| < 1.5; dist [cm]", 100, 0., 2.);
  meEndcapPCASVdiff_ = ibook.book1D("EndcapPCASVdiff", "PCA - SV distance, |eta| > 1.6; dist [cm]", 100, 0., 2.);

  meBarrelNoPIDtype_ = ibook.book1D("BarrelNoPIDtype", "Barrel PID failure category", 4, 0., 4.);
  meEndcapNoPIDtype_ = ibook.book1D("EndcapNoPIDtype", "Endcap PID failure category", 4, 0., 4.);

  meBarrelRecoPiNoPID_ =
      ibook.book1D("BarrelRecoPiNoPID", "Reco pi NoPID momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
  meBarrelRecoKNoPID_ =
      ibook.book1D("BarrelRecoKNoPID", "Reco k NoPID momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
  meEndcapRecoPiNoPID_ =
      ibook.book1D("EndcapRecoPiNoPID", "Reco pi NoPID momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
  meEndcapRecoKNoPID_ =
      ibook.book1D("EndcapRecoKNoPID", "Reco k NoPID momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);

  meBarrelRecoPiAsPi_ =
      ibook.book1D("BarrelRecoPiAsPi", "Reco pi as pi momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
  meBarrelRecoPiAsOth_ =
      ibook.book1D("BarrelRecoPiAsOth", "Reco pi as other momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
  meEndcapRecoPiAsPi_ =
      ibook.book1D("EndcapRecoPiAsPi", "Reco pi as pi momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
  meEndcapRecoPiAsOth_ =
      ibook.book1D("EndcapRecoPiAsOth", "Reco pi as other momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);

  meBarrelRecoKAsOth_ =
      ibook.book1D("BarrelRecoKAsOth", "Reco k as other momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
  meBarrelRecoKAsK_ = ibook.book1D("BarrelRecoKAsK", "Reco k as k momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
  meEndcapRecoKAsOth_ =
      ibook.book1D("EndcapRecoKAsOth", "Reco k as other momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
  meEndcapRecoKAsK_ = ibook.book1D("EndcapRecoKAsK", "Reco k as k momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdSecondaryPvValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/SecondaryV");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));
  desc.add<edm::InputTag>("pruned", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("probPi", edm::InputTag("tofPID:probPi"));
  desc.add<edm::InputTag>("probK", edm::InputTag("tofPID:probK"));
  desc.add<edm::InputTag>("probP", edm::InputTag("tofPID:probP"));
  desc.add<double>("trackMinimumPt", 0.7);  // [GeV]
  desc.add<double>("trackMaximumBtlEta", 1.5);
  desc.add<double>("trackMinimumEtlEta", 1.6);
  desc.add<double>("trackMaximumEtlEta", 3.);
  desc.add<double>("minProbHeavy", 0.75);
  desc.addUntracked<bool>("printMsg", false);

  descriptions.add("mtdSecondaryPvValid", desc);
}

const bool MtdSecondaryPvValidation::mvaGenSel(const HepMC::GenParticle& gp, const float& charge) {
  bool match = false;
  if (gp.status() != 1) {
    return match;
  }
  match = charge != 0.f && gp.momentum().perp() > pTcut_ && std::abs(gp.momentum().eta()) < etacutGEN_;
  return match;
}

const bool MtdSecondaryPvValidation::mvaRecSel(const reco::TrackBase& trk,
                                               const reco::Vertex& vtx,
                                               const double& t0,
                                               const double& st0) {
  bool match = false;
  match = trk.pt() > pTcut_ && std::abs(trk.eta()) < etacutREC_ && std::abs(trk.vz() - vtx.z()) <= deltaZcut_;
  if (st0 > 0.) {
    match = match && std::abs(t0 - vtx.t()) < 3. * st0;
  }
  return match;
}

const bool MtdSecondaryPvValidation::mvaGenRecMatch(const HepMC::GenParticle& genP,
                                                    const double& zsim,
                                                    const reco::TrackBase& trk) {
  bool match = false;
  double dR = reco::deltaR(genP.momentum(), trk.momentum());
  double genPT = genP.momentum().perp();
  match =
      std::abs(genPT - trk.pt()) < trk.pt() * deltaPTcut_ && dR < deltaDRcut_ && std::abs(trk.vz() - zsim) < deltaZcut_;
  return match;
}

const edm::Ref<std::vector<TrackingParticle>>* MtdSecondaryPvValidation::getMatchedTP(
    const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // no matching or no unique matching
  if (found == r2s_->end()) {
    return nullptr;
  }

  // reco track matched
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0 && tp.first->eventId().event() == 0)
      return &tp.first;
  }

  // no match
  return nullptr;
}

const bool MtdSecondaryPvValidation::genSelBPHkstar(const HepMC::GenParticle& genP) {
  bool match = false;
  if (genP.status() != 1 || (std::abs(genP.pdg_id()) != 211 && std::abs(genP.pdg_id()) != 321)) {
    return match;
  }
  HepMC::GenVertex* orig0 = genP.production_vertex();
  if (orig0->particles_in_size() > 0 && std::abs((*orig0->particles_in_const_begin())->pdg_id()) == 313) {
    HepMC::GenVertex* orig1 = (*orig0->particles_in_const_begin())->production_vertex();
    if (orig1->particles_in_size() > 0 && std::abs((*orig1->particles_in_const_begin())->pdg_id()) == 511) {
      match = true;
      //edm::LogWarning("MtdSecondaryPvValidation") << genP.pdg_id() << " " << orig1->particles_in_size() << " " << std::abs((*orig0->particles_in_const_begin())->pdg_id()) << " " << std::abs((*orig1->particles_in_const_begin())->pdg_id());
    }
  }
  return match;
}

const bool MtdSecondaryPvValidation::genSelBPHmu(const HepMC::GenParticle& genP) {
  bool match = false;
  if (genP.status() != 1 || std::abs(genP.pdg_id()) != 13) {
    return match;
  }
  HepMC::GenVertex* orig0 = genP.production_vertex();
  if (orig0->particles_in_size() > 0 && std::abs((*orig0->particles_in_const_begin())->pdg_id()) == 511) {
    unsigned int count = 0;
    for (HepMC::GenVertex::particles_out_const_iterator outP = orig0->particles_out_const_begin();
         outP != orig0->particles_out_const_end();
         outP++) {
      if (std::abs((*outP)->pdg_id()) == 313 || std::abs((*outP)->pdg_id()) == 13) {
        count++;
      }
    }
    if (count == 3) {
      match = true;
    }
  }
  return match;
}

std::string MtdSecondaryPvValidation::getMuCat(reco::Muon const& muon) {
  std::stringstream muCat;
  muCat.str("");

  if (muon.isGlobalMuon() == true) {
    muCat << " GlobalMuon";
    if (muon::isGoodMuon(muon, muon::GlobalMuonPromptTight) == true)
      muCat << " GlobalMuonPromptTight";
  }
  if (muon.isTrackerMuon() == true) {
    muCat << " TrackerMuon";
    if (muon::isGoodMuon(muon, muon::TrackerMuonArbitrated) == true)
      muCat << " TrackerMuonArbitrated";
    if (muon::isGoodMuon(muon, muon::TMLastStationTight) == true)
      muCat << " TMLastStationTight";
    if (muon::isGoodMuon(muon, muon::TMLastStationLoose) == true)
      muCat << " TMLastStationLoose";
    if (muon::isGoodMuon(muon, muon::TM2DCompatibilityTight) == true)
      muCat << " TM2DCompatibilityTight";
    if (muon::isGoodMuon(muon, muon::TM2DCompatibilityLoose) == true)
      muCat << " TM2DCompatibilityLoose";
    if (muon::isGoodMuon(muon, muon::TMOneStationTight) == true)
      muCat << " TMOneStationTight";
    if (muon::isGoodMuon(muon, muon::TMOneStationLoose) == true)
      muCat << " TMOneStationLoose";
    if (muon::isGoodMuon(muon, muon::TMLastStationAngTight) == true)
      muCat << " TMLastStationAngTight";
    if (muon::isGoodMuon(muon, muon::TMLastStationAngLoose) == true)
      muCat << " TMLastStationAngLoose";
    if (muon::isGoodMuon(muon, muon::TMOneStationAngTight) == true)
      muCat << " TMOneStationAngTight";
    if (muon::isGoodMuon(muon, muon::TMOneStationAngLoose) == true)
      muCat << " TMOneStationAngLoose";
  }
  if (muon.isStandAloneMuon() == true)
    muCat << " StandAloneMuon";
  if (muon.isCaloMuon() == true)
    muCat << " CaloMuon";
  if ((muon.isGlobalMuon() == false) && (muon.isTrackerMuon() == false) && (muon.isStandAloneMuon() == false) &&
      (muon.isCaloMuon() == false))
    muCat << " NotInTable";

  return muCat.str();
}

void MtdSecondaryPvValidation::MonteCarloStudies(const edm::Event& iEvent) {
  if (iEvent.isRealData()) {
    return;
  }

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(prunedGenToken_, pruned);

  //edm::Handle<pat::PackedGenParticleCollection> packed;
  //iEvent.getByToken(packedGenToken_,packed);

  double deltaEtaPhi;

  //     const reco::Candidate* genPsi = NULL;
  const reco::Candidate* genMum = NULL;
  const reco::Candidate* genMup = NULL;
  //     const reco::Candidate* genKst = NULL;
  const reco::Candidate* genTrkm = NULL;
  const reco::Candidate* genTrkp = NULL;
  //  const reco::Candidate* genB0save = NULL;

  bool found_mum = false;
  bool found_mup = false;
  bool found_trkp = false;
  bool found_trkm = false;

  for (const reco::GenParticle& bMeson : *pruned) {
    if (abs(bMeson.pdgId()) == 511) {
      if (skipOscillations(bMeson, pruned))
        continue;
      // select only B0 -> K0* mu mu decays
      if (!genSelB0kstarmumu(bMeson)) {
        continue;
      }
      if (printMsg_)
        std::cout << "PdgID: " << bMeson.pdgId() << " pt " << bMeson.pt() << " eta: " << bMeson.eta()
                  << " phi: " << bMeson.phi() << "mother: " << bMeson.mother(0)->pdgId() << std::endl;

      genMum = NULL;
      genMup = NULL;
      genTrkm = NULL;
      genTrkp = NULL;

      found_mum = false;
      found_mup = false;
      found_trkp = false;
      found_trkm = false;

      //for (const pat::PackedGenParticle &dau : *packed) {
      for (const reco::GenParticle& dau : *pruned) {
        //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
        const reco::Candidate* motherInPrunedCollection = dau.mother(0);
        if (motherInPrunedCollection != nullptr && isAncestor(&bMeson, motherInPrunedCollection)) {
          if (printMsg_)
            std::cout << "     PdgID: " << dau.pdgId() << " pt " << dau.pt() << " eta: " << dau.eta()
                      << " phi: " << dau.phi() << std::endl;

          if (dau.pdgId() == 13 && fabs(dau.eta()) < 2.5 && dau.pt() > 2.5) {
            found_mum = true;
            genMum = &dau;
          } else if (dau.pdgId() == -13 && fabs(dau.eta()) < 2.5 && dau.pt() > 2.5) {
            found_mup = true;
            genMup = &dau;
          } else if ((dau.pdgId() == 211 || dau.pdgId() == 321) && fabs(dau.eta()) < 2.5 && dau.pt() > 0.4) {
            found_trkp = true;
            genTrkp = &dau;
          } else if ((dau.pdgId() == -211 || dau.pdgId() == -321) && fabs(dau.eta()) < 2.5 && dau.pt() > 0.4) {
            found_trkm = true;
            genTrkm = &dau;
          }
        }
      }

      if (found_mup && found_mum && found_trkp && found_trkm) {
        NTuple->genMumPx = genMum->px();
        NTuple->genMumPy = genMum->py();
        NTuple->genMumPz = genMum->pz();

        NTuple->genMupPx = genMup->px();
        NTuple->genMupPy = genMup->py();
        NTuple->genMupPz = genMup->pz();

        NTuple->genKstTrkmID = genTrkm->pdgId();
        NTuple->genKstTrkmPx = genTrkm->px();
        NTuple->genKstTrkmPy = genTrkm->py();
        NTuple->genKstTrkmPz = genTrkm->pz();

        NTuple->genKstTrkpID = genTrkp->pdgId();
        NTuple->genKstTrkpPx = genTrkp->px();
        NTuple->genKstTrkpPy = genTrkp->py();
        NTuple->genKstTrkpPz = genTrkp->pz();

        NTuple->genKstPx = genTrkm->px() + genTrkp->px();
        NTuple->genKstPy = genTrkm->py() + genTrkp->py();
        NTuple->genKstPz = genTrkm->pz() + genTrkp->pz();
        NTuple->genKstMass = Utility->computeInvMass(genTrkm->px(),
                                                     genTrkm->py(),
                                                     genTrkm->pz(),
                                                     genTrkm->mass(),
                                                     genTrkp->px(),
                                                     genTrkp->py(),
                                                     genTrkp->pz(),
                                                     genTrkp->mass());

        if (genTrkp->pdgId() == 321)
          NTuple->genSignal = 1;
        else if (genTrkp->pdgId() == 211)
          NTuple->genSignal = 2;

        NTuple->genB0Mass = bMeson.mass();
        NTuple->genB0Px = bMeson.px();
        NTuple->genB0Py = bMeson.py();
        NTuple->genB0Pz = bMeson.pz();
        NTuple->genB0VtxX = bMeson.vx();
        NTuple->genB0VtxY = bMeson.vy();
        NTuple->genB0VtxZ = bMeson.vz();

        //                 genB0save = &bMeson;
      }
    }
  }

  // ####################################
  // # Perform matching with candidates #
  // ####################################
  if (printMsg_) {
    edm::LogPrint("MtdSecondaryPvValidation") << "nB candidates reconstructed # " << NTuple->nB;
  }
  for (unsigned int i = 0; i < NTuple->nB; i++) {
    deltaEtaPhi = Utility->computeEtaPhiDistance(NTuple->genMumPx,
                                                 NTuple->genMumPy,
                                                 NTuple->genMumPz,
                                                 NTuple->mumPx->at(i),
                                                 NTuple->mumPy->at(i),
                                                 NTuple->mumPz->at(i));
    NTuple->mumDeltaRwithMC->push_back(deltaEtaPhi);
    //if (deltaEtaPhi < RCUTMU) {
    if (NTuple->matchMum->at(i)) {
      NTuple->truthMatchMum->push_back(1);
      if (printMsg_)
        std::cout << __LINE__ << " : found matched mu-" << std::endl;
    } else
      NTuple->truthMatchMum->push_back(0);

    deltaEtaPhi = Utility->computeEtaPhiDistance(NTuple->genMupPx,
                                                 NTuple->genMupPy,
                                                 NTuple->genMupPz,
                                                 NTuple->mupPx->at(i),
                                                 NTuple->mupPy->at(i),
                                                 NTuple->mupPz->at(i));
    NTuple->mupDeltaRwithMC->push_back(deltaEtaPhi);
    //if (deltaEtaPhi < RCUTMU) {
    if (NTuple->matchMup->at(i)) {
      NTuple->truthMatchMup->push_back(1);
      if (printMsg_)
        std::cout << __LINE__ << " : found matched mu+" << std::endl;
    } else
      NTuple->truthMatchMup->push_back(0);

    deltaEtaPhi = Utility->computeEtaPhiDistance(NTuple->genKstTrkmPx,
                                                 NTuple->genKstTrkmPy,
                                                 NTuple->genKstTrkmPz,
                                                 NTuple->kstTrkmPx->at(i),
                                                 NTuple->kstTrkmPy->at(i),
                                                 NTuple->kstTrkmPz->at(i));
    NTuple->kstTrkmDeltaRwithMC->push_back(deltaEtaPhi);
    //if (deltaEtaPhi < RCUTTRK) {
    if (NTuple->matchTkm->at(i)) {
      NTuple->truthMatchTrkm->push_back(1);
      if (printMsg_)
        std::cout << __LINE__ << " : found matched track-" << std::endl;
    } else
      NTuple->truthMatchTrkm->push_back(0);

    deltaEtaPhi = Utility->computeEtaPhiDistance(NTuple->genKstTrkpPx,
                                                 NTuple->genKstTrkpPy,
                                                 NTuple->genKstTrkpPz,
                                                 NTuple->kstTrkpPx->at(i),
                                                 NTuple->kstTrkpPy->at(i),
                                                 NTuple->kstTrkpPz->at(i));
    NTuple->kstTrkpDeltaRwithMC->push_back(deltaEtaPhi);
    //if (deltaEtaPhi < RCUTTRK) {
    if (NTuple->matchTkp->at(i)) {
      NTuple->truthMatchTrkp->push_back(1);
      if (printMsg_)
        std::cout << __LINE__ << " : found matched track+" << std::endl;
    } else
      NTuple->truthMatchTrkp->push_back(0);

    // ####################################################
    // # Check matching with B0 --> track+ track- mu+ mu- #
    // ####################################################
    if ((NTuple->truthMatchTrkm->back() == 1) && (NTuple->truthMatchTrkp->back() == 1) &&
        (NTuple->truthMatchMum->back() == 1) && (NTuple->truthMatchMup->back() == 1)) {
      NTuple->truthMatchSignal->push_back(1);
      if (printMsg_)
        std::cout << __LINE__ << " : @@@ Found matched B0 --> track+ track- mu+ mu- @@@" << std::endl;
    } else
      NTuple->truthMatchSignal->push_back(0);
  }
  //   else
  //     for (unsigned int i = 0; i < NTuple->nB; i++)
  //       {
  //     NTuple->mumDeltaRwithMC->push_back(-1.0);
  //     NTuple->mupDeltaRwithMC->push_back(-1.0);
  //     NTuple->kstTrkmDeltaRwithMC->push_back(-1.0);
  //     NTuple->kstTrkpDeltaRwithMC->push_back(-1.0);
  //
  //     NTuple->truthMatchMum->push_back(0);
  //     NTuple->truthMatchMup->push_back(0);
  //     NTuple->truthMatchTrkm->push_back(0);
  //     NTuple->truthMatchTrkp->push_back(0);
  //     NTuple->truthMatchSignal->push_back(0);
}

bool MtdSecondaryPvValidation::genSelB0kstarmumu(const reco::GenParticle& bMeson) {
  bool match = false;
  unsigned int count(0);
  for (unsigned int i = 0; i < bMeson.numberOfDaughters(); i++) {
    if (std::abs(bMeson.daughter(i)->pdgId()) == 313 || std::abs(bMeson.daughter(i)->pdgId()) == 13) {
      count++;
    }
    //if (printMsg_) {
    //std::cout << "B meson daughther ID " << bMeson.daughter(i)->pdgId() << " count = " << count << std::endl;
    //}
  }
  if (count == 3) {
    match = true;
  }
  return match;
}

bool MtdSecondaryPvValidation::skipOscillations(const reco::GenParticle& bMeson,
                                                edm::Handle<reco::GenParticleCollection> pruned) {
  for (unsigned int i = 0; i < bMeson.numberOfDaughters(); i++) {
    if (bMeson.daughter(i)->pdgId() == 511 || bMeson.daughter(i)->pdgId() == 531 || bMeson.daughter(i)->pdgId() == 5122)
      return true;
    // std::cout << "oscillating to:     PdgID: " << bMeson.daughter(i)->pdgId() << " pt " << bMeson.daughter(i)->pt() << " eta: " << bMeson.daughter(i)->eta() << " phi: " << bMeson.daughter(i)->phi() << std::endl;
  }

  for (const reco::GenParticle& bMother : *pruned) {
    if (fabs(bMother.pdgId()) == 511) {
      const reco::Candidate* mother = bMother.mother(0);
      if (mother != nullptr && isAncestor(&bMeson, mother))
        return true;
    }
  }

  //     if ( fabs(bMeson.mother(0)->pdgId()) == 511) return true;
  //   if (abs(Mother->pdgId()) != 521)
  //     for (unsigned int i = 0; i < Mother->numberOfDaughters(); i++)
  //       if ((abs(Mother->daughter(i)->pdgId()) == 511) || (abs(Mother->daughter(i)->pdgId()) == 531) || (abs(Mother->daughter(i)->pdgId()) == 5122))
  //       {
  //         if (printMsg_) std::cout << __LINE__ << " : @@@ Found oscillating B0/B0bar OR Bs/Bsbar OR Lambda_b/Lambda_bbar @@@" << std::endl;
  //         Mother = Mother->daughter(i);
  //       }

  return false;
}

bool MtdSecondaryPvValidation::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle) {
  //particle is already the ancestor
  if (ancestor == particle)
    return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  for (size_t i = 0; i < particle->numberOfMothers(); i++) {
    if (isAncestor(ancestor, particle->mother(i)))
      return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

const void MtdSecondaryPvValidation::pidFromTrackRef(
    double probpi, double probk, double probp, unsigned int& noPIDtype, bool& noPID, bool& isPi, bool& isK, bool& isP) {
  noPIDtype = 0;
  if (probpi == -1) {
    noPIDtype = 1;
  } else if (isnan(probpi)) {
    noPIDtype = 2;
  } else if (probpi == 1 && probk == 0 && probp == 0) {
    noPIDtype = 3;
  }

  noPID = noPIDtype > 0;
  isPi = !noPID && 1. - probpi < minProbHeavy_;
  isK = !noPID && !isPi && probk > probp;
  isP = !noPID && !isPi && !isK;
}

DEFINE_FWK_MODULE(MtdSecondaryPvValidation);
