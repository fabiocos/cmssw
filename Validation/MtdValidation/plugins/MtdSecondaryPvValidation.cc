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

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;
  const float minProbHeavy_;

  bool optionalPlots_;
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

  MonitorElement* meMVATrackEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtMtd_;
  MonitorElement* meMVATrackEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaMtd_;
  MonitorElement* meMVATrackResTot_;
  MonitorElement* meMVATrackPullTot_;
  MonitorElement* meMVATrackZposResTot_;

  MonitorElement* meBarrelPiDBetavsp_;
  MonitorElement* meEndcapPiDBetavsp_;
  MonitorElement* meBarrelKDBetavsp_;
  MonitorElement* meEndcapKDBetavsp_;
  MonitorElement* meBarrelPDBetavsp_;
  MonitorElement* meEndcapPDBetavsp_;

  MonitorElement* meBarrelPIDp_;
  MonitorElement* meEndcapPIDp_;

  MonitorElement* meBarrelPID3dip_;
  MonitorElement* meEndcapPID3dip_;

  MonitorElement* meBarrelNoPIDtype_;
  MonitorElement* meEndcapNoPIDtype_;

  MonitorElement* meBarrelTruePiNoPID_;
  MonitorElement* meBarrelTrueKNoPID_;
  MonitorElement* meBarrelTruePNoPID_;
  MonitorElement* meEndcapTruePiNoPID_;
  MonitorElement* meEndcapTrueKNoPID_;
  MonitorElement* meEndcapTruePNoPID_;

  MonitorElement* meBarrelTruePiAsPi_;
  MonitorElement* meBarrelTruePiAsK_;
  MonitorElement* meBarrelTruePiAsP_;
  MonitorElement* meEndcapTruePiAsPi_;
  MonitorElement* meEndcapTruePiAsK_;
  MonitorElement* meEndcapTruePiAsP_;

  MonitorElement* meBarrelTrueKAsPi_;
  MonitorElement* meBarrelTrueKAsK_;
  MonitorElement* meBarrelTrueKAsP_;
  MonitorElement* meEndcapTrueKAsPi_;
  MonitorElement* meEndcapTrueKAsK_;
  MonitorElement* meEndcapTrueKAsP_;

  MonitorElement* meBarrelTruePAsPi_;
  MonitorElement* meBarrelTruePAsK_;
  MonitorElement* meBarrelTruePAsP_;
  MonitorElement* meEndcapTruePAsPi_;
  MonitorElement* meEndcapTruePAsK_;
  MonitorElement* meEndcapTruePAsP_;

  MonitorElement* meSVpid_;
  MonitorElement* meSV_pi_vs_k_pid_;
  MonitorElement* meSV_pi_vs_k_uncpid_;
  MonitorElement* meSV_pi_vs_k_4dpid_;
  MonitorElement* meSV_pi_vs_k_muconst_pid_;

  MonitorElement* me3GeVSVpid_;
  MonitorElement* me3GeVSV_pi_vs_k_pid_;
  MonitorElement* me3GeVSV_pi_vs_k_uncpid_;
  MonitorElement* me3GeVSV_pi_vs_k_4dpid_;
  MonitorElement* me3GeVSV_pi_vs_k_muconst_pid_;

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
      optionalPlots_(iConfig.getUntrackedParameter<bool>("optionalPlots")),
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

  auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
  const MTDTopology* topology = topologyHandle.product();

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

  unsigned int index = 0;

  const auto& primRecoVtx = *(RecVertexHandle.product()->begin());
  double treco = primRecoVtx.t();

  // generator level information (HepMC format)
  auto GenEventHandle = makeValid(iEvent.getHandle(HepMCProductToken_));
  const HepMC::GenEvent* mc = GenEventHandle->GetEvent();
  double zsim = convertMmToCm((*(mc->vertices_begin()))->position().z());
  double tsim = (*(mc->vertices_begin()))->position().t() * CLHEP::mm / CLHEP::c_light;

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
                edm::LogPrint("MtdSecondaryPvValidation")
                    << "Mu- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                    << genP->momentum().eta() << " phi " << genP->momentum().phi();
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
                edm::LogPrint("MtdSecondaryPvValidation")
                    << "Mu- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                    << genP->momentum().eta() << " phi " << genP->momentum().phi();
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
                edm::LogPrint("MtdSecondaryPvValidation")
                    << "Trk- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                    << genP->momentum().eta() << " phi " << genP->momentum().phi();
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
                edm::LogPrint("MtdSecondaryPvValidation")
                    << "Trk- matched " << genP->pdg_id() << " pt " << genP->momentum().perp() << " eta "
                    << genP->momentum().eta() << " phi " << genP->momentum().phi();
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

  /*
  // select events with reco vertex close to true simulated primary vertex
  if (std::abs(primRecoVtx.z() - zsim) < deltaZcut_) {
    index = 0;
    for (const auto& trackGen : *GenRecTrackHandle) {
      const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
      index++;

      // select the reconstructed track

      if (trackAssoc[trackref] == -1) {
        continue;
      }

      bool noCrack = std::abs(trackGen.eta()) < trackMaxBtlEta_ || std::abs(trackGen.eta()) > trackMinEtlEta_;

      // reco-gen matching used for MVA quality flag
      if (mvaRecSel(trackGen, primRecoVtx, t0Safe[trackref], Sigmat0Safe[trackref])) {
        if (noCrack) {
          meMVATrackEffPtTot_->Fill(trackGen.pt());
        }
        meMVATrackEffEtaTot_->Fill(std::abs(trackGen.eta()));

        double dZ = trackGen.vz() - zsim;
        double dT(-9999.);
        double pullT(-9999.);
        if (Sigmat0Safe[trackref] != -1.) {
          dT = t0Safe[trackref] - tsim;
          pullT = dT / Sigmat0Safe[trackref];
        }
        for (const auto& genP : mc->particle_range()) {
          // select status 1 genParticles and match them to the reconstructed track

          float charge = pdTable->particle(HepPDT::ParticleID(genP->pdg_id())) != nullptr
                             ? pdTable->particle(HepPDT::ParticleID(genP->pdg_id()))->charge()
                             : 0.f;
          if (mvaGenSel(*genP, charge)) {
            if (mvaGenRecMatch(*genP, zsim, trackGen)) {
              meMVATrackZposResTot_->Fill(dZ);
              if (noCrack) {
                meMVATrackMatchedEffPtTot_->Fill(trackGen.pt());
              }
              meMVATrackMatchedEffEtaTot_->Fill(std::abs(trackGen.eta()));
              if (pullT > -9999.) {
                meMVATrackResTot_->Fill(dT);
                meMVATrackPullTot_->Fill(pullT);
                if (noCrack) {
                  meMVATrackMatchedEffPtMtd_->Fill(trackGen.pt());
                }
                meMVATrackMatchedEffEtaMtd_->Fill(std::abs(trackGen.eta()));
              }
              break;
            }
          }
        }
      }

      // for PID study select only high purity tracks with associated time information with good MVA quality

      if (!trackGen.quality(reco::TrackBase::TrackQuality::highPurity) || Sigmat0Safe[trackref] == -1. ||
          mtdQualMVA[trackref] < mvaSel_) {
        continue;
      }

      reco::TrackBaseRef tbrTrk(trackref);
      auto tp_info = getMatchedTP(tbrTrk);
      if (tp_info != nullptr) {
        // select for BPH
        if ((*tp_info)->g4Tracks()[0].genpartIndex() == -1) {
          continue;
        }
        const HepMC::GenParticle* genP = mc->barcode_to_particle((*tp_info)->g4Tracks()[0].genpartIndex());
        if (genSelBPHmu(*genP)) {
          candTrkMu.push_back(trackGen);
          candRefMu.push_back(trackref);
          candGenMu.push_back(genP);
          continue;
        }
        if (!genSelBPHkstar(*genP)) {
          continue;
        }

        candTrk.push_back(trackGen);
        candRef.push_back(trackref);
        candGen.push_back(genP);
        double dx = genP->production_vertex()->position().x() - (*(mc->vertices_begin()))->position().x();
        double dy = genP->production_vertex()->position().y() - (*(mc->vertices_begin()))->position().y();
        double dz = genP->production_vertex()->position().z() - (*(mc->vertices_begin()))->position().z();
        double ip3d = convertMmToCm(std::sqrt(dx * dx + dy * dy + dz * dz));

        double dbetaPi = c_cm_ns * (tMtd[trackref] - treco - tofPi[trackref]) / pathLength[trackref];
        double dbetaK = c_cm_ns * (tMtd[trackref] - treco - tofK[trackref]) / pathLength[trackref];
        double dbetaP = c_cm_ns * (tMtd[trackref] - treco - tofP[trackref]) / pathLength[trackref];

        unsigned int noPIDtype = 0;
        if (probPi[trackref] == -1) {
          noPIDtype = 1;
        } else if (isnan(probPi[trackref])) {
          noPIDtype = 2;
        } else if (probPi[trackref] == 1 && probK[trackref] == 0 && probP[trackref] == 0) {
          noPIDtype = 3;
        }
        bool noPID = noPIDtype > 0;
        bool isPi = !noPID && 1. - probPi[trackref] < minProbHeavy_;
        bool isK = !noPID && !isPi && probK[trackref] > probP[trackref];
        bool isP = !noPID && !isPi && !isK;

        if ((isPi && std::abs(tMtd[trackref] - tofPi[trackref] - t0Pid[trackref]) > tol_) ||
            (isK && std::abs(tMtd[trackref] - tofK[trackref] - t0Pid[trackref]) > tol_) ||
            (isP && std::abs(tMtd[trackref] - tofP[trackref] - t0Pid[trackref]) > tol_)) {
          edm::LogWarning("MtdSecondaryPvValidation")
              << "No match between mass hyp. and time: " << std::abs((*tp_info)->pdgId()) << " mass hyp pi/k/p " << isPi
              << " " << isK << " " << isP << " t0/t0safe " << t0Pid[trackref] << " " << t0Safe[trackref]
              << " tMtd - tof pi/K/p " << tMtd[trackref] - tofPi[trackref] << " " << tMtd[trackref] - tofK[trackref]
              << " " << tMtd[trackref] - tofP[trackref] << " Prob pi/K/p " << probPi[trackref] << " " << probK[trackref]
              << " " << probP[trackref];
        }

        //if ( noPID || isPi ) {
        if (isPi) {
          candPid.push_back(0);
        } else if (isK) {
          candPid.push_back(1);
        } else if (isP) {
          candPid.push_back(2);
        } else {
          candPid.push_back(3);
        }

        if (std::abs(trackGen.eta()) < trackMaxBtlEta_) {
          meBarrelPIDp_->Fill(trackGen.p());
          meBarrelPID3dip_->Fill(ip3d, trackGen.p());
          meBarrelNoPIDtype_->Fill(noPIDtype + 0.5);
          if (optionalPlots_) {
            if (std::abs((*tp_info)->pdgId()) == 211) {
              meBarrelPiDBetavsp_->Fill(trackGen.p(), dbetaPi);
              if (noPID) {
                meBarrelTruePiNoPID_->Fill(trackGen.p());
              } else if (isPi) {
                meBarrelTruePiAsPi_->Fill(trackGen.p());
              } else if (isK) {
                meBarrelTruePiAsK_->Fill(trackGen.p());
              } else if (isP) {
                meBarrelTruePiAsP_->Fill(trackGen.p());
              } else {
                edm::LogWarning("MtdSecondaryPvValidation")
                    << "No PID class: " << std::abs((*tp_info)->pdgId()) << " t0/t0safe " << t0Pid[trackref] << " "
                    << t0Safe[trackref] << " Prob pi/K/p " << probPi[trackref] << " " << probK[trackref] << " "
                    << probP[trackref];
              }
            } else if (std::abs((*tp_info)->pdgId()) == 321) {
              meBarrelKDBetavsp_->Fill(trackGen.p(), dbetaK);
              if (noPID) {
                meBarrelTrueKNoPID_->Fill(trackGen.p());
              } else if (isPi) {
                meBarrelTrueKAsPi_->Fill(trackGen.p());
              } else if (isK) {
                meBarrelTrueKAsK_->Fill(trackGen.p());
              } else if (isP) {
                meBarrelTrueKAsP_->Fill(trackGen.p());
              } else {
                edm::LogWarning("MtdSecondaryPvValidation")
                    << "No PID class: " << std::abs((*tp_info)->pdgId()) << " t0/t0safe " << t0Pid[trackref] << " "
                    << t0Safe[trackref] << " Prob pi/K/p " << probPi[trackref] << " " << probK[trackref] << " "
                    << probP[trackref];
              }
            } else if (std::abs((*tp_info)->pdgId()) == 2212) {
              meBarrelPDBetavsp_->Fill(trackGen.p(), dbetaP);
              if (noPID) {
                meBarrelTruePNoPID_->Fill(trackGen.p());
              } else if (isPi) {
                meBarrelTruePAsPi_->Fill(trackGen.p());
              } else if (isK) {
                meBarrelTruePAsK_->Fill(trackGen.p());
              } else if (isP) {
                meBarrelTruePAsP_->Fill(trackGen.p());
              } else {
                edm::LogWarning("MtdSecondaryPvValidation")
                    << "No PID class: " << std::abs((*tp_info)->pdgId()) << " t0/t0safe " << t0Pid[trackref] << " "
                    << t0Safe[trackref] << " Prob pi/K/p " << probPi[trackref] << " " << probK[trackref] << " "
                    << probP[trackref];
              }
            }
          }
        } else if (std::abs(trackGen.eta()) > trackMinEtlEta_ && std::abs(trackGen.eta()) < trackMaxEtlEta_) {
          meEndcapPIDp_->Fill(trackGen.p());
          meEndcapPID3dip_->Fill(ip3d, trackGen.p());
          meEndcapNoPIDtype_->Fill(noPIDtype + 0.5);
          if (optionalPlots_) {
            if (std::abs((*tp_info)->pdgId()) == 211) {
              meEndcapPiDBetavsp_->Fill(trackGen.p(), dbetaPi);
              if (noPID) {
                meEndcapTruePiNoPID_->Fill(trackGen.p());
              } else if (isPi) {
                meEndcapTruePiAsPi_->Fill(trackGen.p());
              } else if (isK) {
                meEndcapTruePiAsK_->Fill(trackGen.p());
              } else if (isP) {
                meEndcapTruePiAsP_->Fill(trackGen.p());
              } else {
                edm::LogWarning("MtdSecondaryPvValidation")
                    << "No PID class: " << std::abs((*tp_info)->pdgId()) << " t0/t0safe " << t0Pid[trackref] << " "
                    << t0Safe[trackref] << " Prob pi/K/p " << probPi[trackref] << " " << probK[trackref] << " "
                    << probP[trackref];
              }
            } else if (std::abs((*tp_info)->pdgId()) == 321) {
              meEndcapKDBetavsp_->Fill(trackGen.p(), dbetaK);
              if (noPID) {
                meEndcapTrueKNoPID_->Fill(trackGen.p());
              } else if (isPi) {
                meEndcapTrueKAsPi_->Fill(trackGen.p());
              } else if (isK) {
                meEndcapTrueKAsK_->Fill(trackGen.p());
              } else if (isP) {
                meEndcapTrueKAsP_->Fill(trackGen.p());
              } else {
                edm::LogWarning("MtdSecondaryPvValidation")
                    << "No PID class: " << std::abs((*tp_info)->pdgId()) << " t0/t0safe " << t0Pid[trackref] << " "
                    << t0Safe[trackref] << " Prob pi/K/p " << probPi[trackref] << " " << probK[trackref] << " "
                    << probP[trackref];
              }
            } else if (std::abs((*tp_info)->pdgId()) == 2212) {
              meEndcapPDBetavsp_->Fill(trackGen.p(), dbetaP);
              if (noPID) {
                meEndcapTruePNoPID_->Fill(trackGen.p());
              } else if (isPi) {
                meEndcapTruePAsPi_->Fill(trackGen.p());
              } else if (isK) {
                meEndcapTruePAsK_->Fill(trackGen.p());
              } else if (isP) {
                meEndcapTruePAsP_->Fill(trackGen.p());
              } else {
                edm::LogWarning("MtdSecondaryPvValidation")
                    << "No PID class: " << std::abs((*tp_info)->pdgId()) << " t0/t0safe " << t0Pid[trackref] << " "
                    << t0Safe[trackref] << " Prob pi/K/p " << probPi[trackref] << " " << probK[trackref] << " "
                    << probP[trackref];
              }
            }
          }
        }
      }
    }

    if (candRef.size() == 2 && candRefMu.size() == 2) {
      bool less3GeV = candTrk[0].p() < 3. && candTrk[1].p() < 3.;

      // array on candidates, for each array of hyptheses (0 = pi, 1 = k)
      std::array<std::array<double, 3>, 2> hypTof;
      std::array<std::array<double, 3>, 2> hypUncTof;
      std::array<unsigned int, 2> hypOk;

      for (unsigned int index = 0; index < 2; index++) {
        double dx =
            convertMmToCm(candGen[index]->production_vertex()->position().x()) - candTrk[index].referencePoint().x();
        double dy =
            convertMmToCm(candGen[index]->production_vertex()->position().y()) - candTrk[index].referencePoint().y();
        double dz =
            convertMmToCm(candGen[index]->production_vertex()->position().z()) - candTrk[index].referencePoint().z();
        double dscorr = std::sqrt(dx * dx + dy * dy + dz * dz);
        double gammasq_pi = 1. + candTrk[index].p() * candTrk[index].p() * m_pi_inv2;
        double beta_pi = std::sqrt(1. - 1. / gammasq_pi);
        double dtcorr_pi = dscorr / beta_pi * c_inv;
        double gammasq_k = 1. + candTrk[index].p() * candTrk[index].p() * m_k_inv2;
        double beta_k = std::sqrt(1. - 1. / gammasq_k);
        double dtcorr_k = dscorr / beta_k * c_inv;
        double gammasq_p = 1. + candTrk[index].p() * candTrk[index].p() * m_p_inv2;
        double beta_p = std::sqrt(1. - 1. / gammasq_p);
        double dtcorr_p = dscorr / beta_p * c_inv;
        std::array<double, 3> tof;
        tof[0] = tMtd[candRef[index]] - tofPi[candRef[index]] + dtcorr_pi;
        tof[1] = tMtd[candRef[index]] - tofK[candRef[index]] + dtcorr_k;
        tof[2] = tMtd[candRef[index]] - tofP[candRef[index]] + dtcorr_p;
        hypTof[index] = tof;
        tof[0] = tMtd[candRef[index]] - tofPi[candRef[index]];
        tof[1] = tMtd[candRef[index]] - tofK[candRef[index]];
        tof[2] = tMtd[candRef[index]] - tofP[candRef[index]];
        hypUncTof[index] = tof;
        if (std::abs(candGen[index]->pdg_id()) == 211) {
          hypOk[index] = 0;
        } else if (std::abs(candGen[index]->pdg_id()) == 321) {
          hypOk[index] = 1;
        } else if (std::abs(candGen[index]->pdg_id()) == 2212) {
          hypOk[index] = 2;
        } else {
          edm::LogWarning("MtdSecondaryPvValid") << "Not expected particle " << candGen[index]->pdg_id();
        }
      }

      double sigmaDT =
          std::sqrt(SigmatMtd[candRef[0]] * SigmatMtd[candRef[0]] + SigmatMtd[candRef[1]] * SigmatMtd[candRef[1]]);
      double sigdtAtSV = 9999.;
      unsigned int id1(3);
      unsigned int id2(3);
      bool idOK = false;
      for (unsigned int ihyp1 = 0; ihyp1 < 3; ihyp1++) {
        for (unsigned int ihyp2 = 0; ihyp2 < 3; ihyp2++) {
          double sigdtTmp = std::abs((hypTof[0])[ihyp1] - (hypTof[1])[ihyp2]) / sigmaDT;
          bool idOkTmp = (hypOk[0] == ihyp1 && hypOk[1] == ihyp2);
          if (sigdtTmp < sigdtAtSV) {
            sigdtAtSV = sigdtTmp;
            idOK = idOkTmp;
            id1 = ihyp1;
            id2 = ihyp2;
          }
        }
      }

      if (idOK) {
        meSVpid_->Fill(0.5);
        if (less3GeV) {
          me3GeVSVpid_->Fill(0.5);
        }
      } else {
        meSVpid_->Fill(1.5);
        if (less3GeV) {
          me3GeVSVpid_->Fill(1.5);
        }
      }
      if (std::abs(candGen[0]->pdg_id()) == 211) {
        meSV_pi_vs_k_pid_->Fill(id2 + 0.5, id1 + 0.5);
        if (less3GeV) {
          me3GeVSV_pi_vs_k_pid_->Fill(id2 + 0.5, id1 + 0.5);
        }
      } else {
        meSV_pi_vs_k_pid_->Fill(id1 + 0.5, id2 + 0.5);
        if (less3GeV) {
          me3GeVSV_pi_vs_k_pid_->Fill(id1 + 0.5, id2 + 0.5);
        }
      }

      sigdtAtSV = 9999.;
      id1 = 3;
      id2 = 3;
      idOK = false;
      for (unsigned int ihyp1 = 0; ihyp1 < 3; ihyp1++) {
        for (unsigned int ihyp2 = 0; ihyp2 < 3; ihyp2++) {
          double sigdtTmp = std::abs((hypUncTof[0])[ihyp1] - (hypUncTof[1])[ihyp2]) / sigmaDT;
          bool idOkTmp = (hypOk[0] == ihyp1 && hypOk[1] == ihyp2);
          if (sigdtTmp < sigdtAtSV) {
            sigdtAtSV = sigdtTmp;
            idOK = idOkTmp;
            id1 = ihyp1;
            id2 = ihyp2;
          }
        }
      }

      if (std::abs(candGen[0]->pdg_id()) == 211) {
        meSV_pi_vs_k_uncpid_->Fill(id2 + 0.5, id1 + 0.5);
        if (less3GeV) {
          me3GeVSV_pi_vs_k_uncpid_->Fill(id2 + 0.5, id1 + 0.5);
        }
      } else {
        meSV_pi_vs_k_uncpid_->Fill(id1 + 0.5, id2 + 0.5);
        if (less3GeV) {
          me3GeVSV_pi_vs_k_uncpid_->Fill(id1 + 0.5, id2 + 0.5);
        }
      }

      id1 = candPid[0];
      id2 = candPid[1];
      if (std::abs(candGen[0]->pdg_id()) == 211) {
        meSV_pi_vs_k_4dpid_->Fill(id2 + 0.5, id1 + 0.5);
        if (less3GeV) {
          me3GeVSV_pi_vs_k_4dpid_->Fill(id2 + 0.5, id1 + 0.5);
        }
      } else {
        meSV_pi_vs_k_4dpid_->Fill(id1 + 0.5, id2 + 0.5);
        if (less3GeV) {
          meSV_pi_vs_k_4dpid_->Fill(id1 + 0.5, id2 + 0.5);
        }
      }

      sigdtAtSV = 9999.;
      id1 = 3;
      id2 = 3;
      double w1 = 1. / SigmatMtd[candRefMu[0]] / SigmatMtd[candRefMu[0]];
      double t1 = tMtd[candRefMu[0]] - tofPi[candRefMu[0]];
      double w2 = 1. / SigmatMtd[candRefMu[1]] / SigmatMtd[candRefMu[1]];
      double t2 = tMtd[candRefMu[1]] - tofPi[candRefMu[1]];
      for (unsigned int ihyp1 = 0; ihyp1 < 3; ihyp1++) {
        for (unsigned int ihyp2 = 0; ihyp2 < 3; ihyp2++) {
          double w3 = 1. / SigmatMtd[candRef[0]] / SigmatMtd[candRef[0]];
          double t3 = (hypTof[0])[ihyp1];
          double w4 = 1. / SigmatMtd[candRef[1]] / SigmatMtd[candRef[1]];
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

      if (std::abs(candGen[0]->pdg_id()) == 211) {
        meSV_pi_vs_k_muconst_pid_->Fill(id2 + 0.5, id1 + 0.5);
        if (less3GeV) {
          me3GeVSV_pi_vs_k_muconst_pid_->Fill(id2 + 0.5, id1 + 0.5);
        }
      } else {
        meSV_pi_vs_k_muconst_pid_->Fill(id1 + 0.5, id2 + 0.5);
        if (less3GeV) {
          me3GeVSV_pi_vs_k_muconst_pid_->Fill(id1 + 0.5, id2 + 0.5);
        }
      }
    }
  }
  */

  NTuple->ClearNTuple();
}

// ------------ method for histogram booking ------------
void MtdSecondaryPvValidation::bookHistograms(DQMStore::IBooker& ibook,
                                              edm::Run const& run,
                                              edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking
  meMVATrackEffPtTot_ = ibook.book1D("MVAEffPtTot", "Pt of tracks associated to LV; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtTot_ =
      ibook.book1D("MVAMatchedEffPtTot", "Pt of tracks associated to LV matched to GEN; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtMtd_ = ibook.book1D(
      "MVAMatchedEffPtMtd", "Pt of tracks associated to LV matched to GEN with time; track pt [GeV] ", 110, 0., 11.);
  meMVATrackEffEtaTot_ = ibook.book1D("MVAEffEtaTot", "Pt of tracks associated to LV; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaTot_ =
      ibook.book1D("MVAMatchedEffEtaTot", "Pt of tracks associated to LV matched to GEN; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaMtd_ = ibook.book1D(
      "MVAMatchedEffEtaMtd", "Pt of tracks associated to LV matched to GEN with time; track eta ", 66, 0., 3.3);
  meMVATrackResTot_ = ibook.book1D(
      "MVATrackRes", "t_{rec} - t_{sim} for LV associated tracks; t_{rec} - t_{sim} [ns] ", 120, -0.15, 0.15);
  meMVATrackPullTot_ =
      ibook.book1D("MVATrackPull", "Pull for associated tracks; (t_{rec}-t_{sim})/#sigma_{t}", 50, -5., 5.);
  meMVATrackZposResTot_ = ibook.book1D(
      "MVATrackZposResTot", "Z_{PCA} - Z_{sim} for associated tracks;Z_{PCA} - Z_{sim} [cm] ", 100, -0.1, 0.1);

  meBarrelPIDp_ = ibook.book1D("BarrelPIDp", "PID track with MTD momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
  meEndcapPIDp_ = ibook.book1D("EndcapPIDp", "PID track with MTD momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);

  meBarrelPID3dip_ = ibook.book2D(
      "BarrelPID3dip", "PID track with MTD 3d p vs ip, |eta| < 1.5;ip [cm]; p [GeV]", 60, 0., 2., 25, 0., 10.);
  meEndcapPID3dip_ = ibook.book2D(
      "EndcapPID3dip", "PID track with MTD 3d p vs ip, |eta| > 1.6;ip [cm]; p [GeV]", 60, 0., 2., 25, 0., 10.);

  meBarrelNoPIDtype_ = ibook.book1D("BarrelNoPIDtype", "Barrel PID failure category", 4, 0., 4.);
  meEndcapNoPIDtype_ = ibook.book1D("EndcapNoPIDtype", "Endcap PID failure category", 4, 0., 4.);

  if (optionalPlots_) {
    meBarrelPiDBetavsp_ = ibook.book2D(
        "BarrelPiDBetavsp", "DeltaBeta true pi as pi vs p, |eta| < 1.5;p [GeV]; dBeta", 25, 0., 10., 50, -0.1, 0.1);
    meEndcapPiDBetavsp_ = ibook.book2D(
        "EndcapPiDBetavsp", "DeltaBeta true pi as pi vs p, |eta| > 1.6;p [GeV]; dBeta", 25, 0., 10., 50, -0.1, 0.1);
    meBarrelKDBetavsp_ = ibook.book2D(
        "BarrelKDBetavsp", "DeltaBeta true K as K vs p, |eta| < 1.5;p [GeV]; dBeta", 25, 0., 10., 50, -0.1, 0.1);
    meEndcapKDBetavsp_ = ibook.book2D(
        "EndcapKDBetavsp", "DeltaBeta true K as K vs p, |eta| > 1.6;p [GeV]; dBeta", 25, 0., 10., 50, -0.1, 0.1);
    meBarrelPDBetavsp_ = ibook.book2D(
        "BarrelPDBetavsp", "DeltaBeta true p as p vs p, |eta| < 1.5;p [GeV]; dBeta", 25, 0., 10., 50, -0.1, 0.1);
    meEndcapPDBetavsp_ = ibook.book2D(
        "EndcapPDBetavsp", "DeltaBeta true p as p vs p, |eta| > 1.6;p [GeV]; dBeta", 25, 0., 10., 50, -0.1, 0.1);

    meBarrelTruePiNoPID_ =
        ibook.book1D("BarrelTruePiNoPID", "True pi NoPID momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTrueKNoPID_ =
        ibook.book1D("BarrelTrueKNoPID", "True k NoPID momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTruePNoPID_ =
        ibook.book1D("BarrelTruePNoPID", "True p NoPID momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meEndcapTruePiNoPID_ =
        ibook.book1D("EndcapTruePiNoPID", "True NoPIDpi momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTrueKNoPID_ =
        ibook.book1D("EndcapTrueKNoPID", "True k NoPID momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTruePNoPID_ =
        ibook.book1D("EndcapTruePNoPID", "True p NoPID momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);

    meBarrelTruePiAsPi_ =
        ibook.book1D("BarrelTruePiAsPi", "True pi as pi momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTruePiAsK_ =
        ibook.book1D("BarrelTruePiAsK", "True pi as k momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTruePiAsP_ =
        ibook.book1D("BarrelTruePiAsP", "True pi as p momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meEndcapTruePiAsPi_ =
        ibook.book1D("EndcapTruePiAsPi", "True pi as pi momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTruePiAsK_ =
        ibook.book1D("EndcapTruePiAsK", "True pi as k momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTruePiAsP_ =
        ibook.book1D("EndcapTruePiAsP", "True pi as p momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);

    meBarrelTrueKAsPi_ =
        ibook.book1D("BarrelTrueKAsPi", "True k as pi momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTrueKAsK_ =
        ibook.book1D("BarrelTrueKAsK", "True k as k momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTrueKAsP_ =
        ibook.book1D("BarrelTrueKAsP", "True k as p momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meEndcapTrueKAsPi_ =
        ibook.book1D("EndcapTrueKAsPi", "True k as pi momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTrueKAsK_ =
        ibook.book1D("EndcapTrueKAsK", "True k as k momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTrueKAsP_ =
        ibook.book1D("EndcapTrueKAsP", "True k as p momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);

    meBarrelTruePAsPi_ =
        ibook.book1D("BarrelTruePAsPi", "True p as pi momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTruePAsK_ =
        ibook.book1D("BarrelTruePAsK", "True p as k momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meBarrelTruePAsP_ =
        ibook.book1D("BarrelTruePAsP", "True p as p momentum spectrum, |eta| < 1.5;p [GeV]", 25, 0., 10.);
    meEndcapTruePAsPi_ =
        ibook.book1D("EndcapTruePAsPi", "True p as pi momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTruePAsK_ =
        ibook.book1D("EndcapTruePAsK", "True p as k momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);
    meEndcapTruePAsP_ =
        ibook.book1D("EndcapTruePAsP", "True p as p momentum spectrum, |eta| > 1.6;p [GeV]", 25, 0., 10.);

    meSVpid_ = ibook.book1D("SVpid", "SV correct identification, 0.5 = OK", 2, 0., 2.);
    meSV_pi_vs_k_pid_ =
        ibook.book2D("SV_pi_vs_k_pid", "SV pi vs k correct identification, pi/k/p", 3, 0., 3., 3, 0., 3.);
    meSV_pi_vs_k_uncpid_ =
        ibook.book2D("SV_pi_vs_k_uncpid", "SV pi vs k correct identification no Sv corr, pi/k/p", 3, 0., 3., 3, 0., 3.);
    meSV_pi_vs_k_4dpid_ =
        ibook.book2D("SV_pi_vs_k_4dpid", "SV pi vs k correct identification 4D corr, pi/k/p", 3, 0., 3., 3, 0., 3.);
    meSV_pi_vs_k_muconst_pid_ = ibook.book2D(
        "SV_pi_vs_k_muconst_pid", "SV pi vs k correct identification mu constraint, pi/k/p", 3, 0., 3., 3, 0., 3.);

    me3GeVSVpid_ = ibook.book1D("3GeVSVpid", "Sv correct identification, tracks p < 3 GeV, 0.5 = OK", 2, 0., 2.);
    me3GeVSV_pi_vs_k_pid_ = ibook.book2D(
        "3GeVSV_pi_vs_k_pid", "Sv pi vs k correct identification, tracks p < 3 GeV, pi/k/p", 3, 0., 3., 3, 0., 3.);
    me3GeVSV_pi_vs_k_uncpid_ = ibook.book2D("3GeVSV_pi_vs_k_uncpid",
                                            "Sv pi vs k correct identification no Sv corr, tracks p < 3 GeV, pi/k/p",
                                            3,
                                            0.,
                                            3.,
                                            3,
                                            0.,
                                            3.);
    me3GeVSV_pi_vs_k_4dpid_ = ibook.book2D("3GeVSV_pi_vs_k_4dpid",
                                           "Sv pi vs k correct identification 4D corr, tracks p < 3 GeV, pi/k/p",
                                           3,
                                           0.,
                                           3.,
                                           3,
                                           0.,
                                           3.);
    me3GeVSV_pi_vs_k_muconst_pid_ =
        ibook.book2D("3GeVSV_pi_vs_k_muconst_pid",
                     "Sv pi vs k correct identification mu constraint, tracks p < 3 GeV, pi/k/p",
                     3,
                     0.,
                     3.,
                     3,
                     0.,
                     3.);
  }
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
  desc.addUntracked<bool>("optionalPlots", false);
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

DEFINE_FWK_MODULE(MtdSecondaryPvValidation);
