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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

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
  const bool genSelBPH(const HepMC::GenParticle&);

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;
  const float minProbHeavy_;

  bool optionalPlots_;

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

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;

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

  MonitorElement* me3GeVSVpid_;
  MonitorElement* me3GeVSV_pi_vs_k_pid_;
  MonitorElement* me3GeVSV_pi_vs_k_uncpid_;
  MonitorElement* me3GeVSV_pi_vs_k_4dpid_;
};

// ------------ constructor and destructor --------------
MtdSecondaryPvValidation::MtdSecondaryPvValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMaxBtlEta_(iConfig.getParameter<double>("trackMaximumBtlEta")),
      trackMinEtlEta_(iConfig.getParameter<double>("trackMinimumEtlEta")),
      trackMaxEtlEta_(iConfig.getParameter<double>("trackMaximumEtlEta")),
      minProbHeavy_(iConfig.getParameter<double>("minProbHeavy")),
      optionalPlots_(iConfig.getUntrackedParameter<bool>("optionalPlots")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));
  HepMCProductToken_ = consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("inputTagH"));
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

  std::vector<reco::Track> candTrk;
  std::vector<reco::TrackRef> candRef;
  std::vector<const HepMC::GenParticle*> candGen;
  std::vector<int> candPid;

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
        if (!genSelBPH(*genP)) {
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

    if (candRef.size() == 2) {
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
    }
  }
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
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdSecondaryPvValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/SecondaryV");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));
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

const bool MtdSecondaryPvValidation::genSelBPH(const HepMC::GenParticle& genP) {
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

DEFINE_FWK_MODULE(MtdSecondaryPvValidation);
