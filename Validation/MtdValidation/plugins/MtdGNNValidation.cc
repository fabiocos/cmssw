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
#include "DataFormats/Math/interface/angle_units.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToRecoClusterAssociationMap.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "MTDHit.h"

class MtdGNNValidation : public DQMEDAnalyzer {
public:
  explicit MtdGNNValidation(const edm::ParameterSet&);
  ~MtdGNNValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const edm::Ref<std::vector<TrackingParticle>>* getMatchedTP(const reco::TrackBaseRef&);

  // ------------ member data ------------

  const std::string folder_;
  const bool optionalPlots_;
  const float trackMaxPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;

  static constexpr double simUnit_ = 1e9;                // sim time in s while reco time in ns
  static constexpr double etacutGEN_ = 4.;               // |eta| < 4;
  static constexpr double etacutREC_ = 3.;               // |eta| < 3;
  static constexpr double pTcutBTL_ = 0.7;               // PT > 0.7 GeV
  static constexpr double pTcutETL_ = 0.2;               // PT > 0.2 GeV
  static constexpr double depositBTLthreshold_ = 1;      // threshold for energy deposit in BTL cell [MeV]
  static constexpr double depositETLthreshold_ = 0.001;  // threshold for energy deposit in ETL cell [MeV]
  static constexpr double rBTL_ = 110.0;
  static constexpr double zETL_ = 290.0;
  static constexpr double etaMatchCut_ = 0.05;
  static constexpr double cluDRradius_ = 0.05;  // to cluster rechits around extrapolated track

  const reco::RecoToSimCollection* r2s_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::TPToSimCollectionMtd> tp2SimAssociationMapToken_;
  edm::EDGetTokenT<reco::SimToTPCollectionMtd> Sim2tpAssociationMapToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;

  edm::EDGetTokenT<FTLRecHitCollection> btlRecHitsToken_;
  edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;
  edm::EDGetTokenT<FTLClusterCollection> btlRecCluToken_;
  edm::EDGetTokenT<FTLClusterCollection> etlRecCluToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> outermostHitPositionToken_;

  edm::EDGetTokenT<edm::ValueMap<float>>  betaTok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  phiTok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  logitPiTok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  logitKTok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  logitPTok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  emb0Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  emb1Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  emb2Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  pca0Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  pca1Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>>  pca2Tok_;

  // histogram declaration
  //
  MonitorElement* meVtxVsZ_;
  MonitorElement* meVtxSpreadVsZ_;
  MonitorElement* meVtxVsPC0_;
  MonitorElement* meVtxSpreadVsPC0_;
};

// ------------ constructor and destructor --------------
MtdGNNValidation::MtdGNNValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      optionalPlots_(iConfig.getParameter<bool>("optionalPlots")),
      trackMaxPt_(iConfig.getParameter<double>("trackMaximumPt")),
      trackMaxBtlEta_(iConfig.getParameter<double>("trackMaximumBtlEta")),
      trackMinEtlEta_(iConfig.getParameter<double>("trackMinimumEtlEta")),
      trackMaxEtlEta_(iConfig.getParameter<double>("trackMaximumEtlEta")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));

  trackingVertexCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  tp2SimAssociationMapToken_ =
      consumes<reco::TPToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMapTag"));
  Sim2tpAssociationMapToken_ =
      consumes<reco::SimToTPCollectionMtd>(iConfig.getParameter<edm::InputTag>("Sim2tpAssociationMapTag"));
  r2sAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(
      iConfig.getParameter<edm::InputTag>("r2sAssociationMapTag"));
  btlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("btlRecHits"));
  etlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("etlRecHits"));
  btlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagBTL"));
  etlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagETL"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  btlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2"));
  etlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2"));
  btlMatchChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchChi2"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  SigmaTofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofPi"));
  SigmaTofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofK"));
  SigmaTofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofP"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  outermostHitPositionToken_ =
      consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("outermostHitPositionSrc"));

  auto mkTag = [&](const std::string& instance)->edm::InputTag{
    return edm::InputTag(pvModule_.label(), instance, pvModule_.process());
  };

  betaTok_     = consumes<edm::ValueMap<float>>( mkTag("gnnBeta") );
  phiTok_      = consumes<edm::ValueMap<float>>( mkTag("gnnPhi") );
  logitPiTok_  = consumes<edm::ValueMap<float>>( mkTag("gnnPidLogitPi") );
  logitKTok_   = consumes<edm::ValueMap<float>>( mkTag("gnnPidLogitK") );
  logitPTok_   = consumes<edm::ValueMap<float>>( mkTag("gnnPidLogitP") );
  emb0Tok_     = consumes<edm::ValueMap<float>>( mkTag("gnnEmb0") );
  emb1Tok_     = consumes<edm::ValueMap<float>>( mkTag("gnnEmb1") );
  emb2Tok_     = consumes<edm::ValueMap<float>>( mkTag("gnnEmb2") );
  pca0Tok_     = consumes<edm::ValueMap<float>>( mkTag("gnnPCA0") );
  pca1Tok_     = consumes<edm::ValueMap<float>>( mkTag("gnnPCA1") );
  pca2Tok_     = consumes<edm::ValueMap<float>>( mkTag("gnnPCA2") );
}

MtdGNNValidation::~MtdGNNValidation() {}

// ------------ method called for each event  ------------
void MtdGNNValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));

  const auto& tp2SimAssociationMap = iEvent.get(tp2SimAssociationMapToken_);
  const auto& Sim2tpAssociationMap = iEvent.get(Sim2tpAssociationMapToken_);
  const auto& r2sAssociationMap = iEvent.get(r2sAssociationMapToken_);

  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  const auto& SigmaTofPi = iEvent.get(SigmaTofPiToken_);
  const auto& SigmaTofK = iEvent.get(SigmaTofKToken_);
  const auto& SigmaTofP = iEvent.get(SigmaTofPToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);
  const auto& btlMatchTimeChi2 = iEvent.get(btlMatchTimeChi2Token_);
  const auto& etlMatchTimeChi2 = iEvent.get(etlMatchTimeChi2Token_);
  const auto& btlMatchChi2 = iEvent.get(btlMatchChi2Token_);
  const auto& outermostHitPosition = iEvent.get(outermostHitPositionToken_);

    const auto& betaVM  = iEvent.get(betaTok_);
    const auto& phiVM   = iEvent.get(phiTok_);
    const auto& logPiVM = iEvent.get(logitPiTok_);
    const auto& logKVM  = iEvent.get(logitKTok_);
    const auto& logPVM  = iEvent.get(logitPTok_);
    const auto& emb0VM  = iEvent.get(emb0Tok_);
    const auto& emb1VM  = iEvent.get(emb1Tok_);
    const auto& emb2VM  = iEvent.get(emb2Tok_);
    const auto& pca0VM  = iEvent.get(pca0Tok_);
    const auto& pca1VM  = iEvent.get(pca1Tok_);
    const auto& pca2VM  = iEvent.get(pca2Tok_);

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  r2s_ = recoToSimH.product();

  unsigned int index = 0;

  std::map<reco::TrackRef,TrackingVertexRef> trkToTV;

  // --- Loop over all RECO tracks ---
  for (const auto& trackGen : *GenRecTrackHandle) {
    const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
    index++;

    if (trackAssoc[trackref] == -1) {
      LogWarning("mtdTracks") << "Extended track not associated";
      continue;
    }

    const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecTrackToken_), trackAssoc[trackref]);
    const reco::Track& track = *mtdTrackref;
      // == TrackingParticle based matching
      const reco::TrackBaseRef trkrefb(trackref);
      auto tp_info = getMatchedTP(trkrefb);
      if (tp_info != nullptr) {
        if (tp_info->parentVertex() != nullptr) {
	   trkToTV[trackref] = tp_info->parentVertex();
	}
      }  // TP matching

  }  // RECO tracks loop
}

// ------------ method for histogram booking ------------
void MtdGNNValidation::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking
  //
  meVtxVsZ_= ibook.book1D("VtxVsZ", "True vtx rec center vs z", 300, -15., 15.);
  meVtxSpreadVsZ_= ibook.Profile("VtxSpreadVsZ", "True vtx rec spread vs z", 300, -15., 15.,100,0., 10.);
  meVtxVsPC0_= ibook.book1D("VtxVsPC0", "True vtx rec center vs PC0", 300, -15., 15.);
  meVtxSpreadVsPC0_= ibook.Profile("VtxSpreadVsPC0", "True vtx rec spread vs PC0", 300, -15., 15.,100,0., 10.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdGNNValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks");
  desc.add<bool>("optionalPlots", false);
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("tp2SimAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  desc.add<edm::InputTag>("Sim2tpAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  desc.add<edm::InputTag>("r2sAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("btlRecHits", edm::InputTag("mtdRecHits", "FTLBarrel"));
  desc.add<edm::InputTag>("etlRecHits", edm::InputTag("mtdRecHits", "FTLEndcap"));
  desc.add<edm::InputTag>("recCluTagBTL", edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("recCluTagETL", edm::InputTag("mtdClusters", "FTLEndcap"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("btlMatchTimeChi2", edm::InputTag("trackExtenderWithMTD:btlMatchTimeChi2"));
  desc.add<edm::InputTag>("etlMatchTimeChi2", edm::InputTag("trackExtenderWithMTD:etlMatchTimeChi2"));
  desc.add<edm::InputTag>("btlMatchChi2", edm::InputTag("trackExtenderWithMTD:btlMatchChi2"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("sigmaTofPi", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("sigmaTofK", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("sigmaTofP", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("outermostHitPositionSrc",
                          edm::InputTag("trackExtenderWithMTD:generalTrackOutermostHitPosition"));
  desc.add<double>("trackMaximumPt", 12.);  // [GeV]
  desc.add<double>("trackMaximumBtlEta", 1.5);
  desc.add<double>("trackMinimumEtlEta", 1.6);
  desc.add<double>("trackMaximumEtlEta", 3.);

  descriptions.add("mtdTracksValid", desc);
}

const edm::Ref<std::vector<TrackingParticle>>* MtdGNNValidation::getMatchedTP(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return nullptr;

  //matched TP equal to any TP associated to in time events
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0)
      return &tp.first;
  }

  // reco track not matched to any TP from vertex
  return nullptr;
}

DEFINE_FWK_MODULE(MtdGNNValidation);
