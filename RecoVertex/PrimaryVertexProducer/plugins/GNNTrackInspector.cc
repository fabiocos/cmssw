#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

class GNNTrackInspector : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit GNNTrackInspector(const edm::ParameterSet& iConfig);
  ~GNNTrackInspector() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event& iEvent, const edm::EventSetup&) override;

  const edm::InputTag trackSrc_;
  const edm::InputTag pvModule_;
  const uint32_t printFirstN_;
  const bool dropNaNs_;

  edm::EDGetTokenT<reco::TrackCollection> tkTok_;
  edm::EDGetTokenT<edm::ValueMap<float>> betaTok_;
  edm::EDGetTokenT<edm::ValueMap<float>> phiTok_;
  edm::EDGetTokenT<edm::ValueMap<float>> logitPiTok_;
  edm::EDGetTokenT<edm::ValueMap<float>> logitKTok_;
  edm::EDGetTokenT<edm::ValueMap<float>> logitPTok_;
  edm::EDGetTokenT<edm::ValueMap<float>> emb0Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>> emb1Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>> emb2Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>> pca0Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>> pca1Tok_;
  edm::EDGetTokenT<edm::ValueMap<float>> pca2Tok_;

  TH1F *h_beta_, *h_phi_, *h_logitPi_, *h_logitK_, *h_logitP_;
  TH1F *h_emb0_, *h_emb1_, *h_emb2_, *h_pca0_, *h_pca1_, *h_pca2_;
};

GNNTrackInspector::GNNTrackInspector(const edm::ParameterSet& iConfig)
    : trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
      pvModule_(iConfig.getParameter<edm::InputTag>("pvModule")),
      printFirstN_(iConfig.getParameter<uint32_t>("printFirstN")),
      dropNaNs_(iConfig.getParameter<bool>("dropNaNs")) {
  usesResource("TFileService");

  tkTok_ = consumes<reco::TrackCollection>(trackSrc_);

  auto mkTag = [&](const std::string& instance) -> edm::InputTag {
    return edm::InputTag(pvModule_.label(), instance, pvModule_.process());
  };

  betaTok_ = consumes<edm::ValueMap<float>>(mkTag("gnnBeta"));
  phiTok_ = consumes<edm::ValueMap<float>>(mkTag("gnnPhi"));
  logitPiTok_ = consumes<edm::ValueMap<float>>(mkTag("gnnPidLogitPi"));
  logitKTok_ = consumes<edm::ValueMap<float>>(mkTag("gnnPidLogitK"));
  logitPTok_ = consumes<edm::ValueMap<float>>(mkTag("gnnPidLogitP"));
  emb0Tok_ = consumes<edm::ValueMap<float>>(mkTag("gnnEmb0"));
  emb1Tok_ = consumes<edm::ValueMap<float>>(mkTag("gnnEmb1"));
  emb2Tok_ = consumes<edm::ValueMap<float>>(mkTag("gnnEmb2"));
  pca0Tok_ = consumes<edm::ValueMap<float>>(mkTag("gnnPCA0"));
  pca1Tok_ = consumes<edm::ValueMap<float>>(mkTag("gnnPCA1"));
  pca2Tok_ = consumes<edm::ValueMap<float>>(mkTag("gnnPCA2"));
}

void GNNTrackInspector::beginJob() {
  edm::Service<TFileService> fs;
  h_beta_ = fs->make<TH1F>("beta", "GNN beta;beta;tracks", 100, 0.0, 1.0);
  h_phi_ = fs->make<TH1F>("phi", "GNN phi;phi;tracks", 100, 0.0, 1.0);
  h_logitPi_ = fs->make<TH1F>("logitPi", "PID logit #pi;logit;tracks", 100, -15, 15);
  h_logitK_ = fs->make<TH1F>("logitK", "PID logit K;logit;tracks", 100, -15, 15);
  h_logitP_ = fs->make<TH1F>("logitP", "PID logit p;logit;tracks", 100, -15, 15);

  h_emb0_ = fs->make<TH1F>("emb0", "Embedding[0];value;tracks", 100, -5, 5);
  h_emb1_ = fs->make<TH1F>("emb1", "Embedding[1];value;tracks", 100, -5, 5);
  h_emb2_ = fs->make<TH1F>("emb2", "Embedding[2];value;tracks", 100, -5, 5);

  h_pca0_ = fs->make<TH1F>("pca0", "PCA[0];value;tracks", 100, -5, 5);
  h_pca1_ = fs->make<TH1F>("pca1", "PCA[1];value;tracks", 100, -5, 5);
  h_pca2_ = fs->make<TH1F>("pca2", "PCA[2];value;tracks", 100, -5, 5);
}

void GNNTrackInspector::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<reco::TrackCollection> hTracks;
  iEvent.getByToken(tkTok_, hTracks);
  auto const& tracks = *hTracks;
  auto const& betaVM = iEvent.get(betaTok_);
  auto const& phiVM = iEvent.get(phiTok_);
  auto const& logPiVM = iEvent.get(logitPiTok_);
  auto const& logKVM = iEvent.get(logitKTok_);
  auto const& logPVM = iEvent.get(logitPTok_);
  auto const& emb0VM = iEvent.get(emb0Tok_);
  auto const& emb1VM = iEvent.get(emb1Tok_);
  auto const& emb2VM = iEvent.get(emb2Tok_);
  auto const& pca0VM = iEvent.get(pca0Tok_);
  auto const& pca1VM = iEvent.get(pca1Tok_);
  auto const& pca2VM = iEvent.get(pca2Tok_);

  const size_t nTk = tracks.size();

  auto get = [&](const edm::ValueMap<float>& vm, size_t i) -> float {
    reco::TrackRef tref(hTracks, i);
    float v = vm[tref];
    return (dropNaNs_ && !std::isfinite(v)) ? 0.f : v;
  };

  for (size_t i = 0; i < std::min<size_t>(nTk, printFirstN_); ++i) {
    float beta = get(betaVM, i), phi = get(phiVM, i);
    float lpi = get(logPiVM, i), lk = get(logKVM, i), lp = get(logPVM, i);
    float e0 = get(emb0VM, i), e1 = get(emb1VM, i), e2 = get(emb2VM, i);
    float p0 = get(pca0VM, i), p1 = get(pca1VM, i), p2 = get(pca2VM, i);

    auto const& trk = tracks[i];
    edm::LogInfo("GNNTrackInspector") << "track[" << i << "] pt=" << trk.pt() << " eta=" << trk.eta()
                                      << " beta=" << beta << " phi=" << phi << " log(pi,k,p)=(" << lpi << "," << lk
                                      << "," << lp << ")"
                                      << " emb[0..2]=(" << e0 << "," << e1 << "," << e2 << ")"
                                      << " pca[0..2]=(" << p0 << "," << p1 << "," << p2 << ")";
  }

  auto fillIf = [&](TH1F* h, const edm::ValueMap<float>& vm, size_t i) {
    reco::TrackRef tref(hTracks, i);
    float v = vm[tref];
    if (dropNaNs_ && !std::isfinite(v))
      return;
    h->Fill(v);
  };
  for (size_t i = 0; i < nTk; ++i) {
    fillIf(h_beta_, betaVM, i);
    fillIf(h_phi_, phiVM, i);
    fillIf(h_logitPi_, logPiVM, i);
    fillIf(h_logitK_, logKVM, i);
    fillIf(h_logitP_, logPVM, i);
    fillIf(h_emb0_, emb0VM, i);
    fillIf(h_emb1_, emb1VM, i);
    fillIf(h_emb2_, emb2VM, i);
    fillIf(h_pca0_, pca0VM, i);
    fillIf(h_pca1_, pca1VM, i);
    fillIf(h_pca2_, pca2VM, i);
  }
}

void GNNTrackInspector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("generalTracks"))
      ->setComment("The TrackCollection these ValueMaps are keyed against.");
  desc.add<edm::InputTag>("pvModule", edm::InputTag("unsortedOfflinePrimaryVertices4D"))
      ->setComment("Module label of the PrimaryVertexProducer that produced gnn* ValueMaps.");
  desc.add<uint32_t>("printFirstN", 10)->setComment("Print first N tracks per event to the Log.");
  desc.add<bool>("dropNaNs", true)->setComment("If true, skip NaNs when filling histos.");
  descriptions.add("GNNTrackInspector", desc);
}
DEFINE_FWK_MODULE(GNNTrackInspector);