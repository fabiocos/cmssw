#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <CLHEP/Units/GlobalPhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

using namespace std;
using namespace edm;

class TOFPIDProducerMCTruth : public edm::stream::EDProducer<> {
public:
  TOFPIDProducerMCTruth(const ParameterSet& pset);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  template <class H, class T>
  void fillValueMap(edm::Event& iEvent,
                    const edm::Handle<H>& handle,
                    const std::vector<T>& vec,
                    const std::string& name) const;

  void produce(edm::Event& ev, const edm::EventSetup& es) final;
  const edm::Ref<std::vector<TrackingParticle>>* getAnyMatchedTP(const reco::TrackBaseRef&);

private:
  static constexpr char t0Name[] = "t0";
  static constexpr char sigmat0Name[] = "sigmat0";
  static constexpr char t0safeName[] = "t0safe";
  static constexpr char sigmat0safeName[] = "sigmat0safe";
  static constexpr char probPiName[] = "probPi";
  static constexpr char probKName[] = "probK";
  static constexpr char probPName[] = "probP";

  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofkToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofpToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofkToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  const double fixedT0Error_;
  const double probPion_;
  const double probKaon_;
  const double probProton_;

  const reco::RecoToSimCollection* r2s_;
};

TOFPIDProducerMCTruth::TOFPIDProducerMCTruth(const ParameterSet& iConfig)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"))),
      t0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"))),
      tmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtdSrc"))),
      sigmat0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"))),
      sigmatmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtdSrc"))),
      tofkToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofkSrc"))),
      tofpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofpSrc"))),
      sigmatofpiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpiSrc"))),
      sigmatofkToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofkSrc"))),
      sigmatofpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpSrc"))),
      trackingParticleCollectionToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"))),
      recoToSimAssociationToken_(consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"))),
      fixedT0Error_(iConfig.getParameter<double>("fixedT0Error")),	
      probPion_(iConfig.getParameter<double>("probPion")),
      probKaon_(iConfig.getParameter<double>("probKaon")),
      probProton_(iConfig.getParameter<double>("probProton")){
  produces<edm::ValueMap<float>>(t0Name);
  produces<edm::ValueMap<float>>(sigmat0Name);
  produces<edm::ValueMap<float>>(t0safeName);
  produces<edm::ValueMap<float>>(sigmat0safeName);
  produces<edm::ValueMap<float>>(probPiName);
  produces<edm::ValueMap<float>>(probKName);
  produces<edm::ValueMap<float>>(probPName);
}

// Configuration descriptions
void TOFPIDProducerMCTruth::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"))->setComment("Input tracks collection");
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"))
      ->setComment("Input ValueMap for track time at beamline");
  desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"))
      ->setComment("Input ValueMap for track time at MTD");
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"))
      ->setComment("Input ValueMap for track time uncertainty at beamline");
  desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"))
      ->setComment("Input ValueMap for track time uncertainty at MTD");
  desc.add<edm::InputTag>("tofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"))
      ->setComment("Input ValueMap for track tof as kaon");
  desc.add<edm::InputTag>("tofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"))
      ->setComment("Input ValueMap for track tof as proton");
  desc.add<edm::InputTag>("sigmatofpiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"))
      ->setComment("Input ValueMap for track sigma(tof) as pion");
  desc.add<edm::InputTag>("sigmatofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"))
      ->setComment("Input ValueMap for track sigma(tof) as kaon");
  desc.add<edm::InputTag>("sigmatofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"))
      ->setComment("Input ValueMap for track sigma(tof) as proton");
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<double>("fixedT0Error", 0.)->setComment("Use a fixed T0 uncertainty [ns]");
  desc.add<double>("probPion", 1.)->setComment("A priori probability pions");
  desc.add<double>("probKaon", 1.)->setComment("A priori probability kaons");
  desc.add<double>("probProton", 1.)->setComment("A priori probability for protons");

  descriptions.add("tofPIDProducerMCTruth", desc);
}

template <class H, class T>
void TOFPIDProducerMCTruth::fillValueMap(edm::Event& iEvent,
                                  const edm::Handle<H>& handle,
                                  const std::vector<T>& vec,
                                  const std::string& name) const {
  auto out = std::make_unique<edm::ValueMap<T>>();
  typename edm::ValueMap<T>::Filler filler(*out);
  filler.insert(handle, vec.begin(), vec.end());
  filler.fill();
  iEvent.put(std::move(out), name);
}

const edm::Ref<std::vector<TrackingParticle>>* TOFPIDProducerMCTruth::getAnyMatchedTP(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return nullptr;

  //matched TP equal to any TP
  for (const auto& tp : found->val) {
    return &tp.first;
  }

  // reco track not matched to any TP from vertex
  return nullptr;
}

void TOFPIDProducerMCTruth::produce(edm::Event& ev, const edm::EventSetup& es) {
  edm::Handle<reco::TrackCollection> tracksH;
  ev.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  const auto& t0In = ev.get(t0Token_);

  const auto& tmtdIn = ev.get(tmtdToken_);

  const auto& sigmat0In = ev.get(sigmat0Token_);

  const auto& sigmatmtdIn = ev.get(sigmatmtdToken_);

  const auto& tofkIn = ev.get(tofkToken_);

  const auto& tofpIn = ev.get(tofpToken_);

  const auto& sigmatofpiIn = ev.get(sigmatofpiToken_);

  const auto& sigmatofkIn = ev.get(sigmatofkToken_);

  const auto& sigmatofpIn = ev.get(sigmatofpToken_);

  edm::Handle<TrackingParticleCollection> TPCollectionH;
  ev.getByToken(trackingParticleCollectionToken_, TPCollectionH);

  edm::Handle<reco::RecoToSimCollection> recoToSimH;
  ev.getByToken(recoToSimAssociationToken_, recoToSimH);
  r2s_ = recoToSimH.product();

  //output value maps (PID probabilities and recalculated time at beamline)
  std::vector<float> t0OutRaw;
  std::vector<float> sigmat0OutRaw;
  std::vector<float> t0safeOutRaw;
  std::vector<float> sigmat0safeOutRaw;
  std::vector<float> probPiOutRaw;
  std::vector<float> probKOutRaw;
  std::vector<float> probPOutRaw;

  //Do work here
  for (unsigned int itrack = 0; itrack < tracks.size(); ++itrack) {
    const reco::TrackRef trackref(tracksH, itrack);
    const reco::TrackBaseRef trackBaseRef(trackref);
    float t0 = t0In[trackref];
    float t0safe = t0;
    float sigmat0safe = sigmat0In[trackref];
    float sigmatmtd = (sigmatmtdIn[trackref] > 0. && fixedT0Error_ > 0.) ? fixedT0Error_ : sigmatmtdIn[trackref];
    float sigmat0 = sigmatmtd;
    float sigmatofpi = sigmatofpiIn[trackref];
    float sigmatofk = sigmatofkIn[trackref];
    float sigmatofp = sigmatofpIn[trackref];

    float prob_pi = -1.;
    float prob_k = -1.;
    float prob_p = -1.;
    //If track has time measurement
    if (sigmat0 > 0.){
       //recompute t0 for alternate mass hypotheses
       double tmtd = tmtdIn[trackref];
       double t0_k = tmtd - tofkIn[trackref];
       double t0_p = tmtd - tofpIn[trackref];

       //match and define t0, t0safe, sigmat0, sigmat0safe, corresponding prob 
       auto anytp_info = getAnyMatchedTP(trackBaseRef);
       if (anytp_info != nullptr) {
         int tp_pdgId = std::abs((*anytp_info)->pdgId());
         if(tp_pdgId == 211 || tp_pdgId == 11 || tp_pdgId == 13){
           t0safe = t0;
           sigmat0 = std::sqrt(sigmatmtd * sigmatmtd + sigmatofpi * sigmatofpi);
           sigmat0safe = sigmat0;
           prob_pi = 1.;
           prob_k = 0.;
           prob_p = 0.;	  
         }else if (tp_pdgId == 321){
           t0 = t0_k;
           t0safe = t0;
           sigmat0 = std::sqrt(sigmatmtd * sigmatmtd + sigmatofk * sigmatofk);
           sigmat0safe = sigmat0;
           prob_pi = 0.;
           prob_k = 1.;
           prob_p = 0.;
         }else if (tp_pdgId == 2212 || tp_pdgId == 3112 || tp_pdgId == 3222 || tp_pdgId == 3312){
           t0 = t0_p;
           t0safe = t0;
           sigmat0 = std::sqrt(sigmatmtd * sigmatmtd + sigmatofp * sigmatofp);
           sigmat0safe = sigmat0;
           prob_pi = 0.;
           prob_k = 0.;
           prob_p = 1.;	
         }else{
           t0 = 0.;
           t0safe = 0.;
           sigmat0 = 0.2; 
           sigmat0safe = 0.2;
         }
       }else{
         t0 = 0.;
         t0safe = 0.;
         sigmat0 = 0.2;
         sigmat0safe = 0.2;
       }
    }else{
      t0 = 0.;
      t0safe = 0.;
      sigmat0 = 0.2;
      sigmat0safe = 0.2;
    }
    t0OutRaw.push_back(t0);
    sigmat0OutRaw.push_back(sigmat0);
    t0safeOutRaw.push_back(t0safe);
    sigmat0safeOutRaw.push_back(sigmat0safe);
    probPiOutRaw.push_back(prob_pi);
    probKOutRaw.push_back(prob_k);
    probPOutRaw.push_back(prob_p);
  }

  fillValueMap(ev, tracksH, t0OutRaw, t0Name);
  fillValueMap(ev, tracksH, sigmat0OutRaw, sigmat0Name);
  fillValueMap(ev, tracksH, t0safeOutRaw, t0safeName);
  fillValueMap(ev, tracksH, sigmat0safeOutRaw, sigmat0safeName);
  fillValueMap(ev, tracksH, probPiOutRaw, probPiName);
  fillValueMap(ev, tracksH, probKOutRaw, probKName);
  fillValueMap(ev, tracksH, probPOutRaw, probPName);
}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TOFPIDProducerMCTruth);
