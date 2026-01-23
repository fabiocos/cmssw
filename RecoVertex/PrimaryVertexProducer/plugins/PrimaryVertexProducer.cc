#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"

#define cputime
#ifdef cputime
#include <chrono>
typedef std::chrono::duration<int, std::micro> microseconds_type;
#endif

#define checkrss
#ifdef checkrss
int parseLine(char* line) {
  // This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p < '0' || *p > '9')
    p++;
  line[i - 3] = '\0';
  i = atoi(p);
  return i;
}

int getValue() {  //Note: this value is in KB!
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL) {
    if (strncmp(line, "VmRSS:", 6) == 0) {
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}
#endif

PrimaryVertexProducer::PrimaryVertexProducer(const edm::ParameterSet& conf, const ONNXRuntime* cache)
    : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf) {
  fVerbose = conf.getUntrackedParameter<bool>("verbose", false);
  useMVASelection_ = conf.getParameter<bool>("useMVACut");
  trkToken = consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("TrackLabel"));
  bsToken = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  useTransientTrackTime_ = false;

  // select and configure the track selection
  std::string trackSelectionAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
  if (trackSelectionAlgorithm == "filter") {
    theTrackFilter = new TrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else if (trackSelectionAlgorithm == "filterWithThreshold") {
    theTrackFilter = new HITrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else {
    throw VertexException("PrimaryVertexProducer: unknown track selection algorithm: " + trackSelectionAlgorithm);
  }

  // select and configure the track clusterizer
  std::string clusteringAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
  if (clusteringAlgorithm == "gap") {
    theTrackClusterizer = new GapClusterizerInZ(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  } else if (clusteringAlgorithm == "DA_vect") {
    theTrackClusterizer = new DAClusterizerInZ_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  } else if (clusteringAlgorithm == "DA2D_vect") {
    theTrackClusterizer = new DAClusterizerInZT_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
    useTransientTrackTime_ = true;
  } else if (clusteringAlgorithm == "GNN2D_vect") {
    const ONNXRuntime* onnxRuntime = cache;
    theTrackClusterizer = new GNNClusterizer(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"),
        onnxRuntime);
    useTransientTrackTime_ = true;
  } else if (clusteringAlgorithm == "GNN2D_alpaka") {
    // Alpaka backend: consume pre-computed SoA from upstream Alpaka producers
    // Inference happens in GNNVertexProducerAlpaka, we only build vertices here
    useAlpakaGNN_ = true;
    useTransientTrackTime_ = true;
    theTrackClusterizer = nullptr;  // Not used for Alpaka path
    
    const auto& clusParams = conf.getParameter<edm::ParameterSet>("TkClusParameters")
                                .getParameter<edm::ParameterSet>("TkDAClusParameters");
    slotPredictionsToken_ = consumes<vertexgnn::SlotPredictionsHostCollection>(
        clusParams.getParameter<edm::InputTag>("slotPredictions"));
    assignmentsToken_ = consumes<vertexgnn::AssignmentHostCollection>(
        clusParams.getParameter<edm::InputTag>("assignments"));
    alpakaClusterizer_ = std::make_unique<vertexgnn::GNNClusterizerFromAlpaka>(clusParams);
  } else {
    throw VertexException("PrimaryVertexProducer: unknown clustering algorithm: " + clusteringAlgorithm);
  }

  if (useTransientTrackTime_) {
    trkTimesToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("TrackTimeResosLabel"));
    trackMTDTimeQualityToken =
        consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"));
    minTrackTimeQuality_ = conf.getParameter<double>("minTrackTimeQuality");
    trkMTDAssocToken = consumes<edm::ValueMap<int>>(conf.getParameter<edm::InputTag>("trackAssocSrc"));
    MTDtimeToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("tmtdSrc"));
    sigmaMTDtimeToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("sigmatmtdSrc"));
    pathLengthToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("pathmtd"));
    btlMatchChi2Token = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("btlMatchChi2Src"));
    btlMatchTime_Chi2Token = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("btlMatchTimeChi2Src"));
    etlMatchChi2Token = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("etlMatchChi2Src"));
    etlMatchTime_Chi2Token = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("etlMatchTimeChi2Src"));
    trkTimePiToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("tofPi"));
    trkTimeKToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("tofK"));
    trkTimePToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("tofP"));
    sigmaTrkTimePiToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("sigmatofpiSrc"));
    sigmaTrkTimeKToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("sigmatofkSrc"));
    sigmaTrkTimePToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("sigmatofpSrc"));
    npixBarrelToken = consumes<edm::ValueMap<int>>(conf.getParameter<edm::InputTag>("npixBarrelSrc"));
    npixEndcapToken = consumes<edm::ValueMap<int>>(conf.getParameter<edm::InputTag>("npixEndcapSrc"));
  }
  // VertexSlotModel outputs
  produces<edm::ValueMap<float>>("gnnSlotAssignment");
  produces<edm::ValueMap<float>>("gnnMaxProb");
  produces<edm::ValueMap<float>>("gnnPiWeight0");
  produces<edm::ValueMap<float>>("gnnPiWeight1");
  produces<edm::ValueMap<float>>("gnnPiWeight2");

  // select and configure the vertex fitters
  std::vector<edm::ParameterSet> vertexCollections =
      conf.getParameter<std::vector<edm::ParameterSet>>("vertexCollections");

  for (std::vector<edm::ParameterSet>::const_iterator algoconf = vertexCollections.begin();
       algoconf != vertexCollections.end();
       algoconf++) {
    algo algorithm;

    algorithm.label = algoconf->getParameter<std::string>("label");

    // configure the fitter and selector
    std::string fitterAlgorithm = algoconf->getParameter<std::string>("algorithm");
    if (fitterAlgorithm == "KalmanVertexFitter") {
      algorithm.pv_fitter = new SequentialPrimaryVertexFitterAdapter(new KalmanVertexFitter());
    } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
      auto fitter = new AdaptiveVertexFitter(GeometricAnnealing(algoconf->getParameter<double>("chi2cutoff")));
      algorithm.pv_fitter = new SequentialPrimaryVertexFitterAdapter(fitter);
    } else if (fitterAlgorithm.empty()) {
      algorithm.pv_fitter = nullptr;
    } else if (fitterAlgorithm == "AdaptiveChisquareVertexFitter") {
      algorithm.pv_fitter = new AdaptiveChisquarePrimaryVertexFitter(algoconf->getParameter<double>("chi2cutoff"),
                                                                     algoconf->getParameter<double>("zcutoff"),
                                                                     algoconf->getParameter<double>("mintrkweight"),
                                                                     false);
    } else if (fitterAlgorithm == "MultiPrimaryVertexFitter") {
      algorithm.pv_fitter = new AdaptiveChisquarePrimaryVertexFitter(algoconf->getParameter<double>("chi2cutoff"),
                                                                     algoconf->getParameter<double>("zcutoff"),
                                                                     algoconf->getParameter<double>("mintrkweight"),
                                                                     true);
    } else if (fitterAlgorithm == "WeightedMeanFitter") {
      algorithm.pv_fitter = new WeightedMeanPrimaryVertexEstimator();
    } else {
      throw VertexException("PrimaryVertexProducer: unknown algorithm: " + fitterAlgorithm);
    }
    algorithm.minNdof = algoconf->getParameter<double>("minNdof");
    algorithm.useBeamConstraint = algoconf->getParameter<bool>("useBeamConstraint");
    algorithm.vertexSelector =
        new VertexCompatibleWithBeam(VertexDistanceXY(), algoconf->getParameter<double>("maxDistanceToBeam"));

    // configure separate vertex time reconstruction if applicable
    // note that the vertex time could, in principle, also come from the clusterizer or the vertex fit

    const auto& pv_time_conf = algoconf->getParameter<edm::ParameterSet>("vertexTimeParameters");
    const std::string vertexTimeAlgorithm = pv_time_conf.getParameter<std::string>("algorithm");
    edm::ConsumesCollector&& collector = consumesCollector();

    if (vertexTimeAlgorithm.empty()) {
      algorithm.pv_time_estimator = nullptr;
    } else if (vertexTimeAlgorithm == "legacy4D") {
      useTransientTrackTime_ = true;
      algorithm.pv_time_estimator =
          new VertexTimeAlgorithmLegacy4D(pv_time_conf.getParameter<edm::ParameterSet>("legacy4D"), collector);
    } else if (vertexTimeAlgorithm == "fromTracksPID") {
      algorithm.pv_time_estimator = new VertexTimeAlgorithmFromTracksPID(
          pv_time_conf.getParameter<edm::ParameterSet>("fromTracksPID"), collector);
    } else {
      edm::LogWarning("MisConfiguration") << "unknown vertexTimeParameters.algorithm" << vertexTimeAlgorithm;
    }
    algorithms.push_back(algorithm);

    produces<reco::VertexCollection>(algorithm.label);
  }

  //check if this is a recovery iteration
  fRecoveryIteration = conf.getParameter<bool>("isRecoveryIteration");
  if (fRecoveryIteration) {
    if (algorithms.empty()) {
      throw VertexException("PrimaryVertexProducer: No algorithm specified. ");
    } else if (algorithms.size() > 1) {
      throw VertexException(
          "PrimaryVertexProducer: Running in Recovery mode and more than one algorithm specified.  Please "
          "only one algorithm.");
    }
    recoveryVtxToken = consumes<reco::VertexCollection>(conf.getParameter<edm::InputTag>("recoveryVtxCollection"));
  }
}

PrimaryVertexProducer::~PrimaryVertexProducer() {
  if (theTrackFilter)
    delete theTrackFilter;
  if (theTrackClusterizer)
    delete theTrackClusterizer;
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    if (algorithm->pv_fitter)
      delete algorithm->pv_fitter;
    if (algorithm->pv_time_estimator)
      delete algorithm->pv_time_estimator;
    if (algorithm->vertexSelector)
      delete algorithm->vertexSelector;
  }
}
std::unique_ptr<ONNXRuntime> PrimaryVertexProducer::initializeGlobalCache(const edm::ParameterSet& conf) {
  // Only initialize ONNX runtime for GNN clustering algorithm
  const auto& clusParams = conf.getParameter<edm::ParameterSet>("TkClusParameters");
  std::string algorithm = clusParams.getParameter<std::string>("algorithm");
  
  if (algorithm == "GNN2D_vect") {
    // Delegate to GNNClusterizer for GPU backend support
    return GNNClusterizer::initializeGlobalCache(
        clusParams.getParameter<edm::ParameterSet>("TkDAClusParameters"));
  }
  
  // For non-GNN algorithms, use default
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
}
void PrimaryVertexProducer::globalEndJob(const ONNXRuntime* cache) {}

void PrimaryVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get the BeamSpot, it will always be needed, even when not used as a constraint
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }

  bool validBS = true;
  VertexState beamVertexState(beamSpot);
  if ((beamVertexState.error().cxx() <= 0.) || (beamVertexState.error().cyy() <= 0.) ||
      (beamVertexState.error().czz() <= 0.)) {
    edm::LogError("UnusableBeamSpot") << "Beamspot with invalid errors " << beamVertexState.error().matrix();
    validBS = false;
  }

  //if this is a recovery iteration, check if we already have a valid PV
  if (fRecoveryIteration) {
    auto const& oldVertices = iEvent.get(recoveryVtxToken);
    //look for the first valid (not-BeamSpot) vertex
    for (auto const& old : oldVertices) {
      if (!(old.isFake())) {
        //found a valid vertex, write the first one to the collection and return
        //otherwise continue with regular vertexing procedure
        auto result = std::make_unique<reco::VertexCollection>();
        result->push_back(old);
        iEvent.put(std::move(result), algorithms.begin()->label);
        return;
      }
    }
  }

  // get RECO tracks from the event
  // `tks` can be used as a ptr to a reco::TrackCollection
  edm::Handle<reco::TrackCollection> tks;
  iEvent.getByToken(trkToken, tks);

  // mechanism to put the beamspot if the track collection is empty
  if (!tks.isValid()) {
    for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
      auto result = std::make_unique<reco::VertexCollection>();
      reco::VertexCollection& vColl = (*result);

      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl.push_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << "Beamspot with invalid errors " << bse.matrix() << std::endl;
          std::cout << "Will put Vertex derived from dummy-fake BeamSpot into Event.\n";
        }
      } else {
        vColl.push_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << " will put Vertex derived from BeamSpot into Event.\n";
        }
      }
      iEvent.put(std::move(result), algorithm->label);
    }

    return;  // early return
  }

  for (auto& algo : algorithms) {
    if (algo.pv_time_estimator) {
      algo.pv_time_estimator->setEvent(iEvent, iSetup);
    }
  }

  // interface RECO tracks to vertex reconstruction
  const auto& theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;

  if (useTransientTrackTime_) {
    auto const& trackTimeResos_ = iEvent.get(trkTimeResosToken);
    auto trackTimes_ = iEvent.get(trkTimesToken);
    auto const& trackMTDAssoc_ = iEvent.get(trkMTDAssocToken);
    auto trackMTDTimes_ = iEvent.get(MTDtimeToken);
    auto trackMTDTimesRes_ = iEvent.get(sigmaMTDtimeToken);
    auto pathlength_ = iEvent.get(pathLengthToken);
    auto btlMatch_ = iEvent.get(btlMatchChi2Token);
    auto btlMatchTime_ = iEvent.get(btlMatchTime_Chi2Token);
    auto etlMatch_ = iEvent.get(etlMatchChi2Token);
    auto etlMatchTime_ = iEvent.get(etlMatchTime_Chi2Token);
    auto trkPiTime_ = iEvent.get(trkTimePiToken);
    auto trkKTime_ = iEvent.get(trkTimeKToken);
    auto trkPTime_ = iEvent.get(trkTimePToken);
    auto sigmatrkPiTime_ = iEvent.get(sigmaTrkTimePiToken);
    auto sigmatrkKTime_ = iEvent.get(sigmaTrkTimeKToken);
    auto sigmatrkPTime_ = iEvent.get(sigmaTrkTimePToken);
    auto npixbarrel_ = iEvent.get(npixBarrelToken);
    auto npixendcap_ = iEvent.get(npixEndcapToken);
    auto trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken);

    if (useMVASelection_) {
      //   trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken);

      for (unsigned int i = 0; i < (*tks).size(); i++) {
        const reco::TrackRef ref(tks, i);
        auto const trkTimeQuality = trackMTDTimeQualities_[ref];
        if (trkTimeQuality < minTrackTimeQuality_) {
          trackTimes_[ref] = std::numeric_limits<double>::max();
        }
      }
      //  t_tks = (*theB).build(tks, beamSpot, trackTimes_, trackTimeResos_);
      t_tks = (*theB).build(tks,
                            beamSpot,
                            trackTimes_,
                            trackTimeResos_,
                            trackMTDAssoc_,
                            trackMTDTimes_,
                            trackMTDTimesRes_,
                            trackMTDTimeQualities_,
                            pathlength_,
                            btlMatch_,
                            btlMatchTime_,
                            etlMatch_,
                            etlMatchTime_,
                            trkPiTime_,
                            trkKTime_,
                            trkPTime_,
                            sigmatrkPiTime_,
                            sigmatrkKTime_,
                            sigmatrkPTime_,
                            npixbarrel_,
                            npixendcap_);
    } else {
      t_tks = (*theB).build(tks,
                            beamSpot,
                            trackTimes_,
                            trackTimeResos_,
                            trackMTDAssoc_,
                            trackMTDTimes_,
                            trackMTDTimesRes_,
                            trackMTDTimeQualities_,
                            pathlength_,
                            btlMatch_,
                            btlMatchTime_,
                            etlMatch_,
                            etlMatchTime_,
                            trkPiTime_,
                            trkKTime_,
                            trkPTime_,
                            sigmatrkPiTime_,
                            sigmatrkKTime_,
                            sigmatrkPTime_,
                            npixbarrel_,
                            npixendcap_);
    }
  } else {
    t_tks = (*theB).build(tks, beamSpot);
  }

  // select tracks
  std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);

#ifdef checkrss
  edm::LogInfo("PrimaryVertexProducer") << "###RSS before clusterizer (kB) " << getValue();
#endif

#ifdef cputime
  std::chrono::duration<int, std::micro> tcpu_clustering(0), tcpu_vtxfit(0), tcpu_vtxtime(0);
  auto start_clustering = std::chrono::high_resolution_clock::now();
#endif

  // clusterize tracks in Z
  std::vector<TransientVertex> clusters;
  
  // Variables for ValueMap production (used by both paths)
  std::vector<float> gnn_A_flat;  // N*K assignment probs
  std::vector<float> gnn_pi_0, gnn_pi_1, gnn_pi_2, gnn_pi_3;  // K PID weights
  int gnn_N = 0, gnn_K = 0;
  bool hasGNNOutputs = false;
  
  if (useAlpakaGNN_) {
    // Alpaka path: consume pre-computed SoA from upstream GNNVertexProducerAlpaka
    const auto& slotPredictions = iEvent.get(slotPredictionsToken_);
    const auto& assignments = iEvent.get(assignmentsToken_);
    
    gnn_N = seltks.size();
    gnn_K = slotPredictions.const_view().metadata().size();
    
    auto slotView = slotPredictions.const_view();
    auto assignView = assignments.const_view();
    
    // Extract data from SoA to arrays for clusterizer
    gnn_A_flat.resize(gnn_N * gnn_K);
    std::vector<float> z_hat(gnn_K), t_hat(gnn_K), p(gnn_K);
    gnn_pi_0.resize(gnn_K); gnn_pi_1.resize(gnn_K); gnn_pi_2.resize(gnn_K); gnn_pi_3.resize(gnn_K);
    
    for (int k = 0; k < gnn_K; ++k) {
      z_hat[k] = slotView[k].z_hat();
      t_hat[k] = slotView[k].t_hat();
      p[k] = slotView[k].p();
      gnn_pi_0[k] = slotView[k].pi_0();
      gnn_pi_1[k] = slotView[k].pi_1();
      gnn_pi_2[k] = slotView[k].pi_2();
      gnn_pi_3[k] = slotView[k].pi_3();
    }
    
    for (int i = 0; i < gnn_N; ++i) {
      for (int k = 0; k < gnn_K; ++k) {
        gnn_A_flat[i * gnn_K + k] = assignView[i * gnn_K + k].prob();
      }
    }
    
    // Build vertices using Alpaka clusterizer
    clusters = alpakaClusterizer_->vertices(
        seltks, gnn_A_flat.data(), z_hat.data(), t_hat.data(), p.data(),
        gnn_pi_0.data(), gnn_pi_1.data(), gnn_pi_2.data(), gnn_pi_3.data(), gnn_N, gnn_K);
    
    hasGNNOutputs = true;
        
    if (fVerbose) {
      edm::LogInfo("PrimaryVertexProducer") << "Alpaka GNN: built " << clusters.size() << " vertices from SoA";
    }
  } else {
    // ONNX path: use theTrackClusterizer (GNNClusterizer or DA)
    clusters = theTrackClusterizer->vertices(seltks);
  }
  
  // per-track GNN outputs as ValueMaps for the full input TrackCollection
  // Alpaka path: produce ValueMaps from the SoA data we extracted
  if (hasGNNOutputs && useAlpakaGNN_) {
    const auto& trkHandle = iEvent.getHandle(trkToken);
    const size_t Nall = trkHandle->size();
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    
    std::vector<float> vmSlotAssign(Nall, NaN);
    std::vector<float> vmMaxProb(Nall, NaN);
    std::vector<float> vmPi0(Nall, NaN), vmPi1(Nall, NaN), vmPi2(Nall, NaN);
    
    for (size_t i = 0; i < seltks.size() && static_cast<int>(i) < gnn_N; ++i) {
      const auto& tt = seltks[i];
      reco::TrackRef tref = tt.trackBaseRef().castTo<reco::TrackRef>();
      if (tref.isNull()) continue;
      const size_t idx = tref.key();
      if (idx >= Nall) continue;
      
      // Find best slot assignment and max probability
      if (gnn_K > 0) {
        float maxProb = -1.0f;
        int bestSlot = -1;
        for (int k = 0; k < gnn_K; ++k) {
          float prob = gnn_A_flat[i * gnn_K + k];
          if (prob > maxProb) {
            maxProb = prob;
            bestSlot = k;
          }
        }
        vmSlotAssign[idx] = static_cast<float>(bestSlot);
        vmMaxProb[idx] = maxProb;
        
        // PID weights - average over slots (simplified)
        vmPi0[idx] = gnn_pi_0[bestSlot];
        vmPi1[idx] = gnn_pi_1[bestSlot];
        vmPi2[idx] = gnn_pi_2[bestSlot];
      }
    }
    
    auto putVM = [&](const std::vector<float>& vals, const std::string& label) {
      auto out = std::make_unique<edm::ValueMap<float>>();
      edm::ValueMap<float>::Filler filler(*out);
      filler.insert(trkHandle, vals.begin(), vals.end());
      filler.fill();
      iEvent.put(std::move(out), label);
    };
    
    putVM(vmSlotAssign, "gnnSlotAssignment");
    putVM(vmMaxProb, "gnnMaxProb");
    putVM(vmPi0, "gnnPiWeight0");
    putVM(vmPi1, "gnnPiWeight1");
    putVM(vmPi2, "gnnPiWeight2");
  }
  
  // ONNX path ValueMaps
  if (auto gnn = dynamic_cast<GNNClusterizer*>(theTrackClusterizer)) {
    if (gnn->hasLastOutputs()) {
      const auto& trkHandle = iEvent.getHandle(trkToken);
      const size_t Nall = trkHandle->size();

      // VertexSlotModel outputs
      const auto& A = gnn->lastAssignments();  // N*K
      const auto& pi = gnn->lastPi();          // N*3
      const int Nsel = gnn->lastTrackCount();
      const int K = gnn->lastNumSlots();

      const float NaN = std::numeric_limits<float>::quiet_NaN();

      // Initialize with NaN, then backfill for selected tracks
      std::vector<float> vmSlotAssign(Nall, NaN);
      std::vector<float> vmMaxProb(Nall, NaN);
      std::vector<float> vmPi0(Nall, NaN), vmPi1(Nall, NaN), vmPi2(Nall, NaN);

      for (size_t i = 0; i < seltks.size(); ++i) {
        if (static_cast<int>(i) >= Nsel)
          break;
        const auto& tt = seltks[i];
        reco::TrackRef tref = tt.trackBaseRef().castTo<reco::TrackRef>();
        if (tref.isNull()) {
          edm::LogError("PrimaryVertexProducer") << "seltks # " << i << " no trackRef";
          continue;
        }
        const size_t idx = tref.key();
        if (idx >= Nall) {
          edm::LogError("PrimaryVertexProducer") << "seltks trk # " << idx << " outside track collection range";
          continue;
        }

        // Find best slot assignment and max probability for this track
        if (K > 0 && i * K + K - 1 < A.size()) {
          float maxProb = -1.0f;
          int bestSlot = -1;
          for (int k = 0; k < K; ++k) {
            float prob = A[i * K + k];
            if (prob > maxProb) {
              maxProb = prob;
              bestSlot = k;
            }
          }
          vmSlotAssign[idx] = static_cast<float>(bestSlot);
          vmMaxProb[idx] = maxProb;
        }

        // PID weights (already softmax normalized in GNNClusterizer)
        if (3 * i + 2 < pi.size()) {
          vmPi0[idx] = pi[3 * i + 0];
          vmPi1[idx] = pi[3 * i + 1];
          vmPi2[idx] = pi[3 * i + 2];
        }
      }

#ifdef cputime
      auto stop_clustering = std::chrono::high_resolution_clock::now();
      tcpu_clustering = std::chrono::duration_cast<std::chrono::microseconds>(stop_clustering - start_clustering);
      edm::LogInfo("PrimaryVertexProducer") << "###TIME clustering " << tcpu_clustering;
#endif

#ifdef checkrss
      edm::LogInfo("PrimaryVertexProducer") << "###RSS after clusterizer (kB) " << getValue();
#endif

      // Put float ValueMaps
      auto putVM = [&](const std::vector<float>& vals, const std::string& label) {
        auto out = std::make_unique<edm::ValueMap<float>>();
        edm::ValueMap<float>::Filler filler(*out);
        filler.insert(trkHandle, vals.begin(), vals.end());
        filler.fill();
        iEvent.put(std::move(out), label);
      };

      putVM(vmSlotAssign, "gnnSlotAssignment");
      putVM(vmMaxProb, "gnnMaxProb");
      putVM(vmPi0, "gnnPiWeight0");
      putVM(vmPi1, "gnnPiWeight1");
      putVM(vmPi2, "gnnPiWeight2");
    }
  }
  if (fVerbose) {
    edm::LogPrint("PrimaryVertexProducer")
        << " bool useTransientTrackTime_ " << useTransientTrackTime_ << " Clustering returned " << clusters.size()
        << " clusters from " << seltks.size() << " selected tracks";
  }

  // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);
    std::vector<TransientVertex> pvs;

#ifdef cputime
    auto start_fit = std::chrono::high_resolution_clock::now();
#endif

    if (algorithm->pv_fitter == nullptr) {
      pvs = clusters;
    } else {
      pvs = algorithm->pv_fitter->fit(seltks, clusters, beamSpot, algorithm->useBeamConstraint);
    }

#ifdef cputime
    auto stop_fit = std::chrono::high_resolution_clock::now();
    tcpu_vtxfit = std::chrono::duration_cast<std::chrono::microseconds>(stop_fit - start_fit);
    edm::LogInfo("PrimaryVertexProducer") << "###TIME fit " << tcpu_vtxfit;
    auto start_vtxtime = std::chrono::high_resolution_clock::now();
#endif

    if (algorithm->pv_time_estimator != nullptr) {
      algorithm->pv_time_estimator->fill_vertex_times(pvs);
    }

#ifdef cputime
    auto stop_vtxtime = std::chrono::high_resolution_clock::now();
    tcpu_vtxtime = std::chrono::duration_cast<std::chrono::microseconds>(stop_vtxtime - start_vtxtime);
    edm::LogInfo("PrimaryVertexProducer") << "###TIME vtxtime " << tcpu_vtxtime;
#endif

    // sort vertices by pt**2  vertex
    if (pvs.size() > 1) {
      sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }

    // select and convert transient vertices to (reco) vertices
    for (std::vector<TransientVertex>::const_iterator iv = pvs.begin(); iv != pvs.end(); iv++) {
      if (iv->isValid() && (iv->degreesOfFreedom() >= algorithm->minNdof)) {
        reco::Vertex v = *iv;
        if (!validBS || ((*(algorithm->vertexSelector))(v, beamVertexState))) {
          vColl.push_back(v);
        }
      }
    }

    if (fVerbose) {
      edm::LogPrint("PrimaryVertexProducer") << "PrimaryVertexProducer \"" << algorithm->label << "\" contains "
                                             << pvs.size() << " reco::Vertex candidates";
    }

    if (clusters.size() > 2 && clusters.size() > 2 * pvs.size()) {
      edm::LogWarning("PrimaryVertexProducer")
          << "More than 50% of candidate vertices lost (" << pvs.size() << " out of " << clusters.size() << ")";
    }

    if (pvs.empty() && seltks.size() > 5) {
      edm::LogWarning("PrimaryVertexProducer")
          << "No vertex found with " << seltks.size() << " tracks and " << clusters.size() << " vertex candidates";
    }

    if (vColl.empty()) {
      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl.push_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        edm::LogWarning("PrimaryVertexProducer") << "Zero recostructed vertices, will put reco::Vertex derived from "
                                                    "dummy/fake BeamSpot into Event, BeamSpot has invalid errors: "
                                                 << bse.matrix();
      } else {
        vColl.push_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        if (fVerbose) {
          edm::LogWarning("PrimaryVertexProducer")
              << "Zero recostructed vertices, will put reco::Vertex derived from BeamSpot into Event.";
        }
      }
    }
    if (fVerbose) {
      int ivtx = 0;
      for (reco::VertexCollection::const_iterator v = vColl.begin(); v != vColl.end(); ++v) {
        edm::LogPrint("PrimaryVertexProducer")
            << "recvtx " << std::setw(3) << std::fixed << ivtx++ << " #trk " << std::setw(3) << v->tracksSize()
            << " chi2 " << std::setw(5) << std::setprecision(1) << v->chi2() << " ndof " << std::setw(5)
            << std::setprecision(1) << v->ndof() << " x " << std::setw(7) << std::setprecision(4) << v->position().x()
            << " dx " << std::setw(6) << std::setprecision(4) << v->xError() << " y " << std::setw(7)
            << std::setprecision(4) << v->position().y() << " dy " << std::setw(6) << std::setprecision(4)
            << v->yError() << " z " << std::setw(8) << std::setprecision(4) << v->position().z() << " dz "
            << std::setw(6) << std::setprecision(4) << v->zError();
        if (v->tError() > 0) {
          edm::LogPrint("PrimaryVertexProducer") << " t " << std::setw(6) << std::setprecision(3) << v->t() << " dt "
                                                 << std::setw(6) << std::setprecision(3) << v->tError();
        }
      }
    }
    iEvent.put(std::move(result), algorithm->label);
  }
}

void PrimaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription psd_pv_time;
  {
    edm::ParameterSetDescription psd1;
    VertexTimeAlgorithmFromTracksPID::fillPSetDescription(psd1);
    psd_pv_time.add<edm::ParameterSetDescription>("fromTracksPID", psd1);

    edm::ParameterSetDescription psd2;
    VertexTimeAlgorithmLegacy4D::fillPSetDescription(psd2);
    psd_pv_time.add<edm::ParameterSetDescription>("legacy4D", psd2);
  }
  psd_pv_time.add<std::string>("algorithm", "");  // default = none

  // vertex collections
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.add<double>("maxDistanceToBeam", 1.0);
    vpsd1.add<std::string>("algorithm", "AdaptiveVertexFitter");
    vpsd1.add<bool>("useBeamConstraint", false);
    vpsd1.add<std::string>("label", "");
    vpsd1.add<double>("chi2cutoff", 2.5);
    vpsd1.add<double>("zcutoff", 1.0);
    vpsd1.add<double>("mintrkweight", 0.0);
    vpsd1.add<double>("minNdof", 0.0);
    vpsd1.add<edm::ParameterSetDescription>("vertexTimeParameters", psd_pv_time);

    // two default values : with- and without beam constraint
    std::vector<edm::ParameterSet> temp1;
    temp1.reserve(2);
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", false);
      temp2.addParameter<std::string>("label", "");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("zcutoff", 1.0);
      temp2.addParameter<double>("mintrkweight", 0.);
      temp2.addParameter<double>("minNdof", 0.0);
      edm::ParameterSet temp_vertexTime;
      temp_vertexTime.addParameter<std::string>("algorithm", "");
      temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
      temp1.push_back(temp2);
    }
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", true);
      temp2.addParameter<std::string>("label", "WithBS");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("zcutoff", 1.0);
      temp2.addParameter<double>("mintrkweight", 0.);
      temp2.addParameter<double>("minNdof", 2.0);
      edm::ParameterSet temp_vertexTime;
      temp_vertexTime.addParameter<std::string>("algorithm", "");
      temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
      temp1.push_back(temp2);
    }
    desc.addVPSet("vertexCollections", vpsd1, temp1);
  }
  desc.addUntracked<bool>("verbose", false);
  {
    edm::ParameterSetDescription psd0;
    HITrackFilterForPVFinding::fillPSetDescription(psd0);  // extension of TrackFilterForPVFinding
    desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
  }
  desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("TrackLabel", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("TrackTimeResosLabel", edm::InputTag("dummy_default"));                         // 4D only
  desc.add<edm::InputTag>("TrackTimesLabel", edm::InputTag("dummy_default"));                             // 4D only
  desc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));  // 4D only
  desc.add<edm::InputTag>("timingSoA", edm::InputTag("mtdSoA"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"));
  desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("pathmtd", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("btlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchChi2"));
  desc.add<edm::InputTag>("btlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchTimeChi2"));
  desc.add<edm::InputTag>("etlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchChi2"));
  desc.add<edm::InputTag>("etlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchTimeChi2"));
  desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("sigmatofpiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("sigmatofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("sigmatofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
  desc.add<edm::InputTag>("npixBarrelSrc", edm::InputTag("trackExtenderWithMTD", "npixBarrel"));
  desc.add<edm::InputTag>("npixEndcapSrc", edm::InputTag("trackExtenderWithMTD", "npixEndcap"));
  {
    edm::ParameterSetDescription psd0;
    {
      edm::ParameterSetDescription psd1;
      DAClusterizerInZ_vect::fillPSetDescription(psd1);

      edm::ParameterSetDescription psd2;
      DAClusterizerInZT_vect::fillPSetDescription(psd2);

      edm::ParameterSetDescription psd3;
      GapClusterizerInZ::fillPSetDescription(psd3);

      edm::ParameterSetDescription psd4;
      GNNClusterizer::fillPSetDescription(psd4);

      // Alpaka GNN backend - consumes SoA from upstream Alpaka producers
      edm::ParameterSetDescription psd5;
      vertexgnn::GNNClusterizerFromAlpaka::fillPSetDescription(psd5);
      psd5.add<edm::InputTag>("slotPredictions", edm::InputTag("gnnVertexProducer"))
          ->setComment("Input tag for SlotPredictionsHostCollection from GNNVertexProducerAlpaka");
      psd5.add<edm::InputTag>("assignments", edm::InputTag("gnnVertexProducer"))
          ->setComment("Input tag for AssignmentHostCollection from GNNVertexProducerAlpaka");

      psd0.ifValue(
          edm::ParameterDescription<std::string>("algorithm", "DA_vect", true),
          "DA_vect" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd1, true) or
              "DA2D_vect" >>
                  edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd2, true) or
              "GNN2D_vect" >>
                  edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd4, true) or
              "GNN2D_alpaka" >>
                  edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd5, true) or
              "gap" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkGapClusParameters", psd3, true));
    }
    desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
  }

  desc.add<bool>("isRecoveryIteration", false);
  desc.add<edm::InputTag>("recoveryVtxCollection", {""});
  desc.add<bool>("useMVACut", false);
  desc.add<double>("minTrackTimeQuality", 0.8);
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/vertex_slot_model.onnx"))
      ->setComment("Path to VertexSlotModel ONNX export");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryVertexProducer);
