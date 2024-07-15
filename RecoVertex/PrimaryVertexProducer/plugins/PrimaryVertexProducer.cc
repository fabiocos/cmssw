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

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"

PrimaryVertexProducer::PrimaryVertexProducer(const edm::ParameterSet& conf, const ONNXRuntime *cache) 
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
	const ONNXRuntime* onnxRuntime =   cache; 
       theTrackClusterizer = new ClusterizerinGNN(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"),onnxRuntime);
       useTransientTrackTime_ = true;
  } else {
    throw VertexException("PrimaryVertexProducer: unknown clustering algorithm: " + clusteringAlgorithm);
  }

  if (useTransientTrackTime_) {
    trkTimesToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimeResosLabel"));
    trackMTDTimeQualityToken =
        consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"));
    minTrackTimeQuality_ = conf.getParameter<double>("minTrackTimeQuality");
    trkMTDAssocToken = consumes<edm::ValueMap<int> >(conf.getParameter<edm::InputTag>("trackAssocSrc")); 
    MTDtimeToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tmtdSrc"));
    sigmaMTDtimeToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatmtdSrc"));
    pathLengthToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("pathmtd"));
    btlMatchChi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("btlMatchChi2Src"));
    btlMatchTime_Chi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("btlMatchTimeChi2Src"));
    etlMatchChi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("etlMatchChi2Src"));
    etlMatchTime_Chi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("etlMatchTimeChi2Src"
));
    trkTimePiToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tofPi"));
    trkTimeKToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tofK"));
    trkTimePToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tofP"));
    sigmaTrkTimePiToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatofpiSrc"));
    sigmaTrkTimeKToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatofkSrc"));
    sigmaTrkTimePToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatofpSrc"));
    npixBarrelToken = consumes<edm::ValueMap<int> >(conf.getParameter<edm::InputTag>("npixBarrelSrc"));
    npixEndcapToken = consumes<edm::ValueMap<int> >(conf.getParameter<edm::InputTag>("npixEndcapSrc"));
  }

  // select and configure the vertex fitters
  std::vector<edm::ParameterSet> vertexCollections =
      conf.getParameter<std::vector<edm::ParameterSet> >("vertexCollections");

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
std::unique_ptr<ONNXRuntime> PrimaryVertexProducer::initializeGlobalCache(const edm::ParameterSet &conf) {
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
}
void PrimaryVertexProducer::globalEndJob(const ONNXRuntime *cache) {}

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
    auto  trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken); 


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
    t_tks = (*theB).build(tks, beamSpot, trackTimes_, trackTimeResos_, trackMTDAssoc_, trackMTDTimes_, trackMTDTimesRes_, trackMTDTimeQualities_, pathlength_, btlMatch_, btlMatchTime_, etlMatch_, etlMatchTime_, trkPiTime_, trkKTime_, trkPTime_, sigmatrkPiTime_, sigmatrkKTime_, sigmatrkPTime_, npixbarrel_, npixendcap_);
    } else {
      t_tks = (*theB).build(tks, beamSpot, trackTimes_, trackTimeResos_, trackMTDAssoc_, trackMTDTimes_, trackMTDTimesRes_, trackMTDTimeQualities_, pathlength_, btlMatch_, btlMatchTime_, etlMatch_, etlMatchTime_, trkPiTime_, trkKTime_, trkPTime_, sigmatrkPiTime_, sigmatrkKTime_, sigmatrkPTime_, npixbarrel_, npixendcap_);
    }
  } else {
    t_tks = (*theB).build(tks, beamSpot);
  }

  // select tracks
  std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);
 
  // clusterize tracks in Z
  std::vector<TransientVertex>&&  clusters = theTrackClusterizer->vertices(seltks);
 
  if (fVerbose) {
    edm::LogPrint("PrimaryVertexProducer")
        << " bool useTransientTrackTime_ "<<useTransientTrackTime_<<" Clustering returned " << clusters.size() << " clusters from " << seltks.size() << " selected tracks";
  }

  // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);
    std::vector<TransientVertex> pvs;


    if (algorithm->pv_fitter == nullptr) {
      pvs = clusters;
    } else {
          pvs = algorithm->pv_fitter->fit(seltks, clusters, beamSpot, algorithm->useBeamConstraint);
      }

    if (algorithm->pv_time_estimator != nullptr) {
      algorithm->pv_time_estimator->fill_vertex_times(pvs);
    }

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
      ClusterizerinGNN::fillPSetDescription(psd4);

      psd0.ifValue(
          edm::ParameterDescription<std::string>("algorithm", "DA_vect", true),
          "DA_vect" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd1, true) or
              "DA2D_vect" >>
                  edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd2, true) or
	      "GNN2D_vect" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd4, true) or	  
              "gap" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkGapClusParameters", psd3, true));
    }
    desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
  }

  desc.add<bool>("isRecoveryIteration", false);
  desc.add<edm::InputTag>("recoveryVtxCollection", {""});
  desc.add<bool>("useMVACut", false);
  desc.add<double>("minTrackTimeQuality", 0.8);
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/model_v6version24June.onnx"))->setComment("Path to GNN (as ONNX model)");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryVertexProducer);
