#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/ClusterizerinGNN.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath>
#include <cassert>
#include <limits>
#include <iomanip>
#include "FWCore/Utilities/interface/isFinite.h"
#include "vdt/vdtMath.h"

using namespace std;

const int FEATURE_SHAPE_GNN_V1 = 21;
const int FEATURE_SHAPE_MLP = 30;
const int NUM_EDGE_FEATURES = 12;
enum class DNNVersion {
    GNN_V1,
    MLP_EDGE_FEATURES,
    MLP_NO_EDGE_FEATURES
};

// Define constants for DNN versions
std::unordered_map<std::string, DNNVersion> DNN_VERSION_MAP = {
    {"gnn_v1", DNNVersion::GNN_V1},
    {"mlp_edge_features", DNNVersion::MLP_EDGE_FEATURES},
    {"mlp_no_edge_features", DNNVersion::MLP_NO_EDGE_FEATURES}
};
// Define a map to store DNN input configurations
std::unordered_map<DNNVersion, std::pair<std::vector<std::string>, int>> dnnInputConfigurations = {
    {DNNVersion::GNN_V1, {{"features", "edge_index", "edge_features"}, FEATURE_SHAPE_GNN_V1}},
    {DNNVersion::MLP_EDGE_FEATURES, {{"features", "edge_index", "edge_features"}, FEATURE_SHAPE_MLP}},
    {DNNVersion::MLP_NO_EDGE_FEATURES, {{"features", "edge_index"}, FEATURE_SHAPE_MLP}}
};


//ClusterizerinGNN::ClusterizerinGNN(const edm::ParameterSet& conf, const ONNXRuntime *cache) 
ClusterizerinGNN::ClusterizerinGNN(const edm::ParameterSet& conf, const ONNXRuntime* onnxRuntime)
:     onnxRuntime_(onnxRuntime),
      nnVersion_(conf.getParameter<std::string>("nnVersion")),
      nnWorkingPoint_(conf.getParameter<double>("nnWorkingPoint")) {
  // some defaults to avoid uninitialized variables
  verbose_ = conf.getUntrackedParameter<bool>("verbose", false);
  zSep = conf.getParameter<double>("zSeparation");
  d0CutOff_ = conf.getParameter<double>("d0CutOff");
  vertexSize_ = conf.getParameter<double>("vertexSize");
  if (verbose_) {
    std::cout << "TrackClusterizerInZ:  algorithm=gap, zSeparation=" << std::endl;
  }
}

void ClusterizerinGNN::globalEndJob(const ONNXRuntime *cache) {}

namespace {
  inline double local_exp(double const& inp) { return vdt::fast_exp(inp); }
}
std::unique_ptr<ONNXRuntime> ClusterizerinGNN::initializeGlobalCache(const edm::ParameterSet &conf) {
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
}
std::unique_ptr<TrackGraph>  ClusterizerinGNN::produce_tracks_graph(const std::vector<reco::TransientTrack>& transientTracks) const {
    std::vector<Node> allNodes;
    for (size_t i = 0; i < transientTracks.size(); ++i) {
    Node tNode(i);
    const reco::TransientTrack& ttrack_node = transientTracks[i];
   for (size_t j = i + 1; j < transientTracks.size(); ++j) {
    const reco::TransientTrack& ttrack = transientTracks[j];
        double neighbour_cut = std::abs(ttrack_node.track().vz() - ttrack.track().vz());
        if (neighbour_cut < 0.3 ) {
          tNode.addInner(j);
    }
  }
    allNodes.push_back(tNode);
  }
     auto resultGraph = std::make_unique<TrackGraph>(allNodes);

    return resultGraph;
}


vector<TransientVertex> ClusterizerinGNN::vertices(const vector<reco::TransientTrack>& tracks) const{
   cms::Ort::FloatArrays data;
   std::vector<std::vector<unsigned int>> linkedTrackIdToInputTrackId;
   std::vector<std::vector<reco::TransientTrack>> seltks;
   std::vector<TransientVertex> clusters;
   std::vector<float> features;
    std::vector<float> edge_features;
    std::vector<std::vector<int64_t>> input_shapes;
    std::vector<float> edges_src;
    std::vector<float> edges_dst;

// *************************************
   DNNVersion versionEnum = DNN_VERSION_MAP[nnVersion_];

    // Access the input configuration based on the DNN version
    std::vector<std::string> input_names;
    int shapeFeatures = 0;
    auto dnnConfig = dnnInputConfigurations.find(versionEnum);
    if (dnnConfig != dnnInputConfigurations.end()) {
        // Valid DNN version found, retrieve input names and shape features
        input_names = dnnConfig->second.first;
        shapeFeatures = dnnConfig->second.second;
    } else {
        // Invalid DNN version
        edm::LogError("Vertex Producer") << "Architecture not defined: " << nnVersion_ ;
    }

    std::cout << "Input names:";
    for (const auto& name : input_names) {
        std::cout << " " << name;
    }
    std::cout << std::endl;
    std::cout << "Shape of features: " << shapeFeatures << std::endl;
    int N = tracks.size();
    // Network input shapes.
    auto trackgraph = produce_tracks_graph(tracks);
    TrackGraph *trkgrp = trackgraph.get();

   for (size_t i = 0; i < tracks.size() ; ++i) {
        const reco::TransientTrack& ttrack1 = tracks[i];
        float trackpt = ttrack1.track().pt();
        float tracketa = ttrack1.track().eta();
        float trackphi = ttrack1.track().phi();
        float trackvz = ttrack1.track().vz();
        float track_dzError = ttrack1.track().dzError();
        float track_ndof = ttrack1.track().ndof();
        float track_chi2 = ttrack1.track().chi2();
        float mtdtime = ttrack1.MTDtime();
        float mtdtimeErr = ttrack1.MTDtimeErr();
        float mvaquality = ttrack1.MVAquality();
        float pathlength = ttrack1.pathLength();
        float btlMatchchi2 = ttrack1.btlMatch_chi2();
        float btlMatchTimechi2 = ttrack1.btlMatchTime_chi2();
        float etlMatchchi2 = ttrack1.etlMatch_chi2();
        float etlMatchTimechi2 = ttrack1.etlMatchTime_chi2();
        float Time_pi_hypo = ttrack1.trackTime_pi();
        float Time_k_hypo = ttrack1.trackTime_k();
        float Time_p_hypo = ttrack1.trackTime_p();
        float npixBarrel = ttrack1.nPixBarrel();
        float npixEndcap = ttrack1.nPixEndcap();
        float sigmaTime_p = ttrack1.sigma_time_p();

	//Variables we have in model : pt, eta, phi, z_pca, t_pi, t_k, t_p, mva, btl*4, path lenght, npixBarrel, npix Ebdcap, mtdtime, sigma mtdtime, dz, sigma tof_p, trk_ndof, trk_chi2
        features.push_back(trackpt);
        features.push_back(tracketa);
        features.push_back(trackphi);
        features.push_back(trackvz);
        features.push_back(Time_pi_hypo);
        features.push_back(Time_k_hypo);
        features.push_back(Time_p_hypo);
        features.push_back(mvaquality);
        features.push_back(btlMatchchi2);
        features.push_back(btlMatchTimechi2);
        features.push_back(etlMatchchi2);
        features.push_back(etlMatchTimechi2);
        features.push_back(pathlength);
        features.push_back(npixBarrel);
        features.push_back(npixEndcap);
        features.push_back(mtdtime);
        features.push_back(track_dzError);
        features.push_back(mtdtimeErr);
        features.push_back(sigmaTime_p);
        features.push_back(track_ndof);
        features.push_back(track_chi2);

    for (auto &j : trkgrp->getNode(i).getInner()){
        const reco::TransientTrack& ttrack2 = tracks[j];


        double neighbour_cut = std::abs(ttrack1.track().vz() - ttrack2.track().vz());
        if (neighbour_cut < 0.3 ) {
        edges_src.push_back(static_cast<float>(j)); // this should be neighbour of the track 
        edges_dst.push_back(static_cast<float>(i));

	double pt_diff = std::abs(ttrack1.track().pt() - ttrack2.track().pt());
        double zpca_diff = std::abs(ttrack1.track().vz() - ttrack2.track().vz());
        double dz_sign = std::abs(ttrack1.track().vz() - ttrack2.track().vz())/sqrt(((ttrack1.track().dzError())*(ttrack1.track().dzError()))+((ttrack2.track().dzError())*(ttrack2.track().dzError())));
        float time_pi_diff = std::abs(ttrack1.trackTime_pi() -  ttrack2.trackTime_pi());
        float time_k_diff = std::abs(ttrack1.trackTime_k() -  ttrack2.trackTime_k());
        float time_p_diff = std::abs(ttrack1.trackTime_p() -  ttrack2.trackTime_p());
        float time_pi_k_diff = std::abs(ttrack1.trackTime_pi() - ttrack2.trackTime_k());
        float time_pi_p_diff = std::abs(ttrack1.trackTime_pi() - ttrack2.trackTime_p());
        float time_k_pi_diff = std::abs(ttrack1.trackTime_k() - ttrack2.trackTime_pi());
        float time_k_p_diff = std::abs(ttrack1.trackTime_k() - ttrack2.trackTime_p());
        float time_p_pi_diff = std::abs(ttrack1.trackTime_p() - ttrack2.trackTime_pi());
        float time_p_k_diff = std::abs(ttrack1.trackTime_p() - ttrack2.trackTime_k());	

        edge_features.push_back(pt_diff);
         edge_features.push_back(zpca_diff);
         edge_features.push_back(dz_sign);
         edge_features.push_back(time_pi_diff);
         edge_features.push_back(time_k_diff);
         edge_features.push_back(time_p_diff);
         edge_features.push_back(time_pi_k_diff);
         edge_features.push_back(time_pi_p_diff);
         edge_features.push_back(time_k_pi_diff);
         edge_features.push_back(time_k_p_diff);
         edge_features.push_back(time_p_k_diff);
         edge_features.push_back(time_p_pi_diff);
    } // if neighbour loop 
   } // loop for edges 
  } // loop for transient tracks
   auto numEdges = static_cast<int>(edges_src.size());
   if (numEdges < 1) {
       edm::LogPrint("GNN Clusterizer")  << "No edges for the event - no linking is done." ;
        linkedTrackIdToInputTrackId.resize(N);
        std::vector<unsigned int> linkedTracks;

        for (int track_id = 0; track_id < N; track_id++)
        {
            linkedTrackIdToInputTrackId[track_id].push_back(track_id);
            linkedTracks.push_back(seltks.size());
        }
    }
    data.clear();
    input_shapes.push_back({1, N, shapeFeatures});
    data.emplace_back(features);
    input_shapes.push_back({1, 2, numEdges});
    data.emplace_back(edges_src);
    for (auto &dst : edges_dst)
    {
            data.back().push_back(dst);
    }
    input_shapes.push_back({1, numEdges, NUM_EDGE_FEATURES});
    data.emplace_back(edge_features);
// run prediction
 if (numEdges >=1 ) {
    auto edge_predictions = onnxRuntime_->run(input_names, data, input_shapes)[0];

  // Update graph weights
    for (int i = 0; i < static_cast<int>(edge_predictions.size()); i++){

            trkgrp->setEdgeWeight(data[1][i], data[1][numEdges + i], edge_predictions[i]);
    }
    auto connected_components = trkgrp->findSubComponents(nnWorkingPoint_);
    int id = 0;
    std::vector<TrackCollection> outTracks; 
    linkedTrackIdToInputTrackId.resize(connected_components.size());
    for (auto &component : connected_components)
    {
        std::vector<unsigned int> linkedTracks;

	TrackCollection outTrack;
	std::vector<reco::TransientTrack> tracks1;
        for (auto &track_id : component)
        {
            linkedTrackIdToInputTrackId[id].push_back(track_id);
            outTrack.addItem(tracks[track_id]);

	    tracks1.push_back(tracks[track_id]);
        }
        id++;
	outTracks.push_back(outTrack);
	seltks.push_back(tracks1);
    }	


  unsigned int nv = connected_components.size(); 

  std::vector<double> zvtx(connected_components.size());
    std::vector<double> tvtx(connected_components.size()); 
    std::vector<double> dz2;
    std::vector<double> zpca;
    std::vector<double> t_tkwt;
    std::vector<double> t_dz2;
 for (unsigned int k =0; k < nv; k++) {
    double sumwz = 0;
    double sumwt = 0;
    double sumw_z = 0;
    double sumw_t = 0;
    const auto &tracks_sel = seltks[k];  
 
    dz2.resize(tracks_sel.size());
    zpca.resize(tracks_sel.size());
    t_tkwt.resize(tracks_sel.size());
    t_dz2.resize(tracks_sel.size());
     for (size_t index = 0; index < tracks_sel.size(); ++index)  {
     const auto &track = tracks_sel[index];
      if (d0CutOff_ > 0) {
      Measurement1D atIP = track.stateAtBeamLine().transverseImpactParameter();  // error contains beamspot
      t_tkwt[index] = 1. / (1. + local_exp(std::pow(atIP.value() / atIP.error(), 2) -
                                    std::pow(d0CutOff_, 2)));  // reduce weight for high ip tracks
      if (edm::isNotFinite(t_tkwt[index]) || t_tkwt[index] < std::numeric_limits<double>::epsilon()) {
        edm::LogWarning("GNN Clusterizer") << "rejected track t_tkwt " << t_tkwt[index];
        continue;  // usually is > 0.99
      }
    } // d0 cutoff 
      auto const& t_mom = track.stateAtBeamLine().trackStateAtPCA().momentum();
      reco::BeamSpot beamspot = track.stateAtBeamLine().beamSpot();
      t_dz2[index] = std::pow(track.track().dzError(), 2)  // track errror
                   + (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
                         std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
                   + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
      t_dz2[index] = 1. / t_dz2[index];
      if (edm::isNotFinite(t_dz2[index]) || t_dz2[index] < std::numeric_limits<double>::min()) {
      edm::LogWarning("GNN Clusterizer") << "rejected track t_dz2 " << t_dz2[index];
      continue;
    }

      zpca[index] = track.track().vz(); 

      double w_z = t_tkwt[index] * t_dz2[index];
      double w_t = t_tkwt[index] * t_dz2[index];
      sumwz +=  w_z * zpca[index];
      sumwt += w_t * zpca[index];
      sumw_z += w_z;
      sumw_t += w_t;
    }
    zvtx[k] = sumwz / sumw_z;
 }


  GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);
    for (unsigned int k = 0; k < nv; k++) {

     if (!linkedTrackIdToInputTrackId[k].empty()) {
      GlobalPoint pos(0, 0, zvtx[k]);
      std::vector<reco::TransientTrack> vertexTracks;
      for (const auto& track : outTracks[k].tracks) {
        // Here, 'track' is each TransientTrack in the k-th TrackCollection
        vertexTracks.push_back(track);
    }

      TransientVertex v(pos, dummyError, vertexTracks, 0);
      clusters.push_back(v);
    }
  }
     
  }// loop for numEdges 

  return clusters;
}



std::vector<std::vector<reco::TransientTrack>> ClusterizerinGNN::clusterize(const std::vector<reco::TransientTrack>& tracks) const {
    // Implementation of clusterize function
    std::vector<std::vector<reco::TransientTrack>> clusters;
    // Example implementation:
    // clusters.push_back(tracks); // For demonstration, replace with actual logic
    return clusters;
}

void ClusterizerinGNN::fillPSetDescription(edm::ParameterSetDescription& desc) {
  DAClusterizerInZ_vect::fillPSetDescription(desc);	
  desc.add<double>("zSeparation", 1.0);
  desc.addUntracked<bool>("verbose", false);
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/model_v6version24June.onnx"))->setComment("Path to GNN (as ONNX model)");
  desc.add<std::string>("nnVersion", "gnn_v1") // gnn_v1, mlp_no_edge_features, mlp_edge_features
        ->setComment("GNN version tag.");
   desc.add<double>("nnWorkingPoint", 0.80)
        ->setComment("Working point of the GNN (in [0, 1]). GNN score above WP will attempt to supercluster.");
}
