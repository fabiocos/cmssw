/**
 * TrackFeatureSource - Alpaka producer to create TrackFeaturesDeviceCollection
 *
 * This follows the DataSource.cc pattern from PyTorchAlpakaTest:
 * - Uses ALPAKA_ACCELERATOR_NAMESPACE::stream::EDProducer<>
 * - Produces TrackFeaturesDeviceCollection on device
 * 
 * For testing, generates random features. In production, would
 * consume tracks and extract real features.
 */

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/VertexGNNReco/interface/VertexGNNSoA.h"
#include "DataFormats/VertexGNNReco/interface/alpaka/VertexGNNDeviceCollection.h"

#include <random>

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  class TrackFeatureSource : public stream::EDProducer<> {
  public:
    explicit TrackFeatureSource(const edm::ParameterSet& params)
        : EDProducer<>(params),
          trackFeaturesToken_{produces()},
          numTracks_(params.getParameter<int>("numTracks")),
          seed_(params.getParameter<int>("seed")),
          verbose_(params.getUntrackedParameter<bool>("verbose", false)) {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<int>("numTracks", 50);
      desc.add<int>("seed", 42);
      desc.addUntracked<bool>("verbose", false);
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, const device::EventSetup& eventSetup) override {
      const auto N = numTracks_;

      if (verbose_) {
        edm::LogInfo("TrackFeatureSource") << "Generating " << N << " track features";
      }

      // Allocate output collection on device
      auto trackFeatures = TrackFeaturesDeviceCollection(N, event.queue());
      auto view = trackFeatures.view();

      // Generate random features (for testing)
      // In production, this would consume actual tracks
      std::mt19937 gen(seed_ + event.id().event());
      std::normal_distribution<float> vz_dist(0.0f, 5.0f);
      std::normal_distribution<float> dz_dist(0.1f, 0.05f);
      std::uniform_real_distribution<float> pt_dist(0.5f, 50.0f);
      std::uniform_real_distribution<float> eta_dist(-2.5f, 2.5f);
      std::uniform_real_distribution<float> mva_dist(0.5f, 1.0f);
      std::uniform_real_distribution<float> time_dist(0.0f, 1.0f);

      for (int i = 0; i < N; ++i) {
        float vz = vz_dist(gen);
        float dz = std::abs(dz_dist(gen));
        float pt = pt_dist(gen);
        float eta = eta_dist(gen);
        float mva = mva_dist(gen);
        float pl = std::abs(vz / (dz > 0 ? dz : 0.01f));

        // Timing features (some tracks have MTD)
        bool has_mtd = time_dist(gen) > 0.3f;
        float t_pi = has_mtd ? std::normal_distribution<float>(0.0f, 0.05f)(gen) : 0.0f;
        float t_k = has_mtd ? std::normal_distribution<float>(0.0f, 0.05f)(gen) : 0.0f;
        float t_p = has_mtd ? std::normal_distribution<float>(0.0f, 0.05f)(gen) : 0.0f;
        float s_pi = has_mtd ? 0.03f : 0.2f;
        float s_k = has_mtd ? 0.03f : 0.2f;
        float s_p = has_mtd ? 0.03f : 0.2f;
        float has_time = has_mtd ? 1.0f : 0.0f;

        // Store features in SoA
        view[i].vz() = vz;
        view[i].dz() = dz;
        view[i].pt() = pt;
        view[i].eta() = eta;
        view[i].mva() = mva;
        view[i].pl() = pl;
        view[i].t_pi() = t_pi;
        view[i].t_k() = t_k;
        view[i].t_p() = t_p;
        view[i].s_pi() = s_pi;
        view[i].s_k() = s_k;
        view[i].s_p() = s_p;
        view[i].has_time() = has_time;
      }

      // Put output into event
      event.emplace(trackFeaturesToken_, std::move(trackFeatures));
    }

  private:
    const device::EDPutToken<TrackFeaturesDeviceCollection> trackFeaturesToken_;
    const int numTracks_;
    const int seed_;
    const bool verbose_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

DEFINE_FWK_ALPAKA_MODULE(vertexgnn::TrackFeatureSource);
