#include <utility>

#include "DataFormats/FTLRecHitSoA/interface/alpaka/BTLBaseRecHitDeviceCollection.h"
#include "DataFormats/FTLRecHitSoA/interface/alpaka/BTLRecHitDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoLocalFastTime/FTLCommonAlgos/interface/MTDTimeCalib.h"

#include "BTLRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit {

  using namespace ::btlrechit;

  class BTLRecHitSoAProducer : public stream::EDProducer<> {
  public:
    // constructor
    BTLRecHitSoAProducer(edm::ParameterSet const& config)
        : EDProducer<>(config),
          baserh_{consumes(config.getParameter<edm::InputTag>("baserh"))},
          rh_{produces()},
          invLightSpeedLYSO_(config.getParameter<double>("invLightSpeedLYSO")),
          c_LYSO_(1. / invLightSpeedLYSO_),
          thresholdToKeep_(config.getParameter<double>("thresholdToKeep")),
          calibration_(config.getParameter<double>("calibrationConstant")) {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("baserh");
      desc.add<double>("invLightSpeedLYSO");
      desc.add<double>("thresholdToKeep");
      desc.add<double>("calibrationConstant");
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, device::EventSetup const& setup) override {
      // NB should be inserted a method to retrieve calibrations, now they are fixed to default values
      // Get the base from the Event.
      BTLBaseRecHitDeviceCollection const& baserh = event.get(baserh_);

      // Allocate a new SoA for the rechit.
      BTLRecHitDeviceCollection rh(event.queue(), baserh.view().metadata().size());

      // Apply the corrections and fill the new SoA. // these launch the kernel, and will run on gpu async
      BTLRecHitSoAProducerAlgo::fromBaseToReco(
          event.queue(), baserh.view(), rh.view(), c_LYSO_, thresholdToKeep_, calibration_);

      // Move the SoA with the rh into the Event.
      event.emplace(rh_, std::move(rh));
    }

  private:
    const device::EDGetToken<BTLBaseRecHitDeviceCollection> baserh_;
    const device::EDPutToken<BTLRecHitDeviceCollection> rh_;
    const double invLightSpeedLYSO_;
    const double c_LYSO_;
    const double thresholdToKeep_;
    const double calibration_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(btlrechit::BTLRecHitSoAProducer);
