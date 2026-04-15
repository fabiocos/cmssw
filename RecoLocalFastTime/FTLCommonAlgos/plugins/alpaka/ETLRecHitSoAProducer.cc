#include <utility>

#include "CommonTools/Utils/interface/FormulaEvaluator.h"
#include "DataFormats/FTLRecHitSoA/interface/alpaka/ETLBaseRecHitDeviceCollection.h"
#include "DataFormats/FTLRecHitSoA/interface/alpaka/ETLRecHitDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "ETLRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::etlrechit {

  using namespace ::etlrechit;

  class ETLRecHitSoAProducer : public global::EDProducer<> {
  public:
    // constructor
    ETLRecHitSoAProducer(edm::ParameterSet const& config)
        : EDProducer<>(config),
          baserh_{consumes(config.getParameter<edm::InputTag>("baserh"))},
          rh_{produces()},
          thresholdToKeep_(config.getParameter<double>("thresholdToKeep")),
          calibration_(config.getParameter<double>("calibrationConstant")) {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("baserh");
      desc.add<double>("thresholdToKeep");
      desc.add<double>("calibrationConstant");
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(edm::StreamID sid, device::Event& event, device::EventSetup const& setup) const override {
      // NB should be inserted a method to retrieve calibrations, now they are fixed to default values
      // Get the base from the Event.
      ETLBaseRecHitDeviceCollection const& baserh = event.get(baserh_);

      // Allocate a new SoA for the rechit.
      ETLRecHitDeviceCollection rh(event.queue(), baserh.view().metadata().size());

      // Apply the corrections and fill the new SoA. // these launch the kernel, and will run on gpu async
      ETLRecHitSoAProducerAlgo::fromBaseToReco(
          event.queue(), baserh.view(), rh.view(), thresholdToKeep_, calibration_);

      // Move the SoA with the rh into the Event.
      event.emplace(rh_, std::move(rh));
    }

  private:
    const device::EDGetToken<ETLBaseRecHitDeviceCollection> baserh_;
    const device::EDPutToken<ETLRecHitDeviceCollection> rh_;
    //edm::ESGetToken<MTDTimeCalib, MTDTimeCalibRecord> tcToken_;
    const double thresholdToKeep_;
    const double calibration_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::etlrechit

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(etlrechit::ETLRecHitSoAProducer);
