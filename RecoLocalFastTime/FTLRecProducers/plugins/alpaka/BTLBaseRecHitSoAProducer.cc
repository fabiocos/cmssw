#include <utility>

#include "DataFormats/FTLRecHitSoA/interface/alpaka/BTLBaseRecHitDeviceCollection.h"
#include "DataFormats/FTLDigiSoA/interface/BTLDigiHostCollection.h"
#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "BTLBaseRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit {

  using namespace ::btlrechit;

  class BTLBaseRecHitSoAProducer : public stream::EDProducer<> {
  public:
    // constructor
    BTLBaseRecHitSoAProducer(edm::ParameterSet const& config)
        : EDProducer<>(config),
          digi_(consumes<::btldigi::BTLDigiHostCollection>(config.getParameter<edm::InputTag>("digi"))),
          uncalibrh_{produces()},
          npeToADC0_(config.getParameter<double>("npeToADC0")),
          npeToADC1_(config.getParameter<double>("npeToADC1")),
          npeSaturationCorr0_(config.getParameter<double>("npeSaturationCorr0")),
          npeSaturationCorr1_(config.getParameter<double>("npeSaturationCorr1")),
          npePerMeV_(config.getParameter<double>("npePerMeV")) {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("digi");
      desc.add<double>("npeToADC0");
      desc.add<double>("npeToADC1");
      desc.add<double>("npeSaturationCorr0");
      desc.add<double>("npeSaturationCorr1");
      desc.add<double>("npePerMeV");
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, device::EventSetup const& setup) override {
      // NB should be inserted a method to retrieve calibrations, now they are fixed to default values
      // Get the digi from the Event.
      auto const& hostDigi = event.get(digi_);  // SoA DIGI stored in the event
      auto const N = hostDigi.const_view().metadata().size();

      // Copy input to device for GPU inference
      btldigi::BTLDigiDeviceCollection deviceDigi(event.queue(), N);
      alpaka::memcpy(event.queue(), deviceDigi.buffer(), hostDigi.buffer());

      // Allocate a new SoA for the uncalibrh jets. // same number of elements we have in input
      BTLBaseRecHitDeviceCollection uncalibrh(event.queue(), N);

      // Apply the corrections and fill the new SoA. // these launch the kernel, and will run on gpu async
      BTLBaseRecHitSoAProducerAlgo::fromDigiToBase(event.queue(),
                                                   deviceDigi.view(),
                                                   uncalibrh.view(),
                                                   npeToADC0_,
                                                   npeToADC1_,
                                                   npeSaturationCorr0_,
                                                   npeSaturationCorr1_,
                                                   npePerMeV_);

      // Move the SoA with the uncalibrh jets into the Event.
      event.emplace(uncalibrh_, std::move(uncalibrh));
    }

  private:
    const edm::EDGetTokenT<::btldigi::BTLDigiHostCollection> digi_;
    const device::EDPutToken<BTLBaseRecHitDeviceCollection> uncalibrh_;
    const double npeToADC0_;
    const double npeToADC1_;
    const double npeSaturationCorr0_;
    const double npeSaturationCorr1_;
    const double npePerMeV_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(btlrechit::BTLBaseRecHitSoAProducer);
