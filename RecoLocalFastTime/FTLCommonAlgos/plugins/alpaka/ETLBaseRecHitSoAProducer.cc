#include <utility>

#include "DataFormats/FTLRecHitSoA/interface/alpaka/ETLBaseRecHitDeviceCollection.h"
#include "DataFormats/FTLDigiSoA/interface/alpaka/ETLDigiDeviceCollection.h"
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

#include "ETLBaseRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::etlrechit {

  using namespace ::etlrechit;

  class ETLBaseRecHitSoAProducer : public global::EDProducer<> {
  public:
    // constructor
    ETLBaseRecHitSoAProducer(edm::ParameterSet const& config)
        : EDProducer<>(config),
          digi_{consumes(config.getParameter<edm::InputTag>("digi"))},
          uncalibrh_{produces()},
          adcNBits_(config.getParameter<uint32_t>("adcNbits")),
          adcSaturation_(config.getParameter<double>("adcSaturation")),
          adcLSB_(adcSaturation_ / (1 << adcNBits_)),
          toaLSB_ns_(config.getParameter<double>("toaLSB_ns")),
          tdcWindowStart_(config.getParameter<double>("tdcWindowStart")),
          timeCorr_p0_(config.getParameter<double>("timeCorr_p0")),
          timeCorr_p1_(config.getParameter<double>("timeCorr_p1")),
          timeCorr_p2_(config.getParameter<double>("timeCorr_p2")),
          timeCorr_p3_(config.getParameter<double>("timeCorr_p3")) {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("digi");
      desc.add<uint32_t>("adcNbits");
      desc.add<double>("adcSaturation");
      desc.add<double>("toaLSB_ns");
      desc.add<double>("tdcWindowStart");
      desc.add<double>("timeCorr_p0");
      desc.add<double>("timeCorr_p1");
      desc.add<double>("timeCorr_p2");
      desc.add<double>("timeCorr_p3");
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(edm::StreamID sid, device::Event& event, device::EventSetup const& setup) const override {
      // NB should be inserted a method to retrieve calibrations, now they are fixed to default values
      // Get the digi from the Event.
      etldigi::ETLDigiDeviceCollection const& digi = event.get(digi_);  // this should match Claudio Class name

      // Allocate a new SoA for the uncalibrh jets. // same number of elements we have in input
      ETLBaseRecHitDeviceCollection uncalibrh(event.queue(), digi.view().metadata().size());

      // Apply the corrections and fill the new SoA. // these launch the kernel, and will run on gpu async
      ETLBaseRecHitSoAProducerAlgo::fromDigiToBase(
          event.queue(), digi.view(), uncalibrh.view(), adcNBits_, adcSaturation_, adcLSB_,
          toaLSB_ns_, tdcWindowStart_, timeCorr_p0_, timeCorr_p2_, timeCorr_p1_, timeCorr_p3_);

      // Move the SoA with the uncalibrh jets into the Event.
      event.emplace(uncalibrh_, std::move(uncalibrh));
    }

  private:
    const device::EDGetToken<etldigi::ETLDigiDeviceCollection> digi_;
    const device::EDPutToken<ETLBaseRecHitDeviceCollection> uncalibrh_;
    const uint32_t adcNBits_;
    const double adcSaturation_;
    const double adcLSB_;
    const double toaLSB_ns_;
    const double tdcWindowStart_;
    const double timeCorr_p0_;
    const double timeCorr_p1_;
    const double timeCorr_p2_;
    const double timeCorr_p3_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::etlrechit

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(etlrechit::ETLBaseRecHitSoAProducer);
