#include <cstdio>

#include <alpaka/alpaka.hpp>

#include "DataFormats/FTLRecHitSoA/interface/alpaka/BTLBaseRecHitDeviceCollection.h"
#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "BTLBaseRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit {

  using namespace ::btlrechit;

  ALPAKA_FN_ACC uint8_t rowFromId(uint32_t rawId) {  // NB working only with new geometry
    int crys = ((rawId >> BTLDetId::kBTLCrystalOffset) & BTLDetId::kBTLCrystalMask);
    uint8_t row = crys % BTLDetId::kCrystalsPerModuleV2;
    return row;
  }

  ALPAKA_FN_ACC float TcoarseTfineToTime(
      uint32_t rawId, uint8_t chID, uint8_t TACID, uint16_t tcoarse, uint16_t tfine, bool isT1) {
    // tdc calibration parameters
    // (to be modified: these parameters are evaluated by channel and stored in parquet files)
    static constexpr float a0 = 57.244545;
    static constexpr float a1 = 511.27832;
    static constexpr float a2 = -7.8838577;
    static constexpr float t0 = -0.048343264;

    float const qT = (-a1 + sqrt(a1 * a1 - 4.0 * (a0 - float(tfine)) * a2)) / (2.0 * a2);
    float const time = tcoarse - qT - t0;
    return time;
  }

  ALPAKA_FN_ACC uint32_t
  QfineToADC(uint32_t rawId, uint8_t chID, uint8_t TACID, uint16_t qfine, float time1, uint16_t timeEndQ) {
    // qdc calibration parameters
    // (to be modified: these parameters are evaluated by channel and stored in parquet files)
    static constexpr float p0 = 49.542229;
    static constexpr float p1 = -0.323424;
    static constexpr float p2 = 0.062578;
    static constexpr float p3 = -0.002484;
    static constexpr float p4 = 0.0;
    static constexpr float p5 = 0.0;
    static constexpr float p6 = 0.0;
    static constexpr float p7 = 0.0;
    static constexpr float p8 = 0.0;
    static constexpr float p9 = 0.0;
    float const ti = float(timeEndQ) - time1;

    uint32_t pedestal = (  // check the type
        p0 + p1 * ti + p2 * ti * ti + p3 * ti * ti * ti + p4 * ti * ti * ti * ti + p5 * ti * ti * ti * ti * ti +
        p6 * ti * ti * ti * ti * ti * ti + p7 * ti * ti * ti * ti * ti * ti * ti +
        p8 * ti * ti * ti * ti * ti * ti * ti * ti + p9 * ti * ti * ti * ti * ti * ti * ti * ti * ti);

    const uint32_t adc = qfine - pedestal;

    return adc;
  }

  ALPAKA_FN_ACC float timeWalkCorr(float amp) {
    float tdcLSB_ns = 0.020;
    float corr = 1.9e6 / 0.020 * pow(9.389e5 / 0.0348 * (amp + 22.5), -0.663) - 7.5e-4 * amp - 3.5e-3;
    return tdcLSB_ns * corr;
  }

  class BTLdigiToBaseKernel {
  public:
    ALPAKA_FN_ACC void operator()(
        Acc1D const& acc,
        ::btldigi::BTLDigiSoA::ConstView input,
        BTLBaseRecHitSoA::View output,
        const double npeToADC0_,
        const double npeToADC1_,
        const double npeSaturationCorr0_,
        const double npeSaturationCorr1_,
        const double npePerMeV_) const {  // when condformat for calib ready, add also tdc and qdc in inputs

      static constexpr uint32_t adcBitSaturation_ = 1023;
      static constexpr float tclock = 6.25;
      // make a strided loop over the kernel grid, covering up to "size" elements
      for (int32_t i : cms::alpakatools::uniform_elements(acc, input.metadata().size())) {
        auto entry = input[i];
        // here you should call your functions to apply TDC and QDC, reading timecoarse, fine, ... etc from digi

        // for the times at first and second th, still in clock units
        // atm tdc and qdc calibs are fixed to dummy values for each channel, hence rawId, ch, and the bool to select branch 1 or 2 are not used.
        auto time1R =
            TcoarseTfineToTime(entry.rawId(), entry.chIDR(), entry.TACIDR(), entry.T1coarseR(), entry.T1fineR(), true);
        auto time1L =
            TcoarseTfineToTime(entry.rawId(), entry.chIDL(), entry.TACIDL(), entry.T1coarseL(), entry.T1fineL(), true);

        auto time2R =
            TcoarseTfineToTime(entry.rawId(), entry.chIDR(), entry.TACIDR(), entry.T2coarseR(), entry.T2fineR(), false);
        auto time2L =
            TcoarseTfineToTime(entry.rawId(), entry.chIDL(), entry.TACIDL(), entry.T2coarseL(), entry.T2fineL(), false);

        // from qfine to energy in adc, NB you need to pass calibrated time
        auto ampL =
            QfineToADC(entry.rawId(), entry.chIDL(), entry.TACIDL(), entry.ChargeL(), time1L, entry.EOIcoarseL());
        auto ampR =
            QfineToADC(entry.rawId(), entry.chIDR(), entry.TACIDR(), entry.ChargeR(), time1R, entry.EOIcoarseR());

        uint8_t row = rowFromId(entry.rawId());

        // flags for the usability of the channel uint_8: atm 2 bit are used
        //  first bit is channel has signal (1) or not (0)
        //  second bit channel was saturated (1) or not (0)
        uint8_t flagsL = 0;
        uint8_t flagsR = 0;

        if (ampL > 0)
          flagsL |= 0x1;
        if (ampL == adcBitSaturation_)
          flagsL |= (0x1 << 1);
        if (ampR > 0)
          flagsR |= 0x1;
        if (ampR == adcBitSaturation_)
          flagsR |= (0x1 << 1);

        // detId from rawId
        DetId detId(entry.rawId());

        // convert from clock units to ps
        time1R *= tclock;
        time1L *= tclock;
        time2R *= tclock;
        time2L *= tclock;

        // amp walk correction
        auto corrR = timeWalkCorr(ampR);
        auto corrL = timeWalkCorr(ampL);

        // convert from clock units to ps, apply amp walk corrections
        auto time1Rcorr = time1R - corrR;
        auto time1Lcorr = time1L - corrL;

        auto time2Rcorr = time2R - corrR;
        auto time2Lcorr = time2L - corrL;

        // converting the energy from ADC to energy
        auto energyR = float((float(ampR) - npeToADC0_) / npeToADC1_);
        // Correction for SiPM saturation (just invert the function used to model this effect in BTLElectronicsSim)
        float dR = npeSaturationCorr1_ * npeSaturationCorr1_ + 4. * npeSaturationCorr0_ * energyR;
        energyR = (-npeSaturationCorr1_ + sqrt(dR)) / (2. * (npeSaturationCorr0_));
        energyR /= npePerMeV_;

        auto energyL = float((float(ampL) - npeToADC0_) / npeToADC1_);
        // Correction for SiPM saturation (just invert the function used to model this effect in BTLElectronicsSim)
        float dL = npeSaturationCorr1_ * npeSaturationCorr1_ + 4. * npeSaturationCorr0_ * energyL;
        energyL = (-npeSaturationCorr1_ + sqrt(dL)) / (2. * (npeSaturationCorr0_));
        energyL /= npePerMeV_;

#ifdef EDM_ML_DEBUG
        printf("Base recHit SoA with raw id %i \n", entry.rawId());
        printf("Time 1 before corrections L,R (%f, %f) - ", time1L, time1R);
        printf("Amp Walk Corrections L,R (%f , %f) -->  ", corrL, corrR);
        printf("Time 1 after corrections L,R (%f, %f) \n", time1Lcorr, time1Rcorr);

        printf("Time 2 before corrections L,R (%f, %f) - ", time2L, time2R);
        printf("Amp Walk Corrections L,R (%f , %f) -->  ", corrL, corrR);
        printf("Time 2 after corrections L,R (%f, %f) \n", time2Lcorr, time2Rcorr);

        printf("Amplidute in ADC L,R (%i, %i) - ", ampL, ampR);
        printf("converting to energy L,R (%f, %f) --> ", npeToADC0_, invADCPerMeV_);
        printf("Energy in MeV L,R (%f, %f) \n", energyL, energyR);
#endif

        // fill the base rechit
        output[i] = {
            detId,
            row,
            time1Rcorr,  // in ns
            time2Rcorr,
            energyR,  // energy
            entry.IdleTimeR(),
            flagsR,
            time1Lcorr,  // in ns
            time2Lcorr,
            energyL,  // energy
            entry.IdleTimeL(),
            flagsL,

        };
      }
    }
  };

  void BTLBaseRecHitSoAProducerAlgo::fromDigiToBase(Queue& queue,
                                                    ::btldigi::BTLDigiSoA::ConstView const& input,
                                                    BTLBaseRecHitSoA::View& output,
                                                    const double npeToADC0_,
                                                    const double npeToADC1_,
                                                    const double npeSaturationCorr0_,
                                                    const double npeSaturationCorr1_,
                                                    const double npePerMeV_) {
    //,
    //Table const& tdc,
    //Table const& qdc) {
    // Use 64 items per group.
    // This value is arbitrary, but it's a reasonable starting point.
    uint32_t items = 64;

    // Use as many groups as needed to cover the whole problem.
    // If this value is too large, a smaller number of blocks can give better performance.
    uint32_t groups = cms::alpakatools::divide_up_by(input.metadata().size(), items);

    auto grid = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        BTLdigiToBaseKernel{},
                        input,
                        output,
                        npeToADC0_,
                        npeToADC1_,
                        npeSaturationCorr0_,
                        npeSaturationCorr1_,
                        npePerMeV_);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit
