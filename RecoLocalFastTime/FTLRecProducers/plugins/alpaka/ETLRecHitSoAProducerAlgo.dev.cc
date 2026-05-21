#include <alpaka/alpaka.hpp>

#include "CommonTools/Utils/interface/FormulaEvaluator.h"

#include "RecoLocalFastTime/FTLCommonAlgos/interface/MTDTimeCalib.h"
#include "DataFormats/FTLRecHitSoA/interface/ETLBaseRecHitSoA.h"
#include "DataFormats/FTLRecHitSoA/interface/ETLRecHitSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "RecoLocalFastTime/Records/interface/MTDTimeCalibRecord.h"
#include "RecoLocalFastTime/FTLCommonAlgos/interface/MTDTimeCalib.h"

#include "ETLRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::etlrechit {

  using namespace ::etlrechit;
  ALPAKA_FN_ACC float timeResolutionInNs(float amp) { return 0.0370; }

  class ETLBaseToRecoKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  ETLBaseRecHitSoA::ConstView input,
                                  ETLRecHitSoA::View output,
                                  const double thresholdToKeep_,
                                  const double calibration_) const {
      // make a strided loop over the kernel grid, covering up to "size" elements

      for (int32_t i : cms::alpakatools::uniform_elements(acc, input.metadata().size())) {
        auto entry = input[i];
        float time1 = 0;
        float time2 = 0;
        float position = -1.;        // dummy
        float position_error = -1.;  // dummy
        float time_error = 0;
        float energy = -1;  // dummy
        uint8_t flag = 0;   // assumed to be ok

        //!!!!!!! position error calculation to be added

        // time set
        time1 = entry.toa();
        time2 = entry.tot();

        // --- Energy calibration
        energy = time2;  //for ETL, it is the time_over_threshold
        energy *= calibration_;

        time_error = timeResolutionInNs(energy);

        if (energy > thresholdToKeep_) {
          flag = 1;
        } else {
          flag = 0;
        }

#ifdef EDM_ML_DEBUG

        printf("RecHit SoA with raw id %i \n", entry.detId().rawId());
        printf("Time 1: %f +- %f \n", time1, time_error);
        printf("Time 2: %f \n", time2);
        printf("Energy %f \n", energy);
        printf("Position: %f +- %f \n", position, position_error);

#endif

        // fill the rechit
        output[i] = {
            entry.detId(), entry.row(), entry.column(), time1, time2, energy, position, time_error, position_error, flag};
      }
    }
  };

  void ETLRecHitSoAProducerAlgo::fromBaseToReco(Queue& queue,
                                                ETLBaseRecHitSoA::ConstView const& input,
                                                ETLRecHitSoA::View& output,
                                                const double thresholdToKeep_,
                                                const double calibration_) {
    // Use 64 items per group.
    // This value is arbitrary, but it's a reasonable starting point.
    uint32_t items = 64;

    // Use as many groups as needed to cover the whole problem.
    // If this value is too large, a smaller number of blocks can give better performance.
    uint32_t groups = cms::alpakatools::divide_up_by(input.metadata().size(), items);

    auto grid = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, grid, ETLBaseToRecoKernel{}, input, output, thresholdToKeep_, calibration_);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::etlrechit
