#ifndef DataFormats_FTLRecHitSoA_interface_BTLBaseRecHitSoA_h
#define DataFormats_FTLRecHitSoA_interface_BTLBaseRecHitSoA_h

#include <ostream>

#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

namespace btlrechit {
  GENERATE_SOA_LAYOUT(BTLBaseRecHitSoALayout,
                      SOA_COLUMN(DetId, detId),
                      SOA_COLUMN(uint8_t, row),
                      SOA_COLUMN(float, time1R),
                      SOA_COLUMN(float, time2R),
                      SOA_COLUMN(float, ampR),
                      SOA_COLUMN(uint16_t, idleTimeR),
                      SOA_COLUMN(uint8_t, flagsR),
                      SOA_COLUMN(float, time1L),
                      SOA_COLUMN(float, time2L),
                      SOA_COLUMN(float, ampL),
                      SOA_COLUMN(uint16_t, idleTimeL),
                      SOA_COLUMN(uint8_t, flagsL))

  using BTLBaseRecHitSoA = BTLBaseRecHitSoALayout<>;
  using BTLBaseRecHitSoAView = BTLBaseRecHitSoA::View;
  using BTLBaseRecHitSoAConstView = BTLBaseRecHitSoA::ConstView;

  std::ostream& operator<<(std::ostream& out, BTLBaseRecHitSoA::View::const_element const& btlrh);
}  // namespace btlrechit
#endif  // DataFormats_FTLRecHitSoA_interface_BTLBaseRecHitSoA_h
