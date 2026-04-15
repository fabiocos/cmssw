#include <ostream>
#include <bitset>
#include <fmt/format.h>

#include "DataFormats/FTLRecHitSoA/interface/BTLBaseRecHitSoA.h"

namespace btlrechit {

  std::ostream& operator<<(std::ostream& out, BTLBaseRecHitSoA::View::const_element const& btlrh) {
    out << "BTL uncalib rechit SoA: "
        << " detID: " << btlrh.detId().rawId() << ", row: " << static_cast<int>(btlrh.row())
        << ", time1 R: " << btlrh.time1R() << ", time2 R: " << btlrh.time2R() << ", amplitude R: " << btlrh.ampR()
        << ", idleTime R: " << btlrh.idleTimeR() << ", flags R: " << std::bitset<2>(btlrh.flagsL())
        << ", time1 L: " << btlrh.time1L() << ", time2 L: " << btlrh.time1L() << ", amplitude L: " << btlrh.ampL()
        << ", idleTime L: " << btlrh.idleTimeL() << ", flags L: " << std::bitset<2>(btlrh.flagsL());

    return out;
  }

}  // namespace btlrechit
