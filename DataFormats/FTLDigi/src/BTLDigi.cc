#include "DataFormats/FTLDigi/interface/BTLDigi.h"

namespace btldigi {

  std::ostream& operator<<(std::ostream& out, const BTLDigi& digi) {
    out << "BTL Digi rawId : " << digi.krawId() << ", BC0count = " << digi.kBC0count()
        << ", status = " << digi.kstatus() << ", BCcount = " << digi.kBCcount() << std::endl
        << "\t sample 0 (left side) : chIDL = " << static_cast<int>(digi.kchIDL())
        << ", T1coarseL = " << digi.kT1coarseL() << ", T1fineL = " << digi.kT1fineL()
        << ", T2coarseL = " << digi.kT2coarseL() << ", T2fineL = " << digi.kT2fineL()
        << ", EOIcoarseL = " << digi.kEOIcoarseL() << ", ChargeL = " << digi.kChargeL() << std::endl
        << "\t sample 1 (right side) : chIDR = " << static_cast<int>(digi.kchIDR())
        << ", T1coarseR = " << digi.kT1coarseR() << ", T1fineR = " << digi.kT1fineR()
        << ", T2coarseR = " << digi.kT2coarseR() << ", T2fineR = " << digi.kT2fineR()
        << ", EOIcoarseR = " << digi.kEOIcoarseR() << ", ChargeR = " << digi.kChargeR();
    return out;
  }
}  // namespace btldigi
