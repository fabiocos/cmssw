#include <ostream>

#include "DataFormats/FTLDigiSoA/interface/ETLDigiSoA.h"

namespace etldigi {

  std::ostream& operator<<(std::ostream& out, ETLDigiSoA::View::const_element const& digi) {
    out << "ETL Digi SoA rawId : " << digi.rawId() 
        << ", header = " << digi.header() << ", status = " << digi.status()
        << ", column = " << digi.colID() << ", row = " << digi.rowID()
        << ", ToA = " << digi.ToAdata() << ", ToT = " << digi.ToTdata() << ", CAL = " << digi.CALdata();
    return out;
  }
}  // namespace etldigi
