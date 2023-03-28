#include "DataFormats/FTLDigi/interface/BTLDigi.h"
#include "FWCore/Utilities/interface/Exception.h"

//bool BTLDigi::fillWithSampleArray(const std::vector<BTLSample> in, BTLDigi *out) {
   //if (in.size() != 2) {
     //throw cms::Exception("IncorrectInput") << "Input vector of BTLSample size is " << in.size() << " instead of 2, aborting...";
   //}

   //out->clear();

   //// incorrect code, just to fill with "something" for test

   //out->set(static_cast<uint32_t>(in[0].raw_flag()), static_cast<uint64_t>(in[0].raw_data()), static_cast<uint32_t>(in[1].raw_flag()), static_cast<uint64_t>(in[1].raw_data()));
//}

#include <iomanip>

std::ostream& operator<<(std::ostream& os, const BTLDigi& digi) {
   os << " Lhs SiPM: " << digi.plrhs1() << " " << digi.plrhs2() << " Rhs SiPM: " << digi.pllhs1() << " " << digi.pllhs2() << std::endl;
  return os;
}

