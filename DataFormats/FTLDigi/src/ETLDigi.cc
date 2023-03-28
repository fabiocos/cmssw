#include "DataFormats/FTLDigi/interface/ETLDigi.h"
#include "FWCore/Utilities/interface/Exception.h"

//bool ETLDigi::fillWithSampleArray(const std::vector<ETLSample> in, ETLDigi *out) {
   //if (in.size() != 5) {
     //throw cms::Exception("IncorrectInput") << "Input vector of ETLSample size is " << in.size() << " instead of 5, aborting...";
   //}

   //out->clear();

   //// incorrect code, just to fill with "something" for test

   //out->set(static_cast<uint8_t>(0), static_cast<uint32_t>(in[2].raw()));
//}

#include <iomanip>

std::ostream& operator<<(std::ostream& os, const ETLDigi& digi) {
   os << " row " << digi.row() << " col " << digi.column() << " load " << digi.data() << std::endl;
   return os;
}

