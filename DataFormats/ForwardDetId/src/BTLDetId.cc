#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include <iomanip>

std::ostream& operator<< ( std::ostream& os, const BTLDetId& id ) {
  return os << " BTL " << std::endl
            << " Side        : " << id.mtdSide() << std::endl
            << " Rod         : " << id.mtdRR() << std::endl
            << " Module      : " << id.btlModule() << std::endl
            << " Crystal type: " << id.btlCrysType() << std::endl
            << " Crystal     : " << id.btlCrystal() << std::endl;
}
