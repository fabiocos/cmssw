#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include <iomanip>

std::ostream& operator<<(std::ostream& os, const ETLDetId& id) {
  os << (MTDDetId&)id;
  os << " ETL " << std::endl
     << " Side                  : " << id.mtdSide() << std::endl
     << " Disc                  : " << id.nDisc() << std::endl
     << " Disc Side             : " << id.discSide() << std::endl
     << " Sector                : " << id.sector() << std::endl
     << " Service Hybrid Type   : " << id.servType() << std::endl
     << " Service Hybrid Number : " << id.servCopy() << std::endl
     << " Module Number         : " << id.module() << std::endl
     << " Module Type           : " << id.modType() << std::endl
     << " Sensor                : " << id.sensor() << std::endl;
  return os;
}
