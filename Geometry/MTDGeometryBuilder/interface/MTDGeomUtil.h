#ifndef __RecoLocalFastTime_FTLCommonAlgos_MTDGeomUtil_h__
#define __RecoLocalFastTime_FTLCommonAlgos_MTDGeomUtil_h__

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"

class DetId;
class MTDGeometry;

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace mtd {
  class MTDGeomUtil {
  public:
    MTDGeomUtil()
        : geom_(nullptr),
	  topology_(nullptr) {}
    ~MTDGeomUtil() {}

    void setGeometry(MTDGeometry const* geom);
    void setTopology(MTDTopology const* topo);

    bool isETL(const DetId&) const;
    bool isBTL(const DetId&) const;

    GlobalPoint getPosition(const DetId& id, int row=0, int column=0) const;
    GlobalPoint getPosition(const DetId& id, const LocalPoint& local_point) const;

    // zside returns +/- 1
    int zside(const DetId& id) const;

    unsigned int getLayer(const DetId&) const;
    int getModule(const DetId&) const;
    std::pair<float, float> getPixelInModule(const DetId& id, const int row, const int column) const;
    std::pair<uint8_t, uint8_t> getPixelInModule(const DetId& id, const LocalPoint& local_point) const;
    int getCrystalInModule(const DetId&) const;

    // 4-vector helper functions using GlobalPoint
    float getEta(const GlobalPoint& position, const float& vertex_z = 0.) const;
    float getPhi(const GlobalPoint& position) const;
    float getPt(const GlobalPoint& position, const float& hitEnergy, const float& vertex_z = 0.) const;

    // 4-vector helper functions using DetId
    float getEta(const DetId& id, const LocalPoint& local_point, const float& vertex_z = 0.) const;
    float getPhi(const DetId& id, const LocalPoint& local_point) const;
    float getPt(const DetId& id, const LocalPoint& local_point, const float& hitEnergy, const float& vertex_z = 0.) const;

    inline const MTDGeometry* getGeometry() const { return geom_; };
    inline const MTDTopology* getTopology() const { return topology_; };

  private:
    const MTDGeometry* geom_;
    const MTDTopology* topology_;
  };
}  // namespace mtd

#endif
