#include "Geometry/MTDGeometryBuilder/interface/MTDGeomUtil.h"

#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

using namespace mtd;

void MTDGeomUtil::setGeometry(const MTDGeometry* geom) { geom_ = geom; }

void MTDGeomUtil::setTopology(const MTDTopology* topo) { topology_ = topo; }

bool MTDGeomUtil::isETL(const DetId& id) const {
  MTDDetId hid{id};
  const auto& subDet = hid.mtdSubDetector();
  if (subDet == 0)
    throw cms::Exception("mtd::MTDGeomUtil") << "DetId " << hid.rawId() << " not in MTD!" << std::endl;
  if (subDet == MTDDetId::MTDType::ETL)
    return true;
  return false;
}

bool MTDGeomUtil::isBTL(const DetId& id) const { return !(isETL(id)); }

// row and column set to 0 by default since they are not needed for BTL
GlobalPoint MTDGeomUtil::getPosition(const DetId& id, int row, int column) const {
  GlobalPoint global_point(0., 0., 0.);
  if (isBTL(id)) {
    BTLDetId detId{id};
    DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(getTopology()->getMTDTopologyMode()));
    const MTDGeomDet* thedet = geom_->idToDet(geoId);
    if (thedet == nullptr)
      throw cms::Exception("mtd::MTDGeomUtil") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                               << detId.rawId() << ") is invalid!" << std::dec << std::endl;
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

    Local3DPoint local_point(0., 0., 0.);
    local_point = topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
    global_point = thedet->toGlobal(local_point);
  } else if (isETL(id)) {
    ETLDetId detId{id};
    DetId geoId = detId.geographicalId();
    const MTDGeomDet* thedet = geom_->idToDet(geoId);
    if (thedet == nullptr)
      throw cms::Exception("mtd::MTDGeomUtil") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                               << detId.rawId() << ") is invalid!" << std::dec << std::endl;
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

    Local3DPoint local_point(topo.localX(row), topo.localY(column), 0.);
    global_point = thedet->toGlobal(local_point);
  } else {
    throw cms::Exception("mtd::MTDGeomUtil") << "detId " << id.rawId() << " not in MTD" << std::endl;
  }
  return global_point;
}

GlobalPoint MTDGeomUtil::getPosition(const DetId& id, const LocalPoint& local_point) const {
  auto global_point = GlobalPoint(0., 0., 0.);
  if (isBTL(id)) {
    BTLDetId detId{id};
    DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(getTopology()->getMTDTopologyMode()));
    const MTDGeomDet* thedet = geom_->idToDet(geoId);
    if (thedet == nullptr)
      throw cms::Exception("mtd::MTDGeomUtil") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                               << detId.rawId() << ") is invalid!" << std::dec << std::endl;
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
    auto local_point_sim =
        topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
    global_point = thedet->toGlobal(local_point_sim);
  } else if (isETL(id)) {
    ETLDetId detId{id};
    DetId geoId = detId.geographicalId();
    const MTDGeomDet* thedet = geom_->idToDet(geoId);
    if (thedet == nullptr)
      throw cms::Exception("mtd::MTDGeomUtil") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                               << detId.rawId() << ") is invalid!" << std::dec << std::endl;
    global_point = thedet->toGlobal(local_point);
  } else {
    throw cms::Exception("mtd::MTDGeomUtil") << "detId " << id.rawId() << " not in MTD" << std::endl;
  }
  return global_point;
}

int MTDGeomUtil::zside(const DetId& id) const {
  const MTDDetId hid(id);
  return hid.zside();
}

unsigned int MTDGeomUtil::getLayer(const DetId& id) const {
  unsigned int layer(0);
  if (isETL(id)) {
    ETLDetId hid(id);
    layer = hid.nDisc();
  } else {
    layer = 0;
  }
  return layer;
}

int MTDGeomUtil::getModule(const DetId& id) const {
  int module = -1;
  if (isETL(id)) {
    ETLDetId hid(id);
    module = hid.module();
  } else {
    BTLDetId hid(id);
    module = hid.module();
  }
  return module;
}

// returns the local position as a pair (x, y) - for ETL
std::pair<float, float> MTDGeomUtil::getPixelInModule(const DetId& id, const int row, const int column) const {
  ETLDetId detId(id);
  DetId geoId = detId.geographicalId();
  const MTDGeomDet* thedet = getGeometry()->idToDet(geoId);
  if (thedet == nullptr)
    throw cms::Exception("mtd::MTDGeomUtil") << "GeographicalID: " << std::hex << geoId.rawId() << " (" << detId.rawId()
                                             << ") is invalid!" << std::dec << std::endl;
  const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
  const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
  const Local3DPoint local_point(topo.localX(row), topo.localY(column), 0.);
  return topo.pixel(local_point);
}

// returns row and column as a pair (row, col)
std::pair<uint8_t, uint8_t> MTDGeomUtil::getPixelInModule(const DetId& id, const LocalPoint& local_point) const {
  if (isETL(id)) {
    ETLDetId detId(id);
    DetId geoId = detId.geographicalId();
    const MTDGeomDet* thedet = getGeometry()->idToDet(geoId);
    if (thedet == nullptr)
      throw cms::Exception("mtd::MTDGeomUtil") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                               << detId.rawId() << ") is invalid!" << std::dec << std::endl;
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
    const auto& thepixel = topo.pixel(local_point);
    uint8_t row(thepixel.first), col(thepixel.second);
    return std::pair<uint8_t, uint8_t>(row, col);
  } else {
    BTLDetId detId(id);
    DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(getTopology()->getMTDTopologyMode()));
    const MTDGeomDet* thedet = getGeometry()->idToDet(geoId);
    if (thedet == nullptr)
      throw cms::Exception("mtd::MTDGeomUtil") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                               << detId.rawId() << ") is invalid!" << std::dec << std::endl;
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
    auto topo_point = topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
    const auto& thepixel = topo.pixel(topo_point);
    uint8_t row(thepixel.first), col(thepixel.second);
    return std::pair<uint8_t, uint8_t>(row, col);
  }
}

int MTDGeomUtil::getCrystalInModule(const DetId& id) const {
  BTLDetId hid(id);
  return hid.crystal();
}

float MTDGeomUtil::getEta(const GlobalPoint& position, const float& vertex_z) const {
  GlobalPoint corrected_position = GlobalPoint(position.x(), position.y(), position.z() - vertex_z);
  return corrected_position.eta();
}

float MTDGeomUtil::getEta(const DetId& id, const LocalPoint& local_point, const float& vertex_z) const {
  GlobalPoint position = getPosition(id, local_point);
  float eta = getEta(position, vertex_z);
  return eta;
}

float MTDGeomUtil::getPhi(const GlobalPoint& position) const {
  float phi = atan2(position.y(), position.x());
  return phi;
}

float MTDGeomUtil::getPhi(const DetId& id, const LocalPoint& local_point) const {
  GlobalPoint position = getPosition(id, local_point);
  float phi = atan2(position.y(), position.x());
  return phi;
}

float MTDGeomUtil::getPt(const GlobalPoint& position, const float& hitEnergy, const float& vertex_z) const {
  float eta = getEta(position, vertex_z);
  float pt = hitEnergy / cosh(eta);
  return pt;
}

float MTDGeomUtil::getPt(const DetId& id,
                         const LocalPoint& local_point,
                         const float& hitEnergy,
                         const float& vertex_z) const {
  GlobalPoint position = getPosition(id, local_point);
  float eta = getEta(position, vertex_z);
  float pt = hitEnergy / cosh(eta);
  return pt;
}
