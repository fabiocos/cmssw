#ifndef __RecoLocalFastTime_FTLCommonAlgos_RecHitTools_h__
#define __RecoLocalFastTime_FTLCommonAlgos_RecHitTools_h__

#include <array>
#include <cmath>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"

#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

class DetId;
class MTDGeometry;

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace mtd {
  class RecHitTools {
  public:
    RecHitTools()
        : geom_(nullptr),
	  topology_(nullptr) {}
    ~RecHitTools() {}

    void setGeometry(MTDGeometry const* geom);
    void setTopology(MTDTopology const* topo);
    //const CaloSubdetectorGeometry* getSubdetectorGeometry(const DetId& id) const;
    bool isETL(const DetId&) const;
    bool isBTL(const DetId&) const;

    GlobalPoint getPosition(const DetId& id, int row=0, int column=0) const;
    GlobalPoint getPosition(const DetId& id, const LocalPoint& local_point) const;
    //GlobalPoint getPositionLayer(int layer, bool nose = false) const;
    // zside returns +/- 1
    int zside(const DetId& id) const;

    //unsigned int getLayer(DetId::Detector type, bool nose = false) const;
    //unsigned int getLayer(ForwardSubdetector type) const;
    unsigned int getLayer(const DetId&) const;
    int getModule(const DetId&) const;
    std::pair<float, float> getPixelInModule(const DetId& id, const int row, const int column) const;
    std::pair<uint8_t, uint8_t> getPixelInModule(const DetId& id, const LocalPoint& local_point) const;
    int getCrystalInModule(const DetId&) const;

    // 4-vector helper functions using GlobalPoint
    float getEta(const GlobalPoint& position, const float& vertex_z = 0.) const;
    float getPhi(const GlobalPoint& position) const;
    float getPt(const GlobalPoint& position, const float& hitEnergy, const float& vertex_z = 0.) const;

//TODO
    // 4-vector helper functions using DetId
//    float getEta(const DetId& id, const float& vertex_z = 0.) const;
//    float getPhi(const DetId& id) const;
//    float getPt(const DetId& id, const float& hitEnergy, const float& vertex_z = 0.) const;

    inline const MTDGeometry* getGeometry() const { return geom_; };
    inline const MTDTopology* getTopology() const { return topology_; };
//    unsigned int lastLayerEE(bool nose = false) const { return (nose ? HFNoseDetId::HFNoseLayerEEmax : fhOffset_); }
//    unsigned int lastLayerFH() const { return fhLastLayer_; }
//    unsigned int firstLayerBH() const { return bhFirstLayer_; }
//    unsigned int lastLayerBH() const { return bhLastLayer_; }
//    unsigned int lastLayer(bool nose = false) const { return (nose ? noseLastLayer_ : bhLastLayer_); }
//    std::pair<uint32_t, uint32_t> firstAndLastLayer(DetId::Detector det, int subdet) const;
//    unsigned int maxNumberOfWafersPerLayer(bool nose = false) const {
//      return (nose ? maxNumberOfWafersNose_ : maxNumberOfWafersPerLayer_);
//    }
//    inline int getScintMaxIphi() const { return bhMaxIphi_; }
//    inline int getGeometryType() const { return geometryType_; }
//    bool maskCell(const DetId& id, int corners = 3) const;

  private:
    const MTDGeometry* geom_;
    const MTDTopology* topology_;
  };
}  // namespace mtd

#endif
