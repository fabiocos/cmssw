#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/MTDCommonData/interface/BTLElectronicsMapping.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <ostream>

BTLElectronicsMapping::BTLElectronicsMapping() {}

// BTLElectronicsMapping::BTLElectronicsMapping(uint32_t id) { BTLElectronicsMapping_ = id; }

int BTLElectronicsMapping::SiPMCh(uint32_t zside, uint32_t smodCopy, uint32_t crystal, uint32_t SiPMSide) {

  if (1 > crystal || crystal > BTLDetId::kCrystalsPerModuleV2){
    edm::LogWarning("MTDGeom") << "BTLNumberingScheme::BTLElectronicsMapping(): "
                               << "****************** Bad crystal number = " << crystal;
    return 0;
  }

  if (1 > smodCopy || smodCopy > BTLDetId::kSModulesPerDM) {
    edm::LogWarning("MTDGeom") << "BTLNumberingScheme::getUnitID(): "
                                << "****************** Bad detector module copy = " << smodCopy;
    return 0;
  }

  if (1 < zside) {
    edm::LogWarning("MTDGeom") << "BTLNumberingScheme::getUnitID(): "
                                << "****************** Bad side = " << zside;
    return 0;
  }

  if (zside == 1) {
    if( smodCopy == 1 ) return BTLElectronicsMapping::SiPMChannelMapFWZpos[(crystal-1) + SiPMSide * BTLDetId::kCrystalsPerModuleV2];
    else return BTLElectronicsMapping::SiPMChannelMapBWZpos[(crystal-1) + SiPMSide * BTLDetId::kCrystalsPerModuleV2];
  } else {
    if( smodCopy == 1 ) return BTLElectronicsMapping::SiPMChannelMapFWZneg[(crystal-1) + SiPMSide * BTLDetId::kCrystalsPerModuleV2];
    else return BTLElectronicsMapping::SiPMChannelMapBWZneg[(crystal-1) + SiPMSide * BTLDetId::kCrystalsPerModuleV2];
  }
  
}

int BTLElectronicsMapping::SiPMCh(BTLDetId det, uint32_t SiPMSide) {
  uint32_t zside = det.mtdSide();
  uint32_t smodCopy = det.smodule();
  uint32_t crystal = det.crystal();

  return BTLElectronicsMapping::SiPMCh(zside, smodCopy, crystal, SiPMSide);
}

int BTLElectronicsMapping::SiPMCh(uint32_t rawId, uint32_t SiPMSide) {
  BTLDetId theId(rawId);
  return BTLElectronicsMapping::SiPMCh(theId, SiPMSide);
}


int BTLElectronicsMapping::TOFHIRCh(uint32_t zside, uint32_t smodCopy, uint32_t crystal, uint32_t SiPMSide) {
  int SiPMCh_ = BTLElectronicsMapping::SiPMCh(zside, smodCopy, crystal, SiPMSide);
  return BTLElectronicsMapping::THChannelMap[SiPMCh_];
}

int BTLElectronicsMapping::TOFHIRCh(BTLDetId det, uint32_t SiPMSide) {
  uint32_t zside = det.mtdSide();
  uint32_t smodCopy = det.smodule();
  uint32_t crystal = det.crystal();

  return BTLElectronicsMapping::TOFHIRCh(zside, smodCopy, crystal, SiPMSide);
}

int BTLElectronicsMapping::TOFHIRCh(uint32_t rawId, uint32_t SiPMSide) {
  BTLDetId theId(rawId);
  return BTLElectronicsMapping::TOFHIRCh(theId, SiPMSide);
}

int BTLElectronicsMapping::THChToXtal(uint32_t zside, uint32_t smodCopy, uint32_t THCh){

  if (1 > smodCopy || BTLDetId::kSModulesPerDM < smodCopy) {
    edm::LogWarning("MTDGeom") << "BTLNumberingScheme::getUnitID(): "
                                << "****************** Bad detector module copy = " << smodCopy;
    return 0;
  }

  if (1 < zside) {
    edm::LogWarning("MTDGeom") << "BTLNumberingScheme::getUnitID(): "
                                << "****************** Bad side = " << zside;
    return 0;
  }

  // std::array<uint32_t, BTLDetId::kCrystalsPerModuleV2 * 2> SiPMChMap;
  // if ( (zside + smodCopy) % 2 == 0 ) SiPMChMap = BTLElectronicsMapping::SiPMChannelMapFW;  
  // else SiPMChMap = BTLElectronicsMapping::SiPMChannelMapBW; 
  // return SiPMChMap.find(THCh) % 16;
  return 0;
}

