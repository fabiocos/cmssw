#ifndef DataFormats_BTLDetId_BTLDetId_h
#define DataFormats_BTLDetId_BTLDetId_h

#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include <iostream>
#include <ostream>
#include <array>

/** 
    @class BTLDetId
    @brief Detector identifier class for the Barrel Timing Layer.
    The crystal count must start from 0, copy number must be scaled by 1 unit.

    // Geometry v1, v2
    bit 15-10: module sequential number
    bit 9-8  : crystal type (1 - 3)
    bit 7-6  : readout unit sequential number within a type ( 1 - 2 )
    bit 5-0  : crystal sequential number within a module ( 0 - 15 )

    // Geometry v3 (all type 1 modules)
    bit 12-10: Readout unit number ( 1 - 6 )
    bit 9-6  : Detector Module ( 1 - 12 )
    bit  5   : Sensor Module inside DM ( 0 - 1 )
    bit 4-0  : Crystal number in a SM ( 1 - 16 )

*/

class BTLDetId : public MTDDetId {
public:
  
  // Commmon to all Geometry Versions
  static constexpr uint32_t kBTLCrystalOffset = 0;
  static constexpr uint32_t kBTLCrystalMask = 0x1F;

  // v1, v2
  static constexpr uint32_t kBTLmoduleOffset = 10;
  static constexpr uint32_t kBTLmoduleMask = 0x3F;
  static constexpr uint32_t kBTLmodTypeOffset = 8;
  static constexpr uint32_t kBTLmodTypeMask = 0x3;
  static constexpr uint32_t kBTLRUOffset = 6;
  static constexpr uint32_t kBTLRUMask = 0x3;

  // v3 (current geometry)
  static constexpr uint32_t kBTLRUallOffset = 10;
  static constexpr uint32_t kBTLRUallMask = 0x3F;
  static constexpr uint32_t kBTLdetectorModOffset = 6;
  static constexpr uint32_t kBTLdetectorModMask = 0xF;
  static constexpr uint32_t kBTLsensorModOffset = 5;
  static constexpr uint32_t kBTLsensorModMask = 0x1;
  uint32_t modulesInRURow = 3;
  uint32_t modulesInRULine = 4;
  uint32_t modulesInDM = 2;

  // Geometry old to new form
  static constexpr uint32_t kBTLoldFieldMask = 0xFFFF; 
  static constexpr uint32_t kBTLnewFormat = 0;
  
  /// range constants, need two sets for the time being (one for tiles and one for bars)
  static constexpr uint32_t HALF_ROD = 36;
  static constexpr uint32_t kModulesPerRODBarPhiFlat = 48;
  static constexpr uint32_t kModulePerTypeBarPhiFlat = 48 / 3;
  static constexpr uint32_t kRUPerTypeV2 = 2;
  static constexpr uint32_t kModulesPerRUV2 = 24;
  static constexpr uint32_t kCrystalsPerModuleV2 = 16;
  static constexpr uint32_t kModulesPerTrkV2 = 3;
  static constexpr uint32_t kCrystalTypes = 3;

  // Number of crystals in BTL according to TDR design, valid also for barphiflat scenario:
  // 16 crystals x 24 modules x 2 readout units/type x 3 types x 36 rods/side x 2 sides
  //
  static constexpr uint32_t kCrystalsBTL =
      kCrystalsPerModuleV2 * kModulesPerRUV2 * kRUPerTypeV2 * kCrystalTypes * HALF_ROD * 2;

  enum class CrysLayout { tile = 1, bar = 2, barzflat = 3, barphiflat = 4, v2 = 5, v3 = 6 };

  // ---------- Constructors, enumerated types ----------

  /** Construct a null id */
  BTLDetId() : MTDDetId(DetId::Forward, ForwardSubdetector::FastTime) {
    id_ |= (MTDType::BTL & kMTDsubdMask) << kMTDsubdOffset;
  }

  // /** Construct from a raw value */
  // BTLDetId(const uint32_t& raw_id) { id_ = MTDDetId(tmpId).rawId(); }

  // /** Construct from generic DetId */
  // BTLDetId(const DetId& det_id) { id_ = MTDDetId(tmpId).rawId(); }

  /** Construct from a raw value */
  BTLDetId(const uint32_t& raw_id) : MTDDetId(raw_id) { ; }

  /** Construct from generic DetId */
  BTLDetId(const DetId& det_id) : MTDDetId(det_id.rawId()) { ; }

  /** Construct from complete geometry information, v1 **/
  BTLDetId(uint32_t zside, uint32_t rod, uint32_t module, uint32_t modtyp, uint32_t crystal)
      : MTDDetId(DetId::Forward, ForwardSubdetector::FastTime) {
    id_ |= (MTDType::BTL & kMTDsubdMask) << kMTDsubdOffset | (zside & kZsideMask) << kZsideOffset |
           (rod & kRodRingMask) << kRodRingOffset | (module & kBTLmoduleMask) << kBTLmoduleOffset |
           (modtyp & kBTLmodTypeMask) << kBTLmodTypeOffset | ((crystal - 1) & kBTLCrystalMask) << kBTLCrystalOffset;
  }

  /** Construct from complete geometry information **/
  BTLDetId(uint32_t zside, uint32_t rod, uint32_t runit, uint32_t dmodule, uint32_t smodule, uint32_t crystal)
      : MTDDetId(DetId::Forward, ForwardSubdetector::FastTime) {
    
    // new BTLDetID format (v3)
    if (kBTLnewFormat == 1){
    id_ |= (MTDType::BTL & kMTDsubdMask) << kMTDsubdOffset | (zside & kZsideMask) << kZsideOffset |
           (rod & kRodRingMask) << kRodRingOffset | (runit & kBTLRUallMask) << kBTLRUallOffset | 
           (dmodule & kBTLdetectorModMask) << kBTLdetectorModOffset |
           (smodule & kBTLsensorModMask) << kBTLsensorModOffset  |
           ((crystal - 1) & kBTLCrystalMask) << kBTLCrystalOffset;
    }

    // old BTLDetID format (v1, v2)
    else if (kBTLnewFormat == 0){
      uint32_t module = dmodule;
      uint32_t modtyp = smodule;
      id_ |= (MTDType::BTL & kMTDsubdMask) << kMTDsubdOffset | (zside & kZsideMask) << kZsideOffset |
            (rod & kRodRingMask) << kRodRingOffset | (module & kBTLmoduleMask) << kBTLmoduleOffset |
            (modtyp & kBTLmodTypeMask) << kBTLmodTypeOffset | (runit & kBTLRUMask) << kBTLRUOffset |
            ((crystal - 1) & kBTLCrystalMask) << kBTLCrystalOffset;
    }

  }

  // ---------- Common methods ----------

  /** Returns BTL sensor module number. */
  inline int smodule() const { return (id_ >> kBTLsensorModOffset) & kBTLsensorModMask; }

  /** Returns BTL detector module number. */
  inline int dmodule() const { return (id_ >> kBTLdetectorModOffset) & kBTLdetectorModMask; }

  /** Returns BTL module number from 1 to 24 **/
  inline int module() const { 
    if (kBTLnewFormat==0){
      return (id_ >> kBTLmoduleOffset) & kBTLmoduleMask; 
    } else if (kBTLnewFormat==1) {
      return (int(dmodule()) * modulesInDM + int(smodule()) - 1);
    }
    return 0;
  }


  /** Returns BTL crystal type number. */
  inline int modType() const { 
    if (kBTLnewFormat==0){
      return (id_ >> kBTLmodTypeOffset) & kBTLmodTypeMask; 
    } else if (kBTLnewFormat==1) {
      return 1;
    }
    return 0;
  }

  /** Returns BTL crystal number. */
  inline int crystal() const { return ((id_ >> kBTLCrystalOffset) & kBTLCrystalMask) + 1; }

  /** Returns BTL readout unit number per type. */
  inline int runit() const { 
    if (kBTLnewFormat==0){  
      return (id_ >> kBTLRUOffset) & kBTLRUMask; 
    } else if (kBTLnewFormat == 1){
      return (id_ >> kBTLRUallOffset) & kBTLRUallMask;
    } 
    return 0;
  }

  /** Returns BTL global readout unit number. */
  inline int globalRunit() const {
    if (runit() == 0) {
      // pre-V2: build a RU identifier from available information
      return (module() - 1) / kModulePerTypeBarPhiFlat / kRUPerTypeV2 + 1;
    } else if (runit() > 0 && modType() > 0) {
      // V2/V3: build global RU identifier from RU per type and type
      return (modType() - 1) * kRUPerTypeV2 + runit();
    }
    return 0;
  }

  /** return the row in GeomDet language **/
  inline int row(unsigned nrows = kCrystalsPerModuleV2) const {
    return (crystal() - 1) % nrows;  // anything else for now
  }

  /** return the column in GeomDetLanguage **/
  inline int column(unsigned nrows = kCrystalsPerModuleV2) const { return (crystal() - 1) / nrows; }

  /** create a Geographical DetId for Tracking **/
  BTLDetId geographicalId(CrysLayout lay) const;

  /** conversion from v1,v2 (old BTLDetID) to v3 geometry **/
  uint32_t newForm(const uint32_t& rawid) {


    uint32_t fixedP = rawid & (0xFFFFFFFF - kBTLoldFieldMask);          // unchanged part of id
    
    // convert old module number into detector module + sensor module numbers
    uint32_t oldModule = (rawid >> kBTLmoduleOffset) & kBTLmoduleMask;  
    uint32_t detModule = int( (oldModule-1)/(modulesInRURow*modulesInDM) ) + 1 + modulesInRULine * ( int( (oldModule%6)/modulesInDM ) );
    uint32_t senModule = (oldModule-1) % modulesInDM;
    // uint32_t detModule = ((int((oldModule - 1)/modulesInRURow)/modulesInDM)*3) + ((oldModule-1)%3 + 1);
    // uint32_t senModule = int((oldModule - 1)/modulesInRURow)%2;

    // convert old RU and type number into new RU number
    uint32_t oldRU   = (rawid >> kBTLRUOffset) & kBTLRUMask;  
    uint32_t oldType = (rawid >> kBTLmodTypeOffset) & kBTLmodTypeMask;
    uint32_t newRU   = ((oldType-1) >> 1) + oldRU;
    
    
    // get crystal number
    uint32_t crystal = (rawid & kBTLCrystalMask) >> kBTLCrystalOffset; 

    // std::cout << " oldRU: " << oldRU << std::endl;
    // std::cout << " oldType: " << oldType << std::endl;
    // std::cout << " newRU " << oldRU << std::endl;
    // std::cout << " Raw ID: " << rawid << std::endl;
    // std::cout << " fixedP: " << fixedP << std::endl;
    // std::cout << " Crystal type: " << oldType << std::endl;
    // std::cout << " Readout unit: " << oldRU << std::endl;
    // std::cout << " Global RU   : " << newRU << std::endl;
    // std::cout << " Module (v2 geom) : " << oldModule << std::endl;
    // std::cout << " detModule      : " << detModule << std::endl;
    // std::cout << " senModule      : " << senModule << std::endl;
    // std::cout << " Crystal     : " << crystal << std::endl;

    // return new BTLDetID for v3 geom
    return (fixedP | (newRU & kBTLRUallMask) << kBTLRUallOffset | 
           (detModule & kBTLdetectorModMask) << kBTLdetectorModOffset |
           (senModule & kBTLsensorModMask) << kBTLsensorModOffset  |
           ((crystal & kBTLCrystalMask) << kBTLCrystalOffset));
  }
};

std::ostream& operator<<(std::ostream&, const BTLDetId&);

#endif  // DataFormats_BTLDetId_BTLDetId_h
