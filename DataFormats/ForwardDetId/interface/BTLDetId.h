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

    bit 15-10: module sequential number
    bit 9-8  : crystal type (1 - 3)
    bit 7-6  : readout unit sequential number within a type ( 1 - 2 )
    bit 5-0  : crystal sequential number within a module ( 0 - 15 )

    // Geometry v3 new DetID (all type 1 modules)
    bit 15: kBTLNewFormat (0 - old BTLDetID, 1 - new BTLDetID)
    bit 12-10: Readout unit number ( 1 - 6 )
    bit 9-6  : Detector Module ( 1 - 12 )
    bit  5   : Sensor Module inside DM ( 0 - 1 )
    bit 4-0  : Crystal number in a SM ( 1 - 16 )
*/

class BTLDetId : public MTDDetId {
public:
  // old BTLDetID RU and module number scheme
  static constexpr uint32_t kBTLoldModuleOffset = 10;
  static constexpr uint32_t kBTLoldModuleMask = 0x3F;
  static constexpr uint32_t kBTLoldModTypeOffset = 8;
  static constexpr uint32_t kBTLoldModTypeMask = 0x3;
  static constexpr uint32_t kBTLoldRUOffset = 6;
  static constexpr uint32_t kBTLoldRUMask = 0x3;

  // New BTLDetID
  static constexpr uint32_t kBTLRUOffset = 10;
  static constexpr uint32_t kBTLRUMask = 0x7;
  static constexpr uint32_t kBTLdetectorModOffset = 6;
  static constexpr uint32_t kBTLdetectorModMask = 0xF;
  static constexpr uint32_t kBTLsensorModOffset = 5;
  static constexpr uint32_t kBTLsensorModMask = 0x1;
  static constexpr uint32_t kBTLCrystalOffset = 0;
  static constexpr uint32_t kBTLCrystalMask = 0x1F;

  /// range constants, need two sets for the time being (one for tiles and one for bars)
  static constexpr uint32_t HALF_ROD = 36;
  static constexpr uint32_t kModulesPerRODBarPhiFlat = 48;
  static constexpr uint32_t kModulePerTypeBarPhiFlat = 48 / 3;
  static constexpr uint32_t kRUPerTypeV2 = 2;
  static constexpr uint32_t kRUPerRod = 6;
  static constexpr uint32_t kModulesPerRUV2 = 24;
  static constexpr uint32_t kDModulesPerRU = 12;
  static constexpr uint32_t kSModulesPerDM = 2;
  static constexpr uint32_t kDModulesInRUCol = 3;
  static constexpr uint32_t kDModulesInRURow = 4;
  static constexpr uint32_t kSModulesInDM = 2;
  static constexpr uint32_t kCrystalsPerModuleV2 = 16;
  static constexpr uint32_t kModulesPerTrkV2 = 3;
  static constexpr uint32_t kCrystalTypes = 3;

  // conversion
  static constexpr uint32_t kBTLoldFieldMask = 0xFFFF;
  static constexpr uint32_t kBTLNewFormat = 1 << 15;

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
    id_ |= kBTLNewFormat;
  }

  /** Construct from a raw value */
  BTLDetId(const uint32_t& raw_id) : MTDDetId(raw_id) {
    uint32_t tmpId = raw_id;
    if ((tmpId & kBTLNewFormat) == 0) {
      tmpId = newForm(tmpId);
    }
    id_ = MTDDetId(tmpId).rawId();
  }

  /** Construct from generic DetId */
  BTLDetId(const DetId& det_id) : MTDDetId(det_id.rawId()) {
    uint32_t tmpId = det_id.rawId();
    if ((tmpId & kBTLNewFormat) == 0) {
      tmpId = newForm(tmpId);
    }
    id_ = MTDDetId(tmpId).rawId();
  }

  /** Construct from complete geometry information, v1 **/
  BTLDetId(uint32_t zside, uint32_t rod, uint32_t module, uint32_t modtyp, uint32_t crystal)
      : MTDDetId(DetId::Forward, ForwardSubdetector::FastTime) {
    id_ |= (MTDType::BTL & kMTDsubdMask) << kMTDsubdOffset | (zside & kZsideMask) << kZsideOffset |
           (rod & kRodRingMask) << kRodRingOffset | (module & kBTLoldModuleMask) << kBTLoldModuleOffset |
           (modtyp & kBTLoldModTypeMask) << kBTLoldModTypeOffset |
           ((crystal - 1) & kBTLCrystalMask) << kBTLCrystalOffset;
    id_ |= kBTLNewFormat;
  }

  /** Construct from complete geometry information, v2, v3 **/
  BTLDetId(uint32_t zside, uint32_t rod, uint32_t runit, uint32_t dmodule, uint32_t smodule, uint32_t crystal)
      : MTDDetId(DetId::Forward, ForwardSubdetector::FastTime) {
    //RU, DM, SM & Xtal numbers start from 0
    id_ |= (MTDType::BTL & kMTDsubdMask) << kMTDsubdOffset | (zside & kZsideMask) << kZsideOffset |
           (rod & kRodRingMask) << kRodRingOffset | ((runit - 1) & kBTLRUMask) << kBTLRUOffset |
           ((dmodule - 1) & kBTLdetectorModMask) << kBTLdetectorModOffset |
           ((smodule - 1) & kBTLsensorModMask) << kBTLsensorModOffset |
           ((crystal - 1) & kBTLCrystalMask) << kBTLCrystalOffset;
    id_ |= kBTLNewFormat;
  }

  // ---------- Common methods ----------

  /** Returns BTL crystal number. */
  inline int crystal() const { return ((id_ >> kBTLCrystalOffset) & kBTLCrystalMask) + 1;}

  /** Returns BTL crystal number in construction database. */
  inline int crystalConsDB() const { 
    if (crystal() == kCrystalsPerModuleV2 + 1) return -1;
    if (smodule() == 1) return kCrystalsPerModuleV2 - crystal();
    else return crystal()- 1;
  }

  /** Returns BTL detector module number. */
  inline int dmodule() const { return ((id_ >> kBTLdetectorModOffset) & kBTLdetectorModMask) + 1; }

  /** Returns BTL sensor module number. */
  inline int smodule() const { return ((id_ >> kBTLsensorModOffset) & kBTLsensorModMask) + 1; }

  /** Returns BTL module number [1-24] (OLD BTL NUMBERING). */
  inline int module() const {
    return (((dmodule() - 1) % kDModulesInRURow) * (kSModulesInDM * kDModulesInRUCol) + 1 +
            int((dmodule() - 1) / kDModulesInRURow) + kDModulesInRUCol * (smodule() - 1));
  }

  /** Returns BTL crystal type number [1-3] (OLD BTL NUMBERING). */
  inline int modType() const {
    int gRU = globalRunit();
    return int((gRU - 1) / kRUPerTypeV2 + 1);
  }

  /** Returns BTL readout unit number per type [1-2], from Global RU number [1-6]. */
  inline int runit() const { return ((globalRunit() - 1) % kRUPerTypeV2 + 1); }

  /** Returns BTL global readout unit number. */
  inline int globalRunit() const { return ((id_ >> kBTLRUOffset) & kBTLRUMask) + 1; }
  // old globalRU function
  // inline int globalRunit() const {
  // if (runit() == 0) {
  //   // pre-V2: build a RU identifier from available information
  //   return (module() - 1) / kModulePerTypeBarPhiFlat / kRUPerTypeV2 + 1;
  // } else if (runit() > 0 && modType() > 0) {
  //   // V2/V3: build global RU identifier from RU per type and type
  //   return (modType() - 1) * kRUPerTypeV2 + runit();
  // }
  // }

  /** return the row in GeomDet language **/
  inline int row(unsigned nrows = kCrystalsPerModuleV2) const {
    return (crystal() - 1) % nrows;  // anything else for now
  }

  /** return the column in GeomDetLanguage **/
  inline int column(unsigned nrows = kCrystalsPerModuleV2) const { return (crystal() - 1) / nrows; }

  /** create a Geographical DetId for Tracking **/
  BTLDetId geographicalId(CrysLayout lay) const;

  /** conversion from old to new BTLDetID**/
  uint32_t newForm(const uint32_t& rawid) {
    uint32_t fixedP = rawid & (0xFFFFFFFF - kBTLoldFieldMask);  // unchanged part of id

    // convert old module number into detector module + sensor module numbers
    uint32_t oldModule = (rawid >> kBTLoldModuleOffset) & kBTLoldModuleMask;
    uint32_t detModule = int((oldModule - 1) / (kDModulesInRUCol * kSModulesInDM)) + 1 +
                         kDModulesInRURow * (int((oldModule % 6) / kSModulesInDM));
    uint32_t senModule = (oldModule - 1) % kSModulesInDM;
    // uint32_t detModule = ((int((oldModule - 1)/kDModulesInRUCol)/kSModulesInDM)*3) + ((oldModule-1)%3 + 1);
    // uint32_t senModule = int((oldModule - 1)/kDModulesInRUCol)%2;

    // convert old RU and type number into new RU number
    uint32_t oldRU = (rawid >> kBTLoldRUOffset) & kBTLoldRUMask;
    uint32_t oldType = (rawid >> kBTLoldModTypeOffset) & kBTLoldModTypeMask;
    uint32_t newRU = ((oldType - 1) >> 1) + oldRU;

    // get crystal number
    uint32_t crystal = (rawid & kBTLCrystalMask) >> kBTLCrystalOffset;

    // return new BTLDetID for v3 geom
    return (fixedP | (newRU & kBTLRUMask) << kBTLRUOffset | (detModule & kBTLdetectorModMask) << kBTLdetectorModOffset |
            (senModule & kBTLsensorModMask) << kBTLsensorModOffset |
            ((crystal & kBTLCrystalMask) << kBTLCrystalOffset)) |
           kBTLNewFormat;
  }
};

std::ostream& operator<<(std::ostream&, const BTLDetId&);

#endif  // DataFormats_BTLDetId_BTLDetId_h
