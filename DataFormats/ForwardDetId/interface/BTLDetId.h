#ifndef DataFormats_BTLDetId_BTLDetId_h
#define DataFormats_BTLDetId_BTLDetId_h

#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include <ostream>

/** 
    @class BTLDetId
    @brief Detector identifier class for the Barrel Timing Layer.
    The crystal count must start from 0, copy number must be scaled by 1 unit.
*/

class BTLDetId : public MTDDetId {
  
 private:
  
  static const uint32_t kBTLmoduleOffset           = 10;
  static const uint32_t kBTLmoduleMask             = 0x7F;
  static const uint32_t kBTLmodTypeOffset          = 8;
  static const uint32_t kBTLmodTypeMask            = 0x3;
  static const uint32_t kBTLCrystalOffset          = 0;
  static const uint32_t kBTLCrystalMask            = 0x3F;
  
 public:
  
  // ---------- Constructors, enumerated types ----------
  
  /** Construct a null id */
 BTLDetId() : MTDDetId( DetId::MTD, MTDDetId::BTL ) {;}

  /** Construct from a raw value */
 BTLDetId( const uint32_t& raw_id ) : MTDDetId( raw_id ) {;}
  
  /** Construct from generic DetId */
 BTLDetId( const DetId& det_id )  : MTDDetId( det_id.rawId() ) {;}
  
  /** Construct and fill only the det and sub-det fields. */
 BTLDetId( uint32_t zside, 
           uint32_t rod, 
           uint32_t module, 
           uint32_t modtyp, 
           uint32_t crystal ) : MTDDetId( DetId::MTD, MTDDetId::BTL ) {
    id_ |= ( zside& kZsideMask ) << kZsideOffset |
      ( rod& kRodRingMask ) << kRodRingOffset |
      ( module& kBTLmoduleMask ) << kBTLmoduleOffset |
      ( modtyp& kBTLmodTypeMask ) << kBTLmodTypeOffset |
      ( (crystal-1)& kBTLCrystalMask ) << kBTLCrystalOffset ; 
}

// ---------- Common methods ----------

/** Returns BTL module number. */
inline int btlModule() const { return (id_>>kBTLmoduleOffset)&kBTLmoduleMask; }

/** Returns BTL crystal type number. */
inline int btlmodType() const { return (id_>>kBTLmodTypeOffset)&kBTLmodTypeMask; }

/** Returns BTL crystal number. */
 inline int btlCrystal() const { return ((id_>>kBTLCrystalOffset)&kBTLCrystalMask) + 1; }

};

std::ostream& operator<< ( std::ostream&, const BTLDetId& );

#endif // DataFormats_BTLDetId_BTLDetId_h

