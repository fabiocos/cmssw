#ifndef DataFormats_MTDDetId_MTDDetId_h
#define DataFormats_MTDDetId_MTDDetId_h

#include "DataFormats/DetId/interface/DetId.h"
#include <ostream>

/** 
    @class MTDDetId
    @brief Detector identifier base class for the MIP Timing Layer.
*/

class MTDDetId : public DetId {
  
 protected:
  
  /** Enumerated type for Forward sub-deteector systems. */
  enum SubDetector { subUNKNOWN=0, BTL=1, ETL=2 };
  
  static const uint32_t kZsideOffset               = 23;
  static const uint32_t kZsideMask                 = 0x3;
  static const uint32_t kRodRingOffset             = 17;
  static const uint32_t kRodRingMask               = 0x3F;
  
 public:
  
  // ---------- Constructors, enumerated types ----------
  
  /** Construct a null id */
 MTDDetId()  : DetId() {;}
  
  /** Construct from a raw value */
 MTDDetId( const uint32_t& raw_id ) : DetId( raw_id ) {;}
  
  /** Construct from generic DetId */
 MTDDetId( const DetId& det_id )  : DetId( det_id.rawId() ) {;}
  
  /** Construct and fill only the det and sub-det fields. */
 MTDDetId( Detector det, int subdet ) : DetId( det, subdet ) {;}
  
  // ---------- Common methods ----------
  
  /** Returns enumerated type specifying MTD sub-detector. */
  inline SubDetector subDetector() const { return static_cast<MTDDetId::SubDetector>(subdetId()); }
  
  /** Returns MTD side, i.e. Z-=1 or Z+=2. */
  inline int mtdSide() const { return (id_>>kZsideOffset)&kZsideMask; }
  
  /** Returns MTD rod/ring number. */
  inline int mtdRR() const { return (id_>>kRodRingOffset)&kRodRingMask; }
  
};

std::ostream& operator<< ( std::ostream&, const MTDDetId& );

#endif // DataFormats_MTDDetId_MTDDetId_h

