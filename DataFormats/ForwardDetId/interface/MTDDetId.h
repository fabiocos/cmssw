#ifndef DataFormats_MTDDetId_MTDDetId_h
#define DataFormats_MTDDetId_MTDDetId_h

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include <ostream>

/** 
    @class MTDDetId
    @brief Detector identifier base class for the MIP Timing Layer.
*/

class MTDDetId : public DetId {
  
 protected:
  
  /** Enumerated type for Forward sub-deteector systems. */
  enum SubDetector { subUNKNOWN=0, FastTime=1 };
  
  /** Enumerated type for MTD sub-deteector systems. */
  enum MTDType { typeUNKNOWN=0, BTL=1, ETL=2 };
  
  static const uint32_t kMTDtypeOffset             = 23;
  static const uint32_t kMTDtypeMask               = 0x3;
  static const uint32_t kZsideOffset               = 22;
  static const uint32_t kZsideMask                 = 0x1;
  static const uint32_t kRodRingOffset             = 15;
  static const uint32_t kRodRingMask               = 0x7F;
  
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
  
  /** Returns enumerated type specifying CMS sub-detector. */
  inline SubDetector subDetector() const { return static_cast<MTDDetId::SubDetector>(subdetId()); }
  
  /** Returns enumerated type specifying MTD sub-detector, i.e. BTL or ETL. */
  inline int mtdType() const { return (id_>>kMTDtypeOffset)&kMTDtypeMask; }
  
  /** Returns MTD side, i.e. Z-=1 or Z+=2. */
  inline int mtdSide() const { return (id_>>kZsideOffset)&kZsideMask; }
  
  /** Returns MTD rod/ring number. */
  inline int mtdRR() const { return (id_>>kRodRingOffset)&kRodRingMask; }
  
};

std::ostream& operator<< ( std::ostream&, const MTDDetId& );

#endif // DataFormats_MTDDetId_MTDDetId_h

