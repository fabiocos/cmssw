#ifndef DataFormats_FTLDigi_interface_BTLDigi_h
#define DataFormats_FTLDigi_interface_BTLDigi_h

#include <cstdint>
#include <ostream>

namespace btldigi {

  class BTLDigi {
  public:
    /**
     @short key to sort the collection
  */
    typedef uint32_t key_type;

    BTLDigi()
        : krawId_(0),
          kBC0count_(0),
          kstatus_(false),
          kBCcount_(0),
          kchIDR_(0),
          kT1coarseR_(0),
          kT2coarseR_(0),
          kEOIcoarseR_(0),
          kChargeR_(0),
          kT1fineR_(0),
          kT2fineR_(0),
          kIdleTimeR_(0),
          kPrevTrigFR_(0),
          kTACIDR_(0),
          kchIDL_(0),
          kT1coarseL_(0),
          kT2coarseL_(0),
          kEOIcoarseL_(0),
          kChargeL_(0),
          kT1fineL_(0),
          kT2fineL_(0),
          kIdleTimeL_(0),
          kPrevTrigFL_(0),
          kTACIDL_(0) {}

    BTLDigi(uint32_t rawId,
            uint16_t BC0count,
            bool status,
            uint32_t BCcount,
            uint8_t chIDR,
            uint16_t T1coarseR,
            uint16_t T2coarseR,
            uint16_t EOIcoarseR,
            uint16_t ChargeR,
            uint16_t T1fineR,
            uint16_t T2fineR,
            uint16_t IdleTimeR,
            uint8_t PrevTrigFR,
            uint8_t TACIDR,
            uint8_t chIDL,
            uint16_t T1coarseL,
            uint16_t T2coarseL,
            uint16_t EOIcoarseL,
            uint16_t ChargeL,
            uint16_t T1fineL,
            uint16_t T2fineL,
            uint16_t IdleTimeL,
            uint8_t PrevTrigFL,
            uint8_t TACIDL)
        : krawId_(rawId),
          kBC0count_(BC0count),
          kstatus_(status),
          kBCcount_(BCcount),
          kchIDR_(chIDR),
          kT1coarseR_(T1coarseR),
          kT2coarseR_(T2coarseR),
          kEOIcoarseR_(EOIcoarseR),
          kChargeR_(ChargeR),
          kT1fineR_(T1fineR),
          kT2fineR_(T2fineR),
          kIdleTimeR_(IdleTimeR),
          kPrevTrigFR_(PrevTrigFR),
          kTACIDR_(TACIDR),
          kchIDL_(chIDL),
          kT1coarseL_(T1coarseL),
          kT2coarseL_(T2coarseL),
          kEOIcoarseL_(EOIcoarseL),
          kChargeL_(ChargeL),
          kT1fineL_(T1fineL),
          kT2fineL_(T2fineL),
          kIdleTimeL_(IdleTimeL),
          kPrevTrigFL_(PrevTrigFL),
          kTACIDL_(TACIDL) {}

    uint32_t krawId() const { return krawId_; }
    uint16_t kBC0count() const { return kBC0count_; }
    bool kstatus() const { return kstatus_; }
    uint32_t kBCcount() const { return kBCcount_; }
    uint8_t kchIDR() const { return kchIDR_; }
    uint16_t kT1coarseR() const { return kT1coarseR_; }
    uint16_t kT2coarseR() const { return kT2coarseR_; }
    uint16_t kEOIcoarseR() const { return kEOIcoarseR_; }
    uint16_t kChargeR() const { return kChargeR_; }
    uint16_t kT1fineR() const { return kT1fineR_; }
    uint16_t kT2fineR() const { return kT2fineR_; }
    uint16_t kIdleTimeR() const { return kIdleTimeR_; }
    uint8_t kPrevTrigFR() const { return kPrevTrigFR_; }
    uint8_t kTACIDR() const { return kTACIDR_; }
    uint8_t kchIDL() const { return kchIDL_; }
    uint16_t kT1coarseL() const { return kT1coarseL_; }
    uint16_t kT2coarseL() const { return kT2coarseL_; }
    uint16_t kEOIcoarseL() const { return kEOIcoarseL_; }
    uint16_t kChargeL() const { return kChargeL_; }
    uint16_t kT1fineL() const { return kT1fineL_; }
    uint16_t kT2fineL() const { return kT2fineL_; }
    uint16_t kIdleTimeL() const { return kIdleTimeL_; }
    uint8_t kPrevTrigFL() const { return kPrevTrigFL_; }
    uint8_t kTACIDL() const { return kTACIDL_; }

    // needed for sorting in SortedCollection
    uint32_t id() const { return krawId_; }

  private:
    uint32_t krawId_;     // Raw ID of the module/TOFHIR
    uint16_t kBC0count_;  // BC0 count (reserved)
    bool kstatus_;        // status of the TOFHIR
    uint32_t kBCcount_;
    uint8_t kchIDR_;       // TOFHIR channel ID, right side of crystal
    uint16_t kT1coarseR_;  // data from crystal right side
    uint16_t kT2coarseR_;
    uint16_t kEOIcoarseR_;
    uint16_t kChargeR_;
    uint16_t kT1fineR_;
    uint16_t kT2fineR_;
    uint16_t kIdleTimeR_;
    uint8_t kPrevTrigFR_;
    uint8_t kTACIDR_;
    uint8_t kchIDL_;       // TOFHIR channel ID, left side of crystal
    uint16_t kT1coarseL_;  // data from crystal left side
    uint16_t kT2coarseL_;
    uint16_t kEOIcoarseL_;
    uint16_t kChargeL_;
    uint16_t kT1fineL_;
    uint16_t kT2fineL_;
    uint16_t kIdleTimeL_;
    uint8_t kPrevTrigFL_;
    uint8_t kTACIDL_;
  };

  std::ostream& operator<<(std::ostream&, const BTLDigi&);

}  // namespace btldigi
#endif  // DataFormats_FTLDigi_interface_BTLDigi_h
