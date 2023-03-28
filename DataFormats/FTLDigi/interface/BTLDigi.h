#ifndef DIGIFTL_BTLDIGI_H
#define DIGIFTL_BTLDIGI_H

#include <ostream>
#include <vector>

#include "DataFormats/FTLDigi/interface/BTLSample.h"

/**
   @class BTLDigi
   @short wrapper for BTL readout channel payload
 */

class BTLDigi {
public:

  /**
     @short CTOR
  */
  BTLDigi() : data1rhs_(0), data2rhs_(0), data1lhs_(0), data2lhs_(0) {}
  BTLDigi(uint32_t data1rhs, uint64_t data2rhs, uint32_t data1lhs, uint64_t data2lhs)
      :  data1rhs_(data1rhs), data2rhs_(data2rhs), data1lhs_(data1lhs), data2lhs_(data2lhs) {}
  BTLDigi(const BTLDigi& o) : data1rhs_(o.data1rhs_), data2rhs_(o.data2rhs_), data1lhs_(o.data1lhs_), data2lhs_(o.data2lhs_) {}

  /**
     @short setters
  */

  inline void set(uint32_t pl1rhs, uint64_t pl2rhs, uint32_t pl1lhs, uint64_t pl2lhs) {
    data1rhs_ = pl1rhs;
    data2rhs_ = pl2rhs;
    data1lhs_ = pl1lhs;
    data2lhs_ = pl2lhs;
  }

  /**
     @short getters
  */

  inline uint32_t plrhs1() const { return data1rhs_; }
  inline uint64_t plrhs2() const { return data2rhs_; }
  inline uint32_t pllhs1() const { return data1lhs_; }
  inline uint64_t pllhs2() const { return data2lhs_; }

  /**
     @general utilities
  */

  inline void clear() {
    data1rhs_ = 0;
    data2rhs_ = 0;
    data1lhs_ = 0;
    data2lhs_ = 0;
  }

  //static bool fillWithSampleArray(const std::vector<BTLSample> in, BTLDigi *out);

private:

  // TOFHIR channel 88 bits payload split into two separate words:
  // uint32_t -> bits 0 -> 22
  // uint64_t -> bits 23 -> 86
  // twofold payload for left hand side and right hand side SiPMs attached to the same crystal (defined BTLDetId)
  uint32_t data1rhs_;
  uint64_t data2rhs_;
  uint32_t data1lhs_;
  uint64_t data2lhs_;
};

std::ostream& operator<<(std::ostream&, const BTLDigi&);

#endif
