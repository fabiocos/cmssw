#ifndef DIGIFTL_ETLDIGI_H
#define DIGIFTL_ETLDIGI_H

#include <ostream>
#include <vector>

#include "DataFormats/FTLDigi/interface/ETLSample.h"

/**
   @class ETLDigi
   @short wrapper for a data word
 */

class ETLDigi {
public:

  /**
     @short CTOR
  */
  ETLDigi() : data1_(0), data2_(0) {}
  ETLDigi(uint8_t data1, uint32_t data2)
      :  data1_(data1), data2_(data2) {}
  ETLDigi(const ETLDigi& o) : data1_(o.data1_), data2_(o.data2_) {}

  /**
     @short setters
  */

  inline void set(uint8_t pl1, uint32_t pl2) {
    data1_ = pl1;
    data2_ = pl2;
  }

  /**
     @short getters
  */

  // to be fixed with the correct field implementation
  inline uint32_t row() const { return data1_; }
  inline uint32_t column() const { return data1_; }
  inline uint32_t data() const { return data2_; }

  /**
     @general utilities
  */

  inline void clear() {
    data1_ = 0;
    data2_ = 0;
  }

  //static bool fillWithSampleArray(const std::vector<ETLSample> in, ETLDigi *out);

private:

  // ETROC channel 40 bits payload split into two separate words:
  // uint8_t -> row,col
  // uint32_t -> remaining part of the payload
  uint8_t data1_;
  uint32_t data2_;
};

std::ostream& operator<<(std::ostream&, const ETLDigi&);

#endif
