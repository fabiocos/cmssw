#ifndef DIGIFTL_NEWFTLDATAFRAMET_H
#define DIGIFTL_NEWFTLDATAFRAMET_H

#include <vector>
#include <ostream>
#include <iostream>

/**
   @class NewFTLDataFrameT
*/

template <class D, class S, class DECODE>
class NewFTLDataFrameT {
public:
  /**
     @short key to sort the collection
  */
  typedef D key_type;

  /**
     @short CTOR
  */
  NewFTLDataFrameT() : id_(0) {}
  NewFTLDataFrameT(const D& id) : id_(id) {}
  NewFTLDataFrameT(const NewFTLDataFrameT& o) : data_(o.data_), id_(o.id_) {}

  /**
    @short det id
  */
  const D& id() const { return id_; }

  /**
   @short row
   */
  const int row() const { return DECODE::row(id_, data_); }

  /**
   @short column
   */
  const int column() const { return DECODE::col(id_, data_); }

  /**
     @short assess/set specific samples
  */
  const S& sample() const { return data_; }
void setSample(const S& sample) {
      data_ = sample;
  }
  void print(std::ostream& out = std::cout) {
    out << data_.print(out);
  }

private:

  //channel payload
  S data_;

  // det id for this data frame
  D id_;

};

#endif
