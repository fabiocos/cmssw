#ifndef DIGIFTL_FTLDIGICOLLECTION_H
#define DIGIFTL_FTLDIGICOLLECTION_H

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/FTLDigi/interface/FTLDataFrameT.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/FTLDigi/interface/BTLSample.h"
#include "DataFormats/FTLDigi/interface/ETLSample.h"

#include "DataFormats/FTLDigi/interface/NewFTLDataFrameT.h"
#include "DataFormats/FTLDigi/interface/BTLDigi.h"
#include "DataFormats/FTLDigi/interface/ETLDigi.h"

namespace mtdhelpers {

  struct BTLRowColDecode {
    static inline int row(const DetId& id, const std::vector<BTLSample>& data) { return data.front().row(); }
    static inline int col(const DetId& id, const std::vector<BTLSample>& data) { return data.front().column(); }
  };

  struct ETLRowColDecode {
    static inline int row(const DetId& id, const std::vector<ETLSample>& data) { return data.front().row(); }
    static inline int col(const DetId& id, const std::vector<ETLSample>& data) { return data.front().column(); }
  };

  struct NewBTLRowColDecode {
    static inline int row(const DetId& id, const BTLDigi& data) { return static_cast<BTLDetId>(id).row(); }
    static inline int col(const DetId& id, const BTLDigi& data) { return static_cast<BTLDetId>(id).column(); }
  };

  struct NewETLRowColDecode {
    static inline int row(const DetId& id, const ETLDigi& data) { return data.row(); }
    static inline int col(const DetId& id, const ETLDigi& data) { return data.column(); }
  };
}  // namespace mtdhelpers

typedef FTLDataFrameT<BTLDetId, BTLSample, mtdhelpers::BTLRowColDecode> BTLDataFrame;
typedef edm::SortedCollection<BTLDataFrame> BTLDigiCollection;

typedef FTLDataFrameT<ETLDetId, ETLSample, mtdhelpers::ETLRowColDecode> ETLDataFrame;
typedef edm::SortedCollection<ETLDataFrame> ETLDigiCollection;

typedef NewFTLDataFrameT<BTLDetId, BTLDigi, mtdhelpers::NewBTLRowColDecode> NewBTLDataFrame;
typedef edm::SortedCollection<NewBTLDataFrame> NewBTLDigiCollection;

typedef NewFTLDataFrameT<ETLDetId, ETLDigi, mtdhelpers::NewETLRowColDecode> NewETLDataFrame;
typedef edm::SortedCollection<NewETLDataFrame> NewETLDigiCollection;

#endif
