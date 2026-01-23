#ifndef RecoVertex_PrimaryVertexProducer_interface_VertexGNNHostCollection_h
#define RecoVertex_PrimaryVertexProducer_interface_VertexGNNHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"

namespace vertexgnn {

  // Input: Track features [N, 13]
  using TrackFeaturesHostCollection = PortableHostCollection<TrackFeaturesSoA>;

  // Output: Combined GNN outputs [N, K] with Eigen columns
  using GNNOutputHostCollection = PortableHostCollection<GNNOutputSoA>;

}  // namespace vertexgnn

#endif  // RecoVertex_PrimaryVertexProducer_interface_VertexGNNHostCollection_h
