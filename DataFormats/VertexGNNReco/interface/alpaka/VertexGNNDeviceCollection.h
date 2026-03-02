#ifndef DataFormats_VertexGNNReco_interface_alpaka_VertexGNNDeviceCollection_h
#define DataFormats_VertexGNNReco_interface_alpaka_VertexGNNDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/VertexGNNReco/interface/VertexGNNSoA.h"
#include "DataFormats/VertexGNNReco/interface/VertexGNNHostCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  // Input: Track features [N, 13]
  using TrackFeaturesDeviceCollection = PortableCollection<::vertexgnn::TrackFeaturesSoA>;

  // Output: Combined GNN outputs [N, K] with Eigen columns
  using GNNOutputDeviceCollection = PortableCollection<::vertexgnn::GNNOutputSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

#endif  // DataFormats_VertexGNNReco_interface_alpaka_VertexGNNDeviceCollection_h
