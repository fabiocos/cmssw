#ifndef RecoVertex_PrimaryVertexProducer_interface_alpaka_VertexGNNDeviceCollection_h
#define RecoVertex_PrimaryVertexProducer_interface_alpaka_VertexGNNDeviceCollection_h

/**
 * Alpaka PortableCollection wrappers for GNN vertex producer SoA types.
 */

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  // Input: 13 features per track
  using TrackFeaturesDeviceCollection = PortableCollection<::vertexgnn::TrackFeaturesSoA>;
  
  // Output: Per-slot predictions (K elements)
  using SlotPredictionsDeviceCollection = PortableCollection<::vertexgnn::SlotPredictionsSoA>;
  
  // Output: Assignment matrix (N*K elements)
  using AssignmentDeviceCollection = PortableCollection<::vertexgnn::AssignmentSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

#endif  // RecoVertex_PrimaryVertexProducer_interface_alpaka_VertexGNNDeviceCollection_h
