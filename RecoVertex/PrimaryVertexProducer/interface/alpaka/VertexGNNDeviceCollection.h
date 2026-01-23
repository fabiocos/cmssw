#ifndef RecoVertex_PrimaryVertexProducer_interface_alpaka_VertexGNNDeviceCollection_h
#define RecoVertex_PrimaryVertexProducer_interface_alpaka_VertexGNNDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  using TrackFeaturesDeviceCollection = PortableCollection<::vertexgnn::TrackFeaturesSoA>;
  using SlotPredictionsDeviceCollection = PortableCollection<::vertexgnn::SlotPredictionsSoA>;
  using AssignmentDeviceCollection = PortableCollection<::vertexgnn::AssignmentSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

#endif  // RecoVertex_PrimaryVertexProducer_interface_alpaka_VertexGNNDeviceCollection_h
