#ifndef RecoVertex_PrimaryVertexProducer_interface_VertexGNNHostCollection_h
#define RecoVertex_PrimaryVertexProducer_interface_VertexGNNHostCollection_h

/**
 * Host-side PortableCollection types for GNN vertex producer.
 * These are needed for ROOT dictionaries and EDM product storage.
 */

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"

namespace vertexgnn {

  using TrackFeaturesHostCollection = PortableHostCollection<TrackFeaturesSoA>;
  using SlotPredictionsHostCollection = PortableHostCollection<SlotPredictionsSoA>;
  using AssignmentHostCollection = PortableHostCollection<AssignmentSoA>;

}  // namespace vertexgnn

#endif  // RecoVertex_PrimaryVertexProducer_interface_VertexGNNHostCollection_h
