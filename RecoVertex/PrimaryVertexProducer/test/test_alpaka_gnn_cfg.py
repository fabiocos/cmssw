"""
Test configuration for PyTorchAlpaka-based GNN vertex producer.

Usage:
  cmsRun RecoVertex/PrimaryVertexProducer/test/test_alpaka_gnn_cfg.py
"""

import FWCore.ParameterSet.Config as cms

process = cms.Process("AlpakaGNNTest")

# Message logger - enable INFO level for our producers
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = "INFO"
# Enable our producer's log category
process.MessageLogger.GNNVertexProducerAlpaka = cms.untracked.PSet()
process.MessageLogger.TrackFeatureSource = cms.untracked.PSet()
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

# Enable Alpaka and GPU support
process.load("Configuration.StandardSequences.Accelerators_cff")
process.PyTorchService = cms.Service("PyTorchService")

# Number of events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2))

# Empty source
process.source = cms.Source("EmptySource")

# =============================================================================
# ALPAKA PRODUCERS
# =============================================================================

# 1. TrackFeatureSource: Generate dummy track features
process.trackFeatureSource = cms.EDProducer("alpaka_serial_sync::vertexgnn::TrackFeatureSource",
    numTracks = cms.int32(50),
    seed = cms.int32(42),
    verbose = cms.untracked.bool(True),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("")
    )
)

# 2. GNNVertexProducerAlpaka: Run PyTorch inference
process.gnnVertexProducer = cms.EDProducer("alpaka_serial_sync::vertexgnn::GNNVertexProducerAlpaka",
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
    trackFeatures = cms.InputTag("trackFeatureSource"),
    numSlots = cms.int32(200),
    verbose = cms.untracked.bool(True),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("")
    )
)

# =============================================================================
# PATH
# =============================================================================
process.path = cms.Path(
    process.trackFeatureSource +
    process.gnnVertexProducer
)
