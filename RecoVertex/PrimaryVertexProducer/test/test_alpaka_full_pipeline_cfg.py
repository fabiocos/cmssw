"""
Test configuration for the COMPLETE Alpaka GNN vertex pipeline.

Pipeline:
  TrackFeatureProducerAlpaka → GNNVertexProducerAlpaka → PrimaryVertexProducer (GNN2D_alpaka)

Usage:
  cmsRun RecoVertex/PrimaryVertexProducer/test/test_alpaka_full_pipeline_cfg.py
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('AlpakaGNNFullPipeline', Phase2C17I13M9)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Accelerators_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/02098e56-4e94-4d1d-afea-d901f9842456.root'
    ),
)

# Enable logging for Alpaka producers
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.GNNVertexProducerAlpaka = cms.untracked.PSet()
process.MessageLogger.TrackFeatureSource = cms.untracked.PSet()
process.MessageLogger.GNNClusterizerFromAlpaka = cms.untracked.PSet()
process.MessageLogger.PrimaryVertexProducer = cms.untracked.PSet()
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# =============================================================================
# ALPAKA PRODUCERS (Inference happens here)
# =============================================================================

# 1. TrackFeatureProducerAlpaka - Generate track features from real tracks
#    For now using TrackFeatureSource which generates dummy features
#    TODO: Create TrackFeatureProducerAlpaka that consumes real generalTracks
process.trackFeatureSource = cms.EDProducer("alpaka_serial_sync::vertexgnn::TrackFeatureSource",
    numTracks = cms.int32(50),  # Will be overridden by actual track count in production
    seed = cms.int32(42),
    verbose = cms.untracked.bool(True),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("")
    )
)

# 2. GNNVertexProducerAlpaka - Run PyTorchAlpaka inference
process.gnnVertexProducer = cms.EDProducer("alpaka_serial_sync::vertexgnn::GNNVertexProducerAlpaka",
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
    trackFeatures = cms.InputTag("trackFeatureSource"),
    numSlots = cms.int32(180),  # Must match model's num_slots
    verbose = cms.untracked.bool(True),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("")
    )
)

# =============================================================================
# Standalone Vertex Builder (Alternative to PrimaryVertexProducer integration)
# This consumes the Alpaka SoA outputs directly
# =============================================================================
process.gnnVertexBuilderAlpaka = cms.EDProducer("GNNVertexBuilderFromAlpaka",
    tracks = cms.InputTag("generalTracks"),
    slotPredictions = cms.InputTag("gnnVertexProducer"),
    assignments = cms.InputTag("gnnVertexProducer"),
    numSlots = cms.int32(180),  # Must match model's num_slots
    existenceThreshold = cms.double(0.5),
    trackAssignmentThreshold = cms.double(0.4),
    verbose = cms.untracked.bool(True)
)

# =============================================================================
# TFileService for outputs
# =============================================================================
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("alpaka_gnn_output.root")
)

# =============================================================================
# PATH
# =============================================================================
process.alpaka_gnn_path = cms.Path(
    process.trackFeatureSource +
    process.gnnVertexProducer +
    process.gnnVertexBuilderAlpaka
)

# Output
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('alpaka_gnn_vertices.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_gnnVertexBuilderAlpaka_*_*',
        'keep *_gnnVertexProducer_*_*',
    )
)
process.output_step = cms.EndPath(process.output)

print("=" * 70)
print("Alpaka GNN Full Pipeline Test")
print("=" * 70)
print("Pipeline:")
print("  TrackFeatureSource → GNNVertexProducerAlpaka → GNNVertexBuilderFromAlpaka")
print("=" * 70)
