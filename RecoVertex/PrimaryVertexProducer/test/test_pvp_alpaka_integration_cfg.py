"""
Integration Test: PrimaryVertexProducer with GNN2D_alpaka Algorithm

This test mirrors the vertexTask_cfg.py flow but uses GNN2D_alpaka instead of GNN2D_vect.
The key difference is that GNN2D_alpaka requires upstream Alpaka producers for inference.

Flow comparison:
  GNN2D_vect (ONNX):
    unsortedOfflinePrimaryVerticesGNN (runs ONNX inference internally)

  GNN2D_alpaka (PyTorchAlpaka):
    TrackFeatureProducer -> GNNVertexProducerAlpaka -> unsortedOfflinePrimaryVerticesGNNAlpaka
                                                          ^
                                                     (consumes SoA from upstream)
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('REVTXAlpaka', Phase2C17I13M9)

# Standard configurations (same as vertexTask_cfg.py)
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

# Logging (same as vertexTask_cfg.py)
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.cerr.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.cerr.DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.PrimaryVertexProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNClusterizer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.TrackFeatureProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNVertexProducerAlpaka = cms.untracked.PSet(limit = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# =============================================================================
# ALPAKA UPSTREAM PRODUCERS (required for GNN2D_alpaka)
# =============================================================================

# Step 1: TrackFeatureProducer (standard EDProducer)
process.trackFeatureProducer = cms.EDProducer("vertexgnn::TrackFeatureProducer",
    tracks = cms.InputTag("generalTracks"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    trackTimesLabel = cms.InputTag("tofPID4DnoPID:t0safe"),
    trackTimeResosLabel = cms.InputTag("tofPID4DnoPID:sigmat0safe"),
    trackMTDTimeQualityVMapTag = cms.InputTag("mtdTrackQualityMVA:mtdQualMVA"),
    trackAssocSrc = cms.InputTag("trackExtenderWithMTD:generalTrackassoc"),
    tmtdSrc = cms.InputTag("trackExtenderWithMTD:generalTracktmtd"),
    sigmatmtdSrc = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd"),
    pathmtd = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength"),
    btlMatchChi2Src = cms.InputTag("trackExtenderWithMTD", "btlMatchChi2"),
    btlMatchTimeChi2Src = cms.InputTag("trackExtenderWithMTD", "btlMatchTimeChi2"),
    etlMatchChi2Src = cms.InputTag("trackExtenderWithMTD", "etlMatchChi2"),
    etlMatchTimeChi2Src = cms.InputTag("trackExtenderWithMTD", "etlMatchTimeChi2"),
    tofPi = cms.InputTag("trackExtenderWithMTD:generalTrackTofPi"),
    tofK = cms.InputTag("trackExtenderWithMTD:generalTrackTofK"),
    tofP = cms.InputTag("trackExtenderWithMTD:generalTrackTofP"),
    sigmatofpiSrc = cms.InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"),
    sigmatofkSrc = cms.InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"),
    sigmatofpSrc = cms.InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"),
    npixBarrelSrc = cms.InputTag("trackExtenderWithMTD", "npixBarrel"),
    npixEndcapSrc = cms.InputTag("trackExtenderWithMTD", "npixEndcap"),
    minTrackTimeQuality = cms.double(0.8),
    useMVACut = cms.bool(False),
    verbose = cms.untracked.bool(False)
)

# Step 2: GNNVertexProducerAlpaka (Alpaka PyTorch inference)
process.gnnVertexProducer = cms.EDProducer("alpaka_serial_sync::vertexgnn::GNNVertexProducerAlpaka",
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
    trackFeatures = cms.InputTag("trackFeatureProducer"),
    numSlots = cms.int32(200),
    verbose = cms.untracked.bool(False),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("")
    )
)

# =============================================================================
# GNN VERTEX PRODUCER - GNN2D_alpaka (mirrors unsortedOfflinePrimaryVerticesGNN)
# =============================================================================
# Clone from the existing GNN producer and modify to use Alpaka backend
process.unsortedOfflinePrimaryVerticesGNNAlpaka = process.unsortedOfflinePrimaryVerticesGNN.clone(
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_alpaka"),
        TkDAClusParameters = cms.PSet(
            # Alpaka-specific: consume SoA from upstream producers
            slotPredictions = cms.InputTag("gnnVertexProducer"),
            assignments = cms.InputTag("gnnVertexProducer"),
            # Same parameters as GNN2D_vect
            existenceThreshold = cms.double(0.5),
            trackAssignmentThreshold = cms.double(0.0),
            verbose = cms.untracked.bool(True),
        )
    ),
)

# =============================================================================
# PATH (same flow as vertexTask_cfg.py, but with upstream Alpaka producers)
# =============================================================================
process.exe = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.vertexreco +
    # Alpaka upstream producers
    process.trackFeatureProducer +
    process.gnnVertexProducer +
    # PrimaryVertexProducer with GNN2D_alpaka
    process.unsortedOfflinePrimaryVerticesGNNAlpaka
)

# Output (same as vertexTask_cfg.py)
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:gnn_alpaka_integrated.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_unsortedOfflinePrimaryVerticesGNNAlpaka_*_*',
        'keep *_gnnVertexProducer_*_*',
        'keep *_trackFeatureProducer_*_*',
    ),
    splitLevel = cms.untracked.int32(0)
)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Performance monitoring (same as vertexTask_cfg.py)
from Validation.Performance.TimeMemorySummary import customiseWithTimeMemorySummary
process = customiseWithTimeMemorySummary(process)

print("=" * 70)
print("vertexTask Alpaka - GNN2D_alpaka for CMSSW_16_1_0_pre1")
print("=" * 70)
print("Mode: PyTorchAlpaka (upstream inference + PVP vertex building)")
print("Flow: TrackFeatureProducer -> GNNVertexProducerAlpaka -> PVP(GNN2D_alpaka)")
print("Parameters:")
print("  - existenceThreshold: 0.5")
print("  - numSlots: 200")
print("=" * 70)
