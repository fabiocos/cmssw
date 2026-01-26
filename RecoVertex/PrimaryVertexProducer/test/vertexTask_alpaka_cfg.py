"""
vertexTask_alpaka_cfg.py - Alpaka GNN Vertex Reconstruction (Drop-in replacement)

Same structure as vertexTask_cfg.py but uses Alpaka backend for GPU inference.
Replaces unsortedOfflinePrimaryVerticesGNN in vertexreco with Alpaka pipeline.
"""

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('REVTX', Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Alpaka + CUDA for GPU inference
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32, cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/02098e56-4e94-4d1d-afea-d901f9842456.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

# Debug logging for GNN modules
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.cerr.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.cerr.DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.PrimaryVertexProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNClusterizerFromAlpaka = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.TrackFeatureProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNVertexProducerAlpaka = cms.untracked.PSet(limit = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(10),
    numberOfStreams = cms.untracked.uint32(10),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# ============================================================================
# ALPAKA GNN Vertex Producer - GPU Inference
# ============================================================================

# Step 1: TrackFeatureProducer (CPU, produces HostCollection)
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
    useTrackFilter = cms.bool(True),
    verbose = cms.untracked.bool(False),
    # Use same TkFilterParameters as unsortedOfflinePrimaryVerticesGNN
    TkFilterParameters = process.unsortedOfflinePrimaryVerticesGNN.TkFilterParameters.clone(),
)

# Step 2: GNNVertexProducerAlpaka (GPU inference)
process.gnnVertexProducer = cms.EDProducer("vertexgnn::GNNVertexProducerAlpaka@alpaka",
    trackFeatures = cms.InputTag("trackFeatureProducer"),
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/vertex_slot_production.pt"),
    verbose = cms.untracked.bool(False),
)

# Step 3: Replace unsortedOfflinePrimaryVerticesGNN with Alpaka version
process.unsortedOfflinePrimaryVerticesGNN = process.unsortedOfflinePrimaryVerticesGNN.clone(
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_alpaka"),
        TkDAClusParameters = cms.PSet(
            existenceThreshold = cms.double(0.01),  # Low threshold for untrained model debugging
            trackAssignmentThreshold = cms.double(0.0),
            gnnOutput = cms.InputTag("gnnVertexProducer"),
            verbose = cms.untracked.bool(False),
        )
    ),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Alpaka GNN Vertex Reconstruction (GPU)'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
outputCommand = process.FEVTDEBUGHLTEventContent.outputCommands
outputCommand.append('drop *_offlinePrimaryVertices_*_RECO')
outputCommand.append('drop *_offlinePrimaryVertices4D_*_RECO')
outputCommand.append('drop *_offlinePrimaryVertices4DWithBS_*_RECO')
outputCommand.append('drop *_offlinePrimaryVerticesWithBS_*_RECO')
outputCommand.append('drop *_TriggerResults_*_RECO')

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:revtx_step3_alpaka.root'),
    outputCommands = outputCommand,
    splitLevel = cms.untracked.int32(0)
)

# TFileService for GNN track inspector output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("gnn_track_inspector_alpaka.root")
)

# ============================================================================
# GNN Track Inspector (same as vertexTask)
# ============================================================================
process.gnnInspector = cms.EDAnalyzer("GNNTrackInspector",
    pvModule = cms.InputTag("unsortedOfflinePrimaryVerticesGNN", "", "REVTX"),
    trackSrc = cms.InputTag("generalTracks"),
    printFirstN = cms.uint32(15),
    dropNaNs = cms.bool(True),
)

# ============================================================================
# Path definitions (same structure as vertexTask)
# ============================================================================
# Add Alpaka producers before vertexreco since we're replacing the GNN in it
process.alpaka_producers = cms.Sequence(
    process.trackFeatureProducer *
    process.gnnVertexProducer
)

process.exe = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.alpaka_producers *
    process.vertexreco
)
process.inspect_path = cms.EndPath(process.gnnInspector)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Performance monitoring
from Validation.Performance.TimeMemorySummary import customiseWithTimeMemorySummary
process = customiseWithTimeMemorySummary(process)

print("=" * 70)
print("vertexTask_alpaka_cfg.py - Alpaka GNN VertexSlotModel (GPU)")
print("=" * 70)
print("Mode: GPU Inference via Alpaka/PyTorchAlpaka")
print("Parameters:")
print("  - existenceThreshold: 0.5")
print("  - model: vertex_slot_model.pt")
print("  - backend: CUDA (GPU)")
print("=" * 70)
