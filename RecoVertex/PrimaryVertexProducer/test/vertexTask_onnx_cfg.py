"""
vertexTask_onnx_cfg.py - ONNX GNN Vertex Reconstruction (GPU)

Same structure as vertexTask_alpaka_cfg.py but uses ONNX backend for GPU inference.
Uses GNN2D_vect algorithm with CUDA backend.
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
process.MessageLogger.cerr.GNNClusterizer = cms.untracked.PSet(limit = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(10),
    numberOfStreams = cms.untracked.uint32(10),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# ============================================================================
# ONNX GNN Vertex Producer - GPU Inference
# ============================================================================

# Replace unsortedOfflinePrimaryVerticesGNN with ONNX CUDA version
process.unsortedOfflinePrimaryVerticesGNN = process.unsortedOfflinePrimaryVerticesGNN.clone(
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_vect"),
        TkDAClusParameters = cms.PSet(
            existenceThreshold = cms.double(0.5),
            trackAssignmentThreshold = cms.double(0.4),
            numSlots = cms.int32(180),  # Must match model's num_slots
            nnVersion = cms.string("dummy_vertex_slot"),
            onnxBackend = cms.string("CUDA"),  # GPU inference
            onnxModelPath = cms.FileInPath('RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.onnx'),
            Tmin = cms.double(4.0),
            Tpurge = cms.double(4.0),
            Tstop = cms.double(2.0),
            vertexSize = cms.double(0.01),
            d0CutOff = cms.double(3.0),
            verbose = cms.untracked.bool(False),
        )
    ),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('ONNX GNN Vertex Reconstruction (GPU)'),
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
    fileName = cms.untracked.string('file:revtx_step3_onnx.root'),
    outputCommands = outputCommand,
    splitLevel = cms.untracked.int32(0)
)

# TFileService for GNN track inspector output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("gnn_track_inspector_onnx.root")
)

# ============================================================================
# GNN Track Inspector
# ============================================================================
process.gnnInspector = cms.EDAnalyzer("GNNTrackInspector",
    pvModule = cms.InputTag("unsortedOfflinePrimaryVerticesGNN", "", "REVTX"),
    trackSrc = cms.InputTag("generalTracks"),
    printFirstN = cms.uint32(15),
    dropNaNs = cms.bool(True),
)

# ============================================================================
# Path definitions
# ============================================================================
process.exe = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.vertexreco
)
process.inspect_path = cms.EndPath(process.gnnInspector)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Performance monitoring
from Validation.Performance.TimeMemorySummary import customiseWithTimeMemorySummary
process = customiseWithTimeMemorySummary(process)

print("=" * 70)
print("vertexTask_onnx_cfg.py - ONNX GNN VertexSlotModel (GPU)")
print("=" * 70)
print("Mode: GPU Inference via ONNX Runtime CUDA")
print("Parameters:")
print("  - existenceThreshold: 0.5")
print("  - model: dummy_vertex_slot.onnx")
print("  - backend: CUDA (GPU)")
print("=" * 70)
