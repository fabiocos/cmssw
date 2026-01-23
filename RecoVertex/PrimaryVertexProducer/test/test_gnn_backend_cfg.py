"""
Test configuration to compare ONNX vs TorchScript GNN inference.
Usage:
  cmsRun test_gnn_backend_cfg.py backend=onnx
  cmsRun test_gnn_backend_cfg.py backend=torchscript
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Setup command line options
options = VarParsing('analysis')
options.register('backend', 'onnx',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Inference backend: 'onnx' or 'torchscript'")
options.register('numTracks', 50,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Number of dummy tracks to generate")
options.register('numSlots', 200,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Number of vertex slots (K)")
options.parseArguments()

process = cms.Process("GNNBackendTest")

# Common services
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

# Empty source - we generate dummy data
process.source = cms.Source("EmptySource")

# Dummy track feature generator (for testing without real data)
# This produces the same output as TrackFeatureProducerAlpaka would
process.dummyTrackFeatures = cms.EDProducer("DummyTrackFeatureProducer",
    numTracks = cms.int32(options.numTracks),
    numFeatures = cms.int32(13),
    seed = cms.int32(42)
)

# Select backend
if options.backend.lower() == 'onnx':
    print("=" * 60)
    print("Using ONNX Runtime backend")
    print("=" * 60)
    
    # Standard ONNX-based GNN clusterizer (existing)
    process.gnnVertexProducer = cms.EDProducer("PrimaryVertexProducer",
        # Use GNN clustering
        TkClusParameters = cms.PSet(
            algorithm = cms.string("GNN2D_vect"),
            gnnOnnxFile = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.onnx"),
            gnnExistenceThreshold = cms.double(0.5),
            gnnNumSlots = cms.int32(options.numSlots),
        ),
        # Minimal config for testing
        verbose = cms.untracked.bool(True),
        TrackLabel = cms.InputTag("generalTracks"),
        beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    )
    
elif options.backend.lower() == 'torchscript':
    print("=" * 60)
    print("Using TorchScript (PyTorchAlpaka) backend")
    print("=" * 60)
    
    # PyTorchAlpaka-based producer
    process.gnnVertexProducerAlpaka = cms.EDProducer("vertexgnn::GNNVertexProducerAlpakaSerialSync@alpaka",
        model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
        trackFeatures = cms.InputTag("trackFeatureProducerAlpaka"),
        numSlots = cms.int32(options.numSlots),
        existenceThreshold = cms.double(0.5),
        verbosity = cms.untracked.int32(1)
    )
    
    # Track feature producer (prepares SoA)
    process.trackFeatureProducerAlpaka = cms.EDProducer("vertexgnn::TrackFeatureProducerAlpakaSerialSync@alpaka",
        tracks = cms.InputTag("generalTracks"),
        verbosity = cms.untracked.int32(1)
    )

else:
    raise ValueError(f"Unknown backend: {options.backend}. Use 'onnx' or 'torchscript'")

# Output
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(f'gnn_test_{options.backend}.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_gnn*_*_*',
    )
)

# Path
if options.backend.lower() == 'onnx':
    process.p = cms.Path(process.gnnVertexProducer)
else:
    process.p = cms.Path(
        process.trackFeatureProducerAlpaka +
        process.gnnVertexProducerAlpaka
    )

process.e = cms.EndPath(process.out)

# Print configuration summary
print(f"\nConfiguration Summary:")
print(f"  Backend:     {options.backend}")
print(f"  Num tracks:  {options.numTracks}")
print(f"  Num slots:   {options.numSlots}")
print(f"  Output:      gnn_test_{options.backend}.root")
