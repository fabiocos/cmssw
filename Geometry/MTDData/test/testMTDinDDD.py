import FWCore.ParameterSet.Config as cms

process = cms.Process("CompareGeometryTest")
process.load('Configuration.Geometry.GeometryExtended2023Dmtd_cff')

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )

process.source = cms.Source("EmptyIOVSource",
                            lastValue = cms.uint64(1),
                            timetype = cms.string('runnumber'),
                            firstValue = cms.uint64(1),
                            interval = cms.uint64(1)
                            )

process.myprint = cms.OutputModule("AsciiOutputModule")

process.testBTL = cms.EDAnalyzer("TestMTDNumbering",
                               label = cms.untracked.string(''),
                               isMagField = cms.untracked.bool(False),
                               outFileName = cms.untracked.string('BTL'),
                               numNodesToDump = cms.untracked.uint32(0),
                               ddTopNodeName = cms.untracked.string('btl:BarrelTimingLayer')
                               )

process.testETL = cms.EDAnalyzer("TestMTDNumbering",
                               label = cms.untracked.string(''),
                               isMagField = cms.untracked.bool(False),
                               outFileName = cms.untracked.string('ETL'),
                               numNodesToDump = cms.untracked.uint32(0),
                               ddTopNodeName = cms.untracked.string('etl:EndcapTimingLayer')
                               )

process.MessageLogger = cms.Service("MessageLogger",
                                    cout = cms.untracked.PSet( INFO = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
                                                               noLineBreaks = cms.untracked.bool(True),
                                                               DEBUG = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
                                                               threshold = cms.untracked.string('DEBUG'),
                                                               ),
                                    # For LogDebug/LogTrace output...
                                    debugModules = cms.untracked.vstring('testBTL','testETL'),
                                    categories = cms.untracked.vstring('TestMTDNumbering','MTDGeom'),
                                    destinations = cms.untracked.vstring('cout')
                                    )

process.Timing = cms.Service("Timing")

process.p1 = cms.Path(process.testBTL+process.testETL)

process.e1 = cms.EndPath(process.myprint)
