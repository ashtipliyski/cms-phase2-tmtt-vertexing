import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.Config as cms

# options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')
process = cms.Process('TMTrackAnalyzer')

options.register('skipEvents',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to skip")
options.register('selBx',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Select readout Bx")
options.register('selAllBx',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run over all Bx")
options.register('histFile',
                 'tmTrackHistos.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Name of output histogram file")
options.register('makeStubs',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Make stubs on the fly")

options.parseArguments()

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
    skipEvents=cms.untracked.uint32(options.skipEvents)
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Output definition
# TTree output file
process.load("CommonTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName = cms.string(options.histFile)

# enable debug message logging for our modules
process.MessageLogger = cms.Service(
    "MessageLogger",
    threshold  = cms.untracked.string('DEBUG'),
    categories = cms.untracked.vstring('*'),
    debugModules = cms.untracked.vstring('*')
)


# histograms
process.load('TMTrackTrigger.TMTrackFinder.tmTrackAnalyzer_cfi')
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "keep *_analyzer_*_*",
        "keep *_TMTrackAnalyzer_*_*"
        )
)


# Path and EndPath definitions
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")
process.p = cms.Path(process.TrackTriggerAssociatorClustersStubs * process.TMTrackAnalyzer)
process.e = cms.EndPath(process.out)

