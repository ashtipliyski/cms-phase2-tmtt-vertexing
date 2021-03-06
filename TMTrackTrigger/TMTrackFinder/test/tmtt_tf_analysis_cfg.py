#################################################################################################
# To run execute do
# cmsRun tmtt_tf_analysis_cfg.py Events=50 inputMC=Samples/Muons/PU0.txt histFile=outputHistFile.root
# where the arguments take default values if you don't specify them. You can change defaults below.
#################################################################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

options = VarParsing.VarParsing ('analysis')

#--- Specify input MC
#options.register('inputMC', '../../../Samples/Muons_fixed_pT_10GeV/PU0.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../Samples/Electrons_FixedPt_10GeV/PU0.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
# options.register('inputMC', '../PU140scratch.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
# options.register('inputMC', '../SamplesCMS/MonoJetPU200.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
options.register('inputMC','./3July_0pileup.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")

#--- Specify number of events to process.
options.register('Events',-1,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"Number of Events to analyze")

#--- Specify name of output histogram file.
options.register('histFile','Hist.root',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output histogram file")
#--- Specify if stubs need to be produced i.e. they are not available in the input file
options.register('makeStubs',0,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"Make stubs on the fly")

options.parseArguments()

#--- input and output

list = FileUtils.loadListFromFile(options.inputMC)
readFiles = cms.untracked.vstring(*list)
# readFiles = cms.untracked.vstring("/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU200_DES23_62_V1-v1/110000/004C20AB-4D9E-E611-AE77-00266CFFBDAC.root")
secFiles = cms.untracked.vstring()

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.histFile)
)

process.source = cms.Source ("PoolSource",
                            fileNames = readFiles,
                            secondaryFileNames = secFiles,
                            # skipEvents = cms.untracked.uint32(500)
                            )


process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

#--- Load code that produces our L1 tracks and makes corresponding histograms.

#--- Either use this one for studies of the final 2025 system.
process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_cff')
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_SLR_cff')
#--- Or this one for studies of the 2015 demonstrator with Thomas Schuh's new "daisy chain" firmware,
#--- that includes MUX to solve busy sector problem & Mark's 18 eta sectors.
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_DaisyChainMux_cff')
#--- Or this one for studies of the 2015 demonstrator with Thomas Schuh's original "daisy chain" firmware.
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_DaisyChain_cff')
#--- Or this one for studies of the 2015 demonstrator with Thomas Schuh's old firmware.
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_Thomas_cff')
#--- Or this one for studies of the 2015 demonstrator with systolic array firmware.
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_Systolic_cff')

#--- Optionally override default configuration parameters here (example given of how).

#process.TMTrackProducer.HTArraySpecRz.EnableRzHT = cms.bool(True)
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")
process.p = cms.Path(process.TrackTriggerAssociatorClustersStubs * process.TMTrackProducer)
#process.e = cms.EndPath(process.out)
