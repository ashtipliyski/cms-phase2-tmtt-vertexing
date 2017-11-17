
#########################################################################################################
# This is a copy of Louise Skinnari's script L1TrackNupleMaker_cfg.py obtained via
# getLouiseAnalysisCode.csh .
# It has been modified to read our MC and to run our L1 track producer.
#
# To run execute do
# cmsRun L1TrackNtupleMaker_cfg.py Events=50 inputMC=Samples/Muons/PU0.txt histFile=outputHistFile.root trkFitAlgo=All
# where the arguments take default values if you don't specify them. You can change defaults below.
#
# The option trkFitAlgo should be equal to the name of the TMT track fitting algorithm that you wish to use,
# which must be one of the fitters specified by option TrackFitters in TMTrackProducer_Defaults_cfi.py.
# Alternatively, you can set it equal to "Tracklet4" or "Tracklet5" to use the Tracklet groups 
# tracks with 4 or 5 helix parameters.  Or to "All" to use all our track fitters at once.
#########################################################################################################


############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("L1TrackNtuple")
 
#################################################################################
# Options such as choice of MC sample (specific to us)
################################################################################

options = VarParsing.VarParsing ('analysis')

#--- Specify input MC
#options.register('inputMC', '../../../SamplesCMSseb/Muon_Pt100/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../SamplesCMSseb/Electron_Pt35/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../SamplesCMSseb/StubFix/TTbar/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
options.register('inputMC', '../../../SamplesCMSseb/StubFix/TTbar/PU200.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../SamplesDecemberReview/TTbar/PU0.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")

#--- Specify number of events to process.
options.register('Events',1000,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"Number of Events to analyze")

#--- Specify name of output histogram file.
options.register('histFile','Hist.root',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output histogram file")

options.register('trkFitAlgo','All',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Nuame of track helix fit algorithm")
# options.register('trkFitAlgo','TrackFitLinearAlgo4',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Nuame of track helix fit algorithm")
# options.register('trkFitAlgo','Tracklet4',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Nuame of track helix fit algorithm")

options.parseArguments()

#--- Decode which track fitting algorithm results are wanted for

# Note name of track fitting algorithm.
trkFitAlgo = options.trkFitAlgo
runOverAllTrackTypes = ( trkFitAlgo.find("All") != -1 )

fit5param = None
# Are results wanted for the 4 or 5 parameter helix fit?
if ( (trkFitAlgo.find("4") == -1 and trkFitAlgo.find("5") != -1) or (trkFitAlgo.find("5") == -1 and trkFitAlgo.find("4") != -1) ) :
    # If this is true, it was a 5 param fit; else 4 param.
    fit5param = (trkFitAlgo.find("5") != -1)
elif not runOverAllTrackTypes:
    print "### ERROR: Unable to figure out if helix fit used 4 or 5 params from its algorithm name"
# Are our TMT tracks to be used? If not, then the Tracklet group's tracks will be used.
useTMTtracks = (trkFitAlgo.find("Tracklet") == -1) and not runOverAllTrackTypes
# Instead of running one one particular set of tracks, run on all know tracks
# i.e. Tracklet group, TMT with different track fitters

if runOverAllTrackTypes and useTMTtracks:
  print '### ERROR: Only one of runOverAllTrackTypes and useTMTtracks should be true'
# Figure out name of our TTTrack collection
if (useTMTtracks):
  trackName = "TML1Tracks" + trkFitAlgo
else:
  trackName = "Level1TTTracks"

if (fit5param):
    print "### Will produce results for TTTrack collection ",trackName," with 5 helix parameters."
else:
    print "### Will produce results for TTTrack collection ",trackName," with 4 helix parameters."

############################################################
# input and output (specific to us, as we use our own MC)
############################################################

list = FileUtils.loadListFromFile(options.inputMC)
readFiles = cms.untracked.vstring(*list)
secFiles = cms.untracked.vstring()

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.histFile)
)

process.source = cms.Source ("PoolSource",
                            fileNames = readFiles,
                            secondaryFileNames = secFiles,
                            # skipEvents = cms.untracked.uint32(28),
                            )

process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))
 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

############################################################
# Path definitions & schedule
############################################################

#run the tracking
BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.TT_step = cms.Path(process.TrackTriggerTTTracks)
process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)

############################################################
# primary vertex producer 
############################################################

process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkPrimaryVertexProducer_cfi")
process.pL1TkPrimaryVertex = cms.Path( process.L1TkPrimaryVertex )

process.L1TkPrimaryVertexMC = process.L1TkPrimaryVertex.clone()
process.L1TkPrimaryVertexMC.MonteCarloVertex = cms.bool( True )
process.pL1TkPrimaryVertexMC = cms.Path( process.L1TkPrimaryVertexMC )

############################################################

############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# Valid options are:
#      single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets from ttbar = 6
#      pions from taus = 15
#      inclusively, store all TPs (also those not from primary interaction, if available in samples you are running on) = 1
############################################################

# Modified this to specify our tracks & associator to be used.
process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(1.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.4),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks" ),               # TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), # MCTruth input 
                                       L1PVInputTag = cms.InputTag("L1TkPrimaryVertex", "" ),
                                       L1MCPVInputTag = cms.InputTag("L1TkPrimaryVertexMC", "" ),
                                       ## isolation stuff 
                                       TrackIsolation = cms.bool(True),
                                       # cuts on the central object (the truth muon & track matched to it)
                                       PTmin = cms.double(20.0),           # central object pt > X GeV, ptmin < 0 means no cut applied
                                       ETAmax = cms.double(2.4),           # central object |eta| < X
                                       TrackPTmin = cms.double(20.0),      # for track matched to central object
                                       TrackETAmax = cms.double(2.4),
                                       TrackChi2max = cms.double(100.0),
                                       TrackNStubmin = cms.int32(4),
                                       # cuts on the tracks (used to determine isolation variable)
                                       IsoTrackZmax = cms.double(25.),     # in cm
                                       IsoTrackChi2max = cms.double(1e10),
                                       IsoTrackNStubmin = cms.int32(4),    
                                       IsoTrackPTmin = cms.double(3.0),    # in GeV
                                       IsoDRmin = cms.double(0.0),
                                       IsoDRmax = cms.double(0.3),
                                       IsoDZmax = cms.double(1e10),        # in cm
                                       ## tracking in jets stuff (--> requires AK4 genjet collection present!)
                                       TrackingInJets = cms.bool(True),
                                       JetPTmin = cms.double(30),
                                       ## save primary vertex information? (--> requires that you ran that above)
                                       PrimaryVertex = cms.bool(True),
                                       )

process.ana = cms.Path(process.L1TrackNtuple)

############################################################
############################################################
############################################################
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
############################################################

# OLD, FOR TTI SAMPLES
# from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI
# process = cust_2023TTI(process)

# # NEW, FOR SCOPE DOCUMENT REL VAL SAMPLES
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D
process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

#=== EVERYTHING BELOW WAS ADDED BY THE TMT GROUP.

# Set number of helix fit parameters that analysis code should get results for.
if (fit5param):
    process.L1TrackNtuple.L1Tk_nPar = cms.int32(5)


# If use wants plots for our TMT L1 tracks instead of Tracklet group's ones, then the following is run.
if (useTMTtracks or runOverAllTrackTypes):

    #--- Load code that produces our L1 tracks and makes corresponding histograms.
    #--- Either use this one for studies of the final 2025 system.
    process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_cff')

    #--- Optionally override default configuration parameters here (example given of how).

    #process.TMTrackProducer.RZfilterOpts.UseSeedFilter = cms.bool(False)

if (useTMTtracks):
    # Change tracking step to run our L1 track reconstruction algorithm instead of the tracklet group's one.
    process.TT_step = cms.Path(process.BeamSpotFromSim+process.TMTrackProducer)

    # Change track association to truth step to use our tracks instead of tracklet group's ones.
    process.TMTrackAssociator = process.TTTrackAssociatorFromPixelDigis.clone(
        TTTracks = cms.VInputTag(cms.InputTag("TMTrackProducer", trackName))
    )
    process.TTAssociator_step = cms.Path(process.TMTrackAssociator)

    # Run primary vertex finder using our tracks
    process.L1TkPrimaryVertex.L1TrackInputTag = cms.InputTag("TMTrackProducer", trackName)
    process.L1TkPrimaryVertexMC.L1TrackInputTag = cms.InputTag("TMTrackProducer", trackName)


    # Run analysis code over our L1 tracks.
    process.L1TrackNtuple.L1TrackInputTag      = cms.InputTag("TMTrackProducer", trackName)
    # And use the associator that matches our tracks to the truth tracks.
    process.L1TrackNtuple.MCTruthTrackInputTag = cms.InputTag("TMTrackAssociator", trackName)

    # And use the correct PV collection
#    process.L1TrackNtuple.L1PVInputTag      = cms.InputTag("L1TkPrimaryVertex", "")

elif (runOverAllTrackTypes):

  allTMTFitters = process.TMTrackProducer.TrackFitSettings.TrackFitters

  # for each type of TMT track (from each fitter)
  # Run the tracking
  # All track fitters (that are uncommented in TMTrackProducer.TrackFitSettings.TrackFitters)
  # are run by this step
  process.TT_step = cms.Path(process.BeamSpotFromSim+process.TMTrackProducer)

  # Get the truth association
  # Have to add a process for each track collection produced by each track fit
  process.TTAssociator_step = cms.Path()
  # Also run the (official) analysis code, again for each track type
  process.ana = cms.Path()

  process.pL1TkPrimaryVertex = cms.Path()
  process.pL1TkPrimaryVertexMC = cms.Path()

  if 'TrackFitLinearAlgo4' in allTMTFitters:
    process.TTAssociatorTrackFitLinearAlgo4 = process.TTTrackAssociatorFromPixelDigis.clone(
            TTTracks = cms.VInputTag(cms.InputTag("TMTrackProducer", 'TML1TracksTrackFitLinearAlgo4'))
        )
    process.TTAssociator_step *= process.TTAssociatorTrackFitLinearAlgo4

    process.L1PVTrackFitLinearAlgo4 = process.L1TkPrimaryVertex.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", "TML1TracksTrackFitLinearAlgo4")
      )
    process.L1PVTrackFitLinearAlgo4MC = process.L1TkPrimaryVertexMC.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", "TML1TracksTrackFitLinearAlgo4")
      )
    process.pL1TkPrimaryVertex *= process.L1PVTrackFitLinearAlgo4
    process.pL1TkPrimaryVertexMC *= process.L1PVTrackFitLinearAlgo4MC

    process.analyzerTrackFitLinearAlgo4 = process.L1TrackNtuple.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", 'TML1TracksTrackFitLinearAlgo4'),
        MCTruthTrackInputTag = cms.InputTag("TTAssociatorTrackFitLinearAlgo4", 'TML1TracksTrackFitLinearAlgo4'),
        L1PVInputTag = cms.InputTag("L1PVTrackFitLinearAlgo4",""),
        L1MCPVInputTag = cms.InputTag("L1PVTrackFitLinearAlgo4MC",""),
      )
    process.ana *= process.analyzerTrackFitLinearAlgo4

  if 'KF4ParamsComb' in allTMTFitters:
    process.TTAssociatorKF4ParamsComb = process.TTTrackAssociatorFromPixelDigis.clone(
            TTTracks = cms.VInputTag(cms.InputTag("TMTrackProducer", 'TML1TracksKF4ParamsComb'))
        )
    process.TTAssociator_step *= process.TTAssociatorKF4ParamsComb

    process.L1PVKF4ParamsComb = process.L1TkPrimaryVertex.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", "TML1TracksKF4ParamsComb")
      )
    process.L1PVKF4ParamsCombMC = process.L1TkPrimaryVertexMC.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", "TML1TracksKF4ParamsComb")
      )
    
    process.pL1TkPrimaryVertex *= process.L1PVKF4ParamsComb
    process.pL1TkPrimaryVertexMC *= process.L1PVKF4ParamsCombMC

    process.analyzerKF4ParamsComb = process.L1TrackNtuple.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", 'TML1TracksKF4ParamsComb'),
        MCTruthTrackInputTag = cms.InputTag("TTAssociatorKF4ParamsComb", 'TML1TracksKF4ParamsComb'),
        L1PVInputTag = cms.InputTag("L1PVKF4ParamsComb",""),
        L1MCPVInputTag = cms.InputTag("L1PVKF4ParamsCombMC",""),
      )
    process.ana *= process.analyzerKF4ParamsComb


  if 'globalLinearRegression' in allTMTFitters:
    process.TTAssociatorglobalLinearRegression = process.TTTrackAssociatorFromPixelDigis.clone(
            TTTracks = cms.VInputTag(cms.InputTag("TMTrackProducer", 'TML1TracksglobalLinearRegression'))
        )
    process.TTAssociator_step *= process.TTAssociatorglobalLinearRegression

    process.L1PVglobalLinearRegression = process.L1TkPrimaryVertex.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", "TML1TracksglobalLinearRegression")
      )
    process.L1PVglobalLinearRegressionMC = process.L1TkPrimaryVertexMC.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", "TML1TracksglobalLinearRegression")
      ) 
    process.pL1TkPrimaryVertex *= process.L1PVglobalLinearRegression
    process.pL1TkPrimaryVertexMC *= process.L1PVglobalLinearRegressionMC

    process.analyzerglobalLinearRegression = process.L1TrackNtuple.clone(
        L1TrackInputTag = cms.InputTag("TMTrackProducer", 'TML1TracksglobalLinearRegression'),
        MCTruthTrackInputTag = cms.InputTag("TTAssociatorglobalLinearRegression", 'TML1TracksglobalLinearRegression'),
        L1PVInputTag = cms.InputTag("L1PVglobalLinearRegression",""),
        L1MCPVInputTag = cms.InputTag("L1PVglobalLinearRegressionMC",""),
      )
    process.ana *= process.analyzerglobalLinearRegression

process.schedule = cms.Schedule(process.TT_step,process.TTAssociator_step,process.pL1TkPrimaryVertex,process.pL1TkPrimaryVertexMC,process.ana)
