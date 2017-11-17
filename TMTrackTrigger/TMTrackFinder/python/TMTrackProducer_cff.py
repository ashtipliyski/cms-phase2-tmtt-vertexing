import FWCore.ParameterSet.Config as cms

#---------------------------------------------------------------------------------------------------------
# This describes the full track reconstruction chain, where: the GP divides the tracker into 18 eta 
# sectors (each sub-divided into 2 virtual eta subsectors);  the HT is the new daisy-chain HT using a 
# 64x32 array, with its outputs multiplexed so as to solve the "busy sector" problem;
# and the duplicate removal is Algo50 run after track fit.
#---------------------------------------------------------------------------------------------------------

#=== Import default values for all parameters.

from TMTrackTrigger.TMTrackFinder.TMTrackProducer_Defaults_cfi import TMTrackProducer

#=== Change below any parameters you don't want to take their default values.

# Don't order stubs by bend in DTC, such that highest Pt stubs are transmitted first.
#TMTrackProducer.StubCuts.OrderStubsByBend = cms.bool(False)

#--- Use Thomas's improved Hough transform (high granularity in cells with Pt > 6 GeV; and 2 eta subsectors;
#--- and outputting +ve and -ve charged tracks on separate links.

#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt       = cms.uint32(64) 
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi      = cms.uint32(64)
#TMTrackProducer.HTArraySpecRphi.EnableMerge2x2     = cms.bool(True)
#TMTrackProducer.HTArraySpecRphi.MaxPtToMerge2x2    = cms.double(6.)
#TMTrackProducer.HTArraySpecRphi.NumSubSecsEta      = cms.uint32(2)
#TMTrackProducer.HTFillingRphi.BusySectorEachCharge = cms.bool(True)

#--- Alternatively use Mark's proposal to reduce resource usage. (9/8/2016 meeting)

TMTrackProducer.StubDigitize.PhiSRange             = cms.double(0.78539816340) # Double digitisation range to cope with larger phi sectors.
TMTrackProducer.PhiSectors.NumPhiSectors           = cms.uint32(16)
TMTrackProducer.EtaSectors.EtaRegions              = cms.vdouble(-2.4,-2.16,-1.95,-1.7,-1.43,-1.16,-0.89,-0.61,-0.31,0.0,0.31,0.61,0.89,1.16,1.43,1.7,1.95,2.16,2.4)
TMTrackProducer.EtaSectors.ChosenRofZ              = cms.double(50.)     
TMTrackProducer.EtaSectors.AllowOver2EtaSecs       = cms.bool(True)
TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi      = cms.uint32(64)
TMTrackProducer.HTArraySpecRphi.NumSubSecsEta      = cms.uint32(2)
TMTrackProducer.HTFillingRphi.BusySectorEachCharge = cms.bool(True)

#--- Kill tracks that HT doesn't have time to output during time-multiplexed period.
#--- (Reduces efficiency, but is more realistic).

TMTrackProducer.HTFillingRphi.BusySectorKill       = cms.bool(True)
# Also killed those that the GP doesn't have time to send to HT.
TMTrackProducer.HTFillingRphi.BusyInputSectorKill  = cms.bool(True)

#--- Reduce efficiency loss from tracks that the HT can't output during time multiplexed period,
#--- by subdividing the HT m-bin (q/Pt) axis & sending each division to a different output link.
#--- (Only relevant if you have set BusySectorKill = True).
TMTrackProducer.HTFillingRphi.MuxOutputsHT = cms.bool(True)

# This is the standard solution for the HTmux.
TMTrackProducer.HTFillingRphi.BusySectorMbinRanges = cms.vuint32(6,5,5,5,5,6)

# Or this is the solution for the new HTmux aimed at giving better load balancing in the KF fit.
#TMTrackProducer.HTFillingRphi.MuxOutputsHTforKF = cms.bool(True)
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges = cms.vuint32(3,3,3,3,2,2,3,3,3,3,2,2)
#TMTrackProducer.HTFillingRphi.BusySectorMbinOrder  = cms.vuint32(16,21,26, 17,22,27, 18,23,28, 19,24,29, 20,30, 25,31, 15,10,5, 14,9,4, 13,8,3, 12,7,2, 11,1, 6,0)

#--- Switch on duplicate track removal algorithm that runs after the track fitter. (Ian's talk in 26/8/2016 meeting)
#--- Requires all other duplicate track removal to be disabled.

TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi  = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRz    = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgFit   = cms.uint32(50)

#--- Reduce requirement on number of layers a track must have stubs in, to improve efficiency with official matching criteria.
# This reduces it globally.
#TMTrackProducer.L1TrackDef.MinStubLayers       = cms.uint32(4)
# This reduces it in the specified eta sectors (barrel-endcap transition & extreme forward regions).
TMTrackProducer.L1TrackDef.EtaSecsReduceLayers = cms.vuint32(5,12)

#--- Official TP to track matching criteria.

TMTrackProducer.GenCuts.GenMinStubLayers               = cms.uint32(4)
TMTrackProducer.TrackMatchDef.MinNumMatchLayers        = cms.uint32(4)

#--- Switch off parts of the track reconstruction chain.

#TMTrackProducer.RZfilterOpts.UseSeedFilter    = cms.bool(True)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi   = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRz     = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg  = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgFit    = cms.uint32(0)
#TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring()

#--- Switch on digitisation, using an extra bit in r & z and one less in phi than we have in hardware,
#--- since given an extra few weeks, we would change the hardware to this cfg
#process.TMTrackProducer.StubDigitize.EnableDigitize  = cms.bool(True)
# Remove one bit.
#process.TMTrackProducer.StubDigitize.PhiSBits        = cms.uint32(13)
#process.TMTrackProducer.StubDigitize.PhiOBits        = cms.uint32(14)
# Add one extra bit.
#process.TMTrackProducer.StubDigitize.RtBits          = cms.uint32(11)
#process.TMTrackProducer.StubDigitize.ZBits           = cms.uint32(13)

#--- Do not reduce bend resolution of stubs
TMTrackProducer.StubCuts.BendResReduced = cms.bool(False)
