import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------------------------------------------------

#
# OUT-OF-DATE FILE: Use TMTrackProducer_new_cff.py instead, unless you are completing old studies.
#
# This describes the full track reconstruction chain, where: the GP divides the tracker into 32 phi X 9 eta sectors,
# & the HT is the original daisy-chain HT using a 32x32 array; and the duplicate removal is Algo8 run before
# the track fit.
# In addition, the TP to track matching criteria are the 5 layers ones.
#--------------------------------------------------------------------------------------------------------------------------


#=== Import default values for all parameters.

from TMTrackTrigger.TMTrackFinder.TMTrackProducer_Defaults_cfi import TMTrackProducer

#=== Change below any parameters you don't want to take their default values.

#--- Use Thomas's improved Hough transform (high granularity in cells with Pt > 6 GeV; and 2 eta subsectors;
#--- and outputting +ve and -ve charged tracks on separate links.

#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt       = cms.uint32(64) 
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi      = cms.uint32(64)
#TMTrackProducer.HTArraySpecRphi.EnableMerge2x2     = cms.bool(True)
#TMTrackProducer.HTArraySpecRphi.MaxPtToMerge2x2    = cms.double(6.)
#TMTrackProducer.HTArraySpecRphi.NumSubSecsEta      = cms.uint32(2)
#TMTrackProducer.HTFillingRphi.BusySectorEachCharge = cms.bool(True)

#--- Alternatively use Mark's proposal to reduce resource usage. (9/8/2016 meeting)

#TMTrackProducer.StubDigitize.PhiSRange             = cms.double(0.78539816340) # Double digitisation range to cope with larger phi sectors.
#TMTrackProducer.PhiSectors.NumPhiSectors           = cms.uint32(16)
#TMTrackProducer.EtaSectors.EtaRegions              = cms.vdouble(-2.4,-2.16,-1.95,-1.7,-1.43,-1.16,-0.89,-0.61,-0.31,0.0,0.31,0.61,0.89,1.16,1.43,1.7,1.95,2.16,2.4)
#TMTrackProducer.EtaSectors.ChosenRofZ              = cms.double(50.)     
#TMTrackProducer.EtaSectors.AllowOver2EtaSecs       = cms.bool(True)
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi      = cms.uint32(64)
#TMTrackProducer.HTArraySpecRphi.NumSubSecsEta      = cms.uint32(2)
#TMTrackProducer.HTFillingRphi.BusySectorEachCharge = cms.bool(True)

#--- Kill tracks that HT doesn't have time to output during time-multiplexed period.
#--- (Reduces efficiency, but is more realistic).

#TMTrackProducer.HTFillingRphi.BusySectorKill       = cms.bool(True)

#--- Reduce efficiency loss from tracks that the HT can't output during time multiplexed period,
#--- by subdividing the HT m-bin (q/Pt) axis & sending each division to a different output link.
#--- (Only relevant if you have set BusySectorKill = True).
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges = cms.vuint32(6,5,5,5,5,6)
#TMTrackProducer.HTFillingRphi.MuxOutputsHT = cms.bool(True)

#--- Switch on duplicate track removal algorithm that runs after the track fitter. (Ian's talk in 26/8/2016 meeting)
#--- Requires all other duplicate track removal to be disabled.

#TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi  = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRz    = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgFit   = cms.uint32(50)

#--- Old TP to track matching criteria.

TMTrackProducer.GenCuts.GenMinStubLayers               = cms.uint32(5)
TMTrackProducer.TrackMatchDef.MinNumMatchLayers        = cms.uint32(5)

#--- Switch off parts of the track reconstruction chain.

#TMTrackProducer.RZfilterOpts.UseSeedFilter    = cms.bool(False)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi   = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRz     = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg  = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgFit    = cms.uint32(0)
#TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring()



