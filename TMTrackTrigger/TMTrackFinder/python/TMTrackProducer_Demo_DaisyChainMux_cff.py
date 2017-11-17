import FWCore.ParameterSet.Config as cms

#---------------------------------------------------------------------------------------------------------
# This describes a chain consisting only of the GP + HT, where the GP divides the tracker into 18 eta 
# sectors (each sub-divided into 2 virtual eta subsectors) &  the HT is the new daisy-chain HT using a 
# 64x32 array, with its outputs multiplexed so as to solve the "busy sector" problem.
#---------------------------------------------------------------------------------------------------------

#=== Import default values for all parameters.

from TMTrackTrigger.TMTrackFinder.TMTrackProducer_Defaults_cfi import TMTrackProducer

#=== These changes result in a configuration like that for our 1st hardware demonstrator,
#=== running with Thomas Schuh's firmware.

#--- Cuts applied to stubs before arriving in L1 track finding board.

# Reduce number of bits used by front-end chips to store stub bend info. (CMS is considering this option proposed by Seb. Viret).
TMTrackProducer.StubCuts.BendResReduced = cms.bool(True)
# Don't use stubs whose measured Pt from bend info is significantly below HTArraySpec.HoughMinPt, where "significantly" means allowing for resolution in q/Pt derived from stub bend resolution specified below.
TMTrackProducer.StubCuts.KillLowPtStubs = cms.bool(True)
# Bend resolution assumed by bend filter in units of strip pitch. Also used when assigning stubs to sectors if EtaPhiSectors.CalcPhiTrkRes=True. And by the bend filter if HTFillingRphi.UseBendFilter=True.
TMTrackProducer.StubCuts.BendResolution = cms.double(1.25)

#--- Stub digitization.

TMTrackProducer.StubDigitize = cms.PSet(
   EnableDigitize  = cms.bool(True),   # Digitize stub coords? If not, use floating point coords.
   FirmwareType    = cms.uint32(1),    # 0 = Old Thomas 2-cbin data format, 1 = new Thomas data format for daisy chain, 2-98 = reserved for demonstrator daisy chain use, 99 = Systolic array data format.
   #--- Parameters available in MP board.
   PhiSectorBits   = cms.uint32(6),    # Bits used to store phi sector number
   PhiSBits        = cms.uint32(14),   # Bits used to store phiS coord.
   PhiSRange       = cms.double(0.78539816340),  # Range phiS coord. covers in radians.
   RtBits          = cms.uint32(10),    # Bits used to store Rt coord.
   RtRange         = cms.double(103.0382), # Range Rt coord. covers in units of cm.
   ZBits           = cms.uint32(12),   # Bits used to store z coord.
   ZRange          = cms.double(640.), # Range z coord. covers in units of cm.
   #--- Parameters available in GP board (excluding any in common with MP specified above).
   PhiOBits        = cms.uint32(15),      # Bits used to store PhiO parameter.
   PhiORange       = cms.double(1.5707963268), # Range PhiO parameter covers.
   BendBits        = cms.uint32(6)        # Bits used to store stub bend.
)

#--- Mark's proposal to reduce resource usage. (9/8/2016 meeting)

TMTrackProducer.PhiSectors.NumPhiSectors           = cms.uint32(16)
TMTrackProducer.EtaSectors.EtaRegions              = cms.vdouble(-2.4,-2.16,-1.95,-1.7,-1.43,-1.16,-0.89,-0.61,-0.31,0.0,0.31,0.61,0.89,1.16,1.43,1.7,1.95,2.16,2.4)
TMTrackProducer.EtaSectors.ChosenRofZ              = cms.double(50.)     
TMTrackProducer.EtaSectors.AllowOver2EtaSecs       = cms.bool(True)
TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi      = cms.uint32(64)
TMTrackProducer.HTArraySpecRphi.NumSubSecsEta      = cms.uint32(2)

#--- Rules governing how stubs are filled into the r-phi Hough Transform array.

# Use filter in each r-phi HT cell, filling it only with stubs that have consistent bend information?
# The assumed bend resolution is specified in StubCuts.BendResolution.
TMTrackProducer.HTFillingRphi.UseBendFilter = cms.bool(True)
# Use filter in each HT cell, preventing more than the specified number of stubs being stored in the cell. (Reflecting memory limit of hardware). N.B. 16 would reflect hardware, but results unreliable as depend on assumed order of stubs.
TMTrackProducer.HTFillingRphi.MaxStubsInCell = cms.uint32(16) # Setting this to anything more than 99 disables this option
# If BusySectorKill = True, and more than BusySectorNumStubs stubs are assigned to tracks by an r-phi HT array, then the excess tracks are killed, with lowest Pt ones killed first. This is because HT hardware has finite readout time.
TMTrackProducer.HTFillingRphi.BusySectorKill     = cms.bool(True)
TMTrackProducer.HTFillingRphi.BusySectorNumStubs = cms.uint32(210)
# If BusyInputSectorKill = True, and more than BusyInputSectorNumStubs are input to the HT array from the GP, then
# the excess stubs are killed. This is because HT hardware has finite readin time.
# Results unreliable as depend on assumed order of stubs.
TMTrackProducer.HTFillingRphi.BusyInputSectorKill  = cms.bool(True)
TMTrackProducer.HTFillingRphi.BusyInputSectorNumStubs = cms.uint32(210)

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

#--- Official TP to track matching criteria.

TMTrackProducer.GenCuts.GenMinStubLayers               = cms.uint32(4)
TMTrackProducer.TrackMatchDef.MinNumMatchLayers        = cms.uint32(4)

#--- Disable r-z Hough transform, disable r-z track filters, disable duplicate track removal & disable track fit.

TMTrackProducer.HTArraySpecRz.EnableRzHT   = cms.bool(False)  
TMTrackProducer.RZfilterOpts.UseSeedFilter = cms.bool(False)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi  = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRz    = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgFit   = cms.uint32(0)
TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring()
