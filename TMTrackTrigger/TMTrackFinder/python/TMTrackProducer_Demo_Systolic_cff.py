import FWCore.ParameterSet.Config as cms

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
TMTrackProducer.StubCuts.BendResolution = cms.double(1.30)

#--- Stub digitization.

TMTrackProducer.StubDigitize = cms.PSet(
   EnableDigitize  = cms.bool(True),   # Digitize stub coords? If not, use floating point coords.
   FirmwareType    = cms.uint32(99),   # 0 = Old Thomas 2-cbin data format, 1 = new Thomas data format for daisy chain, 2-98 = reserved for demonstrator daisy chain use, 99 = Systolic array data format.
   #--- Parameters available in MP board.
   PhiSectorBits   = cms.uint32(6),    # Bits used to store phi sector number
   PhiSBits        = cms.uint32(12),   # Bits used to store phiS coord.
   PhiSRange       = cms.double(1.),   # Range phiS coord. covers in radians.
   RtBits          = cms.uint32(9),    # Bits used to store Rt coord.
   RtRange         = cms.double(256.), # Range Rt coord. covers in units of cm.
   ZBits           = cms.uint32(12),   # Bits used to store z coord.
   ZRange          = cms.double(1024), # Range z coord. covers in units of cm.
   DPhiBits        = cms.untracked.uint32(8),    # Bits used to store Delta(phi) track angle.
   DPhiRange       = cms.untracked.double(1.),   # Range Delta(phi) covers in radians.
   RhoBits         = cms.untracked.uint32(6),    # Bits used to store rho parameter.
   RhoRange        = cms.untracked.double(0.25), # Range rho parameter covers. 
   #--- Parameters available in GP board (excluding any in common with MP specified above).
   PhiOBits        = cms.uint32(15),      # Bits used to store PhiO parameter.
   PhiORange       = cms.double(1.57078), # Range PhiO parameter covers.
   BendBits        = cms.uint32(6)        # Bits used to store stub bend.
)

#--- Phi sector definition

# Number of phi sectors.
TMTrackProducer.PhiSectors.NumPhiSectors = cms.uint32(64)
# Use phi of track at this radius for assignment of stubs to phi sectors & also for one of the axes of the r-phi HT. If ChosenRofPhi=0, then use track phi0.
TMTrackProducer.PhiSectors.ChosenRofPhi  = cms.double(45.)

#--- r-phi Hough transform definition

# Number of rows or columns in HT array in each coordinate (ignored if HoughNcellsRphi > 0).
TMTrackProducer.HTArraySpecRphi.HoughNbinsPt  = cms.uint32(36)
TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi = cms.uint32(14)

#--- Rules governing how stubs are filled into the r-phi Hough Transform array.

# Use filter in each r-phi HT cell, filling it only with stubs that have consistent bend information?
# The assumed bend resolution is specified in StubCuts.BendResolution.
TMTrackProducer.HTFillingRphi.UseBendFilter = cms.bool(True)

#--- Rules for deciding when the Hough Transform has found an L1 track candidate

# Define layers using layer ID (true) or by bins in radius of 5 cm width (false).
TMTrackProducer.L1TrackDef.UseLayerID = cms.bool(False)

#--- Disable r-z Hough transform, disable r-z track filters, disable duplicate track removal & disable track fit.

TMTrackProducer.HTArraySpecRz.EnableRzHT   = cms.bool(False)  
TMTrackProducer.RZfilterOpts.UseSeedFilter = cms.bool(False)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi  = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRz    = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg = cms.uint32(0)
TMTrackProducer.DupTrkRemoval.DupTrkAlgFit   = cms.uint32(0)
TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring()
