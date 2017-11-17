// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TMTrackTrigger/TMTrackFinder/interface/InputData.h>
#include <TMTrackTrigger/TMTrackFinder/interface/Settings.h>   
#include <TMTrackTrigger/TMTrackFinder/interface/Sector.h>

#include <TMTrackTrigger/TMTrackFinder/interface/HTpair.h>
#include <TMTrackTrigger/TMTrackFinder/interface/KillDupFitTrks.h>
#include <TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h>
#include <TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h>
#include <TMTrackTrigger/TMTrackFinder/interface/L1fittedTrk4and5.h>
#include <TMTrackTrigger/TMTrackFinder/interface/ConverterToTTTrack.h>
#include "TMTrackTrigger/TMTrackFinder/interface/HTcell.h"
#include "TMTrackTrigger/TMTrackFinder/interface/MuxHToutputs.h"
#include "TMTrackTrigger/TMTrackFinder/interface/VertexFinder.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/HTrphi.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "boost/numeric/ublas/matrix.hpp"

#include "DataFormats/L1Trigger/interface/BXVector.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <set>

#include <algorithm>
#include <array>
#include <unordered_set>
#include "fstream"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

class Settings;
class TrackFitGeneric;
class L1fittedTrack;


using namespace std;
using  boost::numeric::ublas::matrix;


class TMTrackAnalyzer : public edm::EDAnalyzer {

 public:
  explicit TMTrackAnalyzer(const edm::ParameterSet&);
  ~TMTrackAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;
  map<const TP*, string> diagnoseTracking(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs) const;

  //declare the tokens                                                             
  edm::EDGetToken m_tpInput; //check TMTrackProduction                             
  edm::EDGetToken m_stubInput;
  edm::EDGetToken m_stubTruthInput;
  edm::EDGetToken m_clusterTruthInput;

  bool m_doTpInput;
  bool m_doStubInput;
  bool m_doStubTruthInput;
  bool m_doClusterTruthInput;

  bool doText_;
  bool doHistos_;
  
  //declare the folders                                                            

    
  enum ObjectType{
    InputDt=1,
    CheckSectors=2,
    HTrphiHis=3,
    RZfilters=4,
    BusyEvents=5,
    TrackCands=6,
    SimpleLR=7,
  };

  vector< ObjectType > types_;
  vector< std::string > typeStr_;
    
    
    // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted.
    map<std::string, unsigned int> numFitAlgEff_;
    map<std::string, unsigned int> numFitPerfAlgEff_;
    
    // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted post-cut.
    map<std::string, unsigned int> numFitAlgEffPass_;
    map<std::string, unsigned int> numFitPerfAlgEffPass_;
    
   // float settings_->chi2OverNdfCut();

  //declare all the histogram we use

  map< ObjectType, TFileDirectory> dirs_;

  // 1. InputDt
  map< ObjectType, TProfile*> profNumStubs_;
  map< ObjectType, TH1F*> hisStubsVsEta_;
  map< ObjectType, TH1F*> hisStubsVsR_;
  map< ObjectType, TH2F*> hisStubsVsRVsZ_;
  map< ObjectType, TH2F*> hisStubsModuleVsRVsZ_;
  map< ObjectType, TH2F*> hisStubsVsRVsPhi_;
  map< ObjectType, TH2F*> hisStubsModuleVsRVsPhi_;
  map< ObjectType, TH2F*> hisStubsModuleTiltVsZ_;
  map< ObjectType, TH2F*> hisStubsdPhiCorrectionVsZ_;

  map< ObjectType, TH2F*> hisStubsVsRVsZ_outerModuleAtSmallerR_;
  map< ObjectType, TH2F*> hisStubsVsRVsPhi_outerModuleAtSmallerR_;

  map< ObjectType, TProfile*> profNumTPs_;
  map< ObjectType, TH1F*> hisNumStubsPerTP_;
  map< ObjectType, TH1F*> hisNumPSStubsPerTP_;
  map< ObjectType, TH1F*> hisNum2SStubsPerTP_;

            // Study efficiency of tightened front end-electronics cuts.

  map< ObjectType, TProfile*> hisStubKillFE_;
  map< ObjectType, TProfile*> hisStubIneffiVsInvPt_;
  map< ObjectType, TProfile*> hisStubIneffiVsEta_;
  map< ObjectType, TProfile*> hisStubKillDataCorr_;

            // Study stub resolution.

  map< ObjectType, TH1F*> hisPtStub_;
  map< ObjectType, TH1F*> hisPtResStub_;
  map< ObjectType, TH1F*> hisBendFilterPower_;
  map< ObjectType, TH1F*> hisDelPhiStub_;
  map< ObjectType, TH1F*> hisDelPhiResStub_;
  map< ObjectType, TH1F*> hisDelPhiResStub_tilted_;
  map< ObjectType, TH1F*> hisDelPhiResStub_notTilted_;

  map< ObjectType, TH1F*> hisBendStub_;
  map< ObjectType, TH1F*> hisBendResStub_;
  map< ObjectType, TH1F*> hisNumMergedBend_;
  map< ObjectType, TH2F*> hisBendVsLayerOrRing_;
  map< ObjectType, TH2F*> hisBendFEVsLayerOrRing_;

  map< ObjectType, TH1F*> hisPhiStubVsPhiTP_;
  map< ObjectType, TH1F*> hisPhiStubVsPhi0TP_;
  map< ObjectType, TH1F*> hisPhi0StubVsPhi0TP_;
  map< ObjectType, TH1F*> hisPhi0StubVsPhi0TPres_;
  map< ObjectType, TH1F*> hisPhiStubVsPhi65TP_;
  map< ObjectType, TH1F*> hisPhi65StubVsPhi65TP_;
  map< ObjectType, TH1F*> hisPhi65StubVsPhi65TPres_;

                 // Note ratio of sensor pitch to separation (needed to understand how many bits this can be packed into).

  map< ObjectType, TH1F*> hisPitchOverSep_;
  map< ObjectType, TH1F*> hisRhoParameter_;

                 // Count stubs sharing a common cluster.

  map< ObjectType, TH1F*> hisFracStubsSharingClus0_;
  map< ObjectType, TH1F*> hisFracStubsSharingClus1_;
  map< ObjectType, TH1F*> hisPileUpPt_;
  map< ObjectType, TH1F*> hisPhysicsPt_;

  // 2. CheckSectors                                                                  
  map< ObjectType, TH1F*> hisFracStubsInSec_;
  map< ObjectType, TH1F*> hisFracStubsInEtaSec_;
  map< ObjectType, TH1F*> hisFracStubsInPhiSec_;
  map< ObjectType, TH1F*> hisNumSecsPerStub_;
  map< ObjectType, TH1F*> hisNumEtaSecsPerStub_;
  map< ObjectType, TH1F*> hisNumPhiSecsPerStub_;
  map< ObjectType, TH1F*> hisNumStubsPerSec_;
  map< ObjectType, TProfile*> profNumStubsPerEtaSec_;
  map< ObjectType, TH2F*> hisLayerIDvsEtaSec_;
  map< ObjectType, TH2F*> hisLayerIDreducedvsEtaSec_;

  // 3. HTrphi                                                                        
  map< ObjectType, TH1F*> hisIncStubsPerHT_;
  map< ObjectType, TH1F*> hisExcStubsPerHT_;
  map< ObjectType, TH2F*> hisNumStubsInCellVsEta_;
  map< ObjectType, TH1F*> hisStubsOnRphiTracksPerHT_;

  // 4.  RZfilters                                                                   

  map< ObjectType, TH1F*> hisNumZtrkSeedCombinations_;
  map< ObjectType, TH1F*> hisNumSeedCombinations_;
  map< ObjectType, TH1F*> hisNumGoodSeedCombinations_;
  map< ObjectType, TH1F*> hisCorrelationZTrk_;

  // 5. BusyEvents                                                                              

  map< ObjectType, TH1F*> hisNumBusySecsInPerEvent_;
  map< ObjectType, TH1F*> hisNumBusySecsOutPerEvent_;
  map< ObjectType, TProfile*> profFracBusyInVsEtaReg_;
  map< ObjectType, TProfile*> profFracBusyOutVsEtaReg_;
  map< ObjectType, TProfile*> profFracStubsKilledVsEtaReg_;
  map< ObjectType, TProfile*> profFracTracksKilledVsEtaReg_; 
  map< ObjectType, TProfile*> profFracTracksKilledVsInvPt_;
  map< ObjectType, TProfile*> profFracTPKilledVsEta_;
  map< ObjectType, TProfile*> profFracTPKilledVsInvPt_;
  map< ObjectType, TH1F*> hisNumTPkilledBusySec_;
   
    map<string, map< ObjectType, TH1F*> > hisNumInputStubs_ ;
    map<string, map< ObjectType, TH1F*> > hisQoverPtInputStubs_;
    map<string, map< ObjectType, TH1F*> > hisNumOutputStubs_;
    map<string, map< ObjectType, TH1F*> > hisNumTracks_;
    map<string, map< ObjectType, TH1F*> > hisNumStubsPerTrack_;
    map<string, map< ObjectType, TH1F*> > hisTrackQoverPt_;
    map<string, map< ObjectType, TH1F*> > hisTrackPurity_;
    map<string, map< ObjectType, TH1F*> > hisNumTPphysics_;
    map<string, map< ObjectType, TH1F*> > hisNumTPpileup_;
    map<string, map< ObjectType, TH1F*> > hisSumPtTPphysics_;
    map<string, map< ObjectType, TH1F*> > hisSumPtTPpileup_;

    // 6. TrackCands

  map< ObjectType, TProfile*> profNumTrackCands_;
  map< ObjectType, TProfile*> profNumTracksVsEta_;
  map< ObjectType, TH1F*> hisNumTracksVsQoverPt_;
  map< ObjectType, TH1F*> hisNumTrksPerSect_;
  map< ObjectType, TH1F*> hisNumTrksPerOct_;
    
            // Count stubs per event assigned to tracks (determines HT data output rate)
    
    map< ObjectType, TProfile*> profStubsOnTracks_;
    map< ObjectType, TProfile*> profStubsOnTracksVsEta_;
    map< ObjectType, TH1F*> hisStubsOnTracksPerSect_;
    map< ObjectType, TH1F*> hisStubsOnTracksPerOct_;
    map< ObjectType, TH1F*> hisUniqueStubsOnTrksPerSect_;
    map< ObjectType, TH1F*> hisUniqueStubsOnTrksPerOct_;
    
    map< ObjectType, TH1F*> hisStubsPerTrack_;
    map< ObjectType, TH1F*> hisLayersPerTrack_;
    map< ObjectType, TH1F*> hisPSLayersPerTrack_;
    map< ObjectType, TH1F*> hisLayersPerTrueTrack_;
    map< ObjectType, TH1F*> hisPSLayersPerTrueTrack_;
    
            // Checks if tracks have too many stubs to be stored in memory in each cell.
 
    map< ObjectType, TProfile*> profExcessStubsPerTrackVsPt_;
    map< ObjectType, TH1F*> hisFracMatchStubsOnTracks_;
    
            // See how far stubs lie from true trajectory in r-z plane.

    
    map< ObjectType, TH1F*> hisDeltaPhiRtruePS_;
    map< ObjectType, TH1F*> hisDeltaRorZtruePS_;
    map< ObjectType, TH1F*> hisDeltaPhiRtrue2S_;
    map< ObjectType, TH1F*> hisDeltaRorZtrue2S_;
    map< ObjectType, TH1F*> hisDeltaPhiRfakePS_;
    map< ObjectType, TH1F*> hisDeltaRorZfakePS_;
    map< ObjectType, TH1F*> hisDeltaPhiRfake2S_;
    map< ObjectType, TH1F*> hisDeltaRorZfake2S_;
    map< ObjectType, TProfile*> profNsigmaPhiRvsInvPt_;
    map< ObjectType, TProfile*> profNsigmaPhiRvsFracDist_;
    map< ObjectType, TProfile*> profFracTrueStubsVsLayer_;
    
    
            // Check how much stub bend differs from predicted one.
 
    map< ObjectType, TH1F*> hisDeltaBendTrue_;
    map< ObjectType, TH1F*> hisDeltaBendFake_;
    
            // Study duplication of tracks within HT.

    map< ObjectType, TProfile*> profDupTracksVsTPeta_;
    
            // Histos for tracking efficiency vs. TP kinematics

    map< ObjectType, TH1F*> hisTPinvptForEff_;
    map< ObjectType, TH1F*> hisRecoTPinvptForEff_;
    map< ObjectType, TH1F*> hisTPetaForEff_;
    map< ObjectType, TH1F*> hisRecoTPetaForEff_;
    map< ObjectType, TH1F*> hisTPphiForEff_;
    map< ObjectType, TH1F*> hisRecoTPphiForEff_;
    
            // Histo for efficiency to reconstruct track perfectly (no incorrect hits).

    map< ObjectType, TH1F*> hisPerfRecoTPinvptForEff_;
    map< ObjectType, TH1F*> hisPerfRecoTPetaForEff_;
    
            // Histos for algorithmic tracking efficiency vs. TP kinematics

    map< ObjectType, TH1F*> hisTPinvptForAlgEff_;
    map< ObjectType, TH1F*> hisRecoTPinvptForAlgEff_;
    map< ObjectType, TH1F*> hisTPetaForAlgEff_;
    map< ObjectType, TH1F*> hisRecoTPetaForAlgEff_;
    map< ObjectType, TH1F*> hisTPphiForAlgEff_;
    map< ObjectType, TH1F*> hisRecoTPphiForAlgEff_;
    
            // Histo for efficiency to reconstruct track perfectly (no incorrect hits).

    map< ObjectType, TH1F*> hisPerfRecoTPinvptForAlgEff_;
    map< ObjectType, TH1F*> hisPerfRecoTPetaForAlgEff_;
    
            // Histos for algorithmic tracking efficiency vs. TP production point

    map< ObjectType, TH1F*> hisTPd0ForAlgEff_;
    map< ObjectType, TH1F*> hisRecoTPd0ForAlgEff_;
    map< ObjectType, TH1F*> hisTPz0ForAlgEff_;
    map< ObjectType, TH1F*> hisRecoTPz0ForAlgEff_;
    
            // Histos for algorithmic tracking efficiency vs sector number (to check if looser cuts are needed in certain regions)

    map< ObjectType, TH1F*> hisTPphisecForAlgEff_;
    map< ObjectType, TH1F*> hisRecoTPphisecForAlgEff_;
    map< ObjectType, TH1F*> hisTPetasecForAlgEff_;
    map< ObjectType, TH1F*> hisRecoTPetasecForAlgEff_;
    
            // Histo for efficiency to reconstruct tracks perfectly (no incorrect hits).

    map< ObjectType, TH1F*> hisPerfRecoTPphisecForAlgEff_;
    map< ObjectType, TH1F*> hisPerfRecoTPetasecForAlgEff_;
    
            // Histos of track parameter resolution

    map< ObjectType, TH1F*> hisQoverPtRes_;
    map< ObjectType, TH1F*> hisPhi0Res_;
    map< ObjectType, TH1F*> hisEtaRes_;
    map< ObjectType, TH1F*> hisZ0Res_;
    map< ObjectType, TH1F*> hisTanLambdaRes_;
    
            // For those tracking particles causing the algorithmic efficiency to be below 100%, plot a flag indicating why.

    map< ObjectType, TH1F*> hisRecoFailureReason_;
    map< ObjectType, TH2F*> hisWrongSignStubRZ_pBend_;
    map< ObjectType, TH2F*> hisWrongSignStubRZ_nBend_;
    map< ObjectType, TH1F*> hisNumStubsOnLayer_;
    
    
    // 7. SimpleLR
    
    map<string, map< ObjectType, TH1F*> > hisSeedQinvPt_;
    map<string, map< ObjectType, TH1F*> > hisSeedPhi0_;
    map<string, map< ObjectType, TH1F*> > hisSeedD0_;
    map<string, map< ObjectType, TH1F*> > hisSeedZ0_;
    map<string, map< ObjectType, TH1F*> > hisSeedEta_;
    
    map<string, map< ObjectType, TProfile*> > profNumFittedCands_;
    
    map<string, map< ObjectType, TH1F*> > hisFitQinvPtMatched_;
    map<string, map< ObjectType, TH1F*> > hisFitPhi0Matched_;
    map<string, map< ObjectType, TH1F*> > hisFitD0Matched_;
    map<string, map< ObjectType, TH1F*> > hisFitZ0Matched_;
    map<string, map< ObjectType, TH1F*> > hisFitEtaMatched_;
    
    map<string, map< ObjectType, TH1F*> > hisFitChi2Matched_;
    map<string, map< ObjectType, TH1F*> > hisFitChi2DofMatched_;
    
    map<string, map< ObjectType, TH1F*> > hisFitQinvPtUnmatched_;
    map<string, map< ObjectType, TH1F*> > hisFitPhi0Unmatched_;
    map<string, map< ObjectType, TH1F*> > hisFitD0Unmatched_;
    map<string, map< ObjectType, TH1F*> > hisFitZ0Unmatched_;
    map<string, map< ObjectType, TH1F*> > hisFitEtaUnmatched_;
    
    map<string, map< ObjectType, TH1F*> > hisFitChi2Unmatched_;
    map<string, map< ObjectType, TH1F*> > hisFitChi2DofUnmatched_;
    
    map<string, map< ObjectType, TH2F*> > hisFitVsTrueQinvPtGoodChi2_;
    map<string, map< ObjectType, TH2F*> > hisFitVsTruePhi0GoodChi2_;
    map<string, map< ObjectType, TH2F*> > hisFitVsTrueD0GoodChi2_;
    map<string, map< ObjectType, TH2F*> > hisFitVsTrueZ0GoodChi2_;
    map<string, map< ObjectType, TH2F*> > hisFitVsTrueEtaGoodChi2_;
    
    
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedQinvPtGenCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedPhi0GenCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedD0GenCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedZ0GenCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedEtaGenCand_;
    
    
    map<string, map< ObjectType, TH1F*> > hisFitQinvPtResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisFitPhi0ResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisFitD0ResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisFitZ0ResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisFitEtaResGoodChi2_;
    
    map<string, map< ObjectType, TH1F*> > hisSeedQinvPtResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisSeedPhi0ResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisSeedD0ResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisSeedZ0ResGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisSeedEtaResGoodChi2_;
    
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedQinvPtFakeCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedPhi0FakeCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedD0FakeCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedZ0FakeCand_;
    map<string, map< ObjectType, TH2F*> > hisFitVsSeedEtaFakeCand_;
    
    map<string, map< ObjectType, TProfile*> > hisQoverPtResVsTrueEta_;
    map<string, map< ObjectType, TProfile*> > hisPhi0ResVsTrueEta_;
    map<string, map< ObjectType, TProfile*> > hisEtaResVsTrueEta_;
    map<string, map< ObjectType, TProfile*> > hisZ0ResVsTrueEta_;
    map<string, map< ObjectType, TProfile*> > hisD0ResVsTrueEta_;
    
    map<string, map< ObjectType, TProfile*> > hisQoverPtResVsTrueZ0_;
    map<string, map< ObjectType, TProfile*> > hisPhi0ResVsTrueZ0_;
    map<string, map< ObjectType, TProfile*> > hisEtaResVsTrueZ0_;
    map<string, map< ObjectType, TProfile*> > hisZ0ResVsTrueZ0_;
    map<string, map< ObjectType, TProfile*> > hisD0ResVsTrueZ0_;
    
    map<string, map< ObjectType, TProfile*> > hisQoverPtResVsTrueInvPt_;
    map<string, map< ObjectType, TProfile*> > hisPhi0ResVsTrueInvPt_;
    map<string, map< ObjectType, TProfile*> > hisEtaResVsTrueInvPt_;
    map<string, map< ObjectType, TProfile*> > hisZ0ResVsTrueInvPt_;
    map<string, map< ObjectType, TProfile*> > hisD0ResVsTrueInvPt_;
    
    map<string, map< ObjectType, TH2F*> > hisTrueFittedChiSquaredVsTrueEta_;
    map<string, map< ObjectType, TH2F*> > hisTrueFittedChiSquaredDofVsTrueEta_;
    map<string, map< ObjectType, TH2F*> > hisTrueFittedChiSquaredVsFittedEta_;
    map<string, map< ObjectType, TH2F*> > hisTrueFittedChiSquaredDofVsFittedEta_;
    map<string, map< ObjectType, TH2F*> > hisFittedChiSquaredFunctionOfStubs_;
    map<string, map< ObjectType, TH2F*> > hisFittedChiSquaredDofFunctionOfStubs_;
    
    map<string, map< ObjectType, TH1F*> > hisTrueEtaMatchedGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisTrueEtaMatchedBadChi2_;
    map<string, map< ObjectType, TH1F*> > hisStubPurityMatchedGoodChi2_;
    map<string, map< ObjectType, TH1F*> > hisStubPurityMatchedBadChi2_;
    
    map<string, map< ObjectType, TProfile*> > profChi2DofVsInvPtPERF_;
    map<string, map< ObjectType, TProfile*> > profBigChi2DofVsInvPtPERF_;
    map<string, map< ObjectType, TH1F*> > hisD0TPBigChi2DofPERF_;
    map<string, map< ObjectType, TH1F*> > hisD0TPSmallChi2DofPERF_;
    
    map<string, map< ObjectType, TH2F*> > hisNumMatchedStubsKilledVsKilled_;
    map<string, map< ObjectType, TProfile*> > profTrksKilledByFit_;
    map<string, map< ObjectType, TH2F*> > hisNumStubsVsPurity_;
    
    map<string, map< ObjectType, TH1F*> > hisNumIterations_;
    map<string, map< ObjectType, TH1F*> > hisFailingState_;
    map<string, map< ObjectType, TProfile*> > hisTotalStateCalls_;
    map<string, map< ObjectType, TProfile*> > hisRelativeStateCalls_;
    
    map<string, map< ObjectType, TH1F*> > hisNumStubsFitKills_;
    map<string, map< ObjectType, TH2F*> > hisNumStubsFitKillsVsPurity_;
    map<string, map< ObjectType, TH2F*> > hisNumStubsFitKillsVsPurityMatched_;
    map<string, map< ObjectType, TH2F*> > hisNumStubsFitKillsVsPurityUnmatched_;
    
    map<string, map< ObjectType, TH2F*> > hisFitEfficiencyVsChi2Dof_;
    map<string, map< ObjectType, TH2F*> > hisNumStubsVsChi2Dof_;
    map<string, map< ObjectType, TH2F*> > hisNumLayersVsChi2Dof_;
    map<string, map< ObjectType, TH2F*> > hisAvgNumStubsPerLayerVsChi2Dof_;
    
    map<string, map< ObjectType, TProfile*> > profDupFitTrksVsEta_;
    map<string, map< ObjectType, TProfile*> > profDupFitTrksVsInvPt_;
    map<string, map< ObjectType, TProfile*> > profFakeFitTrksVsEta_;
    map<string, map< ObjectType, TProfile*> > profFakeFitTrksVsZ0_;
    map<string, map< ObjectType, TProfile*> > profFitFracTrueStubsVsLayer_;
    map<string, map< ObjectType, TProfile*> > profFitFracTrueStubsVsEta_;
    
    map<string, map< ObjectType, TH1F*> > hisFitTPinvptForEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPetaForEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPphiForEff_;
    map<string, map< ObjectType, TH1F*> > hisPerfFitTPinvptForEff_;
    map<string, map< ObjectType, TH1F*> > hisPerfFitTPetaForEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPinvptForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPetaForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPphiForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisPerfFitTPinvptForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisPerfFitTPetaForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPd0ForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPz0ForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPphisecForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisFitTPetasecForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisPerfFitTPphisecForAlgEff_;
    map<string, map< ObjectType, TH1F*> > hisPerfFitTPetasecForAlgEff_;

    
  int m_bx = 0;
  bool m_allBx = false;
  Settings *settings_;
  map<string, TrackFitGeneric*> fitterWorkerMap_;
    unsigned int numPerfRecoTPforAlg_;
    
    typedef std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > TTTrackCollection;
    typedef std::vector< const L1fittedTrack* > FitTrackCollection;

  const edm::EDGetTokenT<TrackingParticleCollection> tpInputTag;
  const edm::EDGetTokenT<DetSetVec> stubInputTag;
  const edm::EDGetTokenT<TTStubAssMap> stubTruthInputTag;
  const edm::EDGetTokenT<TTClusterAssMap> clusterTruthInputTag;


};


//constructor
TMTrackAnalyzer::TMTrackAnalyzer(const edm::ParameterSet& iConfig):
  tpInputTag( consumes<TrackingParticleCollection>( iConfig.getParameter<edm::InputTag>("tpInputTag") ) ),
  stubInputTag( consumes<DetSetVec>( iConfig.getParameter<edm::InputTag>("stubInputTag") ) ),
  stubTruthInputTag( consumes<TTStubAssMap>( iConfig.getParameter<edm::InputTag>("stubTruthInputTag") ) ),
  clusterTruthInputTag( consumes<TTClusterAssMap>( iConfig.getParameter<edm::InputTag>("clusterTruthInputTag") ) )
{
    //add the rest...
    settings_ = new Settings(iConfig);
    //now do what ever initialization is needed
      
        
    m_bx  = iConfig.getParameter<int>("selBx");
    m_allBx = iConfig.getParameter<bool>("selAllBx");
        
    // register what you consume and keep token for later access:

    edm::InputTag nullTag("None");
    
    // Tame debug printout.
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(4);
    
    // Create track fitting algorithm (& internal histograms if it uses them)
    for (const string& fitterName : settings_->trackFitters()) {
        fitterWorkerMap_[ fitterName ] = TrackFitGeneric::create(fitterName, settings_);
        fitterWorkerMap_[ fitterName ]->bookHists();
    }
    

    /*    
    edm::InputTag tpInputTag = iConfig.getParameter<edm::InputTag>("tpInputTag");
    m_tpInput                = consumes<TrackingParticleCollection>(tpInputTag);
    m_doTpInput              = !(tpInputTag==nullTag);
        
    edm::InputTag stubInputTag = iConfig.getParameter<edm::InputTag>("stubInputTag");
    m_stubInput                = consumes<DetSetVec>(stubInputTag);
    m_doStubInput              = !(stubInputTag==nullTag);
        
    edm::InputTag stubTruthInputTag = iConfig.getParameter<edm::InputTag>("stubTruthInputTag");
    m_stubTruthInput                = consumes<TTStubAssMap>(stubTruthInputTag);
    m_doStubTruthInput              = !(stubTruthInputTag==nullTag);
        
    edm::InputTag clusterTruthInputTag = iConfig.getParameter<edm::InputTag>("clusterTruthInputTag");
    m_clusterTruthInput                = consumes<TTClusterAssMap>(clusterTruthInputTag);
    m_doClusterTruthInput              = !(clusterTruthInputTag==nullTag);
    */
        
}
    
    
//destructor
TMTrackAnalyzer::~TMTrackAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
        
}   
    
    
    ///////// BEGIN JOB ////////
void TMTrackAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){

  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  float bField = theMagneticField->inTesla(GlobalPoint(0,0,0)).z(); // B field in Tesla.
  cout<<endl<<"--- B field = "<<bField<<" Tesla ---"<<endl<<endl;

  settings_->setBfield(bField);
    
    // Initialize track fitting algorithm at start of run (especially with B-field dependent variables).
    for (const string& fitterName : settings_->trackFitters()) {
        fitterWorkerMap_[ fitterName ]->initRun();
    }
}

void TMTrackAnalyzer::beginJob(){

  edm::Service<TFileService> fs;

  unsigned int nEta = settings_->numEtaRegions();
  unsigned int nPhi = settings_->numPhiSectors();
  float maxAbsQoverPt = 1./settings_->houghMinPt();
  
  //1. InputData Histograms
        
  dirs_.insert(std::pair< ObjectType, TFileDirectory >(InputDt, fs->mkdir("InputData") ) );

  // Count stubs & tracking particles.

  profNumStubs_.insert( pair< ObjectType, TProfile* > (InputDt, dirs_.at(InputDt).make<TProfile>("NumStubs","; Category; No. stubs in tracker",4,0.5,4.5)));
  profNumStubs_.at(InputDt)->GetXaxis()->SetBinLabel(1,"All stubs");
  profNumStubs_.at(InputDt)->GetXaxis()->SetBinLabel(2,"Genuine stubs");
  profNumStubs_.at(InputDt)->GetXaxis()->SetBinLabel(3,"Stubs matched to TP");
  profNumStubs_.at(InputDt)->GetXaxis()->SetBinLabel(4,"Stubs matched to TP for eff");
  profNumStubs_.at(InputDt)->LabelsOption("d"); ////check here !!!!!!
    
  hisStubsVsEta_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("StubsVsEta","; #eta; No. stubs in tracker",30,-3.0,3.0)));
    
  hisStubsVsR_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("StubsVsR","; radius (cm); No. stubs in tracker",1200,0.,120.)));
    
  hisStubsVsRVsZ_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsVsRVsZ","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130)));
    
  hisStubsModuleVsRVsZ_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsModuleVsRVsZ","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130)));
    
    hisStubsVsRVsPhi_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsVsRVsPhi","; x (cm); y (cm); No. stubs in tracker",1000,-130,130,1000,-130,130)));
    
    hisStubsModuleVsRVsPhi_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsModuleVsRVsPhi","; x (cm); y (cm); No. stubs in tracker",1000,-130,130,1000,-130,130)));
    
    hisStubsModuleTiltVsZ_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsModuleTiltVsZ","; z (cm); Tilt; Module tilt vs z",1000,-280,280,1000,-3,3)));
    
    hisStubsdPhiCorrectionVsZ_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsdPhiCorrectionVsZ","; z (cm); Correction; dPhi Correction vs z",1000,-280,280,1000,-1,10)));
    
    
    hisStubsVsRVsZ_outerModuleAtSmallerR_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsVsRVsZ_outerModuleAtSmallerR","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130)));
    
    hisStubsVsRVsPhi_outerModuleAtSmallerR_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("StubsVsRVsPhi_outerModuleAtSmallerR","; x (cm); y (cm); No. stubs in tracker",1000,-130,130,1000,-130,130)));
    
    profNumTPs_.insert( pair< ObjectType, TProfile* > (InputDt, dirs_.at(InputDt).make<TProfile>("NumTPs","; Category; No. of TPs in tracker",3,0.5,3.5)));
    profNumTPs_.at(InputDt)->GetXaxis()->SetBinLabel(1,"All TPs");
    profNumTPs_.at(InputDt)->GetXaxis()->SetBinLabel(2,"TPs for eff.");
    profNumTPs_.at(InputDt)->GetXaxis()->SetBinLabel(3,"TPs for alg. eff.");
    profNumTPs_.at(InputDt)->LabelsOption("d");
    
    
    hisNumStubsPerTP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("NumStubsPerTP","; Number of stubs per TP for alg. eff.",50,-0.5,49.5)));
    
    hisNumPSStubsPerTP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("NumPSStubsPerTP","; Number of PS stubs per TP for alg. eff.",50,-0.5,49.5)));
    
    hisNum2SStubsPerTP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("Num2SStubsPerTP","; Number of 2S stubs per TP for alg. eff.",50,-0.5,49.5)));
    
            // Study efficiency of tightened front end-electronics cuts.

    hisStubKillFE_.insert( pair< ObjectType, TProfile* > (InputDt, dirs_.at(InputDt).make<TProfile>("StubKillFE","; barrelLayer or 10+endcapRing; Stub fraction rejected by readout chip",30,-0.5,29.5)));
    
    hisStubIneffiVsInvPt_.insert( pair< ObjectType, TProfile* > (InputDt, dirs_.at(InputDt).make<TProfile>("StubIneffiVsPt","; 1/Pt; Inefficiency of readout chip for good stubs",30,0.0,1.0)));
    
    hisStubIneffiVsEta_.insert( pair< ObjectType, TProfile* > (InputDt, dirs_.at(InputDt).make<TProfile>("StubIneffiVsEta","; |#eta|; Inefficiency of readout chip for good stubs",30,0.0,3.0)));
    
    hisStubKillDataCorr_.insert( pair< ObjectType, TProfile* > (InputDt, dirs_.at(InputDt).make<TProfile>("StubKillDataCorr","; barrelLayer or 10+endcapRing; Stub fraction killed by DataCorrection.h window cut",30,-0.5,29.5)));
    
            // Study stub resolution.

    
    hisPtStub_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("PtStub","; Stub q/Pt",50,-0.5,0.5)));
    
    hisPtResStub_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("PtResStub","; Stub q/Pt minus TP q/Pt",50,-0.5,0.5)));
    
    hisBendFilterPower_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("BendFilterPower","; Fraction of q/Pt range allowed",102,-0.01,1.01)));
    
    hisDelPhiStub_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("DelPhiStub","; Stub bend angle",50,-0.2,0.2)));
    
    hisDelPhiResStub_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("DelPhiResStub","; Stub bend angle minus TP bend angle",200,-0.2,0.2)));
    
    hisDelPhiResStub_tilted_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("DelPhiResStub_tilted","; Stub bend angle minus TP bend angle",200,-0.2,0.2)));
    
    hisDelPhiResStub_notTilted_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("DelPhiResStub_notTilted","; Stub bend angle minus TP bend angle",200,-0.2,0.2)));
    
    
    hisBendStub_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("BendStub","; Stub bend in units of strips",57,-7.125,7.125)));
    
    hisBendResStub_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("BendResStub","; Stub bend minus TP bend in units of strips",100,-5.,5.)));
    
    hisNumMergedBend_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("NumMergedBend","; No. of bend values merged together by loss of bit",10,-0.5,9.5)));
    
    hisBendVsLayerOrRing_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("BendVsLayerOrRing","; barrelLayer or 10+endcapRing; Stub bend",30,-0.5,29.5,57,-7.125,7.125)));
    
    hisBendFEVsLayerOrRing_.insert( pair< ObjectType, TH2F* > (InputDt, dirs_.at(InputDt).make<TH2F>("BendFEVsLayerOrRing","; barrelLayer or 10+endcapRing; Stub bend in FE chip",30,-0.5,29.5,57,-7.125,7.125)));
    
    
    hisPhiStubVsPhiTP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("PhiStubVsPhiTP","; Stub #phi minus TP #phi at stub radius",100,-0.05,0.05)));
    
    hisPhiStubVsPhi0TP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("PhiStubVsPhi0TP","; Stub #phi minus TP #phi0",100,-0.3,0.3)));
    
    hisPhi0StubVsPhi0TP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("Phi0StubVsPhi0TP","; #phi0 of Stub minus TP",100,-0.2,0.2)));
    
    hisPhi0StubVsPhi0TPres_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("Phi0StubVsPhi0TPres","; #phi0 of Stub minus TP / resolution",100,-5.0,5.0)));
    
    hisPhiStubVsPhi65TP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("PhiStubVsPhi65TP","; Stub #phi minus TP phitrk65",100,-0.2,0.2)));
    
    hisPhi65StubVsPhi65TP_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("Phi65StubVsPhi65TP","; phitrk65 of Stub minus TP",100,-0.2,0.2)));
    
    hisPhi65StubVsPhi65TPres_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("Phi65StubVsPhi65TPres","; phitrk65 of Stub minus TP / resolution",100,-5.0,5.0)));
    
    
    hisPitchOverSep_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("PitchOverSep","; ratio of sensor pitch / separation",100,0.0,0.1)));
    
    hisRhoParameter_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("RhoParameter","; rho parameter",100,0.0,0.2)));
    
    hisFracStubsSharingClus0_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("FracStubsSharingClus0","Fraction of stubs sharing cluster in seed sensor",102,-0.01,1.01)));
    
    hisFracStubsSharingClus1_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("FracStubsSharingClus1","Fraction of stubs sharing cluster in correlation sensor",102,-0.01,1.01)));
    
    hisPileUpPt_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("hisPileUpPt_","Pile-Up tracks transverse momentum [GeV/c]; no. Tracks; p_{T} [GeV/c]", 100, 0, 100)));
    
    hisPhysicsPt_.insert( pair< ObjectType, TH1F* > (InputDt, dirs_.at(InputDt).make<TH1F>("hisPhysicsPt_","Physics Collision tracks transverse momentum [GeV/c]; no. Tracks; p_{T} [GeV/c]", 100, 0, 100)));
    
            // Study stub resolution.


  //2. CheckSectors Histograms
    
  dirs_.insert( pair< ObjectType, TFileDirectory >(CheckSectors, fs->mkdir("CheckSectors") ) );    
    
  hisFracStubsInSec_.insert( pair< ObjectType, TH1F* > (CheckSectors, dirs_.at(CheckSectors).make<TH1F>("FracStubsInSec","; Fraction of stubs on TP in best (#eta,#phi) sector;",102,-0.01,1.01)));
    
  hisFracStubsInEtaSec_.insert( pair< ObjectType, TH1F* > (CheckSectors, dirs_.at(CheckSectors).make<TH1F>("FracStubsInEtaSec","; Fraction of stubs on TP in best #eta sector;",102,-0.01,1.01)));
    
  hisFracStubsInPhiSec_.insert( pair< ObjectType, TH1F* > (CheckSectors, dirs_.at(CheckSectors).make<TH1F>("FracStubsInPhiSec","; Fraction of stubs on TP in best #phi sector;",102,-0.01,1.01)));
    
  hisNumSecsPerStub_.insert( pair< ObjectType, TH1F* > (CheckSectors, dirs_.at(CheckSectors).make<TH1F>("NumSecPerStub","; Number of (#eta,#phi) sectors each stub appears in",20,-0.5,19.5)));
    
  hisNumEtaSecsPerStub_.insert( pair< ObjectType, TH1F* > (CheckSectors, dirs_.at(CheckSectors).make<TH1F>("NumEtaSecPerStub","; Number of #eta sectors each stub appears in",20,-0.5,19.5)));

  hisNumPhiSecsPerStub_.insert( pair< ObjectType, TH1F* > (CheckSectors, dirs_.at(CheckSectors).make<TH1F>("NumPhiSecPerStub","; Number of #phi sectors each stub appears in",20,-0.5,19.5)));

  hisNumStubsPerSec_.insert( pair< ObjectType, TH1F* > (CheckSectors, dirs_.at(CheckSectors).make<TH1F>("NumStubsPerSec","; Number of stubs per sector",250,-0.5,249.5)));

  profNumStubsPerEtaSec_.insert( pair< ObjectType, TProfile* > (CheckSectors, dirs_.at(CheckSectors).make<TProfile>("NumStubsPerEtaSec",";#eta sector; Number of stubs per #eta sector",nEta,-0.5,nEta-0.5)));

  hisLayerIDvsEtaSec_.insert( pair< ObjectType, TH2F* > (CheckSectors, dirs_.at(CheckSectors).make<TH2F>("LayerIDvsEtaSec",";#eta sector; layer ID",nEta,-0.5,nEta-0.5,20,0.5,20.5)));

  hisLayerIDreducedvsEtaSec_.insert( pair< ObjectType, TH2F* > (CheckSectors, dirs_.at(CheckSectors).make<TH2F>("LayerIDreducedvsEtaSec",";#eta sector; reduced layer ID",nEta,-0.5,nEta-0.5,20,0.5,20.5)));

  //3. HTrphi Histograms
    
    dirs_.insert( std::pair< ObjectType, TFileDirectory >(HTrphiHis, fs->mkdir("HTrphi") ) );
    
  hisIncStubsPerHT_.insert( pair< ObjectType, TH1F* > (HTrphiHis, dirs_.at(HTrphiHis).make<TH1F>("IncStubsPerHT","; Number of filtered stubs per r#phi HT array (inc. duplicates)",100,0.,-1.)));

  hisExcStubsPerHT_.insert( pair< ObjectType, TH1F* > (HTrphiHis, dirs_.at(HTrphiHis).make<TH1F>("ExcStubsPerHT","; Number of filtered stubs per r#phi HT array (exc. duplicates)",250,-0.5,249.5)));
    
  hisNumStubsInCellVsEta_.insert( pair< ObjectType, TH2F* > (HTrphiHis, dirs_.at(HTrphiHis).make<TH2F>("NumStubsInCellVsEta","; no. of stubs per HT cell summed over phi sector; #eta region",100,-0.5,499.5, settings_->numEtaRegions(), -0.5, settings_->numEtaRegions() - 0.5)));

  hisStubsOnRphiTracksPerHT_.insert( pair< ObjectType, TH1F* > (HTrphiHis, dirs_.at(HTrphiHis).make<TH1F>("StubsOnRphiTracksPerHT","; Number of stubs assigned to tracks per r#phi HT array",500,-0.5,499.5)));


  //4. RZfilters Histograms
    
    if (settings_->useZTrkFilter() || settings_->useSeedFilter()) {
        
      dirs_.insert( pair< ObjectType, TFileDirectory >(RZfilters, fs->mkdir("RZfilters") ) );
            
      if (settings_->useZTrkFilter()) {
                
        hisNumZtrkSeedCombinations_.insert( pair< ObjectType, TH1F* > (RZfilters, dirs_.at(RZfilters).make<TH1F>("NumZtrkSeedCombinations_","; Number of Ztrk seed combinations per track cand; no. seeds ; ", 50, -0.5 , 49.5)));
      }
            
      if (settings_->useSeedFilter()) {
        hisNumSeedCombinations_.insert( pair< ObjectType, TH1F* > (RZfilters, dirs_.at(RZfilters).make<TH1F>("NumSeedCombinations_","; Number of seed combinations per track cand; no. seeds ; ", 50, -0.5 , 49.5)));
        hisNumGoodSeedCombinations_.insert( pair< ObjectType, TH1F* > (RZfilters, dirs_.at(RZfilters).make<TH1F>("NumGoodSeedCombinations_","; Number of good seed combinations per track cand; ", 30, -0.5 , 29.5)));
      }
            
      if (settings_->useZTrkFilter()){
        hisCorrelationZTrk_.insert( pair< ObjectType, TH1F* > (RZfilters, dirs_.at(RZfilters).make<TH1F>("CorrelationZTrk_","; Correlation factor r between stub in layer 1 and stubs in other layers; r ; ", 50,-1,1)));
      }
    }
    
  //5. BusyEvents Histograms
        
  // Look at (eta, phi) sectors with too many input stubs or too many output (= assigned to tracks) stubs.
        
  dirs_.insert( pair< ObjectType, TFileDirectory >(BusyEvents, fs->mkdir("BusyEvents") ) );
        
  hisNumBusySecsInPerEvent_.insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>("NumBusySecsInPerEvent" ,"; No. sectors with too many input stubs/event" , 20, -0.5, 19.5)));
        
  hisNumBusySecsOutPerEvent_.insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>("NumBusySecsOutPerEvent","; No. sectors with too many output stubs/event", 20, -0.5, 19.5)));
        
  profFracBusyInVsEtaReg_.insert( pair< ObjectType, TProfile* > (BusyEvents, dirs_.at(BusyEvents).make<TProfile>("FracBusyInVsEtaReg" ,"; #eta region; Frac. of sectors with too many input stubs" , nEta, -0.5, nEta-0.5)));
        
  profFracBusyOutVsEtaReg_.insert( pair< ObjectType, TProfile* > (BusyEvents, dirs_.at(BusyEvents).make<TProfile>("FracBusyOutVsEtaReg","; #eta region; Frac. of sectors with too many output stubs", nEta, -0.5, nEta-0.5)));
        
  profFracStubsKilledVsEtaReg_.insert( pair< ObjectType, TProfile* > (BusyEvents, dirs_.at(BusyEvents).make<TProfile>("FracStubsKilledInVsEtaReg" ,"; #eta region; Frac. of input stubs killed" , nEta, -0.5, nEta-0.5)));

  profFracTracksKilledVsEtaReg_.insert( pair< ObjectType, TProfile* > (BusyEvents, dirs_.at(BusyEvents).make<TProfile>("FracTracksKilledInVsEtaReg" ,"; #eta region; Frac. of track killed" , nEta, -0.5, nEta-0.5)));

  profFracTracksKilledVsInvPt_.insert( pair< ObjectType, TProfile* > (BusyEvents, dirs_.at(BusyEvents).make<TProfile>("FracTracksKilledInVsInvPt" ,";1/Pt; Frac. of track killed" , 16, 0.,  1./settings_->houghMinPt())));

  profFracTPKilledVsEta_.insert( pair< ObjectType, TProfile* > (BusyEvents, dirs_.at(BusyEvents).make<TProfile>("FracTPKilledInVsEta" ,";#eta; Efficiency loss due to busy sectors" , 16, 0.,  settings_->maxStubEta())));

  profFracTPKilledVsInvPt_.insert( pair< ObjectType, TProfile* > (BusyEvents, dirs_.at(BusyEvents).make<TProfile>("FracTPKilledInVsInvPt" ,";1/Pt; Efficiency loss due to busy sectors" , 16, 0.,  1./settings_->houghMinPt())));

  hisNumTPkilledBusySec_.insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>("NumTPkilledBusySec","; No. of TP killed in each busy sector",30,-0.5,29.5)));
    
    
    const vector<string> tnames = {"BusyOutSec", "QuietOutSec"};
    const vector<string> enames = {" in busy output sector", " in quiet output sector"};
    
    for (unsigned int i = 0; i <= 1; i++) {
        const string tn = tnames[i];
        const string en = enames[i];
        
        hisNumInputStubs_[tn].insert( pair<ObjectType, TH1F*> (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("NumInputStubs"+(tn)).c_str(), ("; No. input stubs"+(en)).c_str(),   250, -0.5, 249.5)));
            
        hisQoverPtInputStubs_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("QoverPtInputStubs"+(tn)).c_str(), ("; q/Pt of input stubs"+(en)).c_str(),   30, 0., 1./settings_->houghMinPt())));
        
        hisNumOutputStubs_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("NumOutputStubs"+(tn)).c_str(),   ("; No. output stubs"+(en)).c_str(), 1000, -0.5, 999.5)));
        
        hisNumTracks_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("NumTracks"+(tn)).c_str(),         ("; No. tracks"+(en)).c_str(),        200, -0.5, 199.5)));
        
        hisNumStubsPerTrack_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("NumStubsPerTrack"+(tn)).c_str(),  ("; No. stubs/track"+(en)).c_str(),    50, -0.5, 49.5)));
        
        hisTrackQoverPt_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("TrackQoverPt"+(tn)).c_str(),      ("; Track q/pt"+(en)).c_str(),      30, 0., 1./settings_->houghMinPt())));
        
        hisTrackPurity_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("TrackPurity"+(tn)).c_str(),       ("; Track purity"+(en)).c_str(),      102, -0.01, 1.01)));
        
        hisNumTPphysics_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("NumTPphysics"+(tn)).c_str(),      ("; No. physics TP"+(en)).c_str(),     30, -0.5, 29.5)));
        
        hisNumTPpileup_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("NumTPpileup"+(tn)).c_str(),       ("; No. pileup TP"+(en)).c_str(),      30, -0.5, 29.5)));
        
        hisSumPtTPphysics_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("SumPtTPphysics"+(tn)).c_str(),    ("; Sum Pt physics TP"+(en)).c_str(), 100,  0.0, 100.)));
        
        hisSumPtTPpileup_[tn].insert( pair< ObjectType, TH1F* > (BusyEvents, dirs_.at(BusyEvents).make<TH1F>(("SumPtTPpileup"+(tn)).c_str(),     ("; Sum Pt pileup TP"+(en)).c_str(),  100,  0.0, 100.)));
        
    }

    
  //6. TrackCands Histogram
    
  dirs_.insert( pair< ObjectType, TFileDirectory >(TrackCands, fs->mkdir("TrackCands") ) );
    
    profNumTrackCands_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("NumTrackCands","; class; N. of tracks in tracker",8,0.5,8.5)));
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(8,"reco tracks including fakes before rz filter");
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(7,"TP for eff recoed");
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(6,"TP recoed");
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(5,"TP recoed x #eta sector dups");
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(4,"TP recoed x sector dups");
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(3,"TP recoed x sector x r-#phi cell dups");
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(2,"TP recoed x sector x cell dups");
    profNumTrackCands_.at(TrackCands)->GetXaxis()->SetBinLabel(1,"reco tracks including fakes");
    profNumTrackCands_.at(TrackCands)->LabelsOption("d");
                              
        
  profNumTracksVsEta_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("NumTracksVsEta","; #eta region; No. of tracks in tracker", nEta, -0.5, nEta - 0.5)));
                              
        
  hisNumTracksVsQoverPt_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("NumTracksVsQoverPt","; Q/Pt; No. of tracks in tracker",100, -maxAbsQoverPt, maxAbsQoverPt)));
        
  hisNumTrksPerSect_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("NumTrksPerSect","; No. tracks per sector;",100,-0.5,99.5)));
                              
        
  hisNumTrksPerOct_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("NumTrksPerOct", "; No. tracks per octant;",200,-0.5,199.5)));
    
        // Count stubs per event assigned to tracks (determines HT data output rate)
    profStubsOnTracks_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("StubsOnTracks","; ; No. of stubs on tracks per event",1,0.5,1.5)));
    profStubsOnTracksVsEta_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("StubsOnTracksVsEta","; #eta region; No. of stubs on tracks per event", nEta, -0.5, nEta - 0.5)));
    hisStubsOnTracksPerSect_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("StubsOnTracksPerSect","; No. of stubs on tracks per sector", 500,-0.5,499.5)));
    hisStubsOnTracksPerOct_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("StubsOnTracksPerOct","; No. of stubs on tracks per octant", 1000,-0.5,999.5)));
    hisUniqueStubsOnTrksPerSect_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("UniqueStubsOnTrksPerSect","; No. of unique stubs on tracks per sector", 500,-0.5,499.5)));
    hisUniqueStubsOnTrksPerOct_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("UniqueStubsOnTrksPerOct","; No. of unique stubs on tracks per octant", 500,-0.5,499.5)));
    
    hisStubsPerTrack_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("StubsPerTrack",";No. of stubs per track;",50,-0.5,49.5)));
    hisLayersPerTrack_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("LayersPerTrack",";No. of layers with stubs per track;",20,-0.5,19.5)));
    hisPSLayersPerTrack_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PSLayersPerTrack",";No. of PS layers with stubs per track;",20,-0.5,19.5)));
    hisLayersPerTrueTrack_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("LayersPerTrueTrack",";No. of layers with stubs per genuine track;",20,-0.5,19.5)));
    hisPSLayersPerTrueTrack_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PSLayersPerTrueTrack",";No. of PS layers with stubs per genuine track;",20,-0.5,19.5)));
    
        // Checks if tracks have too many stubs to be stored in memory in each cell.

    profExcessStubsPerTrackVsPt_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("ExcessStubsPerTrackVsPt",";q/Pt; Prob. of too many stubs per track",16,0.,maxAbsQoverPt)));
    hisFracMatchStubsOnTracks_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("FracMatchStubsOnTracks","; Fraction of stubs on tracks matching best TP;",101,-0.005,1.005)));
    
        // See how far stubs lie from true trajectory in r-z plane.

    hisDeltaPhiRtruePS_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaPhiRtruePS","PS modules; Dist. of true stubs from true trajectory in r*phi;",100,-0.25,0.25)));
    hisDeltaRorZtruePS_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaRorZtruePS","PS modules; Dist. of true stubs from true trajectory in r-z;",100,-10,10)));
    hisDeltaPhiRtrue2S_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaPhiRtrue2S","2S modules; Dist. of true stubs from true trajectory in r*phi;",100,-0.25,0.25)));
    hisDeltaRorZtrue2S_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaRorZtrue2S","2S modules; Dist. of true stubs from true trajectory in r-z;",100,-10,10)));
    hisDeltaPhiRfakePS_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaPhiRfakePS","PS modules; Dist. of fake stubs from true trajectory in r*phi;",100,-0.25,0.25)));
    hisDeltaRorZfakePS_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaRorZfakePS","PS modules; Dist. of fake stubs from true trajectory in r-z;",100,-10,10)));
    hisDeltaPhiRfake2S_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaPhiRfake2S","2S modules; Dist. of fake stubs from true trajectory in r*phi;",100,-0.25,0.25)));
    hisDeltaRorZfake2S_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaRorZfake2S","2S modules; Dist. of fake stubs from true trajectory in r-z;",100,-10,10)));
    profNsigmaPhiRvsInvPt_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("NsigmaPhiRvsInvPt","; 1/Pt; Num #sigma of true stubs from true trajectory",16,0.,maxAbsQoverPt, 0., 10.)));
    profNsigmaPhiRvsFracDist_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("NsigmaPhiRvsFracDist","; Fractional position in tracker; Num #sigma of true stubs from true trajectory",22,0.,1.1, 0., 10.)));
    profFracTrueStubsVsLayer_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("FracTrueStubsVsLayer",";Layer ID; fraction of true stubs",30,0.5,30.5)));

        // Check how much stub bend differs from predicted one.

    hisDeltaBendTrue_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaBendTrue","True stubs; stub bend minus true bend / resolution;",100,-2.,2.)));
    hisDeltaBendFake_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("DeltaBendFake","Fake stubs; stub bend minus true bend / resolution;",100,-2.,2.)));
    
        // Study duplication of tracks within HT.
    profDupTracksVsTPeta_.insert( pair< ObjectType, TProfile* > (TrackCands, dirs_.at(TrackCands).make<TProfile>("DupTracksVsTPeta" ,"; #eta; Number of duplicate tracks in individual HT array;",30,-3.0,3.0)));
    
        // Histos for tracking efficiency vs. TP kinematics
    
    hisTPinvptForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPinvptForEff" ,"; 1/Pt of TP (used for effi. measurement);",16,0.,maxAbsQoverPt)));
    hisRecoTPinvptForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPinvptForEff" ,"; 1/Pt of TP (used for effi. measurement);",16,0.,maxAbsQoverPt)));
    hisTPetaForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPetaForEff","; #eta of TP (used for effi. measurement);",20,-3.,3.)));
    hisRecoTPetaForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPetaForEff","; #eta of TP (used for effi. measurement);",20,-3.,3.)));
    hisTPphiForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPphiForEff","; #phi of TP (used for effi. measurement);",20,-M_PI,M_PI)));
    hisRecoTPphiForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPphiForEff","; #phi of TP (used for effi. measurement);",20,-M_PI,M_PI)));
    
        // Histo for efficiency to reconstruct track perfectly (no incorrect hits).

    hisPerfRecoTPinvptForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PerfRecoTPinvptForEff" ,"; 1/Pt of TP (used for perf. effi. measurement);",16,0.,maxAbsQoverPt)));
    hisPerfRecoTPetaForEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PerfRecoTPetaForEff","; #eta of TP (used for perf. effi. measurement);",20,-3.,3.)));
    
        // Histos for algorithmic tracking efficiency vs. TP kinematics

    hisTPinvptForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPinvptForAlgEff" ,"; 1/Pt of TP (used for alg. effi. measurement);",16,0.,maxAbsQoverPt)));
    hisRecoTPinvptForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPinvptForAlgEff" ,"; 1/Pt of TP (used for alg. effi. measurement);",16,0.,maxAbsQoverPt)));
    hisTPetaForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPetaForAlgEff","; #eta of TP (used for alg. effi. measurement);",20,-3.,3.)));
    hisRecoTPetaForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPetaForAlgEff","; #eta of TP (used for alg. effi. measurement);",20,-3.,3.)));
    hisTPphiForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPphiForAlgEff","; #phi of TP (used for alg. effi. measurement);",20,-M_PI,M_PI)));
    hisRecoTPphiForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPphiForAlgEff","; #phi of TP (used for alg. effi. measurement);",20,-M_PI,M_PI)));
    
        // Histo for efficiency to reconstruct track perfectly (no incorrect hits).

    hisPerfRecoTPinvptForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PerfRecoTPinvptForAlgEff" ,"; 1/Pt of TP (used for perf. alg. effi. measurement);",16,0.,maxAbsQoverPt)));
    hisPerfRecoTPetaForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PerfRecoTPetaForAlgEff","; #eta of TP (used for perf. alg. effi. measurement);",20,-3.,3.)));
    
        // Histos for algorithmic tracking efficiency vs. TP production point
    hisTPd0ForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPd0ForAlgEff" ,"; d0 of TP (used for alg. effi. measurement);",50,0.,1.)));
    hisRecoTPd0ForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPd0ForAlgEff" ,"; d0 of TP (used for alg. effi. measurement);",50,0.,1.)));
    hisTPz0ForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPz0ForAlgEff" ,"; z0 of TP (used for alg. effi. measurement);",50,0.,25.)));
    hisRecoTPz0ForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPz0ForAlgEff" ,"; z0 of TP (used for alg. effi. measurement);",50,0.,25.)));
    
        // Histos for algorithmic tracking efficiency vs sector number (to check if looser cuts are needed in certain regions)

    hisTPphisecForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPphisecForAlgEff" ,"; #phi sectorof TP (used for alg. effi. measurement);",nPhi,-0.5,nPhi-0.5)));
    hisRecoTPphisecForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPphisecForAlgEff" ,"; #phi sector of TP (used for alg. effi. measurement);",nPhi,-0.5,nPhi-0.5)));
    hisTPetasecForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TPetasecForAlgEff" ,"; #eta sector of TP (used for alg. effi. measurement);",nEta,-0.5,nEta-0.5)));
    hisRecoTPetasecForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoTPetasecForAlgEff" ,"; #eta sector of TP (used for alg. effi. measurement);",nEta,-0.5,nEta-0.5)));
        // Histo for efficiency to reconstruct tracks perfectly (no incorrect hits).

    hisPerfRecoTPphisecForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PerfRecoTPphisecForAlgEff" ,"; #phi sector of TP (used for perf. alg. effi. measurement);",nPhi,-0.5,nPhi-0.5)));
    hisPerfRecoTPetasecForAlgEff_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("PerfRecoTPetasecForAlgEff" ,"; #eta sector of TP (used for perf. alg. effi. measurement);",nEta,-0.5,nEta-0.5)));
    
        // Histos of track parameter resolution

    hisQoverPtRes_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("QoverPtRes","; track resolution in q/Pt", 100,-0.06,0.06)));
    hisPhi0Res_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("Phi0Res",   "; track resolution in #phi0",100,-0.04,0.04)));
    hisEtaRes_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("EtaRes",    "; track resolution in #eta", 100,-1.0,1.0)));
    hisZ0Res_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("Z0Res",     "; track resolution in z0",   100,-10.0,10.0)));
    hisTanLambdaRes_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("TanLambdaRes",    "; track resolution in tan #lambda", 100,-1.0,1.0)));
    
        // For those tracking particles causing the algorithmic efficiency to be below 100%, plot a flag indicating why.

    hisRecoFailureReason_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("RecoFailureReason","; Reason TP (used for alg. effi.) not reconstructed;",1,-0.5,0.5)));
    hisWrongSignStubRZ_pBend_.insert( pair< ObjectType, TH2F* > (TrackCands, dirs_.at(TrackCands).make<TH2F>("RZ of stubs with positive bend, but with wrong sign","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130)));
    hisWrongSignStubRZ_nBend_.insert( pair< ObjectType, TH2F* > (TrackCands, dirs_.at(TrackCands).make<TH2F>("RZ of stubs with negative bend, but with wrong sign","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130)));
    hisNumStubsOnLayer_.insert( pair< ObjectType, TH1F* > (TrackCands, dirs_.at(TrackCands).make<TH1F>("NumStubsOnLayer","; Layer occupancy;",16,1,17)));
    
    
    // 7. SimpleRL Histogram
    for(auto &fitName : settings_->trackFitters() ){
        
        dirs_.insert( pair< ObjectType, TFileDirectory >(SimpleLR, fs->mkdir((fitName)) ) );
        //dirs_.insert( pair< ObjectType, TFileDirectory >([fitName], fs->mkdir((fitName)) ) );
        
        
        
        
        hisSeedQinvPt_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedQinvPt_"+(fitName)).c_str(), "; seed q/p_{T}" , 100, -0.5, 0.5 )));
        hisSeedPhi0_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedPhi0_"+(fitName)).c_str(), "; seed #phi_{0}", 70, -3.5, 3.5 )));
        hisSeedD0_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedD0_"+(fitName)).c_str(), "; seed d_{0}"   , 100, -1., 1. )));
        hisSeedZ0_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedZ0_"+(fitName)).c_str(), "; seed z_{0}"   , 100, -25., 25. )));
        hisSeedEta_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedEta_"+(fitName)).c_str(), "; seed #eta"    , 70, -3.5, 3.5 )));
        
        profNumFittedCands_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("NumFittedCands_"+(fitName)).c_str(), "; class; # of fitted tracks", 11, 0.5, 11.5, -0.5, 9.9e6)));
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(11, "Num tracks exc dups passing cut");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(10, "Num tracks exc dups");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(9, "Num rejected fake tracks");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(8, "Num rejected tracks");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(7, "Num TP tracks");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(6, "Num TP track killed by cut");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(5, "TP fitted passing cut");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(4, "TP fitted");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(3, "Accepted track");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(2, "Number of Stubs");
        profNumFittedCands_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(1, "Fitted tracks including fakes");
        profNumFittedCands_[fitName].at(SimpleLR)->LabelsOption("d");
        
        hisFitQinvPtMatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitQinvPtMatched_"+(fitName)).c_str(),"Fitted q/p_{T} for matched tracks", 100, -0.5, 0.5 )));
        hisFitPhi0Matched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitPhi0Matched_"+(fitName)).c_str(), "Fitted #phi_{0} for matched tracks", 70, -3.5, 3.5 )));
        hisFitD0Matched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitD0Matched_"+(fitName)).c_str(), "Fitted d_{0} for matched tracks", 100, -1., 1. )));
        hisFitZ0Matched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitZ0Matched_"+(fitName)).c_str(), "Fitted z_{0} for matched tracks", 100, -25., 25. )));
        hisFitEtaMatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitEtaMatched_"+(fitName)).c_str(), "Fitted #eta for matched tracks", 70, -3.5, 3.5 )));
        
        hisFitChi2Matched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("Chi2Matched_"+(fitName)).c_str(), "#Chi^{2}", 100, 0., 100. )));
        hisFitChi2DofMatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("Chi2DofMatched_"+(fitName)).c_str(), "#Chi^{2}/DOF", 100, 0., 10. )));
        
        hisFitQinvPtUnmatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitQinvPtUnmatched_"+(fitName)).c_str(), "Fitted q/p_{T} for unmatched tracks", 100, -0.5, 0.5 )));
        hisFitPhi0Unmatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitPhi0Unmatched_"+(fitName)).c_str(), "Fitted #phi_{0} for unmatched tracks", 70, -3.5, 3.5 )));
        hisFitD0Unmatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitD0Unmatched_"+(fitName)).c_str(), "Fitted d_{0} for unmatched tracks", 100, -1., 1. )));
        hisFitZ0Unmatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitZ0Unmatched_"+(fitName)).c_str(), "Fitted z_{0} for unmatched tracks", 100, -25., 25. )));
        hisFitEtaUnmatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitEtaUnmatched_"+(fitName)).c_str(), "Fitted #eta for unmatched tracks", 70, -3.5, 3.5 )));
        
        hisFitChi2Unmatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitChi2Unmatched_"+(fitName)).c_str(), "#Chi^{2}", 100, 0., 100. )));
        hisFitChi2DofUnmatched_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitChi2DofUnmatched_"+(fitName)).c_str(), "#Chi^{2}/DOF", 100, 0., 10. )));
        
        hisFitVsTrueQinvPtGoodChi2_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsTrueQinvPtGoodChi2_"+(fitName)).c_str(), "; TP q/p_{T}; Fitted q/p_{T} (good #chi^2)", 100, -0.5, 0.5, 100, -0.5, 0.5 )));
        hisFitVsTruePhi0GoodChi2_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsTruePhi0GoodChi2_"+(fitName)).c_str(), "; TP #phi_{0}; Fitted #phi_{0} (good #chi^2)", 70, -3.5, 3.5, 70, -3.5, 3.5 )));
        hisFitVsTrueD0GoodChi2_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsTrueD0GoodChi2_"+(fitName)).c_str(), "; TP d_{0}; Fitted d_{0} (good #chi^2)", 100, -1., 1., 100, -1., 1. )));
        hisFitVsTrueZ0GoodChi2_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsTrueZ0GoodChi2_"+(fitName)).c_str(), "; TP z_{0}; Fitted z_{0} (good #chi^2)" , 100, -25., 25., 100, -25., 25. )));
        hisFitVsTrueEtaGoodChi2_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsTrueEtaGoodChi2_"+(fitName)).c_str(), "; TP #eta; Fitted #eta (good #chi^2)", 70, -3.5, 3.5, 70, -3.5, 3.5 )));
        hisFitVsSeedQinvPtGenCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedQinvPtGenCand_"+(fitName)).c_str(), "; Seed q/p_{T} (Genuine Cand); Fitted q/p_{T}", 100, -0.5, 0.5, 100, -0.5, 0.5 )));
        
        
        hisFitVsSeedPhi0GenCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedPhi0GenCand_"+(fitName)).c_str(), "; Seed #phi_{0} (Genuine Cand); Fitted #phi_{0}", 70, -3.5, 3.5, 70, -3.5, 3.5 )));
        hisFitVsSeedD0GenCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedD0GenCand_"+(fitName)).c_str(), "; Seed d_{0} (Genuine Cand); Fitted d_{0}", 100, -1., 1., 100, -1., 1. )));
        hisFitVsSeedZ0GenCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedZ0GenCand_"+(fitName)).c_str(), "; Seed z_{0} (Genuine Cand); Fitted z_{0}", 100, -25., 25., 100, -25., 25. )));
        hisFitVsSeedEtaGenCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedEtaGenCand_"+(fitName)).c_str(), "; Seed #eta (Genuine Cand); Fitted #eta", 70, -3.5, 3.5, 70, -3.5, 3.5 )));
        
        
        hisFitQinvPtResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitQinvPtResGoodChi2_"+(fitName)).c_str(), "Fitted minus true q/p_{T} (good #chi^2)", 100, -0.5, 0.5 )));
        hisFitPhi0ResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitPhi0ResGoodChi2_"+(fitName)).c_str(), "Fitted minus true #phi_{0} (good #chi^2)", 100, -0.1, 0.1 )));
        hisFitD0ResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitD0ResGoodChi2_"+(fitName)).c_str(), "Fitted minus true d_{0} (good #chi^2)", 100,  -1., 1. )));
        hisFitZ0ResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitZ0ResGoodChi2_"+(fitName)).c_str(), "Fitted minus true z_{0} (good #chi^2)", 100, -10., 10. )));
        hisFitEtaResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitEtaResGoodChi2_"+(fitName)).c_str(), "Fitted minus true #eta (good #chi^2)", 100, -0.1, 0.1 )));
        
        hisSeedQinvPtResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedQinvPtResGoodChi2_"+(fitName)).c_str(), "True minus seed q/p_{T} (good #chi^2)", 100, -0.5, 0.5 )));
        hisSeedPhi0ResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedPhi0ResGoodChi2_"+(fitName)).c_str(), "True minus seed #phi_{0} (good #chi^2)", 100, -0.1, 0.1 )));
        hisSeedD0ResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedD0ResGoodChi2_"+(fitName)).c_str(), "True minus seed d_{0} (good #chi^2)", 100,  -1., 1. )));
        hisSeedZ0ResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedZ0ResGoodChi2_"+(fitName)).c_str(), "True minus seed z_{0} (good #chi^2)", 100, -10., 10. )));
        hisSeedEtaResGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("SeedEtaResGoodChi2_"+(fitName)).c_str(), "True minus seed #eta (good #chi^2)", 100, -0.1, 0.1 )));
        
        hisFitVsSeedQinvPtFakeCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedQinvPtFakeCand_"+(fitName)).c_str(), "; Seed q/p_{T} (Fake Cand); Fitted q/p_{T}", 100, -0.5, 0.5, 100, -0.5, 0.5 )));
        hisFitVsSeedPhi0FakeCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedPhi0FakeCand_"+(fitName)).c_str(), "; Seed #phi_{0} (Fake Cand); Fitted #phi_{0}", 70, -3.5, 3.5, 70, -3.5, 3.5 )));
        hisFitVsSeedD0FakeCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedD0FakeCand_"+(fitName)).c_str(), "; Seed d_{0} (Fake Cand); Fitted d_{0}", 100, -1., 1., 100, -1., 1. )));
        hisFitVsSeedZ0FakeCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedZ0FakeCand_"+(fitName)).c_str(), "; Seed z_{0} (Fake Cand); Fitted z_{0}", 100, -25., 25., 100, -25., 25. )));
        hisFitVsSeedEtaFakeCand_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FitVsSeedEtaFakeCand_"+(fitName)).c_str(), "; Seed #eta (Fake Cand); Fitted #eta", 70, -3.5, 3.5, 70, -3.5, 3.5 )));
        
        float maxEta = settings_->maxStubEta();
        hisQoverPtResVsTrueEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("QoverPtResVsTrueEta_"+(fitName)).c_str(), "q/p_{T} resolution; |#eta|; q/p_{T} resolution", 24, 0.0, maxEta)));
        hisPhi0ResVsTrueEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("PhiResVsTrueEta_"+(fitName)).c_str(), "#phi_{0} resolution; |#eta|; #phi_{0} resolution", 24, 0.0, maxEta)));
        hisEtaResVsTrueEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("EtaResVsTrueEta_"+(fitName)).c_str(), "#eta resolution; |#eta|; #eta resolution", 24, 0.0, maxEta)));
        hisZ0ResVsTrueEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("Z0ResVsTrueEta_"+(fitName)).c_str(), "z_{0} resolution; |#eta|; z_{0} resolution", 24, 0.0, maxEta)));
        hisD0ResVsTrueEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("D0ResVsTrueEta_"+(fitName)).c_str(), "d_{0} resolution; |#eta|; d_{0} resolution", 24, 0.0, maxEta)));
        
        hisQoverPtResVsTrueZ0_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("QoverPtResVsTrueZ0_"+(fitName)).c_str(), "q/p_{T} resolution; z_{0}; q/ p_{T} resolution", 50, -25.,25.)));
        hisPhi0ResVsTrueZ0_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("PhiResVsTrueZ0_"+(fitName)).c_str(), "#phi_{0} resolution; z_{0}; #phi_{0} resolution", 50, -25.,25.)));
        hisEtaResVsTrueZ0_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("EtaResVsTrueZ0_"+(fitName)).c_str(), "#eta resolution; z_{0}; #eta resolution", 50, -25.,25.)));
        hisZ0ResVsTrueZ0_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("Z0ResVsTrueZ0_"+(fitName)).c_str(), "z_{0} resolution; z_{0}; z_{0} resolution", 50, -25.,25.)));
        hisD0ResVsTrueZ0_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("D0ResVsTrueZ0_"+(fitName)).c_str(), "d_{0} resolution; z_{0}; d_{0} resolution", 50, -25.,25.)));
        
        hisQoverPtResVsTrueInvPt_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("QoverPtResVsTrueInvPt_"+(fitName)).c_str(), "q/p_{T} resolution; 1/p_{T}; q/p_{T} resolution", 25, 0.0, maxAbsQoverPt)));
        hisPhi0ResVsTrueInvPt_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("PhiResVsTrueInvPt_"+(fitName)).c_str(), "#phi_{0} resolution; 1/p_{T}; #phi_{0} resolution", 25, 0.0, maxAbsQoverPt)));
        hisEtaResVsTrueInvPt_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("EtaResVsTrueInvPt_"+(fitName)).c_str(), "#eta resolution; 1/p_{T}; #eta resolution", 25, 0.0, maxAbsQoverPt)));
        hisZ0ResVsTrueInvPt_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("Z0ResVsTrueInvPt_"+(fitName)).c_str(), "z_{0} resolution; 1/p_{T}; z_{0} resolution", 25, 0.0, maxAbsQoverPt)));
        hisD0ResVsTrueInvPt_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("D0ResVsTrueInvPt_"+(fitName)).c_str(), "d_{0} resolution; 1/p_{T}; d_{0} resolution", 25, 0.0, maxAbsQoverPt)));
        
        hisTrueFittedChiSquaredVsTrueEta_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("TrueFittedChiSqauredVsEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} values and TP's #eta; #Chi^{2}; #eta", 50000, -0.5, 50000.5, 250, -2.50, 2.50 )));
        hisTrueFittedChiSquaredDofVsTrueEta_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("TrueFittedChiSqauredDofVsEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} Degrees of Freedom and TP's #eta; #Chi^{2}; #eta", 5001, -0.5, 5000.5, 250, -2.50, 2.50 )));
        hisTrueFittedChiSquaredVsFittedEta_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("TrueFittedChiSqauredVsFittedEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} values and fitted #eta; #Chi^{2}; #eta", 50000, -0.5, 50000.5, 250, -2.50, 2.50 )));
        hisTrueFittedChiSquaredDofVsFittedEta_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("TrueFittedChiSqauredDofVsFittedEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} Degrees of Freedom and fitted  #eta; #Chi^{2}; #eta", 5001, -0.5, 5000.5, 250, -2.50, 2.50 )));
        hisFittedChiSquaredFunctionOfStubs_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FittedChiSquaredFunctionOfStubs_"+(fitName)).c_str(), "Comparison of #Chi^{2} values as a function of stubs produced by TP; # of stubs; #Chi^{2}", 21, -0.5, 20.5, 500000, -0.5, 500000.5 )));
        hisFittedChiSquaredDofFunctionOfStubs_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("FittedChiSquaredDofFunctionOfStubs_"+(fitName)).c_str(), "Comparison of #Chi^{2} Degrees of Freedom as a function of stubs produced by TP; # of stubs; #Chi^{2} Degrees of Freedom", 21, -0.5, 20.5, 5001, -0.5, 5000.5 )));
        
        hisTrueEtaMatchedGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("TrueEtaMatchedGoodChi2_"+(fitName)).c_str(), "True #eta for matched tracks (good #chi^2)", 70, -3.5, 3.5 )));
        hisTrueEtaMatchedBadChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("TrueEtaMatchedBadChi2_"+(fitName)).c_str(), "True #eta for matched tracks (bad #chi^2)", 70, -3.5, 3.5 )));
        hisStubPurityMatchedGoodChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FracMatchedStubsMatchedGoodChi2_"+(fitName)).c_str(), "Purity of stubs on matched tracks (good #chi^2)", 102, -0.01, 1.01 )));
        hisStubPurityMatchedBadChi2_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FracMatchedStubsMatchedBadChi2_"+(fitName)).c_str(), "Purity of stubs on matched tracks (bad #chi^2)", 102, -0.01, 1.01 )));
        
        profChi2DofVsInvPtPERF_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("Chi2DofVsInvPtPERF_"+(fitName)).c_str(), "; q/Pt; Mean sqrt(#Chi^{2}/DOF) of perfect tracks", 16, 0.,  1./settings_->houghMinPt(), 0., 25.)));
        profBigChi2DofVsInvPtPERF_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("BigChi2DofVsInvPtPERF_"+(fitName)).c_str(), "; q/Pt; Frac perfect tracks with #Chi^{2}/DOF > 10", 16, 0.,  1./settings_->houghMinPt())));
        hisD0TPBigChi2DofPERF_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("D0TPBigChi2DofPERF"+(fitName)).c_str(),"; True d0 for high Pt tracks with bad #chi2/DOF", 100,0.0,0.5)));
        hisD0TPSmallChi2DofPERF_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("D0TPSmallChi2DofPERF"+(fitName)).c_str(),"; True d0 for high Pt tracks with good #chi2/DOF",100,0.0,0.5)));
        
        hisNumMatchedStubsKilledVsKilled_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("NumMatchedStubsKilledVsKilled_"+(fitName)).c_str(), "; All stubs killed by fit; Good stubs killed by fit", 10,-0.5,9.5,10,-0.5,9.5)));
        profTrksKilledByFit_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("TrksKilledByFit_"+(fitName)).c_str(), "Track category; Fraction of tracks killed by fit", 2,0.5,2.5,-999.,999.)));
        profTrksKilledByFit_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(1, "matched");
        profTrksKilledByFit_[fitName].at(SimpleLR)->GetXaxis()->SetBinLabel(2, "unmatched");
        hisNumStubsVsPurity_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("NumStubsVsPurity_"+(fitName)).c_str(), "; Number of stubs; Purity", 30, 0.0, 30.0, 100, 0.0, 1.0)));
        
        
        // Histograms specific to Linear Regression track fitter.
        if (fitName.find("LinearRegression") != string::npos) {
            hisNumIterations_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("NumIterations_"+(fitName)).c_str(), "; Number of Iterations", 30, -0.5, 29.5)));
            hisFailingState_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FailingState_"+(fitName)).c_str(), "; State which lost matching", 10, 0.0, 10)));
            hisTotalStateCalls_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("TotalStateCalls_"+(fitName)).c_str(), "; state; Total number of State calls", 10, 0.0, 10)));
            hisRelativeStateCalls_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("RelativeStateCalls_"+(fitName)).c_str(), "; state; Number of State calls relative to total State calls", 10, 0.0, 10)));
        
        }
        
        hisNumStubsFitKills_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("NumStubsFitKills_"+(fitName)).c_str(), "; Stubs per track killed by fit", 30, -0.5, 29.5)));
        hisNumStubsFitKillsVsPurity_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("NumStubsFitKillsVsPurity_"+(fitName)).c_str(), "; Stubs per track killed by fit; Purity", 30, 0.0, 30, 100, -0.01, 1.01 )));
        hisNumStubsFitKillsVsPurityMatched_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("NumStubsFitKillsVsPurityMatched_"+(fitName)).c_str(), "; Stubs per track killed bu fit; Purity", 30, 0.0, 30.0, 100, -0.01, 1.01 )));
        hisNumStubsFitKillsVsPurityUnmatched_[fitName].insert( pair< ObjectType, TH2F* > (SimpleLR, dirs_.at(SimpleLR).make<TH2F>(("NumStubsFitKillsVsPurityUnmatched_"+(fitName)).c_str(), "; Stubs per track killed by fit; Purity", 30, 0.0, 30.0, 100, -0.01, 1.01 )));
        
        // Duplicate track histos.

        profDupFitTrksVsEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("DupFitTrksVsEta_"+(fitName)).c_str() ,"; #eta; Fraction of duplicate tracks",12,0.,3.)));
        profDupFitTrksVsInvPt_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("DupFitTrksVsInvPt_"+(fitName)).c_str() ,"; 1/Pt; Fraction of duplicate tracks",16,0.,maxAbsQoverPt)));
        // Fake track histos.
        profFakeFitTrksVsEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("FakeFitTrksVsEta_"+(fitName)).c_str() ,"; #eta; Fraction of fake tracks",12,0.,3.)));
        profFakeFitTrksVsZ0_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("FakeFitTrksVsZ0_"+(fitName)).c_str() ,"; z_{0}; Fraction of fake tracks",50,-25.,25.)));
        profFitFracTrueStubsVsLayer_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("FitFracTrueStubsVsLayer_"+(fitName)).c_str() ,";Layer ID; fraction of true stubs",30,0.5,30.5)));
        profFitFracTrueStubsVsEta_[fitName].insert( pair< ObjectType, TProfile* > (SimpleLR, dirs_.at(SimpleLR).make<TProfile>(("FitFracTrueStubsVsEta_"+(fitName)).c_str() ,";#eta; fraction of true stubs",24,0.,3.)));
        
        // Histos for tracking efficiency vs. TP kinematics. (Binning must match similar histos in bookSimpleLR()).

        hisFitTPinvptForEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPinvptForEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for effi. measurement);",16,0.,maxAbsQoverPt)));
        hisFitTPetaForEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPetaForEff_"+(fitName)).c_str(),"; #eta of TP (used for effi. measurement);",20,-3.,3.)));
        hisFitTPphiForEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPphiForEff_"+(fitName)).c_str(),"; #phi of TP (used for effi. measurement);",20,-M_PI,M_PI)));
        
        // Histo for efficiency to reconstruct track perfectly (no incorrect hits). (Binning must match similar histos in bookTrackCands()).

        hisPerfFitTPinvptForEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("PerfFitTPinvptForEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for perf. effi. measurement);",16,0.,maxAbsQoverPt)));
        hisPerfFitTPetaForEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("PerfFitTPetaForEff_"+(fitName)).c_str(),"; #eta of TP (used for perfect effi. measurement);",20,-3.,3.)));
        
        // Histos for algorithmic tracking efficiency vs. TP kinematics. (Binning must match similar histos in bookSimpleLR()).

        hisFitTPinvptForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPinvptForAlgEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for alg. effi. measurement);",16,0.,maxAbsQoverPt)));
        hisFitTPetaForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPetaForAlgEff_"+(fitName)).c_str(),"; #eta of TP (used for alg. effi. measurement);",20,-3.,3.)));
        hisFitTPphiForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPphiForAlgEff_"+(fitName)).c_str(),"; #phi of TP (used for alg. effi. measurement);",20,-M_PI,M_PI)));
        
        // Histo for efficiency to reconstruct track perfectly (no incorrect hits). (Binning must match similar histos in bookSimpleLR()).

        hisPerfFitTPinvptForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("PerfFitTPinvptForAlgEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for perf. alg. effi. measurement);",16,0.,maxAbsQoverPt)));
        hisPerfFitTPetaForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("PerfFitTPetaForAlgEff_"+(fitName)).c_str(),"; #eta of TP (used for perf. alg. effi. measurement);",20,-3.,3.)));
        
        // Histos for algorithmic tracking efficiency vs. TP production point. (Binning must match similar histos in bookSimpleLR()).

        hisFitTPd0ForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPd0ForAlgEff_"+(fitName)).c_str() ,"; Production point x of TP (used for alg. effi. measurement);",50,0.,1.)));
        hisFitTPz0ForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPz0ForAlgEff_"+(fitName)).c_str() ,"; Production point y of TP (used for alg. effi. measurement);",50,0.,25.)));
        
        // Histo for algorithmic tracking efficiency vs sector number (to check if looser cuts are needed in certain regions)

        hisFitTPphisecForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPphisecForAlgEff_"+(fitName)).c_str() ,"; #phi sector of TP (used for alg. effi. measurement);",nPhi,-0.5,nPhi-0.5)));
        hisFitTPetasecForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("FitTPetasecForAlgEff_"+(fitName)).c_str() ,"; #eta sector of TP (used for alg. effi. measurement);",nEta,-0.5,nEta-0.5)));
        hisPerfFitTPphisecForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("PerfFitTPphisecForAlgEff_"+(fitName)).c_str() ,"; #phi sector of TP (used for perf. alg. effi. measurement);",nPhi,-0.5,nPhi-0.5)));
        hisPerfFitTPetasecForAlgEff_[fitName].insert( pair< ObjectType, TH1F* > (SimpleLR, dirs_.at(SimpleLR).make<TH1F>(("PerfFitTPetasecForAlgEff_"+(fitName)).c_str() ,"; #eta sector of TP (used for perf. alg. effi. measurement);",nEta,-0.5,nEta-0.5)));
        
    }
    
}



void TMTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    using namespace edm;
    //////// Fill Input Data ////////
  InputData inputData(iEvent, iSetup, settings_, tpInputTag, stubInputTag, stubTruthInputTag, clusterTruthInputTag );
                       
  const vector<const Stub*>& vStubs = inputData.getStubs();
  const vector<TP>&          vTPs   = inputData.getTPs();
    
    unsigned int nStubsGenuine = 0;
    unsigned int nStubsWithTP = 0;
    unsigned int nStubsWithTPforEff = 0;
    for (const Stub* stub : vStubs) {
        if (stub->genuine()) {
            nStubsGenuine++;
            if (stub->assocTP() != nullptr) {
                nStubsWithTP++;
                if (stub->assocTP()->useForEff()) nStubsWithTPforEff++;
            }
        }
    }
    profNumStubs_.at(InputDt)->Fill(1, vStubs.size());
    profNumStubs_.at(InputDt)->Fill(2, nStubsGenuine);
    profNumStubs_.at(InputDt)->Fill(3, nStubsWithTP);
    profNumStubs_.at(InputDt)->Fill(4, nStubsWithTPforEff);
    

  for (const Stub* stub : vStubs) {
    hisStubsVsEta_.at(InputDt)->Fill(stub->eta());
    hisStubsVsR_.at(InputDt)->Fill(stub->r());
    hisStubsVsRVsZ_.at(InputDt)->Fill(stub->z(),stub->r());
    hisStubsModuleVsRVsZ_.at(InputDt)->Fill(stub->minZ(),stub->minR());
  }
    
    
    for (const Stub* stub : vStubs) {
        hisStubsVsEta_.at(InputDt)->Fill(stub->eta());
        hisStubsVsR_.at(InputDt)->Fill(stub->r());
        hisStubsVsRVsZ_.at(InputDt)->Fill( stub->z(), stub->r() );
        hisStubsModuleVsRVsZ_.at(InputDt)->Fill( stub->minZ(), stub->minR() );
        hisStubsModuleVsRVsZ_.at(InputDt)->Fill( stub->maxZ(), stub->maxR() );
        
        hisStubsModuleTiltVsZ_.at(InputDt)->Fill( stub->minZ(), stub->moduleTilt() );
        hisStubsModuleTiltVsZ_.at(InputDt)->Fill( stub->maxZ(), stub->moduleTilt() );
        
        if ( stub->outerModuleAtSmallerR() ) {
            hisStubsVsRVsZ_outerModuleAtSmallerR_.at(InputDt)->Fill( stub->z(), stub->r() );
        }
        
        hisStubsdPhiCorrectionVsZ_.at(InputDt)->Fill(  stub->minZ(), stub->dphiOverBendCorrection() );
        
        hisStubsVsRVsPhi_.at(InputDt)->Fill( stub->r() * sin( stub->phi() ), stub->r() * cos( stub->phi() ) );
        hisStubsModuleVsRVsPhi_.at(InputDt)->Fill( stub->minR() * sin( stub->minPhi() ), stub->minR() * cos( stub->minPhi() ) );
        hisStubsModuleVsRVsPhi_.at(InputDt)->Fill( stub->maxR() * sin( stub->maxPhi() ), stub->maxR() * cos( stub->maxPhi() ) );
        
        if ( stub->outerModuleAtSmallerR() ) {
            hisStubsVsRVsPhi_outerModuleAtSmallerR_.at(InputDt)->Fill( stub->r() * sin( stub->phi() ), stub->r() * cos( stub->phi() ) );
        }
        
    }
    
    // Count tracking particles.
    unsigned int nTPforEff = 0;
    unsigned int nTPforAlgEff = 0;
    for (const TP& tp: vTPs) {
        if (tp.useForEff())  nTPforEff++;
        if (tp.useForAlgEff()) nTPforAlgEff++;
    }
    profNumTPs_.at(InputDt)->Fill(1, vTPs.size());
    profNumTPs_.at(InputDt)->Fill(2, nTPforEff);
    profNumTPs_.at(InputDt)->Fill(3, nTPforAlgEff);
    
    // Study efficiency of stubs to pass front-end electronics cuts.
    
    const vector<Stub>& vAllStubs = inputData.getAllStubs(); // Get all stubs prior to FE cuts to do this.
    for (const Stub s : vAllStubs) {
        unsigned int layerOrTenPlusRing = s.barrel()  ?  s.layerId()  :  10 + s.endcapRing();
        // Fraction of all stubs (good and bad) failing tightened front-end electronics cuts.
        hisStubKillFE_.at(InputDt)->Fill(layerOrTenPlusRing, (! s.frontendPass()));
        // Fraction of stubs rejected by window cut in DataCorrection.h
        // If it is non-zero, then encoding in DataCorrection.h should ideally be changed to make it zero.
        hisStubKillDataCorr_.at(InputDt)->Fill(layerOrTenPlusRing, s.stubFailedDataCorrWindow());
    }
    
    // Study efficiency for good stubs of tightened front end-electronics cuts.
    for (const TP& tp : vTPs) {
        if (tp.useForAlgEff()) {// Only bother for stubs that are on TP that we have a chance of reconstructing.
            hisPhysicsPt_.at(InputDt)->Fill(tp.pt());
            const vector<const Stub*> stubs = tp.assocStubs();
            for (const Stub* s : stubs) {
                hisStubIneffiVsInvPt_.at(InputDt)->Fill(1./tp.pt()    , (! s->frontendPass()) );
                hisStubIneffiVsEta_.at(InputDt)->Fill  (fabs(tp.eta()), (! s->frontendPass()) );
            }
        }
        if(tp.useForVertexReco() and !tp.physicsCollision() ){
            hisPileUpPt_.at(InputDt)->Fill(tp.pt());
        }
    }
    
    // Plot stub bend-derived information.
    for (const Stub* stub : vStubs) {
        hisPtStub_.at(InputDt)->Fill(stub->qOverPt());
        hisDelPhiStub_.at(InputDt)->Fill(stub->dphi());
        hisBendStub_.at(InputDt)->Fill(stub->dphi() / stub->dphiOverBend());
        // Number of bend values merged together by loss of a bit.
        hisNumMergedBend_.at(InputDt)->Fill(stub->numMergedBend());
        // Min. & max allowed q/Pt obtained from stub bend.
        float minQoverPt = max(float(-1./(settings_->houghMinPt())), stub->qOverPt() - stub->qOverPtres());
        float maxQoverPt = min(float(1./(settings_->houghMinPt())), stub->qOverPt() + stub->qOverPtres());
        // Frac. of full q/Pt range allowed by stub bend.
        float fracAllowed = (maxQoverPt - minQoverPt)/(2./(settings_->houghMinPt()));
        hisBendFilterPower_.at(InputDt)->Fill(fracAllowed);
        unsigned int layerOrTenPlusRing = stub->barrel()  ?  stub->layerId()  :  10 + stub->endcapRing();
        hisBendVsLayerOrRing_.at(InputDt)->Fill(layerOrTenPlusRing, stub->bend());
        // Also plot bend prior to degradation.
        hisBendFEVsLayerOrRing_.at(InputDt)->Fill(layerOrTenPlusRing, stub->bendInFrontend());
    }
    
    // Look at stub resolution.
    for (const TP& tp: vTPs) {
        if (tp.useForAlgEff()) {
            const vector<const Stub*>& assStubs= tp.assocStubs();
            hisNumStubsPerTP_.at(InputDt)->Fill( assStubs.size() );
            
            unsigned int numPSstubs = 0;
            unsigned int num2Sstubs = 0;
            
            //cout<<"=== TP === : index="<<tp.index()<<" pt="<<tp.pt()<<" q="<<tp.charge()<<" phi="<<tp.phi0()<<" eta="<<tp.eta()<<" z0="<<tp.z0()<<endl;
            for (const Stub* stub: assStubs) {
                
                if ( stub->psModule() ) ++numPSstubs;
                else ++num2Sstubs;
                
                //cout<<"    stub : index="<<stub->index()<<" barrel="<<stub->barrel()<<" r="<<stub->r()<<" phi="<<stub->phi()<<" z="<<stub->z()<<" bend="<<stub->bend()<<" assocTP="<<stub->assocTP()->index()<<endl;
                hisPtResStub_.at(InputDt)->Fill(stub->qOverPt() - tp.charge()/tp.pt());
                hisDelPhiResStub_.at(InputDt)->Fill(stub->dphi() - tp.dphi(stub->r()));
                
                if ( stub->moduleTilt() > M_PI / 2 - 0.1  || !stub->barrel() ) {
                    hisDelPhiResStub_notTilted_.at(InputDt)->Fill(stub->dphi() - tp.dphi(stub->r()));
                } else {
                    hisDelPhiResStub_tilted_.at(InputDt)->Fill(stub->dphi() - tp.dphi(stub->r()));
                }
                hisBendResStub_.at(InputDt)->Fill( (stub->dphi() - tp.dphi(stub->r())) / stub->dphiOverBend() );
                // This checks if the TP multiple scattered before producing the stub or hit resolution effects.
                hisPhiStubVsPhiTP_.at(InputDt)->Fill( reco::deltaPhi(stub->phi(), tp.trkPhiAtStub( stub )) );
                // This checks how wide overlap must be if using phi0 sectors, with no stub bend info used for assignment.
                hisPhiStubVsPhi0TP_.at(InputDt)->Fill( reco::deltaPhi(stub->phi(), tp.phi0()) );
                // This checks how wide overlap must be if using phi0 sectors, with stub bend info used for assignment
                hisPhi0StubVsPhi0TP_.at(InputDt)->Fill( reco::deltaPhi(stub->trkPhiAtR(0.).first, tp.phi0()) );
                // This normalizes the previous distribution to the predicted resolution to check if the latter is OK.
                hisPhi0StubVsPhi0TPres_.at(InputDt)->Fill( reco::deltaPhi(stub->trkPhiAtR(0.).first, tp.phi0()) / stub->trkPhiAtRres(0.));
                // This checks how wide overlap must be if using phi65 sectors, with no stub bend info used for assignment.
                hisPhiStubVsPhi65TP_.at(InputDt)->Fill( reco::deltaPhi(stub->phi(), tp.trkPhiAtR(65.)) );
                // This checks how wide overlap must be if using phi65 sectors, with stub bend info used for assignment, optionally reducing discrepancy by uncertainty expected from 2S module strip length.
                pair<float, float> phiAndErr = stub->trkPhiAtR(65.);
                double dPhi = reco::deltaPhi( phiAndErr.first, tp.trkPhiAtR(65.));
                hisPhi65StubVsPhi65TP_.at(InputDt)->Fill( dPhi );
                // This normalizes the previous distribution to the predicted resolution to check if the latter is OK.
                hisPhi65StubVsPhi65TPres_.at(InputDt)->Fill( dPhi / stub->trkPhiAtRres(65.));
            }
            
            hisNumPSStubsPerTP_.at(InputDt)->Fill( numPSstubs );
            hisNum2SStubsPerTP_.at(InputDt)->Fill( num2Sstubs );
            
        }
    }
    
    for (const Stub* stub : vStubs) {
        // Note ratio of sensor pitch to separation (needed to understand how many bits this can be packed into).
        hisPitchOverSep_.at(InputDt)->Fill(stub->pitchOverSep());
        // Also note this same quantity times 1.0 in the barrel or z/r in the endcap. This product is known as "rho".
        float rho = stub->pitchOverSep();
        if ( ! stub->barrel() ) rho *= fabs(stub->z())/stub->r();
        hisRhoParameter_.at(InputDt)->Fill(rho);
    }
    
    // Check fraction of stubs sharing a common cluster.
    // Loop over both clusters in each stub, so looking for common clusters in seed (0) or correlation (1) sensor of module.
    typedef pair< unsigned int, pair<float, float> > ClusterLocation;
    for (unsigned int iClus = 0; iClus <= 1; iClus++) {
        map<ClusterLocation, unsigned int> commonClusterMap;
        for (const Stub* stub : vStubs) {
            // Encode detector ID & strip (or pixel) numbers in both dimensions.
            const ClusterLocation loc( stub->idDet(), pair<float, float>(stub->localU_cluster()[iClus], stub->localV_cluster()[iClus]) );
            if (commonClusterMap.find(loc) == commonClusterMap.end()) {
                commonClusterMap[loc] = 1;
            } else {
                commonClusterMap[loc]++;
            }
        }
        unsigned int nShare = 0;
        for (map<ClusterLocation, unsigned int>::const_iterator it = commonClusterMap.begin(); it != commonClusterMap.end(); it++) {
            if (it->second != 1) nShare += it->second; // 2 or more stubs share a cluster at this detid*strip.
        }
        if (iClus == 0) {
            hisFracStubsSharingClus0_.at(InputDt)->Fill(float(nShare)/float(vStubs.size()));
        } else {
            hisFracStubsSharingClus1_.at(InputDt)->Fill(float(nShare)/float(vStubs.size()));
        }
    }


  ////////// Fill CheckSector ///////////
    
    matrix<Sector>  mSectors(settings_->numPhiSectors(), settings_->numEtaRegions());
    matrix<HTpair>  mHtPairs(settings_->numPhiSectors(), settings_->numEtaRegions());
    
    //=== Initialization
    // Create utility for converting L1 tracks from our private format to official CMSSW EDM format.
    const ConverterToTTTrack converter(settings_);
    // Storage for EDM L1 track collection to be produced from Hough transform output (no fit).
    std::unique_ptr<TTTrackCollection>  htTTTracksForOutput(new TTTrackCollection);
    // Storage for EDM L1 track collection to be produced from fitted tracks (one for each fit algorithm being used).
    // auto_ptr cant be stored in std containers, so use C one, together with map noting which element corresponds to which algorithm.
    const unsigned int nFitAlgs = settings_->trackFitters().size();
    std::unique_ptr<TTTrackCollection> allFitTTTracksForOutput[nFitAlgs];
    std::unique_ptr< FitTrackCollection > allTrackFitTracks[nFitAlgs];
    
    map<string, unsigned int> locationInsideArray;
    unsigned int ialg = 0;
    for (const string& fitterName : settings_->trackFitters()) {
        std::unique_ptr<TTTrackCollection> fitTTTracksForOutput(new TTTrackCollection);
        allFitTTTracksForOutput[ialg] =  std::move( fitTTTracksForOutput );
        std::unique_ptr<FitTrackCollection> fitTracks(new FitTrackCollection);
        allTrackFitTracks[ialg] = std::move(fitTracks);
        locationInsideArray[fitterName] = ialg++;
    }
    
    
    
    //=== Do tracking in the r-phi Hough transform within each sector.
    
    unsigned ntracks(0);
    // Fill Hough-Transform arrays with stubs.
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
            
            Sector& sector = mSectors(iPhiSec, iEtaReg);
            HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            
            // Initialize constants for this sector.
            sector.init(settings_, iPhiSec, iEtaReg);
            htPair.init(settings_, iPhiSec, iEtaReg, sector.etaMin(), sector.etaMax(), sector.phiCentre());
            
            for (const Stub* stub: vStubs) {
                // Digitize stub as would be at input to GP. This doesn't need the octant number, since we assumed an integer number of
                // phi digitisation  bins inside an octant. N.B. This changes the coordinates & bend stored in the stub.
                // The cast allows us to ignore the "const".
                if (settings_->enableDigitize()) (const_cast<Stub*>(stub))->digitizeForGPinput(iPhiSec);
                
                // Check if stub is inside this sector
                bool inside = sector.inside( stub );
                
                if (inside) {
                    // Check which eta subsectors within the sector the stub is compatible with (if subsectors being used).
                    const vector<bool> inEtaSubSecs =  sector.insideEtaSubSecs( stub );
                    
                    // Digitize stub if as would be at input to HT, which slightly degrades its coord. & bend resolution, affecting the HT performance.
                    if (settings_->enableDigitize()) (const_cast<Stub*>(stub))->digitizeForHTinput(iPhiSec);
                    
                    // Store stub in Hough transform array for this sector, indicating its compatibility with eta subsectors with sector.
                    htPair.store( stub, inEtaSubSecs );
                }
            }
            
            // Find tracks in r-phi HT array.
            htPair.end(); // Calls htArrayRphi_.end() -> HTBase::end()
            // std::cout << iPhiSec << " " << iEtaReg << " " << htPair.getRphiHT().numTrackCands2D() << std::endl;
            // if ( htPair.getRphiHT().numTrackCands2D() > 0 ) {
            //   std::cout << "Number of tracks after r-phi HT : " << iPhiSec << " " << iEtaReg << " " << htPair.getRphiHT().numTrackCands2D() << std::endl;
            // }
        }
    }
    

    if (settings_->muxOutputsHT() && settings_->busySectorKill()) {
        // Multiplex outputs of several HT onto one pair of output opto-links.
        // This only affects tracking performance if option busySectorKill is enabled, so that tracks that
        // can't be sent down the link within the time-multiplexed period are killed.
        MuxHToutputs muxHT(settings_);
        muxHT.exec(mHtPairs);
    }

    
    //=== Optionally run r-z filters or r-z HT. Make 3D tracks.
    
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
            
            HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            // std::cout << iPhiSec << " " << iEtaReg << " " << htPair.getRphiHT().numTrackCands2D() << std::endl;
            // Convert these to 3D tracks (optionally by running r-z filters etc.)
            htPair.make3Dtracks();
            
            // Convert these tracks to EDM format for output (not used by Histos class).
            const vector<L1track3D>& vecTrk3D = htPair.trackCands3D();
            ntracks += vecTrk3D.size();
            for (const L1track3D& trk : vecTrk3D) {
                TTTrack< Ref_Phase2TrackerDigi_ > htTTTrack = converter.makeTTTrack(trk, iPhiSec, iEtaReg);
                htTTTracksForOutput->push_back( htTTTrack );
            }
        }
    }
    
    
    // Initialize the duplicate track removal algorithm that can optionally be run after the track fit.
    KillDupFitTrks killDupFitTrks;
    killDupFitTrks.init(settings_, settings_->dupTrkAlgFit());
    
    //=== Do a helix fit to all the track candidates.
    FitTrackCollection fitTracks;
    vector<std::pair<std::string, L1fittedTrack>> fittedTracks;
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
            
            HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            
            // In principal, should digitize stubs on track here using digitization relative to this phi sector.
            // However, previously digitized stubs will still be valid if digitization bins in phi align with
            // phi sector boundaries, so couldn't be bothered. If they do not align, the average effect of
            // digitization on the track fit will be correct, but the effect on individual tracks will not.
            
            // Get track candidates found by Hough transform in this sector.
            const vector<L1track3D>& vecTrk3D = htPair.trackCands3D();
            // Loop over all the fitting algorithms we are trying.
            for (const string& fitterName : settings_->trackFitters()) {
                // Fit all tracks in this sector
                vector<L1fittedTrack> fittedTracksInSec;
                for (const L1track3D& trk : vecTrk3D) {
                    L1fittedTrack fitTrack = fitterWorkerMap_[fitterName]->fit(trk, iPhiSec, iEtaReg);
                    // Store fitted tracks, such that there is one fittedTracks corresponding to each HT tracks.
                    // N.B. Tracks rejected by the fit are also stored, but marked.
                    fittedTracksInSec.push_back(fitTrack);
                }
                
                // Run duplicate track removal on the fitted tracks if requested.
                const vector<L1fittedTrack> filtFittedTracksInSec = killDupFitTrks.filter( fittedTracksInSec );
                
                // Store fitted tracks from entire tracker.
                for (const L1fittedTrack& fitTrk : filtFittedTracksInSec) {
                    fittedTracks.push_back(std::make_pair(fitterName, fitTrk));
                    // Convert these fitted tracks to EDM format for output (not used by Histos class).
                    // Only do this for valid fitted tracks, meaning that these EDM tracks do not correspond 1 to 1 with fittedTracks.
                    if (fitTrk.accepted()) {
                        TTTrack< Ref_Phase2TrackerDigi_ > fitTTTrack = converter.makeTTTrack(fitTrk, iPhiSec, iEtaReg);
                        allFitTTTracksForOutput[locationInsideArray[fitterName]]->push_back(fitTTTrack);
                    }
                }
            }
        }
    }
    
    for (const Stub* stub: vStubs) {
        if (settings_->enableDigitize()) (const_cast<Stub*>(stub))->reset_digitize();
    }

     ///////////////////////////////////////////////////////////////////////////////////////////

  for (const TP& tp : vTPs) {
    if (tp.useForAlgEff()) {
      unsigned int nStubs = tp.numAssocStubs(); // no. of stubs in this TP.

      // Number of stubs this TP has in best (eta,phi) sector, and also just dividing sectors in phi or just in eta.
      unsigned int nStubsInBestSec = 0; 
      unsigned int nStubsInBestEtaSec = 0; 
      unsigned int nStubsInBestPhiSec = 0; 
        
      // Loop over (eta, phi) sectors.
      for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
	for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {

	  const Sector& sector = mSectors(iPhiSec, iEtaReg);

	  // Count number of stubs in given tracking particle which are inside this (phi,eta) sector;
	  // or inside it if only the eta cuts are applied; or inside it if only the phi cuts are applied.
	  unsigned int nStubsInSec, nStubsInEtaSec, nStubsInPhiSec;
	  sector.numStubsInside( tp, nStubsInSec, nStubsInEtaSec, nStubsInPhiSec);

	  // Note best results obtained in any sector.
	  nStubsInBestSec    = max( nStubsInBestSec,    nStubsInSec);
	  nStubsInBestEtaSec = max( nStubsInBestEtaSec, nStubsInEtaSec);
	  nStubsInBestPhiSec = max( nStubsInBestPhiSec, nStubsInPhiSec);
        }
      }
        

      // Plot fraction of stubs on each TP in its best sector.
      hisFracStubsInSec_.at(CheckSectors)->Fill   ( float(nStubsInBestSec)    / float(nStubs) );
      hisFracStubsInEtaSec_.at(CheckSectors)->Fill( float(nStubsInBestEtaSec) / float(nStubs) );
      hisFracStubsInPhiSec_.at(CheckSectors)->Fill( float(nStubsInBestPhiSec) / float(nStubs) );
    }
  }

  for (const Stub* stub : vStubs) {

    // Number of (eta,phi), phi & eta sectors containing this stub.
    unsigned int nSecs = 0; 
    unsigned int nEtaSecs = 0; 
    unsigned int nPhiSecs = 0; 

    // Loop over (eta, phi) sectors.
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
      for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {

        const Sector& sector = mSectors(iPhiSec, iEtaReg);

        // Check if sector contains stub stub, and if so count it.
        // Take care to just use one eta (phi) typical region when counting phi (eta) sectors.
        if ( sector.inside   ( stub ) )                 nSecs++;
        if ( iPhiSec == 0 && sector.insideEta( stub ) ) nEtaSecs++;
        if ( iEtaReg == 0 && sector.insidePhi( stub ) ) nPhiSecs++;

	// Also note which tracker layers are present in each eta sector.
        if (iPhiSec == 0 && sector.insideEta( stub)) {
          const TP* assocTP = stub->assocTP();
          if (assocTP != nullptr) {
            if (assocTP->useForAlgEff()) {
              unsigned int lay = stub->layerId();
              if (lay > 20) lay -= 10; // Don't bother distinguishing two endcaps.
              hisLayerIDvsEtaSec_.at(CheckSectors)->Fill(iEtaReg, lay);
              hisLayerIDreducedvsEtaSec_.at(CheckSectors)->Fill(iEtaReg, stub->layerIdReduced()); // Plot also simplified layerID for hardware, which tries to avoid more than 8 ID in any given eta region.
            }
          }
        }
      }
    }

    // Plot number of sectors each stub appears in.
    hisNumSecsPerStub_.at(CheckSectors)->Fill   ( nSecs );
    hisNumEtaSecsPerStub_.at(CheckSectors)->Fill( nEtaSecs );
    hisNumPhiSecsPerStub_.at(CheckSectors)->Fill( nPhiSecs );

    if ( ! settings_->allowOver2EtaSecs()) {
      if (nEtaSecs > 2)  throw cms::Exception("Histos ERROR: Stub assigned to more than 2 eta regions. Please redefine eta regions to avoid this!")<<" stub r="<<stub->r()<<" eta="<<stub->eta()<<endl;
    }
  }

  //=== Loop over all sectors, counting the stubs in each one.
  for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
    unsigned int nStubsInEtaSec = 0; // Also counts stubs in eta sector, summed over all phi.
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
      const Sector& sector = mSectors(iPhiSec, iEtaReg);

      unsigned int nStubs = 0;
      for (const Stub* stub : vStubs) {
	if ( sector.inside( stub ) )  nStubs++;
      }
      hisNumStubsPerSec_.at(CheckSectors)->Fill(nStubs);
      nStubsInEtaSec += nStubs;
    }
    profNumStubsPerEtaSec_.at(CheckSectors)->Fill(iEtaReg, nStubsInEtaSec);
  }


  /////// HTrphi ////////

  for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
      const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
      const HTrphi& htRphi = htPair.getRphiHT();

      // Here, if a stub appears in multiple cells, it is counted multiple times.
      hisIncStubsPerHT_.at(HTrphiHis)->Fill( htRphi.numStubsInc() );
      // Here, if a stub appears in multiple cells, it is counted only once.
      hisExcStubsPerHT_.at(HTrphiHis)->Fill( htRphi.numStubsExc() );
    }
  }

  //--- Count number of stubs in each cell of HT array, summing over all the phi sectors within a given 
  //--- eta region. This determines the buffer size needed to store them in the firmware.

  // Loop over eta regions.
  for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
    // Get dimensions of HT array (assumed same for all phi sectors)
    unsigned int iPhiSecDummy = 0;
    const matrix<HTcell>& rphiHTcellsDummy = mHtPairs(iPhiSecDummy, iEtaReg).getRphiHT().getAllCells();
    const unsigned int nbins1 = rphiHTcellsDummy.size1();
    const unsigned int nbins2 = rphiHTcellsDummy.size2();
    // Loop over cells inside HT array
    for (unsigned int m = 0; m < nbins1; m++) {
      for (unsigned int n = 0; n < nbins2; n++) {
	// Loop over phi sectors
	unsigned int nStubsInCellPhiSum = 0;
        for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
          const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
          const HTrphi& htRphi = htPair.getRphiHT();
          const matrix<HTcell>& rphiHTcells = htRphi.getAllCells(); 
          nStubsInCellPhiSum += rphiHTcells(m,n).numStubs();
        }  
	// Plot total number of stubs in this cell, summed over all phi sectors.
        hisNumStubsInCellVsEta_.at(HTrphiHis)->Fill( nStubsInCellPhiSum, iEtaReg );
      }
    }
  }

  //--- Count number of cells assigned to track candidates by r-phi HT (before any rz filtering
  //--- or rz HT has been run).
  for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
      const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
      const HTrphi& htRphi = htPair.getRphiHT();
      hisStubsOnRphiTracksPerHT_.at(HTrphiHis)->Fill(htRphi.numStubsOnTrackCands2D()); 
    }
  }

  //////// RZfilters ////////
    

    if (settings_->useZTrkFilter() || settings_->useSeedFilter()) {
        
        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
            for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
                
                if (settings_->useZTrkFilter()) {
                    // Check number of track seeds per sector that r-z "Ztrk" filter checked.
                    const vector<unsigned int>  numSeedComb = htPair.getRZfilters().numZtrkSeedCombsPerTrk();
                    for (const unsigned int& num : numSeedComb) {
                        hisNumZtrkSeedCombinations_.at(RZfilters)->Fill(num);
                    }
                }
                
                
                
                
                if (settings_->useSeedFilter()) {
                    // Check number of track seeds per sector that r-z "seed" filter checked.
                    const vector<unsigned int>  numSeedComb = htPair.getRZfilters().numSeedCombsPerTrk();
                    for (const unsigned int& num : numSeedComb) {
                        hisNumSeedCombinations_.at(RZfilters)->Fill(num);
                    }
                    
                    // Same again, but this time only considering seeds the r-z filters defined as "good".
                    const vector<unsigned int>  numGoodSeedComb = htPair.getRZfilters().numGoodSeedCombsPerTrk();
                    for (const unsigned int& num : numGoodSeedComb) {
                        hisNumGoodSeedCombinations_.at(RZfilters)->Fill(num) ;
                    }
                    
                }
                
                 
                 
                 
                // Check correlation factor used inside zTrk filter (only of interest to experts)
                if (settings_->useZTrkFilter()) {
                    const HTrphi& htRphi = htPair.getRphiHT();
                    // Consider all genuine tracks found by r-phi HT (since they are input to Ztrk filter).
                    const vector<L1track2D>& tracksRphi = htRphi.trackCands2D();
                    for (const L1track2D& trk : tracksRphi) {
                        const TP* assocTP = trk.getMatchedTP();
                        if (assocTP != nullptr) {
                            if (assocTP->useForAlgEff()) {
                                const vector < const Stub* > stubs = trk.getStubs();
                                for(const Stub* s : stubs) {
                                    if (s->layerId()==1)  {
                                        for(const Stub* s2 : stubs){
                                            double r = 0.;
                                            if (s2->layerId()!=s->layerId()){
                                                double sum = 0.;
                                                for (int i = 0; i < 100; ++i){
                                                    double zB = -settings_->beamWindowZ() + (0.5+i)*settings_->beamWindowZ()/50; 
                                                    double z1 = s->zTrk() - (settings_->chosenRofZFilter() - s->r())*zB/s->r();
                                                    double z2 = s2->zTrk() - (settings_->chosenRofZFilter() - s2->r())*zB/s2->r();
                                                    sum = sum + z1*z2;
                                                }
                                                sum = (sum/100) - s->zTrk()*s2->zTrk();
                                                r = sum/(s->zTrkRes()*s2->zTrkRes());
                                                hisCorrelationZTrk_.at(RZfilters)->Fill(r);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


  ////// BusyEvents ////////

    const unsigned int numStubsCut = settings_->busySectorNumStubs();   // No. of stubs per HT array the hardware can output.
    const bool         eachCharge  = settings_->busySectorEachCharge(); // +ve & -ve tracks output on separate optical links?

    map<const L1track3D*, bool> trksInEntireTracker;

    unsigned int nBusySecIn  = 0;
    unsigned int nBusySecOut = 0;
    
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
        for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
            const Sector& sector = mSectors(iPhiSec, iEtaReg);
            const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            const HTrphi& htRphi = htPair.getRphiHT();
            const vector<L1track3D>& tracks = htPair.trackCands3D();
            
            //--- Look for too many stubs input to sector.
            
            unsigned int nStubsIn = htRphi.nReceivedStubs();
            // Plot fraction of input stubs that would be killed by 36BX period.
            for (unsigned int j = 0; j < nStubsIn; j++) {
                bool kill = (j >= numStubsCut);
                profFracStubsKilledVsEtaReg_.at(BusyEvents)->Fill(iEtaReg, kill);
            }
            bool tooBusyIn = (nStubsIn > numStubsCut);
            if (tooBusyIn) nBusySecIn++;
            profFracBusyInVsEtaReg_.at(BusyEvents)->Fill(iEtaReg, tooBusyIn); // Sector had too many input stubs.
            
            //--- Look for too many stubs assigned to output tracks.
            
            // Order tracks in increasing order of abs(q/Pt).
            // Use multimap rather than map to do this, as some tracks may have identical q/Pt, and it will store all of them, unlike map.
            multimap<float, const L1track3D*> orderedTrks;
            for (const L1track3D& trk : tracks) {
                orderedTrks.insert( pair<float, const L1track3D*>( fabs(trk.qOverPt()), &trk) );
            }
            
            // Create map containing L1 tracks found in whole of tracker together with flag indicating if the
            // track was killed because it was in a busy sector.
            map<const L1track3D*, bool> trksInSector;
            
            // Check how many tracks would be killed by 36BX period, assuming we kill preferentially low Pt ones.
            bool tooBusyOut = false;
            unsigned int nStubsOut          = 0;
            unsigned int nStubsOutPosCharge = 0;
            unsigned int nStubsOutNegCharge = 0;
            
            for (const auto& oTrk : orderedTrks) {
                float ptInv = oTrk.first;
                const L1track3D* trk = oTrk.second;
                bool kill = false;
                if (eachCharge) { // Are +ve and -ve charged tracks output on separate optical links to increase bandwidth?
                    if (trk->charge() > 0) {
                        nStubsOutPosCharge += trk->getNumStubs();
                        if (nStubsOutPosCharge > numStubsCut) kill = true;
                    } else {
                        nStubsOutNegCharge += trk->getNumStubs();
                        if (nStubsOutNegCharge > numStubsCut) kill = true;
                    }
                } else {
                    nStubsOut += trk->getNumStubs();
                    if (nStubsOut > numStubsCut) kill = true;
                }
                
                if (kill) tooBusyOut = true; // Note that some tracks were killed in this sector.
                
                profFracTracksKilledVsEtaReg_.at(BusyEvents)->Fill(iEtaReg, kill);
                profFracTracksKilledVsInvPt_.at(BusyEvents)->Fill(ptInv, kill);
                
                // Form a map of all tracks in the entire tracker & also just in this sector, with a flag indicating if they were killed as in a busy sector.
                trksInEntireTracker[trk] = kill;
                trksInSector[trk]        = kill;
            }
            
            if (tooBusyOut) nBusySecOut++;
            profFracBusyOutVsEtaReg_.at(BusyEvents)->Fill(iEtaReg, tooBusyOut); // Sector had too many output stubs.
            
            //--- Compare properties of sectors with/without too many output stubs.
            
            const vector<string> tnames = {"BusyOutSec", "QuietOutSec"};
            
            // Loop over sectors with too many/not too many output stubs.
            for (const string& tn : tnames) {
                if ((tn == "BusyOutSec" && tooBusyOut) || (tn == "QuietOutSec" && (! tooBusyOut))) {
                    
                    hisNumInputStubs_[tn].at(BusyEvents)->Fill(nStubsIn);
                    
                    // Check if q/Pt estimated from stub bend differs in busy & quiet sectors.
                    for (const Stub* stub : vStubs) {
                        if ( sector.inside( stub ) ) hisQoverPtInputStubs_[tn].at(BusyEvents)->Fill(abs(stub->qOverPt()));
                    }
                    
                    // Look at reconstructed tracks in this sector.
                    hisNumOutputStubs_[tn].at(BusyEvents)->Fill(nStubsOut);
                    hisNumTracks_[tn].at(BusyEvents)->Fill(tracks.size());
                    for (const L1track3D& trk : tracks) {
                        hisNumStubsPerTrack_[tn].at(BusyEvents)->Fill(trk.getNumStubs());
                        hisTrackQoverPt_[tn].at(BusyEvents)->Fill(trk.qOverPt());
                        hisTrackPurity_[tn].at(BusyEvents)->Fill(trk.getPurity());
                    }
                    
                    // Look at total Pt of truth particles in this sector to understand if it contains a jet.
                    unsigned int num_TP_physics = 0;
                    unsigned int num_TP_pileup  = 0;
                    float sumPt_TP_physics = 0.;
                    float sumPt_TP_pileup  = 0.;
                    for (const TP& tp : vTPs) {
                        bool tpInSector = (fabs(tp.trkPhiAtR(settings_->chosenRofPhi()) - sector.phiCentre()) < sector.sectorHalfWidth() &&
                                           tp.trkZAtR(settings_->chosenRofZ()) > sector.zAtChosenR_Min() &&
                                           tp.trkZAtR(settings_->chosenRofZ()) < sector.zAtChosenR_Max());
                        if (tpInSector) {
                            if (tp.physicsCollision()) { // distinguish truth particles from physics collision vs from pileup.
                                num_TP_physics++;
                                sumPt_TP_physics += tp.pt();
                            } else {
                                num_TP_pileup++;
                                sumPt_TP_pileup  += tp.pt();
                            }
                        }
                    }
                    hisNumTPphysics_[tn].at(BusyEvents)->Fill(num_TP_physics);
                    hisNumTPpileup_[tn].at(BusyEvents)->Fill(num_TP_pileup);
                    hisSumPtTPphysics_[tn].at(BusyEvents)->Fill(sumPt_TP_physics);
                    hisSumPtTPpileup_[tn].at(BusyEvents)->Fill(sumPt_TP_pileup);
                }
            }
            
            //--- Count tracking particles lost by killing tracks in individual busy sectors.
            if (tooBusyOut) {
                unsigned int nTPkilled = 0;
                for (const TP& tp: vTPs) {
                    if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
                        
                        bool tpRecoed = false;
                        bool tpRecoedSurvived = false;
                        for (const auto& trkm : trksInSector) {
                            const L1track3D* trk = trkm.first;
                            bool kill            = trkm.second;
                            if (trk->getMatchedTP() == &tp) {
                                tpRecoed = true;                        // Truth particle was reconstructed
                                if (! kill) tpRecoedSurvived = true;    // Ditto & reconstructed track wasn't killed by busy sector.
                            }
                        }
                        
                        bool tpKilled = tpRecoed && ( ! tpRecoedSurvived );
                        if (tpKilled) nTPkilled++;
                    }
                }
                hisNumTPkilledBusySec_.at(BusyEvents)->Fill(nTPkilled);
            }
        }
    }
    
    hisNumBusySecsInPerEvent_.at(BusyEvents)->Fill(nBusySecIn); // No. of sectors per event with too many input stubs.
    hisNumBusySecsOutPerEvent_.at(BusyEvents)->Fill(nBusySecOut); // No. of sectors per event with too many output stubs.
    
    //--- Check loss in tracking efficiency caused by killing tracks in busy sectors.
    
    for (const TP& tp: vTPs) {
        if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
            
            bool tpRecoed = false;
            bool tpRecoedSurvived = false;
            for (const auto& trkm : trksInEntireTracker) {
                const L1track3D* trk = trkm.first;
                bool kill            = trkm.second;
                if (trk->getMatchedTP() == &tp) {
                    tpRecoed = true;                        // Truth particle was reconstructed
                    if (! kill) tpRecoedSurvived = true;    // Ditto & reconstructed track wasn't killed by busy sector.
                }
            }
            bool tpKilled = tpRecoed && ( ! tpRecoedSurvived );
            profFracTPKilledVsEta_.at(BusyEvents)->Fill(fabs(tp.eta()), tpKilled);
            profFracTPKilledVsInvPt_.at(BusyEvents)->Fill(fabs(tp.qOverPt()), tpKilled);
        }
    }

    
    
  ///////// TrackCand //////////
    
    // Debug histogram for LR track fitter.
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
        const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
        const std::vector< L1track3D >& tracks = htPair.trackCands3D();
        for ( auto t : tracks ) {
            const std::vector< const Stub* > stubs = t.getStubs();
            std::map< unsigned int, unsigned int > layerMap;
            for ( auto s : stubs )
            layerMap[ s->layerIdReduced() ]++;
            for ( auto l : layerMap )
            hisNumStubsOnLayer_.at(TrackCands)->Fill( l.second );
        }
    }
    
    //=== Count track candidates found in the tracker.
    
    unsigned int nTracks = 0;
    unsigned int nTracks_sf = 0;
    const unsigned int numPhiOctants = 8;
    vector<unsigned int> nTracksInOctant(numPhiOctants, 0);
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
        unsigned int nTracksInEtaReg = 0;
        for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
            const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            hisNumTrksPerSect_.at(TrackCands)->Fill(htPair.numTrackCands3D());
            unsigned int iOctant = floor(iPhiSec*numPhiOctants/(settings_->numPhiSectors())); // phi octant number
            nTracksInOctant[iOctant] += htPair.numTrackCands3D();
            // Note number of tracks in this eta region (summed over all phi).
            nTracksInEtaReg += htPair.numTrackCands3D();
            nTracks_sf += htPair.getRphiHT().numTrackCands2D();
            if (settings_->debug() == 1 && htPair.numTrackCands3D() > 0) cout<<"Sector ("<<iPhiSec<<","<<iEtaReg<<") has ntracks = "<<htPair.numTrackCands3D()<<endl;
        }
        nTracks += nTracksInEtaReg;
        profNumTracksVsEta_.at(TrackCands)->Fill(iEtaReg, nTracksInEtaReg);
    }
    profNumTrackCands_.at(TrackCands)->Fill(1.0, nTracks); // Plot mean number of tracks/event.
    profNumTrackCands_.at(TrackCands)->Fill(8.0, nTracks_sf);
    for (unsigned int k = 0; k < numPhiOctants; k++) {
        hisNumTrksPerOct_.at(TrackCands)->Fill(nTracksInOctant[k]); // Plots tracks in each phi octant.
    }
    if(settings_->useSeedFilter()) cout << " number of tracks prior to seed filter = "<< nTracks_sf << endl;
    cout<<"Number of tracks prior to track helix fit = "<<nTracks<<endl;
    
    //=== Count stubs per event assigned to track candidates in the Tracker
    
    unsigned int nStubsOnTracks = 0;
    vector<unsigned int> nStubsOnTracksInOctant(numPhiOctants, 0);
    map< unsigned int, set<const Stub*> > uniqueStubsOnTracksInOctant;
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
        unsigned int nStubsOnTracksInEtaReg = 0;
        for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
            unsigned int iOctant = floor(iPhiSec*numPhiOctants/(settings_->numPhiSectors())); // phi octant number
            const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            unsigned int nStubsOnTrksInSec = htPair.numStubsOnTrackCands3D();
            hisStubsOnTracksPerSect_.at(TrackCands)->Fill(nStubsOnTrksInSec); // Number of stubs assigned to tracks in this sector.
            nStubsOnTracksInOctant[iOctant] += nStubsOnTrksInSec; // Number of stubs assigned to tracks in this octant.
            nStubsOnTracksInEtaReg += nStubsOnTrksInSec;
            set<const Stub*> uniqueStubsOnTracksInSector;
            // Loop over all stubs on all tracks in this sector, and add to std::set(), so each individual stub recorded at most once.
            for (const L1track3D& trk : htPair.trackCands3D() ) {
                const vector<const Stub*>& stubs = trk.getStubs();
                uniqueStubsOnTracksInSector.insert(stubs.begin(), stubs.end());
                uniqueStubsOnTracksInOctant[iOctant].insert(stubs.begin(), stubs.end());
            }
            // Plot number of stubs assigned to tracks per sector, never counting each individual stub more than once.
            hisUniqueStubsOnTrksPerSect_.at(TrackCands)->Fill(uniqueStubsOnTracksInSector.size());
        }
        nStubsOnTracks += nStubsOnTracksInEtaReg;
        profStubsOnTracksVsEta_.at(TrackCands)->Fill(iEtaReg, nStubsOnTracksInEtaReg);
    }
    profStubsOnTracks_.at(TrackCands)->Fill(1.0, nStubsOnTracks);
    for (unsigned int k = 0; k < numPhiOctants; k++) {
        hisStubsOnTracksPerOct_.at(TrackCands)->Fill(nStubsOnTracksInOctant[k]); // Plots stubs on tracks in each phi octant.
        hisUniqueStubsOnTrksPerOct_.at(TrackCands)->Fill(uniqueStubsOnTracksInOctant[k].size());
    }
    
    // Plot q/pt spectrum of track candidates, and number of stubs/track.
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
            const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            // Loop over all reconstructed tracks in this sector
            const vector<L1track3D>& vecTrk3D = htPair.trackCands3D();
            for (const L1track3D& trk : vecTrk3D) {
                hisNumTracksVsQoverPt_.at(TrackCands)->Fill(trk.qOverPt()); // Plot reconstructed q/Pt of track cands.
                hisStubsPerTrack_.at(TrackCands)->Fill(trk.getNumStubs());  // Stubs per track.
                // For genuine tracks, check how often they have too many stubs to be stored in cell memory. (Perhaps worse for high Pt particles in jets?).
                const TP* tp = trk.getMatchedTP();
                if (tp != nullptr) {
                    if (tp->useForAlgEff()) profExcessStubsPerTrackVsPt_.at(TrackCands)->Fill(1./tp->pt(), trk.getNumStubs() > 16);
                }
                hisLayersPerTrack_.at(TrackCands)->Fill(trk.getNumLayers()); // Number of reduced layers with stubs per track.
                hisPSLayersPerTrack_.at(TrackCands)->Fill( Utility::countLayers(settings_, trk.getStubs(), false, true) ); // Number of reduced PS layers with stubs per track.
                // Also plot just for genuine tracks.
                if (tp != nullptr && tp->useForAlgEff()) {
                    hisLayersPerTrueTrack_.at(TrackCands)->Fill(trk.getNumLayers()); // Number of reduced layers with stubs per track.
                    hisPSLayersPerTrueTrack_.at(TrackCands)->Fill( Utility::countLayers(settings_, trk.getStubs(), false, true) ); // Number of reduced PS layers with stubs per track.
                }
            }
        }
    }
    
    // Count fraction of stubs on each track matched to a TP that are from same TP.
    
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
            const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
            // Loop over all reconstructed tracks in this sector
            const vector<L1track3D>& vecTrk3D = htPair.trackCands3D();
            for (const L1track3D& trk : vecTrk3D) {
                // Only consider tracks that match a tracking particle used for the alg. efficiency measurement.
                const TP* tp = trk.getMatchedTP();
                if (tp != nullptr) {
                    if (tp->useForAlgEff()) {
                        hisFracMatchStubsOnTracks_.at(TrackCands)->Fill( trk.getPurity() );
                        const vector<const Stub*> stubs = trk.getStubs();
                        for (const Stub* s : stubs) {
                            // Was this stub produced by correct truth particle?
                            const set<const TP*> stubTPs = s->assocTPs();
                            bool trueStub = (stubTPs.find(tp) != stubTPs.end());
                            
                            // Distance of stub from true trajectory in z (barrel) or r (endcap)
                            float deltaRorZ =  s->barrel()  ?  (s->z() - tp->trkZAtStub( s ))  :  (s->r() - tp->trkRAtStub( s ));
                            // Distance of stub from true trajectory in r*phi.
                            float deltaPhiR  = s->r() * reco::deltaPhi(s->phi(), tp->trkPhiAtStub(s));
                            
                            // Nasty correction to stub phi coordinate to take into account non-radial strips in endcap 2S modules.
                            // Largely taken from TrackFitLinearAlgo.cc
                            float phiCorr = 0.;
                            float stripAngle = 999.;
                            if ( ! (s->barrel() || s->psModule()) ) {
                                float fracPosInModule = (float(2 * s->iphi()) - float(s->nstrip())) / float(s->nstrip());
                                stripAngle = 0.5 * s->width() * fracPosInModule / s->r();
                                if (s->z() > 0.) stripAngle *= -1.;
                                phiCorr = (tp->trkRAtStub(s) - s->r()) * stripAngle;
                            }
                            deltaPhiR += phiCorr;
                            
                            if (trueStub) {
                                if (s->psModule()) {
                                    hisDeltaPhiRtruePS_.at(TrackCands)->Fill(deltaPhiR);
                                    hisDeltaRorZtruePS_.at(TrackCands)->Fill(deltaRorZ);
                                } else {
                                    //if (tp->pt() > 20. && s->assocTP() != nullptr && fabs(tp->d0()) < 0.01) {
                                    //		    if ( ! (s->barrel()) ) cout<<"RATS "<<s->width()<<" "<<s->iphi()<<" "<<s->nstrip()<<" "<<s->r()<<" "<<s->z()<<endl;
                                    //if ( ! (s->barrel()) ) cout<<"DELTAPHI "<<stripAngle<<" "<<(tp->trkRAtStub(s) - s->r())<<" "<<phiCorr<<" : "<<deltaPhiR<<" "<<deltaPhiR1<<" "<<deltaPhiR2<<endl;
                                    //}
                                    hisDeltaPhiRtrue2S_.at(TrackCands)->Fill(deltaPhiR);
                                    hisDeltaRorZtrue2S_.at(TrackCands)->Fill(deltaRorZ);
                                }
                                // More detailed plots for true stubs to study effect of multiple scattering.
                                float sigPerp = s->sigmaPerp(); // detector resolution
                                float ptThresh = 40.; // threshold where scattering dominates detector resolution
                                float relpos = s->barrel()  ?   s->r() / settings_->trackerOuterRadius()  :  fabs(s->z()) / settings_->trackerHalfLength();
                                float sigmaScat = 0.01 * (ptThresh/tp->pt()) * pow(relpos, 1.5);
                                sigPerp += sigmaScat; // Estimated resolution allowing for scattering.
                                profNsigmaPhiRvsInvPt_.at(TrackCands)->Fill(1./tp->pt(), fabs(deltaPhiR)/sigPerp);
                                profNsigmaPhiRvsFracDist_.at(TrackCands)->Fill(relpos,   fabs(deltaPhiR)/sigPerp);
                            } else {
                                if (s->psModule()) {
                                    hisDeltaPhiRfakePS_.at(TrackCands)->Fill(deltaPhiR);
                                    hisDeltaRorZfakePS_.at(TrackCands)->Fill(deltaRorZ);
                                } else {
                                    hisDeltaPhiRfake2S_.at(TrackCands)->Fill(deltaPhiR);
                                    hisDeltaRorZfake2S_.at(TrackCands)->Fill(deltaRorZ);
                                }
                            }
                            
                            // Fraction of wrong stubs vs. tracker layer.
                            profFracTrueStubsVsLayer_.at(TrackCands)->Fill(s->layerId(), trueStub);
                            
                            // Check how much stub bend differs from predicted one, relative to nominal bend resolution.
                            float diffBend = (s->qOverPt() - trk.qOverPt()) / s->qOverPtOverBend();
                            if (trueStub) {
                                hisDeltaBendTrue_.at(TrackCands)->Fill(diffBend/s->bendRes());
                            } else {
                                hisDeltaBendFake_.at(TrackCands)->Fill(diffBend/s->bendRes());
                            }
                            
                            // Debug printout to understand for matched tracks, how far stubs lie from true particle trajectory
                            // Only prints for tracks with huge number of stubs, to also understand why these tracks exist.
                            //if (trk.getNumStubs() > 20) {
                            /*
                             if (trk.pt() > 20) {
                             cout<<"--- Checking how far stubs on matched tracks lie from true particle trajectory. ---"<<endl;
                             cout<<"    Track "<<trk.getPurity()<<" "<<tp->pt()<<" "<<tp->d0()<<endl;
                             float sigPhiR = deltaPhiR/s->sigmaPerp();
                             float sigRorZ = deltaRorZ/s->sigmaPar();
                             string ohoh =  (fabs(sigPhiR) > 5 || fabs(sigRorZ) > 5)  ?  "FAR"  :  "NEAR";
                             if (trueStub) {
                             cout<<"    Real stub "<<ohoh<<" ps="<<s->psModule()<<" bar="<<s->barrel()<<" lay="<<s->layerId()<<" : phi="<<deltaPhiR<<" ("<<sigPhiR<<") rz="<<deltaRorZ<<" ("<<sigRorZ<<")"<<endl;
                             } else {
                             cout<<"    FAKE stub "<<ohoh<<" ps="<<s->psModule()<<" bar="<<s->barrel()<<" lay="<<s->layerId()<<" : phi="<<deltaPhiR<<" ("<<sigPhiR<<") rz="<<deltaRorZ<<" ("<<sigRorZ<<")"<<endl;
                             }
                             cout<<"        coords="<<s->r()<<" "<<s->phi()<<" "<<s->eta()<<" bend="<<s->bend()<<" iphi="<<s->iphi()<<endl;
                             cout<<"        module="<<s->minR()<<" "<<s->minPhi()<<" "<<s->minZ()<<endl;
                             }
                             */
                        }
                    }
                }
            }
        }
    }
    
    // Count total number of tracking particles in the event that were reconstructed,
    // counting also how many of them were reconstructed multiple times (duplicate tracks).
    
    unsigned int nRecoedTPsForEff = 0; // Total no. of TPs used for the efficiency measurement that were reconstructed as at least one track.
    unsigned int nRecoedTPs = 0; // Total no. of TPs that were reconstructed as at least one track.
    unsigned int nEtaSecsMatchingTPs = 0; // Total no. of eta sectors that all TPs were reconstructed in
    unsigned int nSecsMatchingTPs = 0; // Total no. of eta x phi sectors that all TPs were reconstructed in
    unsigned int nTrksMatchingTPs = 0; // Total no. of tracks that all TPs were reconstructed as
    unsigned int nTrksMatchingTPsIgnoringRzDups = 0; // Ditto, but if TP reconstructed in multiple cells of r-z HT, just count them as 1.
    
    for (const TP& tp: vTPs) {
        
        bool tpRecoed = false;
        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
            bool tpRecoedInEtaSec = false;
            for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
                // Get reconstructed tracks in this sector corresponding to this TP.
                const vector<const L1track3D*> trks = htPair.assocTrackCands3D( tp );
                // Count them
                unsigned int nTrk = trks.size();
                // Same again, but if TP reconstructed in multiple cells of r-z HT, just count them as 1.
                // (If no r-z HT used, then nCellsRphi will equal nTrk).
                unsigned int nCellsRphi = htPair.numRphiCells( trks );
                if (nTrk > 0) {
                    tpRecoed = true;            // This TP was reconstructed at least once in tracker.
                    tpRecoedInEtaSec = true;    // This TP was reconstructed at least once in this eta sector.
                    nSecsMatchingTPs += 1;      // Increment sum by no. of sectors this TP was reconstructed in
                    nTrksMatchingTPs += nTrk; // Increment sum by no. of tracks this TP was reconstructed as
                    nTrksMatchingTPsIgnoringRzDups += nCellsRphi; // Ditto, but if TP reconstructed in multiple cells of r-z HT, just count them as 1.
                    profDupTracksVsTPeta_.at(TrackCands)->Fill(tp.eta(), nTrk); // Study duplication of tracks within an individual HT array.
                }
            }
            if (tpRecoedInEtaSec) nEtaSecsMatchingTPs++; // Increment each time TP found in an eta sector.
        }
        
        if (tpRecoed) {
            // Increment sum each time a TP is reconstructed at least once inside Tracker
            if (tp.useForEff()) nRecoedTPsForEff++;
            nRecoedTPs++;
        }
    }
    
    //--- Plot mean number of tracks/event, counting number due to different kinds of duplicates
    
    // Plot number of TPs used for the efficiency measurement that are reconstructed.
    profNumTrackCands_.at(TrackCands)->Fill(7.0, nRecoedTPsForEff);
    // Plot number of TPs that are reconstructed.
    profNumTrackCands_.at(TrackCands)->Fill(6.0, nRecoedTPs);
    // Plot number of TPs that are reconstructed. (Count +1 for each eta sector they are reconstructed in).
    profNumTrackCands_.at(TrackCands)->Fill(5.0, nEtaSecsMatchingTPs);
    // Plot number of TPs that are reconstructed. (Count +1 for each etaxphisector they are reconstructed in).
    profNumTrackCands_.at(TrackCands)->Fill(4.0, nSecsMatchingTPs);
    // Plot number of TP that are reconstructed. (Ditto, but now multiplying by duplicate cells in r-phi HT).
    profNumTrackCands_.at(TrackCands)->Fill(3.0, nTrksMatchingTPsIgnoringRzDups);
    // Plot number of TP that are reconstructed. (Ditto, but now multiplying by duplicate cells in r-phi x r-z HTs).
    profNumTrackCands_.at(TrackCands)->Fill(2.0, nTrksMatchingTPs);
    
    //=== Study tracking efficiency by looping over tracking particles.
    
    for (const TP& tp: vTPs) {
        
        if (tp.useForEff()) { // Check TP is good for efficiency measurement.
            
            // Check which eta and phi sectors this TP is in.
            int iPhiSec_TP = -1;
            int iEtaReg_TP = -1;
            for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                const Sector& sector = mSectors(iPhiSec, 0);
                if (sector.insidePhiSec(tp)) iPhiSec_TP = iPhiSec;
            }
            for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
                const Sector& sector = mSectors(0, iEtaReg);
                if (sector.insideEtaReg(tp)) iEtaReg_TP = iEtaReg;
            }
            
            // Plot kinematics of all good TP.
            hisTPinvptForEff_.at(TrackCands)->Fill(1./tp.pt());
            hisTPetaForEff_.at(TrackCands)->Fill(tp.eta());
            hisTPphiForEff_.at(TrackCands)->Fill(tp.phi0());
            
            if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
                hisTPinvptForAlgEff_.at(TrackCands)->Fill(1./tp.pt());
                hisTPetaForAlgEff_.at(TrackCands)->Fill(tp.eta());
                hisTPphiForAlgEff_.at(TrackCands)->Fill(tp.phi0());
                // Plot also production point of all good TP.
                hisTPd0ForAlgEff_.at(TrackCands)->Fill(fabs(tp.d0()));
                hisTPz0ForAlgEff_.at(TrackCands)->Fill(fabs(tp.z0()));
                // Plot sector nunber.
                hisTPphisecForAlgEff_.at(TrackCands)->Fill(iPhiSec_TP);
                hisTPetasecForAlgEff_.at(TrackCands)->Fill(iEtaReg_TP);
            }
            
            // Check if this TP was reconstructed anywhere in the tracker..
            bool tpRecoed = false;
            bool tpRecoedPerfect = false;
            for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
                    const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
                    if (htPair.numAssocTrackCands3D( tp ) > 0) {
                        tpRecoed = true;
                        // Also note if TP was reconstructed perfectly (no incorrect hits on reco track).
                        const vector<const L1track3D*> assTrkVec = htPair.assocTrackCands3D( tp );
                        for (const L1track3D* assTrk : assTrkVec) {
                            if (assTrk->getPurity() == 1.) tpRecoedPerfect = true; 
                        }
                    }
                }
            }
            
            // Count perfectly reconstructed TP (no incorrect hits on reconstructed track) used for alg. effi. measurement.
            if (tpRecoedPerfect && tp.useForAlgEff()) numPerfRecoTPforAlg_++;
            
            // If TP was reconstucted by HT, then plot its kinematics.
            if (tpRecoed) {
                hisRecoTPinvptForEff_.at(TrackCands)->Fill(1./tp.pt());
                hisRecoTPetaForEff_.at(TrackCands)->Fill(tp.eta());
                hisRecoTPphiForEff_.at(TrackCands)->Fill(tp.phi0());
                // Also plot efficiency to perfectly reconstruct the track (no fake hits)
                if (tpRecoedPerfect) {
                    hisPerfRecoTPinvptForEff_.at(TrackCands)->Fill(1./tp.pt());
                    hisPerfRecoTPetaForEff_.at(TrackCands)->Fill(tp.eta());
                }
                if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
                    hisRecoTPinvptForAlgEff_.at(TrackCands)->Fill(1./tp.pt());
                    hisRecoTPetaForAlgEff_.at(TrackCands)->Fill(tp.eta());
                    hisRecoTPphiForAlgEff_.at(TrackCands)->Fill(tp.phi0());
                    // Plot also production point of all good reconstructed TP.
                    hisRecoTPd0ForAlgEff_.at(TrackCands)->Fill(fabs(tp.d0()));
                    hisRecoTPz0ForAlgEff_.at(TrackCands)->Fill(fabs(tp.z0()));
                    // Plot sector number to understand if looser cuts are needed in certain eta regions.
                    hisRecoTPphisecForAlgEff_.at(TrackCands)->Fill(iPhiSec_TP);
                    hisRecoTPetasecForAlgEff_.at(TrackCands)->Fill(iEtaReg_TP);
                    // Also plot efficiency to perfectly reconstruct the track (no fake hits)
                    if (tpRecoedPerfect) {
                        hisPerfRecoTPinvptForAlgEff_.at(TrackCands)->Fill(1./tp.pt());
                        hisPerfRecoTPetaForAlgEff_.at(TrackCands)->Fill(tp.eta());
                        hisPerfRecoTPphisecForAlgEff_.at(TrackCands)->Fill(iPhiSec_TP);
                        hisPerfRecoTPetasecForAlgEff_.at(TrackCands)->Fill(iEtaReg_TP);
                    }
                }
            }
        }
    }
    
    // Histos of track parameter resolution
    
    for (const TP& tp: vTPs) {
        
        if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
            
            // For each tracking particle, find the corresponding reconstructed track(s).
            for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
                    const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
                    const vector<const L1track3D*> trkVec = htPair.assocTrackCands3D( tp );
                    for (const L1track3D* trk : trkVec) {
                        hisQoverPtRes_.at(TrackCands)->Fill(trk->qOverPt() - tp.qOverPt());
                        hisPhi0Res_.at(TrackCands)->Fill(reco::deltaPhi(trk->phi0(), tp.phi0()));
                        hisEtaRes_.at(TrackCands)->Fill(trk->eta() - tp.eta());
                        hisZ0Res_.at(TrackCands)->Fill(trk->z0() - tp.z0());
                    }
                }
            }
        }
    }
    
    // Diagnose reason why not all viable tracking particles were reconstructed.
    const map<const TP*, string> diagnosis = this->diagnoseTracking(inputData, mSectors, mHtPairs);
    for (const auto& iter: diagnosis) {
        hisRecoFailureReason_.at(TrackCands)->Fill(iter.second.c_str(), 1.); // Stores flag indicating failure reason.
    }
    
    
    
    //CERN
     ////// SimplelLR ////////
//    const std::vector<std::pair<std::string,L1fittedTrack>> fittedTracks;
    
    map<std::string,uint> nFittedTracks;
    map<std::string,uint> nStubsOnTrack;
    map<std::string,uint> nTracksGenuine;
    map<std::string,uint> nTracksGenuineTP;
    map<std::string,uint> nTracksGenuinePass;
    map<std::string,uint> nTracksFakeCut;
    map<std::string,uint> nTracksFake;
    map<std::string,uint> nRejectedTracks;
    map<std::string,uint> nRejectedFake;
    map<std::string,uint> nTracksExcDups;
    map<std::string,uint> nTracksExcDupsPass;
    
    for(auto &j : settings_->trackFitters() ){
        nFittedTracks[j] = 0;
        nStubsOnTrack[j] = 0;
        nTracksGenuine[j] = 0;
        nTracksGenuineTP[j] = 0;
        nTracksGenuinePass[j] = 0;
        nTracksFakeCut[j] = 0;
        nTracksFake[j] = 0;
        nRejectedTracks[j] = 0;
        nRejectedFake[j] = 0;
        nTracksExcDups[j] = 0;
        nTracksExcDupsPass[j] = 0;
    }
    
    
    // Do track fitting algorithmic efficiences
    for (const TP& tp: vTPs) {
        
        // Setup bools re. duplicate checks
        std::map<std::string, bool> tpRecoed;
        std::map<std::string, bool> tpRecoedPass;
        // Setup bools re. tp efficiencies
        std::map<std::string, bool> tpRecoedForAlgEff;
        std::map<std::string, bool> tpRecoedForAlgEffPerfect;
        std::map<std::string, bool> tpRecoedForAlgEffPassCut;
        std::map<std::string, bool> tpRecoedForAlgEffPerfectPassCut;
        
        // Init bools for all track fit algos to false
        for( auto &j : settings_->trackFitters() ){
            tpRecoed[j] = false;
            tpRecoedPass[j] = false;
            tpRecoedForAlgEff[j] = false;
            tpRecoedForAlgEffPerfect[j] = false;
            tpRecoedForAlgEffPassCut[j] = false;
            tpRecoedForAlgEffPerfectPassCut[j] = false;
        }
        
        // Only consider TP useful for algorithmic efficeincy.
        
        for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : fittedTracks ){
            
            std::string algoName(fittedTrackPair.first); // Get fitting algo name
            const L1fittedTrack& fitTrk = fittedTrackPair.second; // Get fitted track
            //IRT
            //      const TP*   assocTP =  fitTrk.getL1track3D().getMatchedTP(); // Get the TP the fitted track matches to, if any.
            const TP*   assocTP =  fitTrk.getMatchedTP(); // Get the TP the fitted track matches to, if any.
            
            if ( !fitTrk.accepted() ) continue;
            
            // Does this match desired TP?
            if (assocTP == &tp)  tpRecoed[algoName] = true;
            if (assocTP == &tp && fitTrk.chi2dof() <= settings_->chi2OverNdfCut()) tpRecoedPass[algoName] = true;
            
            // Does this match desired TP & useForAlgEff
            if (assocTP == &tp && tp.useForAlgEff()) {
                tpRecoedForAlgEff[algoName] = true; // If it does, set tpRecoed to true
                // Also note if TP was reconstructed perfectly (no incorrect hits on reco track).
                if (fitTrk.getPurity() == 1.) tpRecoedForAlgEffPerfect[algoName] = true;
                if (fitTrk.chi2dof() <= settings_->chi2OverNdfCut()) tpRecoedForAlgEffPassCut[algoName] = true; // Does track with assoc. TP pass the cut
                if (fitTrk.getPurity() == 1. && fitTrk.chi2dof() <= settings_->chi2OverNdfCut()) tpRecoedForAlgEffPerfectPassCut[algoName] = true; // Is the track passing the cut perfectly recoed?
            }
        }
        
        // See if any of the TPs had a corresponding fitting track and if it/they was/were pure?
        for(auto &j : settings_->trackFitters() ){
            if( tpRecoed[j] ) nTracksExcDups[j]++;                            // Number of tracks excluding duplicates
            if( tpRecoedPass[j] ) nTracksExcDupsPass[j]++;                    // Number of tracks excluding duplicates that passed cut
            if (tpRecoedForAlgEff[j]) numFitAlgEff_[j]++;			// Algo efficiency before cut
            if (tpRecoedForAlgEffPerfect[j]) numFitPerfAlgEff_[j]++;		// Perfect algo efficiency before cut
            if (tpRecoedForAlgEffPassCut[j]) numFitAlgEffPass_[j]++;		// Algo efficiency after cut
            if (tpRecoedForAlgEffPerfectPassCut[j]) numFitPerfAlgEffPass_[j]++; // Perfect algo efficiency after cut
        }
    }
    
    for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : fittedTracks ){
        
        std::string j (fittedTrackPair.first);
        const L1fittedTrack& fitTrk = fittedTrackPair.second;
        
        // Get original HT track candidate prior to fit for comparison.
        const L1track3D& htTrk = fitTrk.getL1track3D();
        
        // Get matched truth particle, if any.
        const TP* tp = fitTrk.getMatchedTP();
        
        // Increment nFittedTrack and nStubsOnTrack counters
        nFittedTracks[j] += 1;
        if( fitTrk.chi2dof() <= settings_->chi2OverNdfCut() && fitTrk.accepted() ) nStubsOnTrack[j] += fitTrk.getNumStubs();
        
        if ( fitTrk.accepted() ) {
            
            nTracksGenuine[j] += 1;
            if (tp != nullptr) {
                nTracksGenuineTP[j] += 1;
                hisFitVsSeedQinvPtGenCand_[j].at(SimpleLR)->Fill( htTrk.qOverPt(), fitTrk.qOverPt() );
                hisFitVsSeedPhi0GenCand_[j].at(SimpleLR)->Fill( htTrk.phi0(), fitTrk.phi0() );
                hisFitVsSeedD0GenCand_[j].at(SimpleLR)->Fill( htTrk.d0(), fitTrk.d0() );
                hisFitVsSeedZ0GenCand_[j].at(SimpleLR)->Fill( htTrk.z0(), fitTrk.z0() );
                hisFitVsSeedEtaGenCand_[j].at(SimpleLR)->Fill( htTrk.eta(), fitTrk.eta() );
                
                if ( fitTrk.chi2dof() <= settings_->chi2OverNdfCut() ){
                    nTracksGenuinePass[j] += 1;
                }
                
                // Check if chi2/NDF is well behaved for perfectly reconstructed tracks.
                if (fitTrk.getPurity() == 1.) {
                    profChi2DofVsInvPtPERF_[j].at(SimpleLR)->Fill(fabs(tp->qOverPt()), sqrt(fitTrk.chi2dof()));
                    profBigChi2DofVsInvPtPERF_[j].at(SimpleLR)->Fill(fabs(tp->qOverPt()), (fitTrk.chi2dof() > 10));
                    // Are high Pt tracks sensitive to d0 impact parameter?
                    if (tp->pt() > 10.) {
                        if (fitTrk.chi2dof() > 10.) {
                            hisD0TPBigChi2DofPERF_[j].at(SimpleLR)->Fill(fabs(tp->d0()));
                        } else {
                            hisD0TPSmallChi2DofPERF_[j].at(SimpleLR)->Fill(fabs(tp->d0()));
                        }
                    }
                }
            }
            
            if (tp == nullptr){
                nTracksFake[j] += 1;
                if ( fitTrk.chi2dof() > settings_->chi2OverNdfCut() ){
                    nTracksFakeCut[j] += 1;
                }
                hisFitVsSeedQinvPtFakeCand_[j].at(SimpleLR)->Fill( htTrk.qOverPt(), fitTrk.qOverPt() );
                hisFitVsSeedPhi0FakeCand_[j].at(SimpleLR)->Fill( htTrk.phi0(), fitTrk.phi0() );
                hisFitVsSeedD0FakeCand_[j].at(SimpleLR)->Fill( htTrk.d0(), fitTrk.d0() );
                hisFitVsSeedZ0FakeCand_[j].at(SimpleLR)->Fill( htTrk.z0(), fitTrk.z0() );
                hisFitVsSeedEtaFakeCand_[j].at(SimpleLR)->Fill( htTrk.eta(), fitTrk.eta() );
            }
            
        }
        else if ( !fitTrk.accepted() ) {
            nRejectedTracks[j] += 1;
            if ( tp == nullptr ) nRejectedFake[j] += 1;
        }
        
        if ( fitTrk.accepted() && fitTrk.chi2dof() < settings_->chi2OverNdfCut() ){
            // Study fake tracks.
            bool fakeTrk = (tp == nullptr);
            profFakeFitTrksVsEta_[j].at(SimpleLR)->Fill(fabs(fitTrk.eta()), fakeTrk);
            // Study incorrect hits on good tracks.
            if ( ! fakeTrk) {
                const vector<const Stub*> stubs = fitTrk.getStubs();
                for (const Stub* s : stubs) {
                    // Was this stub produced by correct truth particle?
                    const set<const TP*> stubTPs = s->assocTPs();
                    bool trueStub = (stubTPs.find(tp) != stubTPs.end());
                    profFitFracTrueStubsVsLayer_[j].at(SimpleLR)->Fill(s->layerId(), trueStub);
                    profFitFracTrueStubsVsEta_[j].at(SimpleLR)->Fill(fabs(s->eta()), trueStub);
                }
            }
        }
        
        // Get matched truth particle to (pre-fit) HT track, if any.
        const TP* tpHT = fitTrk.getL1track3D().getMatchedTP();
        // Count fraction of matched (unmatched) tracks killed by track fit (because lost too many hits).
        unsigned int ibin = ( tpHT != nullptr) ? 1 : 2; // Set histogram bin to 1 if track matches TP, or 2 if not.
        
        // --- Study effect of the track fitter killing stubs with large residuals.
        // Count stubs per track removed track by fit (because they had large residuals),
        // distinguishing those which were good stubs (matched the best TP).
        hisNumMatchedStubsKilledVsKilled_[j].at(SimpleLR)->Fill( fitTrk.getNumKilledStubs(), fitTrk.getNumKilledMatchedStubs() );
        
        // --- Study purity against number of stubs to augment track fitter killing stubs due to large residual study
        hisNumStubsVsPurity_[j].at(SimpleLR)->Fill( fitTrk.getNumStubs(), fitTrk.getPurity() );
        
        // Old way before accepted Kalman came along.
        profTrksKilledByFit_[j].at(SimpleLR)->Fill(ibin, !fitTrk.accepted());
        
        // Histograms specific to Linear Regression track fitter.
        if (j.find("LinearRegression") != string::npos) {
            
            int numIterations;
            std::string lostMatchingState;
            std::unordered_map< std::string, int > stateCalls;
            fitTrk.getInfoLR( numIterations, lostMatchingState, stateCalls );
            hisNumIterations_[j].at(SimpleLR)->Fill( numIterations );
            if ( tp != tpHT )
            hisFailingState_[j].at(SimpleLR)->Fill( lostMatchingState.c_str(), 1. );
            int totalCalls = 0;
            for ( auto state : stateCalls ) {
                totalCalls += state.second;
                hisTotalStateCalls_[ j ].at(SimpleLR)->Fill( state.first.c_str(), state.second );
            }
            for ( auto state : stateCalls )
            hisRelativeStateCalls_[ j ].at(SimpleLR)->Fill( state.first.c_str(), state.second / float(totalCalls) );
        }
        
        hisNumStubsFitKills_[j].at(SimpleLR)->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs() );
        hisNumStubsFitKillsVsPurity_[j].at(SimpleLR)->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs(), fitTrk.getPurity() );
        if (fitTrk.accepted()) hisNumStubsFitKillsVsPurityMatched_[j].at(SimpleLR)->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs(), fitTrk.getPurity() );
        else hisNumStubsFitKillsVsPurityUnmatched_[j].at(SimpleLR)->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs(), fitTrk.getPurity() );
        
        if ( !fitTrk.accepted() ) continue; // If a rejected track, do not make plots for these.
        
        // Fill fitted parameter histograms and histograms of seed parameters against fitted parameters
        
        // Seed track parameter distributions
        hisSeedQinvPt_[j].at(SimpleLR)->Fill( htTrk.qOverPt() );
        hisSeedPhi0_[j].at(SimpleLR)->Fill( htTrk.phi0() );
        hisSeedD0_[j].at(SimpleLR)->Fill( htTrk.d0() );
        hisSeedZ0_[j].at(SimpleLR)->Fill( htTrk.z0() );
        hisSeedEta_[j].at(SimpleLR)->Fill( htTrk.eta() );
        // Fitted track parameter distributions & chi2, separately for tracks that do/do not match a truth particle
        if ( tp != nullptr){
            
            hisFitQinvPtMatched_[j].at(SimpleLR)->Fill( fitTrk.qOverPt() );
            hisFitPhi0Matched_[j].at(SimpleLR)->Fill( fitTrk.phi0() );
            hisFitD0Matched_[j].at(SimpleLR)->Fill( fitTrk.d0() );
            hisFitZ0Matched_[j].at(SimpleLR)->Fill( fitTrk.z0() );
            hisFitEtaMatched_[j].at(SimpleLR)->Fill( fitTrk.eta() );
            
            hisFitChi2Matched_[j].at(SimpleLR)->Fill( fitTrk.chi2() );
            hisFitChi2DofMatched_[j].at(SimpleLR)->Fill( fitTrk.chi2dof() );
            
        } else {
            
            hisFitQinvPtUnmatched_[j].at(SimpleLR)->Fill( fitTrk.qOverPt() );
            hisFitPhi0Unmatched_[j].at(SimpleLR)->Fill( fitTrk.phi0() );
            hisFitD0Unmatched_[j].at(SimpleLR)->Fill( fitTrk.d0() );
            hisFitZ0Unmatched_[j].at(SimpleLR)->Fill( fitTrk.z0() );
            hisFitEtaUnmatched_[j].at(SimpleLR)->Fill( fitTrk.eta() );
            
            hisFitChi2Unmatched_[j].at(SimpleLR)->Fill( fitTrk.chi2() );
            hisFitChi2DofUnmatched_[j].at(SimpleLR)->Fill( fitTrk.chi2dof() );
            
        }
        
        // If there is an associated tracking particle, fill up histograms of TP's parameters against fitted parameters
        if ( tp != nullptr ){
            
            // Do seperately for those with good/poor chi2.
            if ( fitTrk.chi2dof() <= settings_->chi2OverNdfCut() ){
                // Fitted vs True parameter distribution 2D plots
                hisFitVsTrueQinvPtGoodChi2_[j].at(SimpleLR)->Fill( tp->qOverPt(), fitTrk.qOverPt() );
                hisFitVsTruePhi0GoodChi2_[j].at(SimpleLR)->Fill( tp->phi0(), fitTrk.phi0( ));
                hisFitVsTrueD0GoodChi2_[j].at(SimpleLR)->Fill( tp->d0(), fitTrk.d0() );
                hisFitVsTrueZ0GoodChi2_[j].at(SimpleLR)->Fill( tp->z0(), fitTrk.z0() );
                hisFitVsTrueEtaGoodChi2_[j].at(SimpleLR)->Fill( tp->eta(), fitTrk.eta() );
                // Residuals between fitted and true helix params as 1D plot.
                hisFitQinvPtResGoodChi2_[j].at(SimpleLR)->Fill( fitTrk.qOverPt() - tp->qOverPt());
                hisFitPhi0ResGoodChi2_[j].at(SimpleLR)->Fill( reco::deltaPhi(fitTrk.phi0(), tp->phi0()) );
                hisFitD0ResGoodChi2_[j].at(SimpleLR)->Fill( fitTrk.d0() - tp->d0() );
                hisFitZ0ResGoodChi2_[j].at(SimpleLR)->Fill( fitTrk.z0() - tp->z0() );
                hisFitEtaResGoodChi2_[j].at(SimpleLR)->Fill( fitTrk.eta() - tp->eta() );
                // Residuals between true and seed helix params as 1D plot.
                hisSeedQinvPtResGoodChi2_[j].at(SimpleLR)->Fill( tp->qOverPt() - htTrk.qOverPt());
                hisSeedPhi0ResGoodChi2_[j].at(SimpleLR)->Fill( reco::deltaPhi(tp->phi0(), htTrk.phi0()) );
                hisSeedD0ResGoodChi2_[j].at(SimpleLR)->Fill( tp->d0() - htTrk.d0() );
                hisSeedZ0ResGoodChi2_[j].at(SimpleLR)->Fill( tp->z0() - htTrk.z0() );
                hisSeedEtaResGoodChi2_[j].at(SimpleLR)->Fill( tp->eta() - htTrk.eta() );
                
                // Understand which matched tracks have good/bad chi2.
                hisTrueEtaMatchedGoodChi2_[j].at(SimpleLR)->Fill( tp->eta() );
                hisStubPurityMatchedGoodChi2_[j].at(SimpleLR)->Fill( fitTrk.getPurity() );
            } else {
                
                // Plot rapidity of matched tracks with bad chi2.
                hisTrueEtaMatchedBadChi2_[j].at(SimpleLR)->Fill( tp->eta() );
                hisStubPurityMatchedBadChi2_[j].at(SimpleLR)->Fill( fitTrk.getPurity() );
                if (fitTrk.getPurity() > 0.99) {

                }
            }
            
            // Plot helix parameter resolution against eta or Pt.
            
            if (fitTrk.chi2dof() <= settings_->chi2OverNdfCut()) {
                hisQoverPtResVsTrueEta_[j].at(SimpleLR)->Fill( std::abs(tp->eta()), std::abs( fitTrk.qOverPt() - tp->qOverPt() ) );
                hisPhi0ResVsTrueEta_[j].at(SimpleLR)->Fill( std::abs(tp->eta()), std::abs(reco::deltaPhi(fitTrk.phi0(), tp->phi0()) ) );
                hisEtaResVsTrueEta_[j].at(SimpleLR)->Fill( std::abs(tp->eta()), std::abs( fitTrk.eta() - tp->eta() ) );
                hisZ0ResVsTrueEta_[j].at(SimpleLR)->Fill( std::abs(tp->eta()), std::abs( fitTrk.z0() - tp->z0() ) );
                hisD0ResVsTrueEta_[j].at(SimpleLR)->Fill( std::abs(tp->eta()), std::abs( fitTrk.d0() - tp->d0() ) );
                
                hisQoverPtResVsTrueInvPt_[j].at(SimpleLR)->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.qOverPt() - tp->qOverPt() ) );
                hisPhi0ResVsTrueInvPt_[j].at(SimpleLR)->Fill( std::abs(tp->qOverPt()), std::abs(reco::deltaPhi(fitTrk.phi0(), tp->phi0()) ) );
                hisEtaResVsTrueInvPt_[j].at(SimpleLR)->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.eta() - tp->eta() ) );
                hisZ0ResVsTrueInvPt_[j].at(SimpleLR)->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.z0() - tp->z0() ) );
                hisD0ResVsTrueInvPt_[j].at(SimpleLR)->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.d0() - tp->d0() ) );
                
                hisQoverPtResVsTrueZ0_[j].at(SimpleLR)->Fill( std::abs(tp->z0()), std::abs( fitTrk.qOverPt()
                                                                              - tp->z0() ) );
                hisPhi0ResVsTrueZ0_[j].at(SimpleLR)->Fill( std::abs(tp->z0()), std::abs(reco::deltaPhi(
                                                                                          fitTrk.phi0(), tp->phi0()) ) );
                hisEtaResVsTrueZ0_[j].at(SimpleLR)->Fill( std::abs(tp->z0()), std::abs( fitTrk.eta() - tp->
                                                                          eta() ) );
                hisZ0ResVsTrueZ0_[j].at(SimpleLR)->Fill( std::abs(tp->z0()), std::abs( fitTrk.z0() - tp->
                                                                         z0() ) );
                hisD0ResVsTrueZ0_[j].at(SimpleLR)->Fill( std::abs(tp->z0()), std::abs( fitTrk.d0() - tp->
                                                                         d0() ) );
                
            }
            
            // Plot chi^2 vs eta, and # stubs vs eta.
            hisTrueFittedChiSquaredVsTrueEta_[j].at(SimpleLR)->Fill( fitTrk.chi2(), fitTrk.getL1track3D().eta() );
            hisTrueFittedChiSquaredDofVsTrueEta_[j].at(SimpleLR)->Fill( fitTrk.chi2dof(), fitTrk.getL1track3D().eta() );
            hisTrueFittedChiSquaredVsFittedEta_[j].at(SimpleLR)->Fill( fitTrk.chi2(), fitTrk.eta() );
            hisTrueFittedChiSquaredDofVsFittedEta_[j].at(SimpleLR)->Fill( fitTrk.chi2dof(), fitTrk.eta() );
            
            hisFittedChiSquaredFunctionOfStubs_[j].at(SimpleLR)->Fill( fitTrk.getStubs().size(), fitTrk.chi2() );
            hisFittedChiSquaredDofFunctionOfStubs_[j].at(SimpleLR)->Fill( fitTrk.getStubs().size(), fitTrk.chi2dof() );
            
            
            
        }
    }
    
    for ( const string& fitterName : settings_->trackFitters() ){
        
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(1.0, nFittedTracks[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(2.0, nStubsOnTrack[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(3.0, nTracksGenuine[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(4.0, nTracksGenuineTP[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(5.0, nTracksGenuinePass[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(6.0, nTracksFakeCut[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(7.0, nTracksFake[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(8.0, nRejectedTracks[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(9.0, nRejectedFake[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(10.0, nTracksExcDups[fitterName]);
        profNumFittedCands_[fitterName].at(SimpleLR)->Fill(11.0, nTracksExcDupsPass[fitterName]);
    }
    
    //=== Study duplicate tracks. IRT
    
    for (const TP& tp: vTPs) {
        
        for (auto &fitName : settings_->trackFitters()) {
            
            unsigned int nMatch = 0;
            vector<const L1fittedTrack*> fitTrkStore;
            
            for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : fittedTracks ){
                
                std::string j (fittedTrackPair.first);
                const L1fittedTrack& fitTrk = fittedTrackPair.second;
                
                if (j == fitName) {
                    
                    if (fitTrk.accepted() && fitTrk.chi2dof() < settings_->chi2OverNdfCut()) {
                        
                        // Does this fitted track match the truth particle?
                        const TP* assocTP = fitTrk.getMatchedTP();
                        if (assocTP != nullptr) {
                            if (assocTP->index() == tp.index()) {
                                nMatch++; // Count matches to this truth particle.
                                fitTrkStore.push_back(&fitTrk); // List of fitted tracks matching this truth particle.
                            }
                        }
                    }
                }
            }
            
            if (nMatch > 0) {
                bool dup = (nMatch > 1);
                profDupFitTrksVsEta_[fitName].at(SimpleLR)->Fill(fabs(tp.eta()), dup);
                profDupFitTrksVsInvPt_[fitName].at(SimpleLR)->Fill(fabs(tp.qOverPt()), dup);

            }
        }
    }
    
    //=== Study tracking efficiency by looping over tracking particles.
    
    for (auto &fitName : settings_->trackFitters()) {
        
        for (const TP& tp: vTPs) {
            
            if (tp.useForEff()) { // Check TP is good for efficiency measurement.
                
                // Check which phi & eta sectors this TP is in.
                int iEtaReg_TP = -1;
                int iPhiSec_TP = -1;
                for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                    Sector secTmp;
                    secTmp.init(settings_, iPhiSec, 0);
                    if (secTmp.insidePhiSec(tp)) iPhiSec_TP = iPhiSec;
                }
                for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
                    Sector secTmp;
                    secTmp.init(settings_, 0, iEtaReg);
                    if (secTmp.insideEtaReg(tp)) iEtaReg_TP = iEtaReg;
                }
                
                // Check if this TP was reconstructed anywhere in the tracker..
                bool tpRecoed = false;
                bool tpRecoedPerfect = false;
                
                for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : fittedTracks ){
                    
                    const std::string j (fittedTrackPair.first);
                    const L1fittedTrack& fitTrk = fittedTrackPair.second;
                    
                    if (j == fitName) {
                        
                        if (fitTrk.accepted() && fitTrk.chi2dof() < settings_->chi2OverNdfCut()) {
                            
                            const TP* assocTP = fitTrk.getMatchedTP();
                            if (assocTP != nullptr) {
                                if (assocTP->index() == tp.index()) {
                                    tpRecoed = true;
                                    if (fitTrk.getPurity() == 1.) tpRecoedPerfect = true;
                                }
                            }
                        }
                    }
                }
                
                // If TP was reconstucted by HT, then plot its kinematics.
                if (tpRecoed) {
                    hisFitTPinvptForEff_[fitName].at(SimpleLR)->Fill(1./tp.pt());
                    hisFitTPetaForEff_[fitName].at(SimpleLR)->Fill(tp.eta());
                    hisFitTPphiForEff_[fitName].at(SimpleLR)->Fill(tp.phi0());
                    // Also plot efficiency to perfectly reconstruct the track (no fake hits)
                    if (tpRecoedPerfect) {
                        hisPerfFitTPinvptForEff_[fitName].at(SimpleLR)->Fill(1./tp.pt());
                        hisPerfFitTPetaForEff_[fitName].at(SimpleLR)->Fill(tp.eta());
                    }
                    if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
                        hisFitTPinvptForAlgEff_[fitName].at(SimpleLR)->Fill(1./tp.pt());
                        hisFitTPetaForAlgEff_[fitName].at(SimpleLR)->Fill(tp.eta());
                        hisFitTPphiForAlgEff_[fitName].at(SimpleLR)->Fill(tp.phi0());
                        // Plot also production point of all good reconstructed TP.
                        hisFitTPd0ForAlgEff_[fitName].at(SimpleLR)->Fill(fabs(tp.d0()));
                        hisFitTPz0ForAlgEff_[fitName].at(SimpleLR)->Fill(fabs(tp.z0()));
                        // Plot sector number to understand if looser cuts are needed in certain regions.
                        hisFitTPphisecForAlgEff_[fitName].at(SimpleLR)->Fill(iPhiSec_TP);
                        hisFitTPetasecForAlgEff_[fitName].at(SimpleLR)->Fill(iEtaReg_TP);
                        // Also plot efficiency to perfectly reconstruct the track (no fake hits)
                        if (tpRecoedPerfect) {
                            hisPerfFitTPinvptForAlgEff_[fitName].at(SimpleLR)->Fill(1./tp.pt());
                            hisPerfFitTPetaForAlgEff_[fitName].at(SimpleLR)->Fill(tp.eta());
                            hisPerfFitTPphisecForAlgEff_[fitName].at(SimpleLR)->Fill(iPhiSec_TP);
                            hisPerfFitTPetasecForAlgEff_[fitName].at(SimpleLR)->Fill(iEtaReg_TP);
                        }
                    }
                }
            }
        }
    }


    /////// Effi_SimpleLR ///////


}

// ------------ method called once each job just after ending the event loop  ------------
void
TMTrackAnalyzer::endJob()
{
}

map<const TP*, string> TMTrackAnalyzer::diagnoseTracking(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs) const {
    
    const vector<TP>&  vTPs = inputData.getTPs();
    
    map<const TP*, string> diagnosis;
    
    for (const TP& tp: vTPs) {
        
        string recoFlag = "unknown";
        
        if ( tp.useForAlgEff()) { //--- Only consider TP that are reconstructable.
            
            //--- Check if this TP was reconstructed anywhere in the tracker..
            bool tpRecoed = false;
            for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
                    const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
                    if (htPair.numAssocTrackCands3D( tp ) > 0) tpRecoed = true;
                }
            }
            
            if ( tpRecoed) {
                
                recoFlag = "success"; // successfully reconstructed so don't bother studying.
                
            } else {
                
                //--- Check if TP was still reconstructable after cuts applied to stubs by front-end electronics.
                vector<const Stub*> fePassStubs;
                for (const Stub* s : tp.assocStubs()) {
                    if (s->frontendPass()) fePassStubs.push_back(s);
                }
                bool fePass = ( Utility::countLayers(settings_, fePassStubs) >= settings_->genMinStubLayers() );
                
                if ( ! fePass) {
                    
                    recoFlag = "FE electronics"; // Tracking failed because of front-end electronics cuts.
                    
                } else {
                    
                    //--- Check if assignment to (eta,phi) sectors prevented this TP being reconstruted.
                    bool insideSecPass = false;
                    bool insidePhiSecPass = false;
                    bool insideEtaRegPass = false;
                    unsigned int nLayers = 0;
                    // The next to variables are vectors in case track could be recontructed in more than one sector.
                    vector< vector<const Stub*> > insideSecStubs;
                    vector<const Sector*> sectorBest;
                    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
                        for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
                            
                            const Sector& sector = mSectors(iPhiSec, iEtaReg);
                            
                            // Get stubs on given tracking particle which are inside this (phi,eta) sector;
                            vector<const Stub*> insideSecStubsTmp;
                            vector<const Stub*> insidePhiSecStubsTmp;
                            vector<const Stub*> insideEtaRegStubsTmp;
                            for (const Stub* s: fePassStubs) {
                                if (sector.inside(s))    insideSecStubsTmp.push_back(s);
                                if (sector.insidePhi(s)) insidePhiSecStubsTmp.push_back(s);
                                if (sector.insideEta(s)) insideEtaRegStubsTmp.push_back(s);
                            }
                            // Check if TP could be reconstructed in this (phi,eta) sector.
                            unsigned int nLayersTmp = Utility::countLayers(settings_, insideSecStubsTmp);
                            if ( nLayersTmp >= settings_->genMinStubLayers() ) {
                                insideSecPass = true;
                                if (nLayers <= nLayersTmp) {
                                    if (nLayers < nLayersTmp) {
                                        nLayers = nLayersTmp;
                                        insideSecStubs.clear();
                                        sectorBest.clear();
                                    }
                                    insideSecStubs.push_back( insideSecStubsTmp );
                                    sectorBest.push_back( &sector );
                                }
                            }
                            // Check if TP could be reconstructed in this (phi) sector.
                            unsigned int nLayersPhiTmp = Utility::countLayers(settings_, insidePhiSecStubsTmp);
                            if ( nLayersPhiTmp >= settings_->genMinStubLayers() ) insidePhiSecPass = true;
                            // Check if TP could be reconstructed in this (eta) region.
                            unsigned int nLayersEtaTmp = Utility::countLayers(settings_, insideEtaRegStubsTmp);
                            if ( nLayersEtaTmp >= settings_->genMinStubLayers() ) insideEtaRegPass = true;
                        }
                    }
                    
                    if ( ! insideSecPass) {
                        
                        // Tracking failed because of stub to sector assignment.
                        if ( ! insideEtaRegPass) {
                            recoFlag = "#eta sector"; // failed because of stub assignment to eta region.
                        } else if ( ! insidePhiSecPass) {
                            recoFlag = "#phi sector"; // failed because of stub assignment to phi sector.
                        } else {
                            recoFlag = "sector";      // failed because of stub assignment to (eta,phi) sector.
                        }
                        
                    } else {
                        
                        //--- Check if TP was reconstructed by r-phi Hough transform with its bend filted turned off.
                        
                        // Consider all sectors in which the track might be reconstructed.
                        bool rphiHTunfilteredPass = false;
                        for (unsigned int iSec = 0; iSec < sectorBest.size(); iSec++) {
                            HTrphi htRphiUnfiltered;
                            htRphiUnfiltered.init(settings_, sectorBest[iSec]->iPhiSec(), sectorBest[iSec]->iEtaReg(),
                                                  sectorBest[iSec]->etaMin(), sectorBest[iSec]->etaMax(), sectorBest[iSec]->phiCentre());
                            htRphiUnfiltered.disableBendFilter(); // Switch off bend filter
                            for (const Stub* s: insideSecStubs[iSec]) {
                                // Check which eta subsectors within the sector the stub is compatible with (if subsectors being used).
                                const vector<bool> inEtaSubSecs =  sectorBest[iSec]->insideEtaSubSecs( s );
                                htRphiUnfiltered.store(s, inEtaSubSecs);
                            }
                            htRphiUnfiltered.end();
                            // Check if  r-phi HT with its filters switched off found the track
                            if (htRphiUnfiltered.numTrackCands2D() > 0) rphiHTunfilteredPass = true;
                        }
                        
                        if ( ! rphiHTunfilteredPass ) {
                            
                            recoFlag = "r-#phi HT UNfiltered"; // Tracking failed r-phi HT even with its bend filter turned off.
                            
                        } else {
                            
                            //--- Check if TP was reconstructed by filtered r-phi & r-z Hough transforms.
                            
                            // Consider all sectors in which the track might be reconstructed.
                            bool rphiHTpass   = false;
                            bool rzFilterPass = false;
                            bool rzHTpass     = false;
                            for (unsigned int iSec = 0; iSec < sectorBest.size(); iSec++) {
                                HTpair htPair;
                                htPair.init(settings_, sectorBest[iSec]->iPhiSec(), sectorBest[iSec]->iEtaReg(),
                                            sectorBest[iSec]->etaMin(), sectorBest[iSec]->etaMax(), sectorBest[iSec]->phiCentre());
                                for (const Stub* s: insideSecStubs[iSec]) {
                                    // Check which eta subsectors within the sector the stub is compatible with (if subsectors being used).
                                    const vector<bool> inEtaSubSecs =  sectorBest[iSec]->insideEtaSubSecs( s );
                                    htPair.store(s, inEtaSubSecs);
                                }
                                htPair.end();
                                htPair.make3Dtracks();
                                
                                // Check if  r-phi HT found the track
                                if (htPair.getRphiHT().numTrackCands2D() > 0) rphiHTpass = true;
                                // Check if track r-z filters run after r-phi HT kept track.
                                if (rphiHTpass) {
                                    // Do so by getting tracks found by r-phi HT and running them through r-z filter.
                                    const vector<L1track2D>& trksRphi     = htPair.getRphiHT().trackCands2D();
                                    TrkRZfilter rzFilter(htPair.getRZfilters());
                                    const vector<L1track2D>& trksRphiFilt = rzFilter.filterTracks(trksRphi);
                                    if (trksRphiFilt.size() > 0) rzFilterPass = true;
                                }
                                // Check if  r-phi * r-z HTs found the track
                                if (htPair.numTrackCands3D() > 0)             rzHTpass   = true;
                            }
                            
                            if ( ! rphiHTpass) {
                                
                                recoFlag = "r-#phi HT BENDfiltered"; // Tracking failed r-phi HT with its bend filter on.
                                
                                //--- Debug printout to understand stubs failing bend filter.
                                
                                
                                cout<<"TRACK FAILING BEND FILTER: pt="<<tp.pt()<<" eta="<<tp.eta()<<endl;
                                bool okIfBendMinus1 = true;
                                for (unsigned int iSec = 0; iSec < sectorBest.size(); iSec++) {
                                    if (sectorBest.size() > 1) cout<<" SECTOR "<<iSec<<endl;
                                    for (const Stub* s: insideSecStubs[iSec]) {
                                        float bend = s->bend();
                                        float bendRes = s->bendRes();
                                        float theory = tp.qOverPt()/s->qOverPtOverBend();
                                        cout<<" BEND: measured="<<bend<<" theory="<<theory<<" res="<<bendRes<<endl;
                                        cout << s->r() << " " << s->z() << " " << s->layerId()<<" PS="<<s->psModule()<<" Barrel="<<s->barrel() << endl;
                                        
                                        if (fabs(bend - theory) > bendRes) {
                                            bool cluster0_OK = false;
                                            if (s->genuineCluster()[0]) cluster0_OK = (s->assocTPofCluster()[0]->index() == tp.index());
                                            bool cluster1_OK = false;
                                            if (s->genuineCluster()[1]) cluster1_OK = (s->assocTPofCluster()[1]->index() == tp.index());
                                            cout<< "    STUB FAILED: layer="<<s->layerId()<<" PS="<<s->psModule()<<" clusters match="<<cluster0_OK<<" "<<cluster1_OK<<endl;
                                            cout << s->bend() << " " << s->stripPitch() << " " << s->stripPitch() / s->pitchOverSep() << " " << s->dphiOverBend() << " " << s->dphi() << std::endl;
                                            cout << "Min R, Z : " << s->minR() << " " << s->minZ() << std::endl;
                                            
                                            if ( fabs( bend * -1.0 - theory ) > bendRes ) {
                                                okIfBendMinus1 = false;
                                            }
                                            else { 
                                                if ( bend > 0 ) hisWrongSignStubRZ_pBend_.at(TrackCands)->Fill( s->z(), s->r() );
                                                else if ( bend < 0 ) hisWrongSignStubRZ_nBend_.at(TrackCands)->Fill( s->z(), s->r() );
                                            }
                                        }
                                    }
                                }
                                
                                if ( okIfBendMinus1 ) {
                                    recoFlag = "BEND WRONG SIGN"; // Tracking failed r-phi HT with its bend filter on, but would have passed if bend of stubs had opposite sign.
                                }
                                
                                
                            } else {
                                
                                
                                if ( ! rzFilterPass) {
                                    
                                    recoFlag = "r-z filter"; // Tracking failed r-z filter.
                                    
                                } else {
                                    
                                    
                                    if ( ! rzHTpass) {
                                        
                                        recoFlag = "r-z HT"; // Tracking failed r-z HT.
                                        
                                    } else {
                                        
                                        recoFlag = "mystery"; // Mystery: logically this tracking particle should have been reconstructed. This may be a result of a duplicate track removal algorithm (or something else where finding one track candidate prevents another being found).
                                    }	    
                                }
                            }
                        }
                    }
                }
                diagnosis[&tp] = recoFlag;
            }
        }
    }
    return diagnosis;
}


void
TMTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}
    

//define this as a plug-in
DEFINE_FWK_MODULE(TMTrackAnalyzer);
