#ifndef __HISTOS_H__
#define __HISTOS_H__

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TMTrackTrigger/TMTrackFinder/interface/Settings.h>

#include "boost/numeric/ublas/matrix.hpp"
using  boost::numeric::ublas::matrix;

#include <vector>
#include <map>
#include <string>

class InputData;
class TP;
class Sector;
class HTpair;
class L1fittedTrack;
class L1fittedTrk4and5;
class TH1F;
class TH2F;
class TH2Poly;
class TF1;
class TProfile;
class TGraphAsymmErrors;
class TGraphErrors;
class VertexFinder;
class TEfficiency;

class Histos {

public:
  // Store cfg parameters.
  Histos(const Settings* settings) : settings_(settings), numPerfRecoTPforAlg_(0) {}

  ~Histos(){}

  // Book all histograms
  void book();

  // Fill histograms with stubs and tracking particles from input data.
  void fillInputData(const InputData& inputData);
  // Fill histograms that check if choice of (eta,phi) sectors is good.
  void fillEtaPhiSectors(const InputData& inputData, const matrix<Sector>& mSectors);
  // Fill histograms checking filling of r-phi HT array.
  void fillRphiHT(const matrix<HTpair>& mHtPairs);
  // Book histograms about r-z track filters (or other filters applied after r-phi HT array).
  void fillRZfilters(const matrix<HTpair>& mHtPairs);
  // Fill histograms studying track candidates found by Hough Transform.
  void fillTrackCands(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs);
  // Fill histograms studying freak, events with too many stubs..
  void fillStudyBusyEvents(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs);
  // Fill histograms relating to track fitting performance.
  void fillTrackFitting(const InputData& inputData, const vector<std::pair<std::string,L1fittedTrack>>& fittedTracks, float chi2dofCutPlots);
  /// Fill histograms relating to vertex reconstruction performance.
  void fillVertexReconstruction(const InputData& inputData, const VertexFinder& vf);
  
  void endJobAnalysis();

private:

  // Book histograms for specific topics.
  void bookInputData();
  void bookEtaPhiSectors();
  void bookRphiHT();
  void bookRZfilters();
  void bookTrackCands();
  void bookStudyBusyEvents();
  void bookTrackFitting();
  void bookVertexReconstruction();

  // Produce plots of tracking efficiency prior to track fit (run at end of job).
  void plotTrackEfficiency();
  // Produce plots of tracking efficiency after track fit (run at end of job).
  void plotTrackEffAfterFit(string fitName);

  void makeEfficiencyPlot( TFileDirectory &inputDir, TGraphAsymmErrors* outputEfficiency, TH1F* pass, TH1F* all, TString name, TString title );

  // Understand why not all tracking particles were reconstructed.
  // Returns list of tracking particles that were not reconstructed and an integer indicating why.
  // Only considers TP used for algorithmic efficiency measurement.
  map<const TP*, string> diagnoseTracking(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs) const;

 private:

  const Settings *settings_; // Configuration parameters.

  edm::Service<TFileService> fs_;

  // Histograms of input data.
  TProfile* profNumStubs_;
  TH1F* hisStubsVsEta_;
  TH1F* hisStubsVsR_;
  TH2F* hisStubsVsRVsZ_;
  TH2F* hisStubsModuleVsRVsZ_;
  TH2F* hisStubsModuleTiltVsZ_;
  TH2F* hisStubsdPhiCorrectionVsZ_;
  TH2F* hisStubsVsRVsPhi_;
  TH2F* hisStubsModuleVsRVsPhi_;

  TH2F* hisStubsVsRVsZ_outerModuleAtSmallerR_;
  TH2F* hisStubsVsRVsPhi_outerModuleAtSmallerR_;

  TProfile* profNumTPs_;
  TH1F* hisNumStubsPerTP_;
  TH1F* hisNumPSStubsPerTP_;
  TH1F* hisNum2SStubsPerTP_;
  TProfile* hisStubKillFE_;
  TProfile* hisStubIneffiVsInvPt_;
  TProfile* hisStubIneffiVsEta_;
  TProfile* hisStubKillDataCorr_;
  TH1F* hisPtStub_;
  TH1F* hisPtResStub_;
  TH1F* hisBendFilterPower_;
  TH1F* hisDelPhiStub_;
  TH1F* hisDelPhiResStub_;
  TH1F* hisDelPhiResStub_tilted_;
  TH1F* hisDelPhiResStub_notTilted_;
  TH1F* hisBendStub_;
  TH1F* hisBendResStub_;
  TH1F* hisNumMergedBend_;
  TH2F* hisBendVsLayerOrRing_;
  TH2F* hisBendFEVsLayerOrRing_;
  TH1F* hisPhiStubVsPhiTP_;
  TH1F* hisPhiStubVsPhi0TP_;
  TH1F* hisPhi0StubVsPhi0TP_;
  TH1F* hisPhi0StubVsPhi0TPres_;
  TH1F* hisPhiStubVsPhi65TP_;
  TH1F* hisPhi65StubVsPhi65TP_;
  TH1F* hisPhi65StubVsPhi65TPres_;
  TH1F* hisPitchOverSep_;
  TH1F* hisRhoParameter_;
  TH1F* hisFracStubsSharingClus0_;
  TH1F* hisFracStubsSharingClus1_;
  TH1F* hisPileUpPt_;
  TH1F* hisPhysicsPt_;

  // Histograms checking that (eta,phi) sector definition is good.
  TH1F* hisFracStubsInSec_;
  TH1F* hisFracStubsInEtaSec_;
  TH1F* hisFracStubsInPhiSec_;
  TH1F* hisNumSecsPerStub_;
  TH1F* hisNumEtaSecsPerStub_;
  TH1F* hisNumPhiSecsPerStub_;
  TH1F* hisNumStubsPerSec_;
  TProfile* profNumStubsPerEtaSec_;
  TH2F* hisLayerIDvsEtaSec_;
  TH2F* hisLayerIDreducedvsEtaSec_;

  TH2F* hisWrongSignStubRZ_pBend_;
  TH2F* hisWrongSignStubRZ_nBend_;
  // Histograms checking filling of r-phi HT array.
  TH1F* hisIncStubsPerHT_;
  TH1F* hisExcStubsPerHT_;
  TH2F* hisNumStubsInCellVsEta_;
  TH1F* hisStubsOnRphiTracksPerHT_;

  // Histograms about r-z track filters (or other filters applied after r-phi HT array).
  TH1F* hisNumZtrkSeedCombinations_;
  TH1F* hisNumSeedCombinations_;
  TH1F* hisNumGoodSeedCombinations_;
  TH1F* hisCorrelationZTrk_;

  // Histograms studying track candidates found by Hough Transform.
  TProfile* profNumTrackCands_;
  TProfile* profNumTracksVsEta_;
  TH1F*     hisNumTracksVsQoverPt_;
  TH1F*     hisNumTrksPerSect_;
  TH1F*     hisNumTrksPerOct_;
  TProfile* profStubsOnTracks_;
  TProfile* profStubsOnTracksVsEta_;
  TH1F*     hisStubsOnTracksPerSect_;
  TH1F*     hisStubsOnTracksPerOct_;
  TH1F*     hisUniqueStubsOnTrksPerSect_;
  TH1F*     hisUniqueStubsOnTrksPerOct_;
  TH1F*     hisStubsPerTrack_;
  TH1F*     hisLayersPerTrack_;
  TH1F*     hisPSLayersPerTrack_;
  TH1F*     hisLayersPerTrueTrack_;
  TH1F*     hisPSLayersPerTrueTrack_;
  TProfile* profExcessStubsPerTrackVsPt_;
  TH1F*     hisFracMatchStubsOnTracks_;
  TH1F* hisDeltaPhiRtruePS_;
  TH1F* hisDeltaRorZtruePS_;
  TH1F* hisDeltaPhiRtrue2S_;
  TH1F* hisDeltaRorZtrue2S_;
  TH1F* hisDeltaPhiRfakePS_;
  TH1F* hisDeltaRorZfakePS_;
  TH1F* hisDeltaPhiRfake2S_;
  TH1F* hisDeltaRorZfake2S_;
  TProfile* profNsigmaPhiRvsInvPt_;
  TProfile* profNsigmaPhiRvsFracDist_;
  TH1F* hisDeltaBendTrue_;
  TH1F* hisDeltaBendFake_;
  TProfile* profFracTrueStubsVsLayer_;
  TProfile* profDupTracksVsTPeta_;

  // Histograms of track parameter resolution after HT transform.
  TH1F* hisQoverPtRes_;
  TH1F* hisPhi0Res_;
  TH1F* hisEtaRes_;
  TH1F* hisZ0Res_;

  // Diagnosis of failed tracking.
  TH1F* hisRecoFailureReason_;
  TH1F* hisRecoFailureLayer_;

  TH1F* hisNumStubsOnLayer_;

  // Histograms used to make efficiency plots with track candidates prior to fit.
  TH1F* hisTPinvptForEff_;
  TH1F* hisRecoTPinvptForEff_;
  TH1F* hisTanLambdaRes_;
  TH1F* hisTPetaForEff_;
  TH1F* hisRecoTPetaForEff_;
  TH1F* hisTPphiForEff_;
  TH1F* hisRecoTPphiForEff_;
  //
  TH1F* hisPerfRecoTPinvptForEff_;
  TH1F* hisPerfRecoTPetaForEff_;
  //
  TH1F* hisTPinvptForAlgEff_;
  TH1F* hisRecoTPinvptForAlgEff_;
  TH1F* hisTPetaForAlgEff_;
  TH1F* hisRecoTPetaForAlgEff_;
  TH1F* hisTPphiForAlgEff_;
  TH1F* hisRecoTPphiForAlgEff_;
  //
  TH1F* hisPerfRecoTPinvptForAlgEff_;
  TH1F* hisPerfRecoTPetaForAlgEff_;
  //
  TH1F* hisTPd0ForAlgEff_;
  TH1F* hisRecoTPd0ForAlgEff_;
  TH1F* hisTPz0ForAlgEff_;
  TH1F* hisRecoTPz0ForAlgEff_;
  //
  TH1F* hisTPphisecForAlgEff_;
  TH1F* hisRecoTPphisecForAlgEff_;
  TH1F* hisPerfRecoTPphisecForAlgEff_;
  TH1F* hisTPetasecForAlgEff_;
  TH1F* hisRecoTPetasecForAlgEff_;
  TH1F* hisPerfRecoTPetasecForAlgEff_;

  // Histograms for studying freak, large events with too many stubs.
  TH1F*     hisNumBusySecsInPerEvent_;
  TH1F*     hisNumBusySecsOutPerEvent_;
  TProfile* profFracBusyInVsEtaReg_;
  TProfile* profFracBusyOutVsEtaReg_;
  TProfile* profFracStubsKilledVsEtaReg_;
  TProfile* profFracTracksKilledVsEtaReg_;
  TProfile* profFracTracksKilledVsInvPt_;
  TProfile* profFracTPKilledVsEta_;
  TProfile* profFracTPKilledVsInvPt_;
  TH1F*     hisNumTPkilledBusySec_;
  map<string, TH1F*> hisNumInputStubs_;
  map<string, TH1F*> hisQoverPtInputStubs_;
  map<string, TH1F*> hisNumOutputStubs_;
  map<string, TH1F*> hisNumTracks_; 
  map<string, TH1F*> hisNumStubsPerTrack_; 
  map<string, TH1F*> hisTrackQoverPt_; 
  map<string, TH1F*> hisTrackPurity_; 
  map<string, TH1F*> hisNumTPphysics_; 
  map<string, TH1F*> hisNumTPpileup_; 
  map<string, TH1F*> hisSumPtTPphysics_; 
  map<string, TH1F*> hisSumPtTPpileup_; 

  // Histograms for track fitting evaluation, where map index specifies name of track fitting algorithm used.
  map<std::string, TH1F*> hisSeedQinvPt_;
  map<std::string, TH1F*> hisSeedPhi0_;
  map<std::string, TH1F*> hisSeedD0_;
  map<std::string, TH1F*> hisSeedZ0_;
  map<std::string, TH1F*> hisSeedEta_;

  map<std::string, TProfile*> profNumFittedCands_;

  map<std::string, TH1F*> hisFitQinvPtMatched_;
  map<std::string, TH1F*> hisFitPhi0Matched_;
  map<std::string, TH1F*> hisFitD0Matched_;
  map<std::string, TH1F*> hisFitZ0Matched_;
  map<std::string, TH1F*> hisFitEtaMatched_;

  map<std::string, TH1F*> hisFitChi2Matched_;
  map<std::string, TH1F*> hisFitChi2DofMatched_;

  map<std::string, TH1F*> hisFitQinvPtUnmatched_;
  map<std::string, TH1F*> hisFitPhi0Unmatched_;
  map<std::string, TH1F*> hisFitD0Unmatched_;
  map<std::string, TH1F*> hisFitZ0Unmatched_;
  map<std::string, TH1F*> hisFitEtaUnmatched_;

  map<std::string, TH1F*> hisFitChi2Unmatched_;
  map<std::string, TH1F*> hisFitChi2DofUnmatched_;

  map<std::string, TH2F*> hisFitVsTrueQinvPtGoodChi2_;
  map<std::string, TH2F*> hisFitVsTruePhi0GoodChi2_;
  map<std::string, TH2F*> hisFitVsTrueD0GoodChi2_;
  map<std::string, TH2F*> hisFitVsTrueZ0GoodChi2_;
  map<std::string, TH2F*> hisFitVsTrueEtaGoodChi2_;

  map<std::string, TH2F*> hisFitVsSeedQinvPtGenCand_;
  map<std::string, TH2F*> hisFitVsSeedPhi0GenCand_;
  map<std::string, TH2F*> hisFitVsSeedD0GenCand_;
  map<std::string, TH2F*> hisFitVsSeedZ0GenCand_;
  map<std::string, TH2F*> hisFitVsSeedEtaGenCand_;

  map<std::string, TH1F*> hisFitQinvPtResGoodChi2_;
  map<std::string, TH1F*> hisFitPhi0ResGoodChi2_;
  map<std::string, TH1F*> hisFitD0ResGoodChi2_;
  map<std::string, TH1F*> hisFitZ0ResGoodChi2_;
  map<std::string, TH1F*> hisFitEtaResGoodChi2_;  

  map<std::string, TH1F*> hisSeedQinvPtResGoodChi2_;
  map<std::string, TH1F*> hisSeedPhi0ResGoodChi2_;
  map<std::string, TH1F*> hisSeedD0ResGoodChi2_;
  map<std::string, TH1F*> hisSeedZ0ResGoodChi2_;
  map<std::string, TH1F*> hisSeedEtaResGoodChi2_;  

  map<std::string, TH2F*> hisFitVsSeedQinvPtFakeCand_;
  map<std::string, TH2F*> hisFitVsSeedPhi0FakeCand_;
  map<std::string, TH2F*> hisFitVsSeedD0FakeCand_;
  map<std::string, TH2F*> hisFitVsSeedZ0FakeCand_;
  map<std::string, TH2F*> hisFitVsSeedEtaFakeCand_;

  map<std::string, TProfile*> hisQoverPtResVsTrueEta_;
  map<std::string, TProfile*> hisPhi0ResVsTrueEta_;
  map<std::string, TProfile*> hisEtaResVsTrueEta_;
  map<std::string, TProfile*> hisZ0ResVsTrueEta_;
  map<std::string, TProfile*> hisD0ResVsTrueEta_;

  map<std::string, TProfile*> hisQoverPtResVsTrueZ0_;
  map<std::string, TProfile*> hisPhi0ResVsTrueZ0_;
  map<std::string, TProfile*> hisEtaResVsTrueZ0_;
  map<std::string, TProfile*> hisZ0ResVsTrueZ0_;
  map<std::string, TProfile*> hisD0ResVsTrueZ0_;

  map<std::string, TProfile*> hisQoverPtResVsTrueInvPt_;
  map<std::string, TProfile*> hisPhi0ResVsTrueInvPt_;
  map<std::string, TProfile*> hisEtaResVsTrueInvPt_;
  map<std::string, TProfile*> hisZ0ResVsTrueInvPt_;
  map<std::string, TProfile*> hisD0ResVsTrueInvPt_;

  map<std::string, TH2F*> hisTrueFittedChiSquaredVsTrueEta_;
  map<std::string, TH2F*> hisTrueFittedChiSquaredDofVsTrueEta_;
  map<std::string, TH2F*> hisTrueFittedChiSquaredVsFittedEta_;
  map<std::string, TH2F*> hisTrueFittedChiSquaredDofVsFittedEta_;
  map<std::string, TH2F*> hisFittedChiSquaredFunctionOfStubs_;
  map<std::string, TH2F*> hisFittedChiSquaredDofFunctionOfStubs_;

  map<std::string, TH1F*> hisTrueEtaMatchedGoodChi2_;
  map<std::string, TH1F*> hisTrueEtaMatchedBadChi2_;
  map<std::string, TH1F*> hisStubPurityMatchedGoodChi2_;
  map<std::string, TH1F*> hisStubPurityMatchedBadChi2_;

  map<std::string, TProfile*> profChi2DofVsInvPtPERF_;
  map<std::string, TProfile*> profBigChi2DofVsInvPtPERF_;
  map<std::string, TH1F*>     hisD0TPBigChi2DofPERF_;
  map<std::string, TH1F*>     hisD0TPSmallChi2DofPERF_;

  map<std::string, TH2F*>     hisNumMatchedStubsKilledVsKilled_;
  map<std::string, TProfile*> profTrksKilledByFit_;
  map<std::string, TH2F*>     hisNumStubsVsPurity_;

  map<std::string, TH1F*> hisNumIterations_;
  map<std::string, TH1F*> hisFailingState_;
  map<std::string, TProfile*> hisTotalStateCalls_;
  map<std::string, TProfile*> hisRelativeStateCalls_;

  map<std::string, TH1F*> hisNumStubsFitKills_;
  map<std::string, TH2F*> hisNumStubsFitKillsVsPurity_;
  map<std::string, TH2F*> hisNumStubsFitKillsVsPurityMatched_;
  map<std::string, TH2F*> hisNumStubsFitKillsVsPurityUnmatched_;

  map<std::string, TH2F*> hisFitEfficiencyVsChi2Dof_;
  map<std::string, TH2F*> hisNumStubsVsChi2Dof_;
  map<std::string, TH2F*> hisNumLayersVsChi2Dof_;
  map<std::string, TH2F*> hisAvgNumStubsPerLayerVsChi2Dof_;

  map<std::string, TProfile*> profDupFitTrksVsEta_;
  map<std::string, TProfile*> profDupFitTrksVsInvPt_;
  map<std::string, TProfile*> profFakeFitTrksVsEta_;
  map<std::string, TProfile*> profFakeFitTrksVsZ0_;
  map<std::string, TProfile*> profFitFracTrueStubsVsLayer_;
  map<std::string, TProfile*> profFitFracTrueStubsVsEta_;

  // Histograms used for efficiency plots made with fitted tracks.
  map<std::string, TH1F*> hisFitTPinvptForEff_;
  map<std::string, TH1F*> hisFitTPetaForEff_;
  map<std::string, TH1F*> hisFitTPphiForEff_;
  map<std::string, TH1F*> hisPerfFitTPinvptForEff_;
  map<std::string, TH1F*> hisPerfFitTPetaForEff_;
  map<std::string, TH1F*> hisFitTPinvptForAlgEff_;
  map<std::string, TH1F*> hisFitTPetaForAlgEff_;
  map<std::string, TH1F*> hisFitTPphiForAlgEff_;
  map<std::string, TH1F*> hisPerfFitTPinvptForAlgEff_;
  map<std::string, TH1F*> hisPerfFitTPetaForAlgEff_;
  map<std::string, TH1F*> hisFitTPd0ForAlgEff_;
  map<std::string, TH1F*> hisFitTPz0ForAlgEff_;
  map<std::string, TH1F*> hisFitTPphisecForAlgEff_;
  map<std::string, TH1F*> hisFitTPetasecForAlgEff_;
  map<std::string, TH1F*> hisPerfFitTPphisecForAlgEff_;
  map<std::string, TH1F*> hisPerfFitTPetasecForAlgEff_;

  // Histograms of tracking efficiency & fake rate after Hough transform based on tracks prior to track fit.
  TGraphAsymmErrors* graphEffVsInvPt_;
  TGraphAsymmErrors* graphEffVsEta_;
  TGraphAsymmErrors* graphEffVsPhi_;
  //
  TGraphAsymmErrors* graphPerfEffVsInvPt_;
  TGraphAsymmErrors* graphPerfEffVsEta_;
  //
  TGraphAsymmErrors* graphAlgEffVsInvPt_;
  TGraphAsymmErrors* graphAlgEffVsEta_;
  TGraphAsymmErrors* graphAlgEffVsPhi_;
  //
  TGraphAsymmErrors* graphPerfAlgEffVsInvPt_;
  TGraphAsymmErrors* graphPerfAlgEffVsEta_;
  //
  TGraphAsymmErrors* graphAlgEffVsD0_;
  TGraphAsymmErrors* graphAlgEffVsZ0_;
  //
  TGraphAsymmErrors* graphAlgEffVsPhiSec_;
  TGraphAsymmErrors* graphAlgEffVsEtaSec_;
  TGraphAsymmErrors* graphPerfAlgEffVsPhiSec_;
  TGraphAsymmErrors* graphPerfAlgEffVsEtaSec_;

  // Histograms of tracking efficiency & fake rate after Hough transform based on tracks after the track fit.
  map<std::string, TGraphAsymmErrors*> graphEffFitVsInvPt_;
  map<std::string, TGraphAsymmErrors*> graphEffFitVsEta_;
  map<std::string, TGraphAsymmErrors*> graphEffFitVsPhi_;
  //
  map<std::string, TGraphAsymmErrors*> graphPerfEffFitVsInvPt_;
  map<std::string, TGraphAsymmErrors*> graphPerfEffFitVsEta_;
  //
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsInvPt_;
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsEta_;
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsPhi_;
  //
  map<std::string, TGraphAsymmErrors*> graphPerfAlgEffFitVsInvPt_;
  map<std::string, TGraphAsymmErrors*> graphPerfAlgEffFitVsEta_;
  //
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsD0_;
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsZ0_;
  //
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsPhiSec_;
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsEtaSec_;
  map<std::string, TGraphAsymmErrors*> graphPerfAlgEffFitVsPhiSec_;
  map<std::string, TGraphAsymmErrors*> graphPerfAlgEffFitVsEtaSec_;

  // Number of perfectly reconstructed tracks amongst TP used for algorithmic efficiency measurement.
  // Perfectly means that all stubs on track were produced by same TP.
  unsigned int numPerfRecoTPforAlg_;

  // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted.
  map<std::string, unsigned int> numFitAlgEff_;
  map<std::string, unsigned int> numFitPerfAlgEff_;

  // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted post-cut.
  map<std::string, unsigned int> numFitAlgEffPass_;
  map<std::string, unsigned int> numFitPerfAlgEffPass_;

  // Histograms for Vertex Reconstruction
  TH1F* hisNoRecoVertices_;
    
  /////// ADDED BY ME ////////
    
  TH1F* hisNoRealPVTracks_;
  TH1F* hisNoRecoPVTracks_;
  TH1F* hisRealPVPt_;
  TH1F* hisRecoPVPt_;
  TH2F* hisRecoTracksvspT_;
  TH2F* hisRealTracksvspT_;
  TH1F* hisNoRealTracksGoodEta_;
  TH1F* hisNoRecoTracksGoodEta_;
  TH1F* hisNoRealTracksBadEta_;
  TH1F* hisNoRecoTracksBadEta_;
  TH1F* hisChi2dofTracks_;
    
  TH1F* numLayersRealTracks_;
  TH1F* numLayersRecoTracks_;
    
  TH1F* numPSLayersRealTracks_;
  TH1F* numPSLayersRecoTracks_;
    
  TH1F* realTracksLayersId_;
  TH1F* recoTracksLayersId_;
    

  TH1F* hisNoRealTracksPerc50_;
  TH1F* hisNoRecoTracksPerc50_;
  TH1F* hisNoRealTracksPerc90_;
  TH1F* hisNoRecoTracksPerc90_;
    
  TH1F* hisz0RealPerc50_;
  TH1F* hisz0RecoPerc50_;
  TH1F* hisz0RealPerc90_;
  TH1F* hisz0RecoPerc90_;
    
  TH1F* hisRealpTPerc50_;
  TH1F* hisRecopTPerc50_;
  TH1F* hisRealpTPerc90_;
  TH1F* hisRecopTPerc90_;
   
  /////////////////////////////
    
  TH1F* hisNoPileUpVertices_;
  TH1F* hisRecoVertexZ0Resolution_;
  TH1F* hisRecoVertexPTResolution_;
  TH2F* hisNoRecoVsNoTruePileUpVertices_;  
  TH2F* hisRecoVertexMETVsTrueMET_;
  TH2F* hisNoTracksFromPrimaryVertex_;
  TProfile* hisRecoVertexPTResolutionVsTruePt_;
  TH2F* hisNoTrueTracksFromPrimaryVertex_;
  TH1F* hisRecoPrimaryVertexZ0width_;
  TH1F* hisRecoPileUpVertexZ0width_;
  TH1F* hisRecoVertexZ0Spacing_;
  TH1F* hisPrimaryVertexZ0width_;  
  TH1F* hisPileUpVertexZ0_;
  TH1F* hisPileUpVertexZ0width_;
  TH1F* hisPileUpVertexZ0Spacing_;
  TH1F* hisRecoPileUpVertexZ0resolution_;
  TH1F* hisRatioMatchedTracksInPV_;
  TH1F* hisFakeTracksRateInPV_;
  TH1F* hisTrueTracksRateInPV_;
  TH2F* hisRecoVertexPTVsTruePt_;
  TH1F* hisUnmatchZ0distance_;
  TH1F* hisUnmatchZ0MinDistance_;
  TH1F* hisUnmatchPt_      ;
  TH1F* hisUnmatchEta_     ;
  TH1F* hisUnmatchTruePt_  ;
  TH1F* hisUnmatchTrueEta_ ;
  TH1F* hisLostPVtracks_      ;
  TH1F* hisUnmatchedPVtracks_ ;
  TH1F* hisNumVxIterations_;
  TH1F* hisNumVxIterationsPerTrack_;
  TProfile* hisTrkMETvsGenMET_;
  TProfile* hisRecoTrkMETvsGenMET_;
  TProfile* hisTDRTrkMETvsGenMET_;
  TH1F* hisRecoPrimaryVertexVsTrueZ0_;
  TH1F* hisTDRPrimaryVertexVsTrueZ0_;
  TH1F* hisPrimaryVertexTrueZ0_;
  TH1F* hisRecoVertexMET_;
  TH1F* hisRecoVertexPT_;
  TH1F* hisRecoPileUpVertexPT_;
  TH1F* hisRecoVertexOffPT_;
  TProfile* hisRecoPrimaryVertexResolutionVsTrueZ0_;
  TProfile* hisTDRPrimaryVertexResolutionVsTrueZ0_;
  TH1F* hisTDRVertexZ0Resolution_           ;
  TH1F* hisTDRVertexPTResolution_           ;
  TProfile* hisTDRVertexPTResolutionVsTruePt_   ;
  TH2F* hisTDRVertexPTVsTruePt_             ;
  TH2F* hisTDRVertexMETVsTrueMET_           ;
  TH2F* hisTDRNoTracksFromPrimaryVertex_    ;
  TH2F* hisTDRNoTrueTracksFromPrimaryVertex_;
  TH1F* hisTDRPrimaryVertexZ0width_         ;
  TH1F* hisRatioMatchedTracksInTDRPV_       ;
  TH1F* hisFakeTracksRateInTDRPV_           ;
  TH1F* hisTrueTracksRateInTDRPV_           ;
  TH1F* hisTDRUnmatchZ0distance_            ;
  TH1F* hisTDRUnmatchZ0MinDistance_         ;
  TH1F* hisTDRUnmatchPt_                    ;
  TH1F* hisTDRUnmatchEta_                   ;
  TH1F* hisTDRUnmatchTruePt_                ;
  TH1F* hisTDRUnmatchTrueEta_               ;
  TH1F* hisTDRUnmatchedPVtracks_            ;
  TH1F* hisUnmatchedVertexZ0distance_       ;
  TH1F* hisTDRUnmatchedVertexZ0distance_    ;
  TH1F* hisTDRVertexMET_;
  TH1F* hisTDRVertexPT_;
  TH1F* hisTDRPileUpVertexPT_;
  TH1F* hisTDRVertexOffPT_;
  
  std::vector<TH1F*> hisMETevents_;
  std::vector<TH1F*> hisRecoVerticesVsTrueZ0PerDistance_;
  std::vector<TH1F*> hisTDRVerticesVsTrueZ0PerDistance_;


  std::vector<TGraphErrors*> grMET_;
  std::vector<TGraphErrors*> grMET_tdr_;
  TEfficiency* PVefficiencyVsTrueZ0_;
  TEfficiency* tdrPVefficiencyVsTrueZ0_;
    
  TEfficiency* PVefficiencyVsNoTracks_;
  TEfficiency* PVefficiencyVsNoTracksGoodPercentage50_;
  TEfficiency* PVefficiencyVsNoTracksGoodPercentage90_;
  TEfficiency* PVefficiencyz0Perc50;
  TEfficiency* PVefficiencyz0Perc90;
  TEfficiency* PVefficiencypTPerc50;
  TEfficiency* PVefficiencypTPerc90;
  TEfficiency* PVefficiencyNoLayers;
  TEfficiency* PVefficiencyNoPSLayers;
  TEfficiency* PVefficiencyLayerID;

    
  TEfficiency* PVefficiencyVsNoTracksVsZ0_;
  TEfficiency* PVefficiencyVspT_;
  TEfficiency* PVefficiencyVsGoodEtaTracks_;
  TEfficiency* PVefficiencyVsBadEtaTracks_;
  std::vector<unsigned int> noSignalEvents;
  std::vector<unsigned int> noBackgroundEvents;
  std::vector<unsigned int> noSignalEventsTDR;
  std::vector<unsigned int> noBackgroundEventsTDR;
  std::vector<std::vector<unsigned int> > noRecoSignalEvents;
  std::vector<std::vector<unsigned int> > noRecoBackgroundEvents;
  std::vector<std::vector<unsigned int> > noTDRSignalEvents;
  std::vector<std::vector<unsigned int> > noTDRBackgroundEvents;

  unsigned int noEvents;
};

#endif
