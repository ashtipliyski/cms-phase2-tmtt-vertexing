#include "TMTrackTrigger/TMTrackFinder/interface/Histos.h"
#include "TMTrackTrigger/TMTrackFinder/interface/InputData.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Sector.h"
#include "TMTrackTrigger/TMTrackFinder/interface/HTpair.h"
#include "TMTrackTrigger/TMTrackFinder/interface/HTrphi.h"
#include "TMTrackTrigger/TMTrackFinder/interface/HTrz.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TrkRZfilter.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrk4and5.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"
#include "TMTrackTrigger/TMTrackFinder/interface/VertexFinder.h"
#include "TGraphErrors.h"


#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <TProfile.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include "TStyle.h"
#include "TROOT.h"

#include <algorithm>
#include <array>
#include <unordered_set>
#include "fstream"

using namespace std;

//=== Book all histogram

void Histos::book() {
  TH1::SetDefaultSumw2(true);

  // Book histograms about various topics.
  this->bookInputData();
  // Book histograms checking if (eta,phi) sector definition choices are good.
  this->bookEtaPhiSectors();
  // Book histograms checking filling of r-phi HT array.
  this->bookRphiHT();
  // Book histograms about r-z track filters (or other filters applied after r-phi HT array).
  this->bookRZfilters();
  // Book histograms studying track candidates found by Hough Transform.
  this->bookTrackCands();
  // Book histograms studying track fitting performance
  this->bookTrackFitting();
  // Book histograms studying vertex reconstruction performance
  this->bookVertexReconstruction();
}

//=== Book histograms using input stubs and tracking particles.

void Histos::bookInputData() {
  TFileDirectory inputDir = fs_->mkdir("InputData");

  // N.B. Histograms of the kinematics and production vertex of tracking particles
  // are booked in bookTrackCands(), since they are used to study the tracking efficiency.

  // Count stubs & tracking particles.

  profNumStubs_        = inputDir.make<TProfile>("NumStubs","; Category; No. stubs in tracker",4,0.5,4.5);
  profNumStubs_->GetXaxis()->SetBinLabel(1,"All stubs");
  profNumStubs_->GetXaxis()->SetBinLabel(2,"Genuine stubs");
  profNumStubs_->GetXaxis()->SetBinLabel(3,"Stubs matched to TP");
  profNumStubs_->GetXaxis()->SetBinLabel(4,"Stubs matched to TP for eff");
  profNumStubs_->LabelsOption("d");

  hisStubsVsEta_      = inputDir.make<TH1F>("StubsVsEta","; #eta; No. stubs in tracker",30,-3.0,3.0);
  hisStubsVsR_        = inputDir.make<TH1F>("StubsVsR","; radius (cm); No. stubs in tracker",1200,0.,120.);

  hisStubsVsRVsZ_      = inputDir.make<TH2F>("StubsVsRVsZ","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130);
  hisStubsModuleVsRVsZ_      = inputDir.make<TH2F>("StubsModuleVsRVsZ","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130);
  hisStubsVsRVsPhi_      = inputDir.make<TH2F>("StubsVsRVsPhi","; x (cm); y (cm); No. stubs in tracker",1000,-130,130,1000,-130,130);
  hisStubsModuleVsRVsPhi_      = inputDir.make<TH2F>("StubsModuleVsRVsPhi","; x (cm); y (cm); No. stubs in tracker",1000,-130,130,1000,-130,130);
  hisStubsModuleTiltVsZ_      = inputDir.make<TH2F>("StubsModuleTiltVsZ","; z (cm); Tilt; Module tilt vs z",1000,-280,280,1000,-3,3);
  hisStubsdPhiCorrectionVsZ_      = inputDir.make<TH2F>("StubsdPhiCorrectionVsZ","; z (cm); Correction; dPhi Correction vs z",1000,-280,280,1000,-1,10);

  hisStubsVsRVsZ_outerModuleAtSmallerR_      = inputDir.make<TH2F>("StubsVsRVsZ_outerModuleAtSmallerR","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130);
  hisStubsVsRVsPhi_outerModuleAtSmallerR_     = inputDir.make<TH2F>("StubsVsRVsPhi_outerModuleAtSmallerR","; x (cm); y (cm); No. stubs in tracker",1000,-130,130,1000,-130,130);

  profNumTPs_          = inputDir.make<TProfile>("NumTPs","; Category; No. of TPs in tracker",3,0.5,3.5);
  profNumTPs_->GetXaxis()->SetBinLabel(1,"All TPs");
  profNumTPs_->GetXaxis()->SetBinLabel(2,"TPs for eff.");
  profNumTPs_->GetXaxis()->SetBinLabel(3,"TPs for alg. eff.");
  profNumTPs_->LabelsOption("d");

  hisNumStubsPerTP_    = inputDir.make<TH1F>("NumStubsPerTP","; Number of stubs per TP for alg. eff.",50,-0.5,49.5);
  hisNumPSStubsPerTP_    = inputDir.make<TH1F>("NumPSStubsPerTP","; Number of PS stubs per TP for alg. eff.",50,-0.5,49.5);
  hisNum2SStubsPerTP_    = inputDir.make<TH1F>("Num2SStubsPerTP","; Number of 2S stubs per TP for alg. eff.",50,-0.5,49.5);

  // Study efficiency of tightened front end-electronics cuts.
  hisStubKillFE_          = inputDir.make<TProfile>("StubKillFE","; barrelLayer or 10+endcapRing; Stub fraction rejected by readout chip",30,-0.5,29.5);
  hisStubIneffiVsInvPt_   = inputDir.make<TProfile>("StubIneffiVsPt","; 1/Pt; Inefficiency of readout chip for good stubs",30,0.0,1.0);
  hisStubIneffiVsEta_     = inputDir.make<TProfile>("StubIneffiVsEta","; |#eta|; Inefficiency of readout chip for good stubs",30,0.0,3.0);
  hisStubKillDataCorr_    = inputDir.make<TProfile>("StubKillDataCorr","; barrelLayer or 10+endcapRing; Stub fraction killed by DataCorrection.h window cut",30,-0.5,29.5);

  // Study stub resolution.

  hisPtStub_             = inputDir.make<TH1F>("PtStub","; Stub q/Pt",50,-0.5,0.5);
  hisPtResStub_          = inputDir.make<TH1F>("PtResStub","; Stub q/Pt minus TP q/Pt",50,-0.5,0.5);
  hisBendFilterPower_    = inputDir.make<TH1F>("BendFilterPower","; Fraction of q/Pt range allowed",102,-0.01,1.01);
  hisDelPhiStub_         = inputDir.make<TH1F>("DelPhiStub","; Stub bend angle",50,-0.2,0.2);
  hisDelPhiResStub_      = inputDir.make<TH1F>("DelPhiResStub","; Stub bend angle minus TP bend angle",200,-0.2,0.2);

  hisDelPhiResStub_tilted_      = inputDir.make<TH1F>("DelPhiResStub_tilted","; Stub bend angle minus TP bend angle",200,-0.2,0.2);
  hisDelPhiResStub_notTilted_      = inputDir.make<TH1F>("DelPhiResStub_notTilted","; Stub bend angle minus TP bend angle",200,-0.2,0.2);

  hisBendStub_           = inputDir.make<TH1F>("BendStub","; Stub bend in units of strips",57,-7.125,7.125);
  hisBendResStub_        = inputDir.make<TH1F>("BendResStub","; Stub bend minus TP bend in units of strips",100,-5.,5.);
  hisNumMergedBend_      = inputDir.make<TH1F>("NumMergedBend","; No. of bend values merged together by loss of bit",10,-0.5,9.5);
  hisBendVsLayerOrRing_  = inputDir.make<TH2F>("BendVsLayerOrRing","; barrelLayer or 10+endcapRing; Stub bend",30,-0.5,29.5,57,-7.125,7.125);
  hisBendFEVsLayerOrRing_  = inputDir.make<TH2F>("BendFEVsLayerOrRing","; barrelLayer or 10+endcapRing; Stub bend in FE chip",30,-0.5,29.5,57,-7.125,7.125);

  hisPhiStubVsPhiTP_     = inputDir.make<TH1F>("PhiStubVsPhiTP","; Stub #phi minus TP #phi at stub radius",100,-0.05,0.05);
  hisPhiStubVsPhi0TP_    = inputDir.make<TH1F>("PhiStubVsPhi0TP","; Stub #phi minus TP #phi0",100,-0.3,0.3);
  hisPhi0StubVsPhi0TP_   = inputDir.make<TH1F>("Phi0StubVsPhi0TP","; #phi0 of Stub minus TP",100,-0.2,0.2);
  hisPhi0StubVsPhi0TPres_= inputDir.make<TH1F>("Phi0StubVsPhi0TPres","; #phi0 of Stub minus TP / resolution",100,-5.0,5.0);
  hisPhiStubVsPhi65TP_   = inputDir.make<TH1F>("PhiStubVsPhi65TP","; Stub #phi minus TP phitrk65",100,-0.2,0.2);
  hisPhi65StubVsPhi65TP_ = inputDir.make<TH1F>("Phi65StubVsPhi65TP","; phitrk65 of Stub minus TP",100,-0.2,0.2);
  hisPhi65StubVsPhi65TPres_ = inputDir.make<TH1F>("Phi65StubVsPhi65TPres","; phitrk65 of Stub minus TP / resolution",100,-5.0,5.0);

  // Note ratio of sensor pitch to separation (needed to understand how many bits this can be packed into).
  hisPitchOverSep_       = inputDir.make<TH1F>("PitchOverSep","; ratio of sensor pitch / separation",100,0.0,0.1);
  hisRhoParameter_       = inputDir.make<TH1F>("RhoParameter","; rho parameter",100,0.0,0.2);
  // Count stubs sharing a common cluster.
  hisFracStubsSharingClus0_ = inputDir.make<TH1F>("FracStubsSharingClus0","Fraction of stubs sharing cluster in seed sensor",102,-0.01,1.01);
  hisFracStubsSharingClus1_ = inputDir.make<TH1F>("FracStubsSharingClus1","Fraction of stubs sharing cluster in correlation sensor",102,-0.01,1.01);
  hisPileUpPt_           = inputDir.make<TH1F>("hisPileUpPt_","Pile-Up tracks transverse momentum [GeV/c]; no. Tracks; p_{T} [GeV/c]", 100, 0, 100);
  hisPhysicsPt_           = inputDir.make<TH1F>("hisPhysicsPt_","Physics Collision tracks transverse momentum [GeV/c]; no. Tracks; p_{T} [GeV/c]", 100, 0, 100);
}

//=== Fill histograms using input stubs and tracking particles.

void Histos::fillInputData(const InputData& inputData) {
  const vector<const Stub*>& vStubs = inputData.getStubs();
  const vector<TP>&          vTPs   = inputData.getTPs();

  // Count stubs.
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
  profNumStubs_->Fill(1, vStubs.size());
  profNumStubs_->Fill(2, nStubsGenuine);
  profNumStubs_->Fill(3, nStubsWithTP);
  profNumStubs_->Fill(4, nStubsWithTPforEff);

  for (const Stub* stub : vStubs) {
    hisStubsVsEta_->Fill(stub->eta());
    hisStubsVsR_->Fill(stub->r());
    hisStubsVsRVsZ_->Fill( stub->z(), stub->r() );
    hisStubsModuleVsRVsZ_->Fill( stub->minZ(), stub->minR() );
    hisStubsModuleVsRVsZ_->Fill( stub->maxZ(), stub->maxR() );

    hisStubsModuleTiltVsZ_->Fill( stub->minZ(), stub->moduleTilt() );
    hisStubsModuleTiltVsZ_->Fill( stub->maxZ(), stub->moduleTilt() );

    if ( stub->outerModuleAtSmallerR() ) {
      hisStubsVsRVsZ_outerModuleAtSmallerR_->Fill( stub->z(), stub->r() );
    }

    hisStubsdPhiCorrectionVsZ_->Fill(  stub->minZ(), stub->dphiOverBendCorrection() );

    hisStubsVsRVsPhi_->Fill( stub->r() * sin( stub->phi() ), stub->r() * cos( stub->phi() ) );
    hisStubsModuleVsRVsPhi_->Fill( stub->minR() * sin( stub->minPhi() ), stub->minR() * cos( stub->minPhi() ) );
    hisStubsModuleVsRVsPhi_->Fill( stub->maxR() * sin( stub->maxPhi() ), stub->maxR() * cos( stub->maxPhi() ) );

    if ( stub->outerModuleAtSmallerR() ) {
      hisStubsVsRVsPhi_outerModuleAtSmallerR_->Fill( stub->r() * sin( stub->phi() ), stub->r() * cos( stub->phi() ) );
    }

  }

  // Count tracking particles.
  unsigned int nTPforEff = 0;
  unsigned int nTPforAlgEff = 0;
  for (const TP& tp: vTPs) {
    if (tp.useForEff())  nTPforEff++; 
    if (tp.useForAlgEff()) nTPforAlgEff++; 
  }
  profNumTPs_->Fill(1, vTPs.size());
  profNumTPs_->Fill(2, nTPforEff);
  profNumTPs_->Fill(3, nTPforAlgEff);

  // Study efficiency of stubs to pass front-end electronics cuts.

  const vector<Stub>& vAllStubs = inputData.getAllStubs(); // Get all stubs prior to FE cuts to do this.
  for (const Stub s : vAllStubs) {
    unsigned int layerOrTenPlusRing = s.barrel()  ?  s.layerId()  :  10 + s.endcapRing(); 
    // Fraction of all stubs (good and bad) failing tightened front-end electronics cuts.
    hisStubKillFE_->Fill(layerOrTenPlusRing, (! s.frontendPass()));
    // Fraction of stubs rejected by window cut in DataCorrection.h
    // If it is non-zero, then encoding in DataCorrection.h should ideally be changed to make it zero.
    hisStubKillDataCorr_->Fill(layerOrTenPlusRing, s.stubFailedDataCorrWindow());
  }

  // Study efficiency for good stubs of tightened front end-electronics cuts.
  for (const TP& tp : vTPs) {
    if (tp.useForAlgEff()) {// Only bother for stubs that are on TP that we have a chance of reconstructing.
      hisPhysicsPt_->Fill(tp.pt());
      const vector<const Stub*> stubs = tp.assocStubs();
      for (const Stub* s : stubs) {
        hisStubIneffiVsInvPt_->Fill(1./tp.pt()    , (! s->frontendPass()) );
        hisStubIneffiVsEta_->Fill  (fabs(tp.eta()), (! s->frontendPass()) );
      }
    }
    if(tp.useForVertexReco() and !tp.physicsCollision() ){
      hisPileUpPt_->Fill(tp.pt());
    }
  }

  // Plot stub bend-derived information.
  for (const Stub* stub : vStubs) {
    hisPtStub_->Fill(stub->qOverPt()); 
    hisDelPhiStub_->Fill(stub->dphi()); 
    hisBendStub_->Fill(stub->dphi() / stub->dphiOverBend());
    // Number of bend values merged together by loss of a bit.
    hisNumMergedBend_->Fill(stub->numMergedBend()); 
    // Min. & max allowed q/Pt obtained from stub bend.
    float minQoverPt = max(float(-1./(settings_->houghMinPt())), stub->qOverPt() - stub->qOverPtres());  
    float maxQoverPt = min(float(1./(settings_->houghMinPt())), stub->qOverPt() + stub->qOverPtres());  
    // Frac. of full q/Pt range allowed by stub bend.
    float fracAllowed = (maxQoverPt - minQoverPt)/(2./(settings_->houghMinPt()));
    hisBendFilterPower_->Fill(fracAllowed);
    unsigned int layerOrTenPlusRing = stub->barrel()  ?  stub->layerId()  :  10 + stub->endcapRing(); 
    hisBendVsLayerOrRing_->Fill(layerOrTenPlusRing, stub->bend());
    // Also plot bend prior to degradation.
    hisBendFEVsLayerOrRing_->Fill(layerOrTenPlusRing, stub->bendInFrontend());
  }

  // Look at stub resolution.
  for (const TP& tp: vTPs) {
    if (tp.useForAlgEff()) {
      const vector<const Stub*>& assStubs= tp.assocStubs();
      hisNumStubsPerTP_->Fill( assStubs.size() );

      unsigned int numPSstubs = 0;
      unsigned int num2Sstubs = 0;

      //cout<<"=== TP === : index="<<tp.index()<<" pt="<<tp.pt()<<" q="<<tp.charge()<<" phi="<<tp.phi0()<<" eta="<<tp.eta()<<" z0="<<tp.z0()<<endl;
      for (const Stub* stub: assStubs) {

        if ( stub->psModule() ) ++numPSstubs;
        else ++num2Sstubs;

        //cout<<"    stub : index="<<stub->index()<<" barrel="<<stub->barrel()<<" r="<<stub->r()<<" phi="<<stub->phi()<<" z="<<stub->z()<<" bend="<<stub->bend()<<" assocTP="<<stub->assocTP()->index()<<endl; 
        hisPtResStub_->Fill(stub->qOverPt() - tp.charge()/tp.pt()); 
        hisDelPhiResStub_->Fill(stub->dphi() - tp.dphi(stub->r()));

        if ( stub->moduleTilt() > M_PI / 2 - 0.1  || !stub->barrel() ) {
          hisDelPhiResStub_notTilted_->Fill(stub->dphi() - tp.dphi(stub->r()));
        } else {
          hisDelPhiResStub_tilted_->Fill(stub->dphi() - tp.dphi(stub->r()));
        }
        hisBendResStub_->Fill( (stub->dphi() - tp.dphi(stub->r())) / stub->dphiOverBend() ); 
	// This checks if the TP multiple scattered before producing the stub or hit resolution effects.
        hisPhiStubVsPhiTP_->Fill( reco::deltaPhi(stub->phi(), tp.trkPhiAtStub( stub )) );
	// This checks how wide overlap must be if using phi0 sectors, with no stub bend info used for assignment.
        hisPhiStubVsPhi0TP_->Fill( reco::deltaPhi(stub->phi(), tp.phi0()) );
	// This checks how wide overlap must be if using phi0 sectors, with stub bend info used for assignment
        hisPhi0StubVsPhi0TP_->Fill( reco::deltaPhi(stub->trkPhiAtR(0.).first, tp.phi0()) );
	// This normalizes the previous distribution to the predicted resolution to check if the latter is OK.
        hisPhi0StubVsPhi0TPres_->Fill( reco::deltaPhi(stub->trkPhiAtR(0.).first, tp.phi0()) / stub->trkPhiAtRres(0.));
	// This checks how wide overlap must be if using phi65 sectors, with no stub bend info used for assignment.
        hisPhiStubVsPhi65TP_->Fill( reco::deltaPhi(stub->phi(), tp.trkPhiAtR(65.)) );
	// This checks how wide overlap must be if using phi65 sectors, with stub bend info used for assignment, optionally reducing discrepancy by uncertainty expected from 2S module strip length.
	pair<float, float> phiAndErr = stub->trkPhiAtR(65.);
	double dPhi = reco::deltaPhi( phiAndErr.first, tp.trkPhiAtR(65.));
        hisPhi65StubVsPhi65TP_->Fill( dPhi );
	// This normalizes the previous distribution to the predicted resolution to check if the latter is OK.
        hisPhi65StubVsPhi65TPres_->Fill( dPhi / stub->trkPhiAtRres(65.));
      }

      hisNumPSStubsPerTP_->Fill( numPSstubs );
      hisNum2SStubsPerTP_->Fill( num2Sstubs );

    }
  }

  for (const Stub* stub : vStubs) {
    // Note ratio of sensor pitch to separation (needed to understand how many bits this can be packed into).
    hisPitchOverSep_->Fill(stub->pitchOverSep());
    // Also note this same quantity times 1.0 in the barrel or z/r in the endcap. This product is known as "rho".
    float rho = stub->pitchOverSep();
    if ( ! stub->barrel() ) rho *= fabs(stub->z())/stub->r();
    hisRhoParameter_->Fill(rho);
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
      hisFracStubsSharingClus0_->Fill(float(nShare)/float(vStubs.size()));
    } else {
      hisFracStubsSharingClus1_->Fill(float(nShare)/float(vStubs.size()));
    }
  }
}

//=== Book histograms checking if (eta,phi) sector definition choices are good.

void Histos::bookEtaPhiSectors() {
  TFileDirectory inputDir = fs_->mkdir("CheckSectors");

  // Check if TP lose stubs because not all in same sector.

  hisFracStubsInSec_    = inputDir.make<TH1F>("FracStubsInSec","; Fraction of stubs on TP in best (#eta,#phi) sector;",102,-0.01,1.01);
  hisFracStubsInEtaSec_ = inputDir.make<TH1F>("FracStubsInEtaSec","; Fraction of stubs on TP in best #eta sector;",102,-0.01,1.01);
  hisFracStubsInPhiSec_ = inputDir.make<TH1F>("FracStubsInPhiSec","; Fraction of stubs on TP in best #phi sector;",102,-0.01,1.01);

  // Check if stubs excessively duplicated between overlapping sectors.

  hisNumSecsPerStub_    = inputDir.make<TH1F>("NumSecPerStub","; Number of (#eta,#phi) sectors each stub appears in",20,-0.5,19.5);
  hisNumEtaSecsPerStub_ = inputDir.make<TH1F>("NumEtaSecPerStub","; Number of #eta sectors each stub appears in",20,-0.5,19.5);
  hisNumPhiSecsPerStub_ = inputDir.make<TH1F>("NumPhiSecPerStub","; Number of #phi sectors each stub appears in",20,-0.5,19.5);

  // Count stubs per (eta,phi) sector.
  hisNumStubsPerSec_  = inputDir.make<TH1F>("NumStubsPerSec","; Number of stubs per sector",250,-0.5,249.5);
  // Ditto, summed over all phi. This checks if equal stubs go into each eta region, important for latency.
  unsigned int nEta = settings_->numEtaRegions();
  profNumStubsPerEtaSec_ = inputDir.make<TProfile>("NumStubsPerEtaSec",";#eta sector; Number of stubs per #eta sector",nEta,-0.5,nEta-0.5);

  // Check which tracker layers are present in each eta sector.
  hisLayerIDvsEtaSec_ = inputDir.make<TH2F>("LayerIDvsEtaSec",";#eta sector; layer ID",nEta,-0.5,nEta-0.5,20,0.5,20.5);
  hisLayerIDreducedvsEtaSec_ = inputDir.make<TH2F>("LayerIDreducedvsEtaSec",";#eta sector; reduced layer ID",nEta,-0.5,nEta-0.5,20,0.5,20.5);
}

//=== Fill histograms checking if (eta,phi) sector definition choices are good.

void Histos::fillEtaPhiSectors(const InputData& inputData, const matrix<Sector>& mSectors) {

  const vector<const Stub*>& vStubs = inputData.getStubs();
  const vector<TP>&          vTPs   = inputData.getTPs();

  //=== Loop over good tracking particles, looking for the (eta,phi) sector in which each has the most stubs.
  //=== and checking what fraction of its stubs were in this sector.
 
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
      hisFracStubsInSec_->Fill   ( float(nStubsInBestSec)    / float(nStubs) );
      hisFracStubsInEtaSec_->Fill( float(nStubsInBestEtaSec) / float(nStubs) );
      hisFracStubsInPhiSec_->Fill( float(nStubsInBestPhiSec) / float(nStubs) );
    }
  }

  //=== Loop over all stubs, counting how many sectors each one appears in. 
 
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
	      hisLayerIDvsEtaSec_->Fill(iEtaReg, lay);
	      hisLayerIDreducedvsEtaSec_->Fill(iEtaReg, stub->layerIdReduced()); // Plot also simplified layerID for hardware, which tries to avoid more than 8 ID in any given eta region.
	    }
	  }
	}
      }
    }

    // Plot number of sectors each stub appears in.
    hisNumSecsPerStub_->Fill   ( nSecs );
    hisNumEtaSecsPerStub_->Fill( nEtaSecs );
    hisNumPhiSecsPerStub_->Fill( nPhiSecs );

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
      hisNumStubsPerSec_->Fill(nStubs);
      nStubsInEtaSec += nStubs;
    }
    profNumStubsPerEtaSec_->Fill(iEtaReg, nStubsInEtaSec);
  }
}

//=== Book histograms checking filling of r-phi HT array.

void Histos::bookRphiHT() {

  TFileDirectory inputDir = fs_->mkdir("HTrphi");

  hisIncStubsPerHT_ = inputDir.make<TH1F>("IncStubsPerHT","; Number of filtered stubs per r#phi HT array (inc. duplicates)",100,0.,-1.);
  hisExcStubsPerHT_ = inputDir.make<TH1F>("ExcStubsPerHT","; Number of filtered stubs per r#phi HT array (exc. duplicates)",250,-0.5,249.5);

  hisNumStubsInCellVsEta_ = inputDir.make<TH2F>("NumStubsInCellVsEta","; no. of stubs per HT cell summed over phi sector; #eta region",100,-0.5,499.5, settings_->numEtaRegions(), -0.5, settings_->numEtaRegions() - 0.5);

  hisStubsOnRphiTracksPerHT_ = inputDir.make<TH1F>("StubsOnRphiTracksPerHT","; Number of stubs assigned to tracks per r#phi HT array",500,-0.5,499.5);
}

//=== Fill histograms checking filling of r-phi HT array.

void Histos::fillRphiHT(const matrix<HTpair>& mHtPairs) {

  //--- Loop over (eta,phi) sectors, counting the number of stubs in the HT array of each.
 
  for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
      const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
      const HTrphi& htRphi = htPair.getRphiHT();

      // Here, if a stub appears in multiple cells, it is counted multiple times.
      hisIncStubsPerHT_->Fill( htRphi.numStubsInc() );
      // Here, if a stub appears in multiple cells, it is counted only once.
      hisExcStubsPerHT_->Fill( htRphi.numStubsExc() );
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
        hisNumStubsInCellVsEta_->Fill( nStubsInCellPhiSum, iEtaReg );
      }
    }
  }

  //--- Count number of cells assigned to track candidates by r-phi HT (before any rz filtering 
  //--- or rz HT has been run).
  for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
    for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
      const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
      const HTrphi& htRphi = htPair.getRphiHT();
      hisStubsOnRphiTracksPerHT_->Fill(htRphi.numStubsOnTrackCands2D()); 
    }
  }
}

//=== Book histograms about r-z track filters (or other filters applied after r-phi HT array).

void Histos::bookRZfilters() {

  // Only book histograms if one of the r-z filters was in use.
  if (settings_->useZTrkFilter() || settings_->useSeedFilter()) {

    TFileDirectory inputDir = fs_->mkdir("RZfilters");

    // Check number of track seeds that r-z filters must check.

    if (settings_->useZTrkFilter()) {
      hisNumZtrkSeedCombinations_ = inputDir.make<TH1F>("NumZtrkSeedCombinations_","; Number of Ztrk seed combinations per track cand; no. seeds ; ", 50, -0.5 , 49.5);
    }

    if (settings_->useSeedFilter()) {
      hisNumSeedCombinations_ = inputDir.make<TH1F>("NumSeedCombinations_","; Number of seed combinations per track cand; no. seeds ; ", 50, -0.5 , 49.5);
      hisNumGoodSeedCombinations_ = inputDir.make<TH1F>("NumGoodSeedCombinations_","; Number of good seed combinations per track cand; ", 30, -0.5 , 29.5);
    }

    // Look at correlation term used inside Ztrk filter (for experts only)
    if (settings_->useZTrkFilter()) {
      hisCorrelationZTrk_ = inputDir.make<TH1F>("CorrelationZTrk_","; Correlation factor r between stub in layer 1 and stubs in other layers; r ; ", 50,-1,1);
    }

  }
}

//=== Fill histograms about r-z track filters (or other filters applied after r-phi HT array).

void Histos::fillRZfilters(const matrix<HTpair>& mHtPairs) {

  // Only fill histograms if one of the r-z filters was in use.
  if (settings_->useZTrkFilter() || settings_->useSeedFilter()) {

    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
      for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
	const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);

	if (settings_->useZTrkFilter()) {
	  // Check number of track seeds per sector that r-z "Ztrk" filter checked.
	  const vector<unsigned int>  numSeedComb = htPair.getRZfilters().numZtrkSeedCombsPerTrk();
	  for (const unsigned int& num : numSeedComb) {
	    hisNumZtrkSeedCombinations_->Fill(num) ;
	  }
	}

	if (settings_->useSeedFilter()) {
	  // Check number of track seeds per sector that r-z "seed" filter checked.
	  const vector<unsigned int>  numSeedComb = htPair.getRZfilters().numSeedCombsPerTrk();
	  for (const unsigned int& num : numSeedComb) {
	    hisNumSeedCombinations_->Fill(num) ;
	  }
	  // Same again, but this time only considering seeds the r-z filters defined as "good".
	  const vector<unsigned int>  numGoodSeedComb = htPair.getRZfilters().numGoodSeedCombsPerTrk();
	  for (const unsigned int& num : numGoodSeedComb) {
	    hisNumGoodSeedCombinations_->Fill(num) ;
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
			hisCorrelationZTrk_->Fill(r);
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
}

//=== Book histograms studying track candidates found by Hough Transform.

void Histos::bookTrackCands() {

  // Book histograms for studying freak, extra large events.
  this->bookStudyBusyEvents();

  // Now book histograms for studying tracking in general.

  TFileDirectory inputDir = fs_->mkdir("TrackCands");

  // Count tracks in various ways (including/excluding duplicates, excluding fakes ...)
  profNumTrackCands_  = inputDir.make<TProfile>("NumTrackCands","; class; N. of tracks in tracker",8,0.5,8.5);
  profNumTrackCands_->GetXaxis()->SetBinLabel(8,"reco tracks including fakes before rz filter");
  profNumTrackCands_->GetXaxis()->SetBinLabel(7,"TP for eff recoed");
  profNumTrackCands_->GetXaxis()->SetBinLabel(6,"TP recoed");
  profNumTrackCands_->GetXaxis()->SetBinLabel(5,"TP recoed x #eta sector dups");
  profNumTrackCands_->GetXaxis()->SetBinLabel(4,"TP recoed x sector dups");
  profNumTrackCands_->GetXaxis()->SetBinLabel(3,"TP recoed x sector x r-#phi cell dups");
  profNumTrackCands_->GetXaxis()->SetBinLabel(2,"TP recoed x sector x cell dups");
  profNumTrackCands_->GetXaxis()->SetBinLabel(1,"reco tracks including fakes");
  profNumTrackCands_->LabelsOption("d");

  unsigned int nPhi = settings_->numPhiSectors();
  unsigned int nEta = settings_->numEtaRegions();
  profNumTracksVsEta_ = inputDir.make<TProfile>("NumTracksVsEta","; #eta region; No. of tracks in tracker", nEta, -0.5, nEta - 0.5);
  float maxAbsQoverPt = 1./settings_->houghMinPt(); // Max. |q/Pt| covered by  HT array.
  hisNumTracksVsQoverPt_ = inputDir.make<TH1F>("NumTracksVsQoverPt","; Q/Pt; No. of tracks in tracker",100, -maxAbsQoverPt, maxAbsQoverPt);

  hisNumTrksPerSect_ = inputDir.make<TH1F>("NumTrksPerSect","; No. tracks per sector;",100,-0.5,99.5);
  hisNumTrksPerOct_  = inputDir.make<TH1F>("NumTrksPerOct", "; No. tracks per octant;",200,-0.5,199.5);

  // Count stubs per event assigned to tracks (determines HT data output rate)

  profStubsOnTracks_ = inputDir.make<TProfile>("StubsOnTracks","; ; No. of stubs on tracks per event",1,0.5,1.5);
  profStubsOnTracksVsEta_ = inputDir.make<TProfile>("StubsOnTracksVsEta","; #eta region; No. of stubs on tracks per event", nEta, -0.5, nEta - 0.5); 
  hisStubsOnTracksPerSect_ = inputDir.make<TH1F>("StubsOnTracksPerSect","; No. of stubs on tracks per sector", 500,-0.5,499.5); 
  hisStubsOnTracksPerOct_  = inputDir.make<TH1F>("StubsOnTracksPerOct","; No. of stubs on tracks per octant", 1000,-0.5,999.5); 
  hisUniqueStubsOnTrksPerSect_ = inputDir.make<TH1F>("UniqueStubsOnTrksPerSect","; No. of unique stubs on tracks per sector", 500,-0.5,499.5); 
  hisUniqueStubsOnTrksPerOct_  = inputDir.make<TH1F>("UniqueStubsOnTrksPerOct","; No. of unique stubs on tracks per octant", 500,-0.5,499.5); 

  hisStubsPerTrack_ = inputDir.make<TH1F>("StubsPerTrack",";No. of stubs per track;",50,-0.5,49.5);
  hisLayersPerTrack_ = inputDir.make<TH1F>("LayersPerTrack",";No. of layers with stubs per track;",20,-0.5,19.5);
  hisPSLayersPerTrack_ = inputDir.make<TH1F>("PSLayersPerTrack",";No. of PS layers with stubs per track;",20,-0.5,19.5);
  hisLayersPerTrueTrack_ = inputDir.make<TH1F>("LayersPerTrueTrack",";No. of layers with stubs per genuine track;",20,-0.5,19.5);
  hisPSLayersPerTrueTrack_ = inputDir.make<TH1F>("PSLayersPerTrueTrack",";No. of PS layers with stubs per genuine track;",20,-0.5,19.5);

  // Checks if tracks have too many stubs to be stored in memory in each cell.
  profExcessStubsPerTrackVsPt_ = inputDir.make<TProfile>("ExcessStubsPerTrackVsPt",";q/Pt; Prob. of too many stubs per track",16,0.,maxAbsQoverPt);

  hisFracMatchStubsOnTracks_ = inputDir.make<TH1F>("FracMatchStubsOnTracks","; Fraction of stubs on tracks matching best TP;",101,-0.005,1.005);

  // See how far stubs lie from true trajectory in r-z plane.
  hisDeltaPhiRtruePS_ = inputDir.make<TH1F>("DeltaPhiRtruePS","PS modules; Dist. of true stubs from true trajectory in r*phi;",100,-0.25,0.25);
  hisDeltaRorZtruePS_ = inputDir.make<TH1F>("DeltaRorZtruePS","PS modules; Dist. of true stubs from true trajectory in r-z;",100,-10,10);
  hisDeltaPhiRtrue2S_ = inputDir.make<TH1F>("DeltaPhiRtrue2S","2S modules; Dist. of true stubs from true trajectory in r*phi;",100,-0.25,0.25);
  hisDeltaRorZtrue2S_ = inputDir.make<TH1F>("DeltaRorZtrue2S","2S modules; Dist. of true stubs from true trajectory in r-z;",100,-10,10);
  hisDeltaPhiRfakePS_ = inputDir.make<TH1F>("DeltaPhiRfakePS","PS modules; Dist. of fake stubs from true trajectory in r*phi;",100,-0.25,0.25);
  hisDeltaRorZfakePS_ = inputDir.make<TH1F>("DeltaRorZfakePS","PS modules; Dist. of fake stubs from true trajectory in r-z;",100,-10,10);
  hisDeltaPhiRfake2S_ = inputDir.make<TH1F>("DeltaPhiRfake2S","2S modules; Dist. of fake stubs from true trajectory in r*phi;",100,-0.25,0.25);
  hisDeltaRorZfake2S_ = inputDir.make<TH1F>("DeltaRorZfake2S","2S modules; Dist. of fake stubs from true trajectory in r-z;",100,-10,10);
  profNsigmaPhiRvsInvPt_ = inputDir.make<TProfile>("NsigmaPhiRvsInvPt","; 1/Pt; Num #sigma of true stubs from true trajectory",16,0.,maxAbsQoverPt, 0., 10.); 
  profNsigmaPhiRvsFracDist_ = inputDir.make<TProfile>("NsigmaPhiRvsFracDist","; Fractional position in tracker; Num #sigma of true stubs from true trajectory",22,0.,1.1, 0., 10.); 
  profFracTrueStubsVsLayer_ = inputDir.make<TProfile>("FracTrueStubsVsLayer",";Layer ID; fraction of true stubs",30,0.5,30.5);

  // Check how much stub bend differs from predicted one.
  hisDeltaBendTrue_ = inputDir.make<TH1F>("DeltaBendTrue","True stubs; stub bend minus true bend / resolution;",100,-2.,2.);
  hisDeltaBendFake_ = inputDir.make<TH1F>("DeltaBendFake","Fake stubs; stub bend minus true bend / resolution;",100,-2.,2.);

  // Study duplication of tracks within HT.
  profDupTracksVsTPeta_ = inputDir.make<TProfile>("DupTracksVsTPeta" ,"; #eta; Number of duplicate tracks in individual HT array;",30,-3.0,3.0);

  // Histos for tracking efficiency vs. TP kinematics
  hisTPinvptForEff_     = inputDir.make<TH1F>("TPinvptForEff" ,"; 1/Pt of TP (used for effi. measurement);",16,0.,maxAbsQoverPt);
  hisRecoTPinvptForEff_ = inputDir.make<TH1F>("RecoTPinvptForEff" ,"; 1/Pt of TP (used for effi. measurement);",16,0.,maxAbsQoverPt);
  hisTPetaForEff_       = inputDir.make<TH1F>("TPetaForEff","; #eta of TP (used for effi. measurement);",20,-3.,3.);
  hisRecoTPetaForEff_   = inputDir.make<TH1F>("RecoTPetaForEff","; #eta of TP (used for effi. measurement);",20,-3.,3.);
  hisTPphiForEff_       = inputDir.make<TH1F>("TPphiForEff","; #phi of TP (used for effi. measurement);",20,-M_PI,M_PI);
  hisRecoTPphiForEff_   = inputDir.make<TH1F>("RecoTPphiForEff","; #phi of TP (used for effi. measurement);",20,-M_PI,M_PI);

  // Histo for efficiency to reconstruct track perfectly (no incorrect hits).
  hisPerfRecoTPinvptForEff_ = inputDir.make<TH1F>("PerfRecoTPinvptForEff" ,"; 1/Pt of TP (used for perf. effi. measurement);",16,0.,maxAbsQoverPt);
  hisPerfRecoTPetaForEff_   = inputDir.make<TH1F>("PerfRecoTPetaForEff","; #eta of TP (used for perf. effi. measurement);",20,-3.,3.);

  // Histos for algorithmic tracking efficiency vs. TP kinematics
  hisTPinvptForAlgEff_     = inputDir.make<TH1F>("TPinvptForAlgEff" ,"; 1/Pt of TP (used for alg. effi. measurement);",16,0.,maxAbsQoverPt);
  hisRecoTPinvptForAlgEff_ = inputDir.make<TH1F>("RecoTPinvptForAlgEff" ,"; 1/Pt of TP (used for alg. effi. measurement);",16,0.,maxAbsQoverPt);
  hisTPetaForAlgEff_       = inputDir.make<TH1F>("TPetaForAlgEff","; #eta of TP (used for alg. effi. measurement);",20,-3.,3.);
  hisRecoTPetaForAlgEff_   = inputDir.make<TH1F>("RecoTPetaForAlgEff","; #eta of TP (used for alg. effi. measurement);",20,-3.,3.);
  hisTPphiForAlgEff_       = inputDir.make<TH1F>("TPphiForAlgEff","; #phi of TP (used for alg. effi. measurement);",20,-M_PI,M_PI);
  hisRecoTPphiForAlgEff_   = inputDir.make<TH1F>("RecoTPphiForAlgEff","; #phi of TP (used for alg. effi. measurement);",20,-M_PI,M_PI);

  // Histo for efficiency to reconstruct track perfectly (no incorrect hits).
  hisPerfRecoTPinvptForAlgEff_ = inputDir.make<TH1F>("PerfRecoTPinvptForAlgEff" ,"; 1/Pt of TP (used for perf. alg. effi. measurement);",16,0.,maxAbsQoverPt);
  hisPerfRecoTPetaForAlgEff_   = inputDir.make<TH1F>("PerfRecoTPetaForAlgEff","; #eta of TP (used for perf. alg. effi. measurement);",20,-3.,3.);

  // Histos for algorithmic tracking efficiency vs. TP production point
  hisTPd0ForAlgEff_      = inputDir.make<TH1F>("TPd0ForAlgEff" ,"; d0 of TP (used for alg. effi. measurement);",50,0.,1.);
  hisRecoTPd0ForAlgEff_  = inputDir.make<TH1F>("RecoTPd0ForAlgEff" ,"; d0 of TP (used for alg. effi. measurement);",50,0.,1.);
  hisTPz0ForAlgEff_      = inputDir.make<TH1F>("TPz0ForAlgEff" ,"; z0 of TP (used for alg. effi. measurement);",50,0.,25.);
  hisRecoTPz0ForAlgEff_  = inputDir.make<TH1F>("RecoTPz0ForAlgEff" ,"; z0 of TP (used for alg. effi. measurement);",50,0.,25.);

  // Histos for algorithmic tracking efficiency vs sector number (to check if looser cuts are needed in certain regions)

  hisTPphisecForAlgEff_      = inputDir.make<TH1F>("TPphisecForAlgEff" ,"; #phi sectorof TP (used for alg. effi. measurement);",nPhi,-0.5,nPhi-0.5);
  hisRecoTPphisecForAlgEff_  = inputDir.make<TH1F>("RecoTPphisecForAlgEff" ,"; #phi sector of TP (used for alg. effi. measurement);",nPhi,-0.5,nPhi-0.5);
  hisTPetasecForAlgEff_      = inputDir.make<TH1F>("TPetasecForAlgEff" ,"; #eta sector of TP (used for alg. effi. measurement);",nEta,-0.5,nEta-0.5);
  hisRecoTPetasecForAlgEff_  = inputDir.make<TH1F>("RecoTPetasecForAlgEff" ,"; #eta sector of TP (used for alg. effi. measurement);",nEta,-0.5,nEta-0.5);

  // Histo for efficiency to reconstruct tracks perfectly (no incorrect hits).
  hisPerfRecoTPphisecForAlgEff_  = inputDir.make<TH1F>("PerfRecoTPphisecForAlgEff" ,"; #phi sector of TP (used for perf. alg. effi. measurement);",nPhi,-0.5,nPhi-0.5);
  hisPerfRecoTPetasecForAlgEff_  = inputDir.make<TH1F>("PerfRecoTPetasecForAlgEff" ,"; #eta sector of TP (used for perf. alg. effi. measurement);",nEta,-0.5,nEta-0.5);

  // Histos of track parameter resolution
  hisQoverPtRes_ = inputDir.make<TH1F>("QoverPtRes","; track resolution in q/Pt", 100,-0.06,0.06);
  hisPhi0Res_    = inputDir.make<TH1F>("Phi0Res",   "; track resolution in #phi0",100,-0.04,0.04);
  hisEtaRes_     = inputDir.make<TH1F>("EtaRes",    "; track resolution in #eta", 100,-1.0,1.0);
  hisZ0Res_      = inputDir.make<TH1F>("Z0Res",     "; track resolution in z0",   100,-10.0,10.0);
  hisTanLambdaRes_     = inputDir.make<TH1F>("TanLambdaRes",    "; track resolution in tan #lambda", 100,-1.0,1.0);
  // For those tracking particles causing the algorithmic efficiency to be below 100%, plot a flag indicating why.
  hisRecoFailureReason_ = inputDir.make<TH1F>("RecoFailureReason","; Reason TP (used for alg. effi.) not reconstructed;",1,-0.5,0.5); 
  //hisRecoFailureLayer_ = inputDir.make<TH1F>("RecoFailureLayer","; Layer ID of lost stubs on unreconstructed TP;",30,-0.5,29.5);

  hisWrongSignStubRZ_pBend_ = inputDir.make<TH2F>("RZ of stubs with positive bend, but with wrong sign","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130);
  hisWrongSignStubRZ_nBend_ = inputDir.make<TH2F>("RZ of stubs with negative bend, but with wrong sign","; z (cm); radius (cm); No. stubs in tracker",1000,-280,280,1000,0,130);

  hisNumStubsOnLayer_ = inputDir.make<TH1F>("NumStubsOnLayer","; Layer occupancy;",16,1,17); 
}

//=== Fill histograms studying track candidates found by Hough Transform.

void Histos::fillTrackCands(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs) {

  // Fill histograms for studying freak, extra large events.
  this->fillStudyBusyEvents(inputData, mSectors, mHtPairs);

  // Now fill histograms for studying tracking in general.

  const vector<TP>&  vTPs = inputData.getTPs();

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
          hisNumStubsOnLayer_->Fill( l.second );
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
      hisNumTrksPerSect_->Fill(htPair.numTrackCands3D()); 
      unsigned int iOctant = floor(iPhiSec*numPhiOctants/(settings_->numPhiSectors())); // phi octant number
      nTracksInOctant[iOctant] += htPair.numTrackCands3D();
      // Note number of tracks in this eta region (summed over all phi).
      nTracksInEtaReg += htPair.numTrackCands3D();
      nTracks_sf += htPair.getRphiHT().numTrackCands2D();
      if (settings_->debug() == 1 && htPair.numTrackCands3D() > 0) cout<<"Sector ("<<iPhiSec<<","<<iEtaReg<<") has ntracks = "<<htPair.numTrackCands3D()<<endl;
    }
    nTracks += nTracksInEtaReg;
    profNumTracksVsEta_->Fill(iEtaReg, nTracksInEtaReg);
  }
  profNumTrackCands_->Fill(1.0, nTracks); // Plot mean number of tracks/event.
  profNumTrackCands_->Fill(8.0, nTracks_sf);
  for (unsigned int k = 0; k < numPhiOctants; k++) {
    hisNumTrksPerOct_->Fill(nTracksInOctant[k]); // Plots tracks in each phi octant.
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
      hisStubsOnTracksPerSect_->Fill(nStubsOnTrksInSec); // Number of stubs assigned to tracks in this sector.
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
      hisUniqueStubsOnTrksPerSect_->Fill(uniqueStubsOnTracksInSector.size());
    }
    nStubsOnTracks += nStubsOnTracksInEtaReg;
    profStubsOnTracksVsEta_->Fill(iEtaReg, nStubsOnTracksInEtaReg);
  }
  profStubsOnTracks_->Fill(1.0, nStubsOnTracks);
  for (unsigned int k = 0; k < numPhiOctants; k++) {
    hisStubsOnTracksPerOct_->Fill(nStubsOnTracksInOctant[k]); // Plots stubs on tracks in each phi octant.
    hisUniqueStubsOnTrksPerOct_->Fill(uniqueStubsOnTracksInOctant[k].size());
  }

  // Plot q/pt spectrum of track candidates, and number of stubs/track.
  for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {
      const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);
      // Loop over all reconstructed tracks in this sector
      const vector<L1track3D>& vecTrk3D = htPair.trackCands3D();
      for (const L1track3D& trk : vecTrk3D) {
	hisNumTracksVsQoverPt_->Fill(trk.qOverPt()); // Plot reconstructed q/Pt of track cands.
	hisStubsPerTrack_->Fill(trk.getNumStubs());  // Stubs per track.
	// For genuine tracks, check how often they have too many stubs to be stored in cell memory. (Perhaps worse for high Pt particles in jets?).
	const TP* tp = trk.getMatchedTP();
	if (tp != nullptr) {
	  if (tp->useForAlgEff()) profExcessStubsPerTrackVsPt_->Fill(1./tp->pt(), trk.getNumStubs() > 16);
	}
	hisLayersPerTrack_->Fill(trk.getNumLayers()); // Number of reduced layers with stubs per track.
	hisPSLayersPerTrack_->Fill( Utility::countLayers(settings_, trk.getStubs(), false, true) ); // Number of reduced PS layers with stubs per track.
	// Also plot just for genuine tracks.
	if (tp != nullptr && tp->useForAlgEff()) {
  	  hisLayersPerTrueTrack_->Fill(trk.getNumLayers()); // Number of reduced layers with stubs per track.
	  hisPSLayersPerTrueTrack_->Fill( Utility::countLayers(settings_, trk.getStubs(), false, true) ); // Number of reduced PS layers with stubs per track.
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
            hisFracMatchStubsOnTracks_->Fill( trk.getPurity() );
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
		  hisDeltaPhiRtruePS_->Fill(deltaPhiR);
		  hisDeltaRorZtruePS_->Fill(deltaRorZ);
		} else {
		  //if (tp->pt() > 20. && s->assocTP() != nullptr && fabs(tp->d0()) < 0.01) {
		    //		    if ( ! (s->barrel()) ) cout<<"RATS "<<s->width()<<" "<<s->iphi()<<" "<<s->nstrip()<<" "<<s->r()<<" "<<s->z()<<endl;
		    //if ( ! (s->barrel()) ) cout<<"DELTAPHI "<<stripAngle<<" "<<(tp->trkRAtStub(s) - s->r())<<" "<<phiCorr<<" : "<<deltaPhiR<<" "<<deltaPhiR1<<" "<<deltaPhiR2<<endl;
		  //}
		  hisDeltaPhiRtrue2S_->Fill(deltaPhiR);
		  hisDeltaRorZtrue2S_->Fill(deltaRorZ);
		}
		// More detailed plots for true stubs to study effect of multiple scattering.
		float sigPerp = s->sigmaPerp(); // detector resolution
                float ptThresh = 40.; // threshold where scattering dominates detector resolution 
		float relpos = s->barrel()  ?   s->r() / settings_->trackerOuterRadius()  :  fabs(s->z()) / settings_->trackerHalfLength();
                float sigmaScat = 0.01 * (ptThresh/tp->pt()) * pow(relpos, 1.5);
		sigPerp += sigmaScat; // Estimated resolution allowing for scattering.
		profNsigmaPhiRvsInvPt_->Fill(1./tp->pt(), fabs(deltaPhiR)/sigPerp);
		profNsigmaPhiRvsFracDist_->Fill(relpos,   fabs(deltaPhiR)/sigPerp);
	      } else {
		if (s->psModule()) {
		  hisDeltaPhiRfakePS_->Fill(deltaPhiR);
		  hisDeltaRorZfakePS_->Fill(deltaRorZ);
		} else {
		  hisDeltaPhiRfake2S_->Fill(deltaPhiR);
		  hisDeltaRorZfake2S_->Fill(deltaRorZ);
		}
	      }

	      // Fraction of wrong stubs vs. tracker layer.
	      profFracTrueStubsVsLayer_->Fill(s->layerId(), trueStub);

	      // Check how much stub bend differs from predicted one, relative to nominal bend resolution.
	      float diffBend = (s->qOverPt() - trk.qOverPt()) / s->qOverPtOverBend();
	      if (trueStub) {
		hisDeltaBendTrue_->Fill(diffBend/s->bendRes());
	      } else {
		hisDeltaBendFake_->Fill(diffBend/s->bendRes());
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
          profDupTracksVsTPeta_->Fill(tp.eta(), nTrk); // Study duplication of tracks within an individual HT array.
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
  profNumTrackCands_->Fill(7.0, nRecoedTPsForEff);
  // Plot number of TPs that are reconstructed. 
  profNumTrackCands_->Fill(6.0, nRecoedTPs);
  // Plot number of TPs that are reconstructed. (Count +1 for each eta sector they are reconstructed in).
  profNumTrackCands_->Fill(5.0, nEtaSecsMatchingTPs);
  // Plot number of TPs that are reconstructed. (Count +1 for each etaxphisector they are reconstructed in).
  profNumTrackCands_->Fill(4.0, nSecsMatchingTPs);
  // Plot number of TP that are reconstructed. (Ditto, but now multiplying by duplicate cells in r-phi HT).
  profNumTrackCands_->Fill(3.0, nTrksMatchingTPsIgnoringRzDups);
  // Plot number of TP that are reconstructed. (Ditto, but now multiplying by duplicate cells in r-phi x r-z HTs).
  profNumTrackCands_->Fill(2.0, nTrksMatchingTPs);

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
      hisTPinvptForEff_->Fill(1./tp.pt());
      hisTPetaForEff_->Fill(tp.eta());
      hisTPphiForEff_->Fill(tp.phi0());

      if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
        hisTPinvptForAlgEff_->Fill(1./tp.pt());
        hisTPetaForAlgEff_->Fill(tp.eta());
        hisTPphiForAlgEff_->Fill(tp.phi0());
	// Plot also production point of all good TP.
        hisTPd0ForAlgEff_->Fill(fabs(tp.d0()));
        hisTPz0ForAlgEff_->Fill(fabs(tp.z0()));
	// Plot sector nunber.
        hisTPphisecForAlgEff_->Fill(iPhiSec_TP);
        hisTPetasecForAlgEff_->Fill(iEtaReg_TP);
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
	hisRecoTPinvptForEff_->Fill(1./tp.pt());
	hisRecoTPetaForEff_->Fill(tp.eta());
	hisRecoTPphiForEff_->Fill(tp.phi0());
	// Also plot efficiency to perfectly reconstruct the track (no fake hits)
	if (tpRecoedPerfect) {
	  hisPerfRecoTPinvptForEff_->Fill(1./tp.pt());
	  hisPerfRecoTPetaForEff_->Fill(tp.eta());
	}
	if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
	  hisRecoTPinvptForAlgEff_->Fill(1./tp.pt());
	  hisRecoTPetaForAlgEff_->Fill(tp.eta());
	  hisRecoTPphiForAlgEff_->Fill(tp.phi0());
	  // Plot also production point of all good reconstructed TP.
	  hisRecoTPd0ForAlgEff_->Fill(fabs(tp.d0()));
	  hisRecoTPz0ForAlgEff_->Fill(fabs(tp.z0()));
	  // Plot sector number to understand if looser cuts are needed in certain eta regions.
	  hisRecoTPphisecForAlgEff_->Fill(iPhiSec_TP);
	  hisRecoTPetasecForAlgEff_->Fill(iEtaReg_TP);
	  // Also plot efficiency to perfectly reconstruct the track (no fake hits)
	  if (tpRecoedPerfect) {
	    hisPerfRecoTPinvptForAlgEff_->Fill(1./tp.pt());
	    hisPerfRecoTPetaForAlgEff_->Fill(tp.eta());
  	    hisPerfRecoTPphisecForAlgEff_->Fill(iPhiSec_TP);
	    hisPerfRecoTPetasecForAlgEff_->Fill(iEtaReg_TP);
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
	    hisQoverPtRes_->Fill(trk->qOverPt() - tp.qOverPt());
	    hisPhi0Res_->Fill(reco::deltaPhi(trk->phi0(), tp.phi0()));
	    hisEtaRes_->Fill(trk->eta() - tp.eta());
	    hisZ0Res_->Fill(trk->z0() - tp.z0());
          }
	}
      }
    }
  }

  // Diagnose reason why not all viable tracking particles were reconstructed.
  const map<const TP*, string> diagnosis = this->diagnoseTracking(inputData, mSectors, mHtPairs);
  for (const auto& iter: diagnosis) {
    hisRecoFailureReason_->Fill(iter.second.c_str(), 1.); // Stores flag indicating failure reason.
  }
}

void Histos::bookVertexReconstruction(){
  TFileDirectory inputDir = fs_->mkdir("VertexReconstruction");

  // Plot number of reconstructed vertices against number of true vertices
  hisNoRecoVertices_                 = inputDir.make<TH1F>("hisNoRecoVertices_","No. reconstructed Vertices; No. reco vertices; Events",50,0,50);
   
  // Plot number of real tracks originating from the PV
  hisNoRealPVTracks_                 = inputDir.make<TH1F>("hisNoRealPVTracks_", "No. Real Tracks PV; No. real tracks PV; Events",20,0,40);
  // Plot number of reco tracks originating from the PV
  hisNoRecoPVTracks_                 = inputDir.make<TH1F>("hisNoRecoPVTracks_", "No. Reco Tracks PV; No. reco tracks PV; Events",20,0,40);
    
  // Plot Pt of the real PV
  hisRealPVPt_                       = inputDir.make<TH1F>("hisRealPVPt_", "Pt Real PV; Pt real PV; Events",100,0,300);
  // Plot Pt of the reco PV
  hisRecoPVPt_                       = inputDir.make<TH1F>("hisRecoPVPt_", "Pt Reco PV; Pt recvo PV; Events",100,0,300); 
  // Plot 2D hist of real no tracks vs real pT
  hisRealTracksvspT_                 = inputDir.make<TH2F>("hisRealTracksvspT_","No Real Tracks vs Real pT; real pT; No real tracks",60,0,300,8,0,40);
  // Plot 2D hist of reco no tracks vs real pT  
  hisRecoTracksvspT_                 = inputDir.make<TH2F>("hisRecoTracksvspT_","No. Tracks vs pT; pT [GeV]; No. tracks",60,0,300,8,0,40);
  // Plot number of reconstructed vertices against number of true vertices
  hisNoPileUpVertices_               = inputDir.make<TH1F>("hisNoPileUpVertices_","No. pile-up Vertices; No. pile-up vertices; Events",50,0,50);
  hisNoRecoVsNoTruePileUpVertices_   = inputDir.make<TH2F>("hisNoRecoVsNoTruePileUpVertices_","No. reconstructed pile-up vertices vs. no. true pile-up vertices; No. reco pile-up vertices; No. true pile-up vertices",50,0,50,50,0,50);

  //Fill the histogram with the chi^2/dof for the reco tracks
  hisChi2dofTracks_ = inputDir.make<TH1F>("hisChi2dofTracks_", "Chi2dof; chi2dof; Events", 20,0,10);

  //Number of Layers hit by tracks
    
  numLayersRealTracks_ = inputDir.make<TH1F>("numLayersRealTracks_", "No. Real Layers; No. real layers; Events",4,4,8);
  numLayersRecoTracks_ = inputDir.make<TH1F>("numLayersRecoTracks_", "No. Reco Layers; No. reco layers; Events",4,4,8);
  
  //Number of Layers hit by tracks
    
  numPSLayersRealTracks_ = inputDir.make<TH1F>("numPSLayersRealTracks_", "No. Real PS Layers; No. real PS layers; Events",10,0,10);
  numPSLayersRecoTracks_ = inputDir.make<TH1F>("numPSLayersRecoTracks_", "No. Reco PS Layers; No. reco PS layers; Events",10,0,10);
    
  //Id of the layers hit by tracks
  
  realTracksLayersId_ = inputDir.make<TH1F>("realTracksLayersId_", "Real Layer ID; Real layer ID; Events",2,0,2);
  recoTracksLayersId_ = inputDir.make<TH1F>("recoTracksLayersId_", "Reco layer ID; Reco layer ID; Events",2,0,2);
    
    
  //No tracks (50% in barrel)
  hisNoRealTracksPerc50_ = inputDir.make<TH1F>("hisNoRealTracksPerc50_", "No. Real Tracks; No. real tracks; Events", 20, 0, 40);
  hisNoRecoTracksPerc50_ = inputDir.make<TH1F>("hisNoRecoTracksPerc50_", "No. Reco Tracks; No. reco tracks; Events", 20, 0, 40);
    
  //No tracks (90% in barrel) 
  hisNoRealTracksPerc90_ = inputDir.make<TH1F>("hisNoRealTracksPerc90_", "No. Real Tracks; No. real tracks; Events", 20, 0, 40);
  hisNoRecoTracksPerc90_ = inputDir.make<TH1F>("hisNoRecoTracksPerc90_", "No. Reco Tracks; No. reco tracks; Events", 20, 0, 40);
    
  //Z0 (50% in barrel) 
  hisz0RealPerc50_ =inputDir.make<TH1F>("hisz0RealPerc50_", "z0 PV; z0; Events", 50,-25,25);
  hisz0RecoPerc50_ =inputDir.make<TH1F>("hisz0RecoPerc50_", "z0 PV; z0; Events", 50,-25,25);

  //Z0 (90% in barrel) 
  hisz0RecoPerc90_ =inputDir.make<TH1F>("hisz0RecoPerc90_", "z0 PV; z0; Events", 50,-25,25);
  hisz0RealPerc90_ =inputDir.make<TH1F>("hisz0RealPerc90_", "z0 PV; z0; Events", 50,-25,25);

  //pT (50% in barrel) 
  hisRealpTPerc50_ =inputDir.make<TH1F>("hisRealpTPerc50_", "pT PV; pT; Events", 100,0,300);
  hisRecopTPerc50_ =inputDir.make<TH1F>("hisRecopTPerc50_", "pT PV; pT; Events", 100,0,300);

  //pT (90% in barrel) 
  hisRecopTPerc90_ =inputDir.make<TH1F>("hisRecopTPerc90_", "pT PV; pT; Events", 100,0,300);
  hisRealpTPerc90_ =inputDir.make<TH1F>("hisRealpTPerc90_", "pT PV; pT; Events", 100,0,300);
    
    
  //No Tracks with eta below 15
  hisNoRealTracksGoodEta_ =inputDir.make<TH1F>("hisNoRealTracksGoodEta", "No Real Tracks with eta below 15; No. real tracks; Events",40,0,40);
  hisNoRecoTracksGoodEta_ =inputDir.make<TH1F>("hisNoRecoTracksGoodEta", "No Real Tracks with eta below 15; No. real tracks; Events",40,0,40);
    
   
  //No Tracks with eta above 15
  hisNoRealTracksBadEta_ =inputDir.make<TH1F>("hisNoRealTracksBadEta", "No Real Tracks with eta above 15; No. real tracks; Events",40,0,40);
  hisNoRecoTracksBadEta_ =inputDir.make<TH1F>("hisNoRecoTracksBadEta", "No Real Tracks with eta above 15; No. real tracks; Events",40,0,40);
  
  // *** Vertex Reconstruction algorithm Plots ***
  hisRecoVertexZ0Resolution_         = inputDir.make<TH1F>("hisRecoVertexZ0Resolution","Reconstructed primary vertex z_{0} resolution; z_{0} Resolution [cm]; Counts", 100, 0., 1.);
  hisUnmatchedVertexZ0distance_            = inputDir.make<TH1F>("hisUnmatchedVertexZ0distance_"," Unmatched primary vertex z_{0} - true PV z_{0}; |z_{0}^{reco} - z_{0}^{true}| [cm]; Counts", 200, 1., 5.);

  hisRecoVertexPTResolution_         = inputDir.make<TH1F>("hisRecoVertexPTResolution","Reconstructed primary vertex p_{T} relative resolution; p_{T} relative Resolution; Counts", 100, 0, 1.);
  hisRecoVertexPTResolutionVsTruePt_ = inputDir.make<TProfile>("hisRecoVertexPTResolutionVsTruePt","Reconstructed primary vertex relative p_{T} resolution vs. True Pt; True p_{T}; p_{T} Resolution [GeV]", 100, 0, 500);
  hisRecoVertexPTVsTruePt_ = inputDir.make<TH2F>("hisRecoVertexPtVsTruePt_","Reconstructed primary vertex p_{T}  vs. True Pt; p_{T} [GeV]; True p_{T}", 100, 0, 500.,100, 0, 500.);
  hisRecoVertexMETVsTrueMET_ = inputDir.make<TH2F>("hisRecoVertexMETVsTrueMET_","Reconstructed primary vertex MET vs. True MET; MET [GeV]; True MET", 100, 0, 500.,100, 0, 500.);
  hisNoTracksFromPrimaryVertex_      = inputDir.make<TH2F>("hisNoTracksFromPrimaryVertex_","No. of Tracks from Primary Vertex (Reco vs. Truth); no. Reco Tracks in PV; no. Truth Tracks in PV ",50,0,50,50,0,50);
  hisNoTrueTracksFromPrimaryVertex_  = inputDir.make<TH2F>("hisNoTrueTracksFromPrimaryVertex_","No. of Matched Tracks from Primary Vertex (Reco vs. Truth); no. Reco Tracks in PV; no. Truth Tracks in PV ",50,0,50,50,0,50);
  hisRecoPrimaryVertexZ0width_       = inputDir.make<TH1F>("hisRecoPrimaryVertexZ0width_", "Reconstructed primary vertex z_{0} width", 100, 0, 0.5);
  hisRatioMatchedTracksInPV_         = inputDir.make<TH1F>("hisRatioMatchedTracksInPV", "Primary vertex matching ratio ", 100, 0, 1.1);
  hisFakeTracksRateInPV_             = inputDir.make<TH1F>("hisFakeTracksRateInPV", "Percentage of fake tracks in recontructed primary vertex", 100, 0, 1.1);
  hisTrueTracksRateInPV_             = inputDir.make<TH1F>("hisTrueTracksRateInPV", "Percentage of true tracks in recontructed primary vertex", 100, 0, 1.1);
  hisUnmatchZ0distance_              = inputDir.make<TH1F>("hisUnmatchZ0distance_", "z0 distance from reconstructed privary vertex of unmatched tracks; |z_{0}^{track} - z_{0}^{vertex}|; no. L1 Tracks", 100, 0, 5.);
  hisUnmatchZ0MinDistance_              = inputDir.make<TH1F>("hisUnmatchZ0MinDistance_", "z0 distance from the closest track in reconstructed privary vertex of unmatched tracks; |z_{0}^{track} - z_{0}^{PV track}|; no. L1 Tracks", 100, 0, 5.);
  hisUnmatchPt_                      = inputDir.make<TH1F>("hisUnmatchPt_", "Transverse momentum of unmatched PV tracks; p_{T} [GeV/c]; no. L1 Tracks", 100, 0, 100.);
  hisUnmatchEta_                      = inputDir.make<TH1F>("hisUnmatchEta_", "#eta of unmatched PV tracks; #eta; no. L1 Tracks", 100, 2.4, 2.4);
  hisUnmatchTruePt_                      = inputDir.make<TH1F>("hisUnmatchTruePt_", "True transverse momentum of unmatched PV tracks; p_{T} [GeV/c]; no. L1 Tracks", 100, 0, 100.);
  hisUnmatchTrueEta_                      = inputDir.make<TH1F>("hisUnmatchTrueEta_", "True #eta of unmatched PV tracks; #eta; no. L1 Tracks", 100, 2.4, 2.4);
  hisUnmatchedPVtracks_              = inputDir.make<TH1F>("hisUnmatchedPVtracks_", " No. of tracks from primary collision that are misassigned; No. misassigned Tracks, No. Events ", 100, 0, 100);
  hisRecoPrimaryVertexVsTrueZ0_      = inputDir.make<TH1F>("hisRecoPrimaryVertexVsTrueZ0_","No. of reconstructed primary vertices per true z_{0}; true z_{0} [cm]; No. Vertices", 50, -25.,25.);
  hisRecoPrimaryVertexResolutionVsTrueZ0_      = inputDir.make<TProfile>("hisRecoPrimaryVertexResolutionVsTrueZ0_","No. of reconstructed primary vertices per true z_{0}; true z_{0} [cm]; Resolution [mm]", 100, -25.,25.);
  hisRecoVertexMET_ = inputDir.make<TH1F>("hisRecoVertexMET_","; L1TrkMET (GeV) ; Entries", 50, 0, 200);
  hisRecoVertexPT_ = inputDir.make<TH1F>("hisRecoVertexPT_","; #Sigma_{vtx} p_{T} (GeV) ; Entries", 50, 0, 200);
  hisRecoPileUpVertexPT_ = inputDir.make<TH1F>("hisRecoPileUpVertexPT_","; #Sigma_{vtx} p_{T} (GeV) ; Entries", 50, 0, 200);
  hisRecoVertexOffPT_ = inputDir.make<TH1F>("hisRecoVertexOffPT_","; #Sigma_{vtx} p_{T} (GeV) ; Entries", 50, 0, 200);


  // *** TDR vertex reconstruction algorithm ***
  hisTDRVertexZ0Resolution_            = inputDir.make<TH1F>("hisTDRVertexZ0Resolution","TDR Reconstructed primary vertex z_{0} resolution; z_{0} Resolution [cm]; Counts", 100, 0., 1);
  hisTDRUnmatchedVertexZ0distance_            = inputDir.make<TH1F>("hisTDRUnmatchedVertexZ0distance_","TP Unmatched primary vertex z_{0} - true PV z_{0}; |z_{0}^{TP} - z_{0}^{true}| [cm]; Counts", 200, 1., 5.);

  hisTDRVertexPTResolution_            = inputDir.make<TH1F>("hisTDRVertexPTResolution","TDR Reconstructed primary vertex p_{T} relative resolution; p_{T} relative resolution; Counts", 100, 0, 1.);
  hisTDRVertexPTResolutionVsTruePt_    = inputDir.make<TProfile>("hisTDRVertexPTResolutionVsTruePt","TDR Reconstructed primary vertex relative p_{T} resolution vs. True Pt; True p_{T}; p_{T} Resolution [GeV]", 100, 0, 500);
  hisTDRVertexPTVsTruePt_              = inputDir.make<TH2F>("hisTDRVertexPtVsTruePt_","TDR Reconstructed primary vertex p_{T}  vs. True Pt; p_{T} [GeV]; True p_{T}", 100, 0, 500.,100, 0, 500.);
  hisTDRVertexMETVsTrueMET_            = inputDir.make<TH2F>("hisTDRVertexMETVsTrueMET_","TDR Reconstructed primary vertex MET  vs. True MET; MET [GeV]; True MET", 100, 0, 500.,100, 0, 500.);
  hisTDRNoTracksFromPrimaryVertex_     = inputDir.make<TH2F>("hisTDRNoTracksFromPrimaryVertex_","No. of Tracks from Primary Vertex (TDR vs. Truth); no. Reco Tracks in PV; no. Truth Tracks in PV ",50,0,50,50,0,50);
  hisTDRNoTrueTracksFromPrimaryVertex_ = inputDir.make<TH2F>("hisTDRNoTrueTracksFromPrimaryVertex_","No. of Matched Tracks from Primary Vertex (TDR vs. Truth); no. Reco Tracks in PV; no. Truth Tracks in PV ",50,0,50,50,0,50);
  hisTDRPrimaryVertexZ0width_          = inputDir.make<TH1F>("hisTDRPrimaryVertexZ0width_", "TDRReconstructed primary vertex z_{0} width", 100, 0, 0.5);
  hisRatioMatchedTracksInTDRPV_         = inputDir.make<TH1F>("hisRatioMatchedTracksInTDRPV", "TDR Primary vertex matching ratio ", 100, 0, 1.1);
  hisFakeTracksRateInTDRPV_             = inputDir.make<TH1F>("hisFakeTracksRateInTDRPV", "Percentage of fake tracks in TDR recontructed primary vertex", 100, 0, 1.1);
  hisTrueTracksRateInTDRPV_             = inputDir.make<TH1F>("hisTrueTracksRateInTDRPV", "Percentage of true tracks in TDR recontructed primary vertex", 100, 0, 1.1);
  
  hisTDRUnmatchZ0distance_              = inputDir.make<TH1F>("hisTDRUnmatchZ0distance_", "z0 distance from TDR reconstructed privary vertex of unmatched tracks; |z_{0}^{track} - z_{0}^{vertex}|; no. L1 Tracks", 100, 0, 5.);
  hisTDRUnmatchZ0MinDistance_              = inputDir.make<TH1F>("hisTDRUnmatchZ0MinDistance_", "z0 distance from the closest track in TDR reconstructed privary vertex of unmatched tracks; |z_{0}^{track} - z_{0}^{PV track}|; no. L1 Tracks", 100, 0, 5.);
  hisTDRUnmatchPt_                      = inputDir.make<TH1F>("hisTDRUnmatchPt_", "Transverse momentum of unmatched TDR PV tracks; p_{T} [GeV/c]; no. L1 Tracks", 100, 0, 100.);
  hisTDRUnmatchEta_                      = inputDir.make<TH1F>("hisTDRUnmatchEta_", "#eta of unmatched TDR PV tracks; #eta; no. L1 Tracks", 100, 2.4, 2.4);
  hisTDRUnmatchTruePt_                      = inputDir.make<TH1F>("hisTDRUnmatchTruePt_", "True transverse momentum of unmatched TDR PV tracks; p_{T} [GeV/c]; no. L1 Tracks", 100, 0, 100.);
  hisTDRUnmatchTrueEta_                      = inputDir.make<TH1F>("hisTDRUnmatchTrueEta_", "True #eta of unmatched TDR PV tracks; #eta; no. L1 Tracks", 100, 2.4, 2.4);
  hisTDRUnmatchedPVtracks_              = inputDir.make<TH1F>("hisTDRUnmatchedPVtracks_", " No. of tracks from primary collision that are misassigned (TDR); No. misassigned Tracks, No. Events ", 100, 0, 100);
  hisTDRPrimaryVertexVsTrueZ0_      = inputDir.make<TH1F>("hisTDRPrimaryVertexVsTrueZ0_","No. of reconstructed primary vertices per true z_{0}; true z_{0} [cm]; No. Vertices", 50, -25.,25.);
  hisTDRPrimaryVertexResolutionVsTrueZ0_      = inputDir.make<TProfile>("hisTDRPrimaryVertexResolutionVsTrueZ0_","No. of reconstructed primary vertices per true z_{0} (Technical Proposal Algo); true z_{0} [cm]; Resolution [mm]", 100, -25.,25.);
  hisTDRVertexMET_ = inputDir.make<TH1F>("hisTDRVertexMET_","; L1TrkMET (GeV) ; Entries", 50, 0, 200);
  hisTDRVertexPT_ = inputDir.make<TH1F>("hisTDRVertexPT_","; #Sigma_{vtx} p_{T} (GeV) ; Entries", 50, 0, 200);
  hisTDRPileUpVertexPT_ = inputDir.make<TH1F>("hisTDRPileUpVertexPT_","; #Sigma_{vtx} p_{T} (GeV) ; Entries", 50, 0, 200);
  hisTDRVertexOffPT_ = inputDir.make<TH1F>("hisTDRVertexOffPT_","; #Sigma_{vtx} p_{T} (GeV) ; Entries", 50, 0, 200);


  // ** PileUp vertices plot **

  hisRecoPileUpVertexZ0width_        = inputDir.make<TH1F>("hisRecoPileUPVertexZ0width_", "Reconstructed pile-up vertex z_{0} width", 100, 0, 1.);
  hisRecoPileUpVertexZ0resolution_   = inputDir.make<TH1F>("hisRecoPileUPVertexZ0resolution_", "Reconstructed pile-up vertex z_{0} resolution; #sigma_{z_{0}} [cm]", 100, 0, 1.);
  
  hisRecoVertexZ0Spacing_            = inputDir.make<TH1F>("hisRecoVertexZ0Spacing", "Reconstructed intravertex z_{0} distance", 100, 0., 5.);
  
  hisPrimaryVertexTrueZ0_            = inputDir.make<TH1F>("hisPrimaryVertexTrueZ0_","No. of gen primary vertices per true z_{0}; true z_{0} [cm]; No. Vertices", 50, -25.,25.);

  hisPrimaryVertexZ0width_           = inputDir.make<TH1F>("hisPrimaryVertexZ0width_", "Primary vertex z_{0} width", 100, 0, 0.5);
  hisPileUpVertexZ0_                 = inputDir.make<TH1F>("hisPileUpVertexZ0", "Pile Up vertex z_{0} position", 200, -15., 15.);
  hisPileUpVertexZ0Spacing_          = inputDir.make<TH1F>("hisPileUpVertexZ0Spacing", "Pile Up intravertex z_{0} distance", 100, 0., 5.);
  hisPileUpVertexZ0width_            = inputDir.make<TH1F>("hisPileUpVertexZ0width_", "Pile Up vertex z_{0} width", 100, 0, 0.5);
  
  hisLostPVtracks_                   = inputDir.make<TH1F>("hisLostPVtracks_", " No. of tracks from primary collision that are not found by the L1 Track Finder; No. Lost Tracks, No. Events ", 100, 0, 100);

  hisTrkMETvsGenMET_                 = inputDir.make<TProfile>("hisTrkMETvsGenMET_"," Generated L1 Track MET vs Generated MET ", 20, 0, 500);
  hisRecoTrkMETvsGenMET_             = inputDir.make<TProfile>("hisRecoTrkMETvsGenMET_"," Reconstructed L1 Track MET vs Generated MET ", 20, 0, 500);
  hisTDRTrkMETvsGenMET_              = inputDir.make<TProfile>("hisTDRTrkMETvsGenMET_"," Technical Proposal L1 Track MET vs Generated MET ", 20, 0, 500);

  hisNumVxIterations_                = inputDir.make<TH1F>("hisNumVxIterations_","Number of Iterations (Vertex Reconstruction); No. Iterations; Entries",100,0,500);
  hisNumVxIterationsPerTrack_        = inputDir.make<TH1F>("hisNumVxIterationsPerTrack_","Number of Iterations per Track(Vertex Reconstruction); No. Iterations; Entries",100,0,200);


  grMET_.resize(4);
  grMET_tdr_.resize(4);
  hisMETevents_.resize(4);

  hisRecoVerticesVsTrueZ0PerDistance_.resize(4);
  hisTDRVerticesVsTrueZ0PerDistance_.resize(4);


  for(unsigned int i = 0; i < 4; ++i ){
    float genmet = 25. + i*25.;
    grMET_[i] = inputDir.make<TGraphErrors>(10);
    grMET_tdr_[i] = inputDir.make<TGraphErrors>(10);

    ostringstream title, name;
    title << "Signal Efficiency vs. Bkg. rejection ( GenMET > "<< genmet << " GeV); Signal Efficiency; Bkg Rejection power";
    name << "grSignVsBkgEffGenMET"<< genmet;

    grMET_[i]->SetTitle(title.str().c_str());
    grMET_[i]->SetName(name.str().c_str());
    grMET_[i]->SetMarkerStyle(21);
    grMET_[i]->GetXaxis()->SetRangeUser(0,1.1);
    grMET_[i]->GetYaxis()->SetRangeUser(0,1.1);
    // grMET_[i]->SetMarkerSize()

    title.clear();
    title.str("");
    title << "Signal Efficiency vs. Bkg. rejection ( GenMET > "<< genmet << " GeV, Technical Proposal Algo); Signal Efficiency; Bkg Rejection power";
    name.clear();
    name.str("");
    name << "grSignVsBkgEffGenMET"<< genmet<<"_TP";

    grMET_tdr_[i]->SetTitle(title.str().c_str());
    grMET_tdr_[i]->SetName(name.str().c_str());
    grMET_tdr_[i]->SetMarkerStyle(21);
    grMET_tdr_[i]->GetXaxis()->SetRangeUser(0,1.1);
    grMET_tdr_[i]->GetYaxis()->SetRangeUser(0,1.1);

    name.clear();
    name.str("");
    name << "hisMET"<< genmet;    
    hisMETevents_[i] = inputDir.make<TH1F>(name.str().c_str(), name.str().c_str(), 42, 0 ,42);
  }
  hisRecoVerticesVsTrueZ0PerDistance_[0] = inputDir.make<TH1F>("hisRecoVerticesVsTrueZ0PerDistance10_","No. of reconstructed primary vertices per true z_{0} (|z_{0,reco} - z_{0,gen}| < 1. mm) ; true z_{0} [cm]; No. Vertices", 50, -25.,25.);
  hisTDRVerticesVsTrueZ0PerDistance_[0] = inputDir.make<TH1F>("hisTDRVerticesVsTrueZ0PerDistance10_","No. of Technical proosal primary vertices per true z_{0} (|z_{0,tdr} - z_{0,gen}| < 1. mm); true z_{0} [cm]; No. Vertices", 50, -25.,25.);

  hisRecoVerticesVsTrueZ0PerDistance_[1] = inputDir.make<TH1F>("hisRecoVerticesVsTrueZ0PerDistance15_","No. of reconstructed primary vertices per true z_{0} (|z_{0,reco} - z_{0,gen}| < 1.5 mm) ; true z_{0} [cm]; No. Vertices", 50, -25.,25.);
  hisTDRVerticesVsTrueZ0PerDistance_[1] = inputDir.make<TH1F>("hisTDRVerticesVsTrueZ0PerDistance15_","No. of Technical proosal primary vertices per true z_{0} (|z_{0,tdr} - z_{0,gen}| < 1.5 mm); true z_{0} [cm]; No. Vertices", 50, -25.,25.);

  hisRecoVerticesVsTrueZ0PerDistance_[2] = inputDir.make<TH1F>("hisRecoVerticesVsTrueZ0PerDistance25_","No. of reconstructed primary vertices per true z_{0} (|z_{0,reco} - z_{0,gen}| < 2.5 mm) ; true z_{0} [cm]; No. Vertices", 50, -25.,25.);
  hisTDRVerticesVsTrueZ0PerDistance_[2] = inputDir.make<TH1F>("hisTDRVerticesVsTrueZ0PerDistance25_","No. of Technical proosal primary vertices per true z_{0} (|z_{0,tdr} - z_{0,gen}| < 2.5 mm); true z_{0} [cm]; No. Vertices", 50, -25.,25.);

  hisRecoVerticesVsTrueZ0PerDistance_[3] = inputDir.make<TH1F>("hisRecoVerticesVsTrueZ0PerDistance50_","No. of reconstructed primary vertices per true z_{0} (|z_{0,reco} - z_{0,gen}| < 5. mm) ; true z_{0} [cm]; No. Vertices", 50, -25.,25.);
  hisTDRVerticesVsTrueZ0PerDistance_[3] = inputDir.make<TH1F>("hisTDRVerticesVsTrueZ0PerDistance50_","No. of Technical proosal primary vertices per true z_{0} (|z_{0,tdr} - z_{0,gen}| < 5. mm); true z_{0} [cm]; No. Vertices", 50, -25.,25.);


  noSignalEvents.assign(4,0);
  noBackgroundEvents.assign(4,0);
  noSignalEventsTDR.assign(4,0);
  noBackgroundEventsTDR.assign(4,0);

  std::vector<unsigned int> emptyVector;
  emptyVector.assign(10,0);

  noRecoSignalEvents.assign(4,emptyVector);
  noRecoBackgroundEvents.assign(4,emptyVector);
  noTDRSignalEvents.assign(4,emptyVector);
  noTDRBackgroundEvents.assign(4,emptyVector);

}

void Histos::fillVertexReconstruction(const InputData& inputData, const VertexFinder& vf){


  // noEvents++;
  const Vertex     TruePrimaryVertex = inputData.getPrimaryVertex();
  const RecoVertex RecoPrimaryVertex = vf.PrimaryVertex();
  const RecoVertex TDRVertex         = vf.TDRPrimaryVertex();

  if(settings_->debug()>5){
    cout << "True PrimaryVertex z0 "<< TruePrimaryVertex.z0() << " pT "<< TruePrimaryVertex.pT() << " met "<< TruePrimaryVertex.met() << endl;
    cout << "Reco PrimaryVertex z0 "<< RecoPrimaryVertex.z0() << " pT "<< RecoPrimaryVertex.pT() << " met "<< RecoPrimaryVertex.met() << endl;
    cout << "TP PrimaryVertex z0 "<< TDRVertex.z0() << " pT "<< TDRVertex.pT() << " met "<< RecoPrimaryVertex.met() << endl;
  }

  hisNumVxIterations_->Fill(vf.NumIterations());
  hisNumVxIterationsPerTrack_->Fill(vf.IterationsPerTrack());
  
  hisTrkMETvsGenMET_->Fill(inputData.GenMET(), TruePrimaryVertex.met() );
  hisRecoTrkMETvsGenMET_->Fill(inputData.GenMET(), RecoPrimaryVertex.met());
  hisTDRTrkMETvsGenMET_->Fill(inputData.GenMET(), TDRVertex.met());

  hisNoRecoVertices_->Fill(vf.numVertices());
  hisNoPileUpVertices_->Fill(inputData.getRecoPileUpVertices().size());
  hisNoRecoVsNoTruePileUpVertices_->Fill(vf.numVertices(),inputData.getRecoPileUpVertices().size() );
  hisPrimaryVertexTrueZ0_->Fill(TruePrimaryVertex.z0());
  
  hisNoRealPVTracks_->Fill(TruePrimaryVertex.numTracks()); 
 
  float z0res = TruePrimaryVertex.z0() - RecoPrimaryVertex.z0();  
  
  if (TruePrimaryVertex.percentageTracksGoodEta()>=0.5){
      hisNoRealTracksPerc50_->Fill(TruePrimaryVertex.numTracks());
      hisz0RealPerc50_->Fill(TruePrimaryVertex.z0());
      hisRealpTPerc50_->Fill(TruePrimaryVertex.pT());
  }

  if (TruePrimaryVertex.percentageTracksGoodEta()>=0.9){
      hisNoRealTracksPerc90_->Fill(TruePrimaryVertex.numTracks());
      hisz0RealPerc90_->Fill(TruePrimaryVertex.z0());
      hisRealpTPerc90_->Fill(TruePrimaryVertex.pT());
  }
  
  hisNoRealTracksGoodEta_->Fill(TruePrimaryVertex.noTracksBelow15eta()); 
  hisNoRealTracksBadEta_->Fill(TruePrimaryVertex.noTracksAbove15eta()); 
  hisRealPVPt_->Fill(TruePrimaryVertex.pT());
  hisRealTracksvspT_->Fill(TruePrimaryVertex.pT(),TruePrimaryVertex.numTracks());
  gROOT->ForceStyle();
  gStyle->SetPalette(kBird);
  hisRealTracksvspT_->Draw("COLZ");
    
    
    
  //Layers Hist
    
    vector<unsigned int> noLay = TruePrimaryVertex.numLayersTracks();
    for(unsigned int k=0;k<noLay.size();k++)
        numLayersRealTracks_->Fill(noLay[k]);
    
    vector<unsigned int> noPSLay = TruePrimaryVertex.numPSLayersTracks();
    for(unsigned int k=0;k<noPSLay.size();k++)
        numPSLayersRealTracks_->Fill(noPSLay[k]);
    
    vector<int> layID = TruePrimaryVertex.layerIDTracks();
    for(unsigned int k=0;k<layID.size();k++)
        realTracksLayersId_->Fill(layID[k]);
  
    
    
    
    

  //  float z0res = TruePrimaryVertex.z0() - RecoPrimaryVertex.z0();
  float pTres = fabs(TruePrimaryVertex.pT() - RecoPrimaryVertex.pT());
  hisRecoVertexZ0Resolution_->Fill(fabs(z0res));

  // Vertex has been found
  if(fabs(z0res) < 0.1 ){
    hisRecoVerticesVsTrueZ0PerDistance_[0]->Fill(TruePrimaryVertex.z0());
  } 
  if(fabs(z0res) < 0.15){
    hisRecoVerticesVsTrueZ0PerDistance_[1]->Fill(TruePrimaryVertex.z0());
  } 
  if(fabs(z0res) < 0.25){
    hisRecoVerticesVsTrueZ0PerDistance_[2]->Fill(TruePrimaryVertex.z0());
  } 
  if(fabs(z0res) < 0.5){
    hisRecoVerticesVsTrueZ0PerDistance_[3]->Fill(TruePrimaryVertex.z0());
  }

  if(fabs(z0res) < settings_->vx_recodistance()) {

    hisNoRecoPVTracks_->Fill(TruePrimaryVertex.numTracks()); 
      
    
    vector<float> chi2 = RecoPrimaryVertex.tracksChi2dof();
    for(unsigned int k=0;k<chi2.size();k++)
      hisChi2dofTracks_->Fill(chi2[k]);
      
    for(unsigned int k=0;k<noLay.size();k++)
      numLayersRecoTracks_->Fill(noLay[k]);
      
    for(unsigned int k=0;k<noPSLay.size();k++)
      numPSLayersRecoTracks_->Fill(noPSLay[k]);
      
    for(unsigned int k=0;k<layID.size();k++)
      recoTracksLayersId_->Fill(layID[k]);
      

    if (TruePrimaryVertex.percentageTracksGoodEta()>=0.5){
        hisNoRecoTracksPerc50_->Fill(TruePrimaryVertex.numTracks());
        hisz0RecoPerc50_->Fill(TruePrimaryVertex.z0());
        hisRecopTPerc50_->Fill(TruePrimaryVertex.pT());
    }
    if (TruePrimaryVertex.percentageTracksGoodEta()>=0.9){
        hisNoRecoTracksPerc90_->Fill(TruePrimaryVertex.numTracks());
        hisz0RecoPerc90_->Fill(TruePrimaryVertex.z0());
        hisRecopTPerc90_->Fill(TruePrimaryVertex.pT());
    }
      
    hisNoRecoTracksGoodEta_->Fill(TruePrimaryVertex.noTracksBelow15eta()); 
    hisNoRecoTracksBadEta_->Fill(TruePrimaryVertex.noTracksAbove15eta()); 
    hisRecoPVPt_->Fill(TruePrimaryVertex.pT());
    hisRecoTracksvspT_->Fill(TruePrimaryVertex.pT(),TruePrimaryVertex.numTracks());
    gROOT->ForceStyle();
    gStyle->SetPalette(kBird);
    hisRecoTracksvspT_->Draw("COLZ");
    hisRecoPrimaryVertexVsTrueZ0_->Fill(TruePrimaryVertex.z0());
    hisRecoVertexPTResolution_->Fill(pTres/TruePrimaryVertex.pT());
    
    hisRecoVertexPTResolutionVsTruePt_->Fill(TruePrimaryVertex.pT(), pTres/TruePrimaryVertex.pT() );
  
    hisRecoVertexPTVsTruePt_->Fill(RecoPrimaryVertex.pT(), TruePrimaryVertex.pT());
    hisRecoVertexMETVsTrueMET_->Fill(RecoPrimaryVertex.met(), TruePrimaryVertex.met());
    hisRecoVertexMET_->Fill(RecoPrimaryVertex.met());
    hisNoTracksFromPrimaryVertex_->Fill(RecoPrimaryVertex.numTracks(),TruePrimaryVertex.numTracks());
    hisNoTrueTracksFromPrimaryVertex_->Fill(RecoPrimaryVertex.numTrueTracks(),TruePrimaryVertex.numTracks());
    hisRecoPrimaryVertexZ0width_->Fill(RecoPrimaryVertex.z0width());
    hisRecoVertexPT_->Fill(RecoPrimaryVertex.pT());
    
    float matchratio = float(RecoPrimaryVertex.numTrueTracks())/float(TruePrimaryVertex.numTracks());
    if(matchratio > 1.) matchratio = 1.;
    hisRatioMatchedTracksInPV_->Fill(matchratio);
    float trueRate = float(RecoPrimaryVertex.numTrueTracks())/float(RecoPrimaryVertex.numTracks());
    hisTrueTracksRateInPV_->Fill(trueRate);
    float fakeRate = float(RecoPrimaryVertex.numTracks()-RecoPrimaryVertex.numTrueTracks())/float(RecoPrimaryVertex.numTracks());
    hisFakeTracksRateInPV_->Fill(fakeRate);
    hisRecoPrimaryVertexResolutionVsTrueZ0_->Fill(TruePrimaryVertex.z0(),fabs(z0res));

    for(unsigned int i = 0; i < 4; ++i ){
      float genmet = 25. + i*25.;
      float met_steps = (genmet-10.)/10;
      bool signal = false;
      if(TruePrimaryVertex.met() > genmet){
        noSignalEvents[i]++;
        signal = true;
      }
      else {
        noBackgroundEvents[i]++;
      }

      for(unsigned int j = 0; j < 10; ++j ){
        float cutmet = 10. + j*met_steps;

        if(RecoPrimaryVertex.met() > cutmet){
          if(signal) noRecoSignalEvents[i][j]++;
        } else if(!signal){
          noRecoBackgroundEvents[i][j]++;
        }
      }
    }

  }
  else{
    hisRecoVertexOffPT_->Fill(RecoPrimaryVertex.pT());
    hisUnmatchedVertexZ0distance_->Fill(fabs(z0res));
    if(settings_->debug()> 5) {
      cout << "Vertex Reconstruction Algorithm doesn't find the correct the primary vertex (Delta Z = " << fabs(z0res) << ")"<<endl;
      for(const L1fittedTrack* l1track : RecoPrimaryVertex.tracks()){
        if(l1track->getMatchedTP() == nullptr){
        cout << "FAKE track assigned to PV. Track z0: "<< l1track->z0() << " track pT "<< l1track->pt() << " chi2/ndof " << l1track->chi2dof() << endl;
        } else if(l1track->getMatchedTP()->physicsCollision() != 0){
        cout << "Pile-Up track assigned to PV. Track z0: "<< l1track->z0() << " track pT "<< l1track->pt() << endl;
        } else{
          cout << "Physics Collision track assigned to PV. Track z0: "<< l1track->z0() << " track pT "<< l1track->pt() << endl;
        }
      }
    }
  }

  // ** Technical Proposal algorithm
  float z0res_tdr = (TruePrimaryVertex.z0() - TDRVertex.z0());
  float pTres_tdr = fabs(TruePrimaryVertex.pT() - TDRVertex.pT());
  
  hisTDRVertexZ0Resolution_->Fill(fabs(z0res_tdr));

  if(fabs(z0res_tdr) < 0.1 ){
    hisTDRVerticesVsTrueZ0PerDistance_[0]->Fill(TruePrimaryVertex.z0());
  } 
  if(fabs(z0res_tdr) < 0.15){
    hisTDRVerticesVsTrueZ0PerDistance_[1]->Fill(TruePrimaryVertex.z0());
  } 
  if(fabs(z0res_tdr) < 0.25){
    hisTDRVerticesVsTrueZ0PerDistance_[2]->Fill(TruePrimaryVertex.z0());
  } 
  if(fabs(z0res_tdr) < 0.5){
    hisTDRVerticesVsTrueZ0PerDistance_[3]->Fill(TruePrimaryVertex.z0());
  }

  if(fabs(z0res_tdr) < settings_->vx_recodistance()){
    hisTDRPrimaryVertexVsTrueZ0_->Fill(TruePrimaryVertex.z0());
    hisTDRPrimaryVertexResolutionVsTrueZ0_->Fill(TruePrimaryVertex.z0(),fabs(z0res_tdr));
    hisTDRVertexPT_->Fill(TDRVertex.pT());
    hisTDRVertexMET_->Fill(TDRVertex.met());
    hisTDRVertexPTResolution_->Fill(pTres_tdr/TruePrimaryVertex.pT());
    hisTDRVertexPTResolutionVsTruePt_->Fill(TruePrimaryVertex.pT(), pTres_tdr/TruePrimaryVertex.pT() );
    hisTDRVertexPTVsTruePt_->Fill(TDRVertex.pT(), TruePrimaryVertex.pT());
    hisTDRVertexMETVsTrueMET_->Fill(TDRVertex.met(), TruePrimaryVertex.met());
    hisTDRNoTracksFromPrimaryVertex_->Fill(TDRVertex.numTracks(),TruePrimaryVertex.numTracks());
    hisTDRNoTrueTracksFromPrimaryVertex_->Fill(TDRVertex.numTrueTracks(),TruePrimaryVertex.numTracks());
    hisTDRPrimaryVertexZ0width_->Fill(TDRVertex.z0width());
  
    float matchratio_res = float(TDRVertex.numTrueTracks())/float(TruePrimaryVertex.numTracks());
    hisRatioMatchedTracksInTDRPV_->Fill(matchratio_res);
    float trueRate_res = float(TDRVertex.numTrueTracks())/float(TDRVertex.numTracks());
    hisTrueTracksRateInTDRPV_->Fill(trueRate_res);
    float fakeRate_res = float(TDRVertex.numTracks()-TDRVertex.numTrueTracks())/float(TDRVertex.numTracks());
    hisFakeTracksRateInTDRPV_->Fill(fakeRate_res);

    for(unsigned int i = 0; i < 4; ++i ){
      float genmet = 25. + i*25.;
      float met_steps = (genmet-10.)/10;
      bool signal = false;
      if(TruePrimaryVertex.met() > genmet){
        noSignalEventsTDR[i]++;
        signal = true;
      }
      else {
        noBackgroundEventsTDR[i]++;
      }

      for(unsigned int j = 0; j < 10; ++j ){
        float cutmet = 10. + j*met_steps;

        if(TDRVertex.met() > cutmet){
          if(signal) noTDRSignalEvents[i][j]++;
        } else if(!signal){
          noTDRBackgroundEvents[i][j]++;
        }
      }
    }

  } else{
    hisTDRVertexOffPT_->Fill(TDRVertex.pT());
    hisTDRUnmatchedVertexZ0distance_->Fill(fabs(z0res_tdr));
    if(settings_->debug()> 5) {
      cout << "TP Reconstruction Algorithm doesn't find the correct the primary vertex (Delta Z = " << fabs(z0res_tdr) << ")"<<endl;
      for(const L1fittedTrack* l1track : TDRVertex.tracks()){
        if(l1track->getMatchedTP() == nullptr){
        cout << "FAKE track assigned to PV. Track z0: "<< l1track->z0() << " track pT "<< l1track->pt() << " chi2/ndof " << l1track->chi2dof() << endl;
        } else if(l1track->getMatchedTP()->physicsCollision() != 0){
        cout << "Pile-Up track assigned to PV. Track z0: "<< l1track->z0() << " track pT "<< l1track->pt() << endl;
        } else{
          cout << "Physics Collision track assigned to PV. Track z0: "<< l1track->z0() << " track pT "<< l1track->pt() << endl;
        }
      }
    }

  }

  unsigned int lostTracks = 0;
  unsigned int misassignedTracks = 0;
  unsigned int misassignedTracks_tdr = 0;


  if(settings_->debug() > 5) cout << "*** Misassigned primary vertex tracks ***"<< endl;
  for(TP tp : TruePrimaryVertex.tracks()){
    bool found = false;
    // cout << tp.index() << " "<< endl;
    for(const L1fittedTrack* l1track : RecoPrimaryVertex.tracks()){
      if(l1track->getMatchedTP()!= nullptr){
        if(tp.index() == l1track->getMatchedTP()->index() ) {
          found = true;
          break;
        }
      } 
    }
    if(!found){
      bool TrackIsReconstructed = false;
      for(const L1fittedTrack* l1track: vf.FitTracks()){
        if(l1track->getMatchedTP()!= nullptr){
          if(tp.index() == l1track->getMatchedTP()->index() ){
            TrackIsReconstructed = true;
            hisUnmatchZ0distance_->Fill(fabs(l1track->z0()-RecoPrimaryVertex.z0()));
            hisUnmatchPt_->Fill(l1track->pt());
            hisUnmatchEta_->Fill(l1track->eta());
            hisUnmatchTruePt_->Fill(tp.pt());
            hisUnmatchTrueEta_->Fill(tp.eta());

            double mindistance = 999.;
            for(const L1fittedTrack* vertexTrack : RecoPrimaryVertex.tracks()){
              if( fabs(vertexTrack->z0()-l1track->z0()) < mindistance ) mindistance = fabs(vertexTrack->z0()-l1track->z0());
            }
            hisUnmatchZ0MinDistance_->Fill(mindistance);
            
            if(settings_->debug()>5){
              cout << "PV Track assigned to wrong vertex. Track z0: "<< l1track->z0() << " PV z0: "<< RecoPrimaryVertex.z0() << " tp z0 "<< tp.z0() << " track pT "<< l1track->pt() << " tp pT "<< tp.pt() << " tp d0 "<< tp.d0() << endl;
            }
            break;  
          }
        }
      }
      
      if(!TrackIsReconstructed){
        lostTracks++;
      } else{
        misassignedTracks++;
      }
    }
    
    
    
    found = false;

    for(const L1fittedTrack* l1track : TDRVertex.tracks()){
      if(l1track->getMatchedTP()!= nullptr){
        // cout << l1track->getMatchedTP()->index() << " ";
        if(tp.index() == l1track->getMatchedTP()->index() ) {
          found = true;
          break;
        }
      }
    }

    if(!found){
      for(const L1fittedTrack* l1track: vf.FitTracks()){
        if(l1track->getMatchedTP()!= nullptr){
          if(tp.index() == l1track->getMatchedTP()->index() ){
            hisTDRUnmatchZ0distance_->Fill(fabs(l1track->z0()-TDRVertex.z0()));
            hisTDRUnmatchPt_->Fill(l1track->pt());
            hisTDRUnmatchEta_->Fill(l1track->eta());
            hisTDRUnmatchTruePt_->Fill(tp.pt());
            hisTDRUnmatchTrueEta_->Fill(tp.eta());
            misassignedTracks_tdr++;
            double mindistance = 999.;
            for(const L1fittedTrack* vertexTrack : TDRVertex.tracks()){
              if( fabs(vertexTrack->z0()-l1track->z0()) < mindistance ) mindistance = fabs(vertexTrack->z0()-l1track->z0());
            }
            hisTDRUnmatchZ0MinDistance_->Fill(mindistance);
            if(settings_->debug()>5){
              cout << "TP PV Track assigned to wrong vertex. Track z0: "<< l1track->z0() << " PV z0: "<< TDRVertex.z0() << " tp z0 "<< tp.z0() << " track pT "<< l1track->pt() << " tp pT "<< tp.pt() << " tp d0 "<< tp.d0() << endl;
            }
            break;  
          }
        }
      }
    }


  }

  hisLostPVtracks_->Fill(lostTracks);
  hisUnmatchedPVtracks_->Fill(misassignedTracks);
  hisTDRUnmatchedPVtracks_->Fill(misassignedTracks_tdr);

  hisPrimaryVertexZ0width_->Fill(TruePrimaryVertex.z0width());

  float z0distance = 0.;

  for(unsigned int i = 0; i<vf.Vertices().size(); ++i){
    if(i < vf.Vertices().size()-1){
      z0distance = vf.Vertices()[i+1].z0() - vf.Vertices()[i].z0();
      hisRecoVertexZ0Spacing_->Fill(z0distance);
    }
    if(i != vf.PrimaryVertexId()  ){ 
      hisRecoPileUpVertexZ0width_->Fill(vf.Vertices()[i].z0width());
      hisRecoPileUpVertexPT_->Fill(vf.Vertices()[i].pT());
      double PUres = 999.;
      for(unsigned int j = 0; j<inputData.getRecoPileUpVertices().size(); ++j){
        if(fabs(vf.Vertices()[i].z0()-inputData.getRecoPileUpVertices()[j].z0()) < PUres){
          PUres = fabs(vf.Vertices()[i].z0()-inputData.getRecoPileUpVertices()[j].z0());
        }
      }
      hisRecoPileUpVertexZ0resolution_->Fill(PUres);
    }
  }

  for(unsigned int i = 0; i<inputData.getRecoPileUpVertices().size(); ++i){
    if(i < inputData.getRecoPileUpVertices().size()-1){
      z0distance = inputData.getRecoPileUpVertices()[i+1].z0() - inputData.getRecoPileUpVertices()[i].z0();
      hisPileUpVertexZ0Spacing_->Fill(z0distance);
    }
    hisPileUpVertexZ0_->Fill(inputData.getRecoPileUpVertices()[i].z0());
    hisPileUpVertexZ0width_->Fill(inputData.getRecoPileUpVertices()[i].z0width());
  }


  if(settings_->debug()>5) cout << "================ End of Event =============="<< endl;





}


//=== Understand why not all tracking particles were reconstructed.
//=== Returns list of tracking particles that were not reconstructed and an string indicating why.
//=== Only considers TP used for algorithmic efficiency measurement.

// (If string = "mystery", reason for loss unknown. This may be a result of reconstruction of one 
// track candidate preventing reconstruction of another. e.g. Due to duplicate track removal).

map<const TP*, string> Histos::diagnoseTracking(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs) const {

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
              if ( bend > 0 ) hisWrongSignStubRZ_pBend_->Fill( s->z(), s->r() );
              else if ( bend < 0 ) hisWrongSignStubRZ_nBend_->Fill( s->z(), s->r() );
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

//=== Book histograms studying freak, large events with too many stubs.

void Histos::bookStudyBusyEvents() {

  TFileDirectory inputDir = fs_->mkdir("BusyEvents");

  // Look at (eta, phi) sectors with too many input stubs or too many output (= assigned to tracks) stubs.

  unsigned int nEta = settings_->numEtaRegions();

  hisNumBusySecsInPerEvent_  = inputDir.make<TH1F>("NumBusySecsInPerEvent" ,"; No. sectors with too many input stubs/event" , 20, -0.5, 19.5);
  hisNumBusySecsOutPerEvent_ = inputDir.make<TH1F>("NumBusySecsOutPerEvent","; No. sectors with too many output stubs/event", 20, -0.5, 19.5);
  profFracBusyInVsEtaReg_   = inputDir.make<TProfile>("FracBusyInVsEtaReg" ,"; #eta region; Frac. of sectors with too many input stubs" , nEta, -0.5, nEta-0.5);
  profFracBusyOutVsEtaReg_  = inputDir.make<TProfile>("FracBusyOutVsEtaReg","; #eta region; Frac. of sectors with too many output stubs", nEta, -0.5, nEta-0.5);
  profFracStubsKilledVsEtaReg_ = inputDir.make<TProfile>("FracStubsKilledInVsEtaReg" ,"; #eta region; Frac. of input stubs killed" , nEta, -0.5, nEta-0.5);
  profFracTracksKilledVsEtaReg_ = inputDir.make<TProfile>("FracTracksKilledInVsEtaReg" ,"; #eta region; Frac. of track killed" , nEta, -0.5, nEta-0.5);
  profFracTracksKilledVsInvPt_ = inputDir.make<TProfile>("FracTracksKilledInVsInvPt" ,";1/Pt; Frac. of track killed" , 16, 0.,  1./settings_->houghMinPt());
  profFracTPKilledVsEta_   = inputDir.make<TProfile>("FracTPKilledInVsEta" ,";#eta; Efficiency loss due to busy sectors" , 16, 0.,  settings_->maxStubEta());
  profFracTPKilledVsInvPt_ = inputDir.make<TProfile>("FracTPKilledInVsInvPt" ,";1/Pt; Efficiency loss due to busy sectors" , 16, 0.,  1./settings_->houghMinPt());
  hisNumTPkilledBusySec_ = inputDir.make<TH1F>("NumTPkilledBusySec","; No. of TP killed in each busy sector",30,-0.5,29.5);

  // Compare properties of sectors with/without too many output stubs.

  const vector<string> tnames = {"BusyOutSec", "QuietOutSec"};
  const vector<string> enames = {" in busy output sector", " in quiet output sector"};
  for (unsigned int i = 0; i <= 1; i++) {
    const string tn = tnames[i];
    const string en = enames[i];

    hisNumInputStubs_[tn]     = inputDir.make<TH1F>(("NumInputStubs"+(tn)).c_str(),     ("; No. input stubs"+(en)).c_str(),   250, -0.5, 249.5); 
    hisQoverPtInputStubs_[tn] = inputDir.make<TH1F>(("QoverPtInputStubs"+(tn)).c_str(), ("; q/Pt of input stubs"+(en)).c_str(),   30, 0., 1./settings_->houghMinPt()); 
    hisNumOutputStubs_[tn]     = inputDir.make<TH1F>(("NumOutputStubs"+(tn)).c_str(),   ("; No. output stubs"+(en)).c_str(), 1000, -0.5, 999.5); 
    hisNumTracks_[tn]         = inputDir.make<TH1F>(("NumTracks"+(tn)).c_str(),         ("; No. tracks"+(en)).c_str(),        200, -0.5, 199.5);
    hisNumStubsPerTrack_[tn]  = inputDir.make<TH1F>(("NumStubsPerTrack"+(tn)).c_str(),  ("; No. stubs/track"+(en)).c_str(),    50, -0.5, 49.5);
    hisTrackQoverPt_[tn]      = inputDir.make<TH1F>(("TrackQoverPt"+(tn)).c_str(),      ("; Track q/pt"+(en)).c_str(),      30, 0., 1./settings_->houghMinPt());
    hisTrackPurity_[tn]       = inputDir.make<TH1F>(("TrackPurity"+(tn)).c_str(),       ("; Track purity"+(en)).c_str(),      102, -0.01, 1.01);
    hisNumTPphysics_[tn]      = inputDir.make<TH1F>(("NumTPphysics"+(tn)).c_str(),      ("; No. physics TP"+(en)).c_str(),     30, -0.5, 29.5);
    hisNumTPpileup_[tn]       = inputDir.make<TH1F>(("NumTPpileup"+(tn)).c_str(),       ("; No. pileup TP"+(en)).c_str(),      30, -0.5, 29.5);
    hisSumPtTPphysics_[tn]    = inputDir.make<TH1F>(("SumPtTPphysics"+(tn)).c_str(),    ("; Sum Pt physics TP"+(en)).c_str(), 100,  0.0, 100.);
    hisSumPtTPpileup_[tn]     = inputDir.make<TH1F>(("SumPtTPpileup"+(tn)).c_str(),     ("; Sum Pt pileup TP"+(en)).c_str(),  100,  0.0, 100.);
  }
}

//=== Fill histograms studying freak, large events with too many stubs.

void Histos::fillStudyBusyEvents(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs) {

  const unsigned int numStubsCut = settings_->busySectorNumStubs();   // No. of stubs per HT array the hardware can output.
  const bool         eachCharge  = settings_->busySectorEachCharge(); // +ve & -ve tracks output on separate optical links?

  const vector<const Stub*>&  vStubs = inputData.getStubs();
  const vector<TP>&           vTPs   = inputData.getTPs();

  // Create map containing L1 tracks found in whole of tracker together with flag indicating if the
  // track was killed because it was in a busy sector.
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
	profFracStubsKilledVsEtaReg_->Fill(iEtaReg, kill); 
      }
      bool tooBusyIn = (nStubsIn > numStubsCut);      
      if (tooBusyIn) nBusySecIn++;
      profFracBusyInVsEtaReg_->Fill(iEtaReg, tooBusyIn); // Sector had too many input stubs.
 
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

	profFracTracksKilledVsEtaReg_->Fill(iEtaReg, kill); 
	profFracTracksKilledVsInvPt_->Fill(ptInv, kill); 

	// Form a map of all tracks in the entire tracker & also just in this sector, with a flag indicating if they were killed as in a busy sector.
	trksInEntireTracker[trk] = kill;
	trksInSector[trk]        = kill;
      }

      if (tooBusyOut) nBusySecOut++;
      profFracBusyOutVsEtaReg_->Fill(iEtaReg, tooBusyOut); // Sector had too many output stubs.

      //--- Compare properties of sectors with/without too many output stubs.

      const vector<string> tnames = {"BusyOutSec", "QuietOutSec"};

      // Loop over sectors with too many/not too many output stubs.
      for (const string& tn : tnames) {
	if ((tn == "BusyOutSec" && tooBusyOut) || (tn == "QuietOutSec" && (! tooBusyOut))) {

	  hisNumInputStubs_[tn]->Fill(nStubsIn);

	  // Check if q/Pt estimated from stub bend differs in busy & quiet sectors.
          for (const Stub* stub : vStubs) {
    	    if ( sector.inside( stub ) ) hisQoverPtInputStubs_[tn]->Fill(abs(stub->qOverPt()));
          }

	  // Look at reconstructed tracks in this sector.
	  hisNumOutputStubs_[tn]->Fill(nStubsOut);
	  hisNumTracks_[tn]->Fill(tracks.size());
	  for (const L1track3D& trk : tracks) {
	    hisNumStubsPerTrack_[tn]->Fill(trk.getNumStubs());
	    hisTrackQoverPt_[tn]->Fill(trk.qOverPt());
	    hisTrackPurity_[tn]->Fill(trk.getPurity());
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
	  hisNumTPphysics_[tn]->Fill(num_TP_physics);
	  hisNumTPpileup_[tn]->Fill(num_TP_pileup);
	  hisSumPtTPphysics_[tn]->Fill(sumPt_TP_physics);
	  hisSumPtTPpileup_[tn]->Fill(sumPt_TP_pileup);
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
	hisNumTPkilledBusySec_->Fill(nTPkilled);
      }
    }
  }

  hisNumBusySecsInPerEvent_->Fill(nBusySecIn); // No. of sectors per event with too many input stubs.
  hisNumBusySecsOutPerEvent_->Fill(nBusySecOut); // No. of sectors per event with too many output stubs.

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
      profFracTPKilledVsEta_->Fill(fabs(tp.eta()), tpKilled);
      profFracTPKilledVsInvPt_->Fill(fabs(tp.qOverPt()), tpKilled);
    }
  }
}

//=== Book histograms for studying track fitting.

void Histos::bookTrackFitting() {
 
  // Define two inline functions to add number of helix parameters used in track fit to any character string.
  // One is used to modify the histogram name and the other to modify its title.
 
  // Book histograms for 4 and 5 parameter helix fits.
//  for ( std::vector<std::string>::const_iterator lHitIt = fitList_.begin(); lHitIt != fitList_.end(); ++lHitIt ) {
  for(auto &fitName : settings_->trackFitters() ){

    //std::cout << "Booking histograms for " << fitName << std::endl;
    TFileDirectory inputDir = fs_->mkdir( (fitName)  );
 
    hisSeedQinvPt_[fitName] = inputDir.make<TH1F>(("SeedQinvPt_"+(fitName)).c_str(), "; seed q/p_{T}" , 100, -0.5, 0.5 );
    hisSeedPhi0_[fitName]   = inputDir.make<TH1F>(("SeedPhi0_"+(fitName)).c_str(), "; seed #phi_{0}", 70, -3.5, 3.5 );
    hisSeedD0_[fitName]     = inputDir.make<TH1F>(("SeedD0_"+(fitName)).c_str(), "; seed d_{0}"   , 100, -1., 1. );
    hisSeedZ0_[fitName]     = inputDir.make<TH1F>(("SeedZ0_"+(fitName)).c_str(), "; seed z_{0}"   , 100, -25., 25. );
    hisSeedEta_[fitName]    = inputDir.make<TH1F>(("SeedEta_"+(fitName)).c_str(), "; seed #eta"    , 70, -3.5, 3.5 );
 
    profNumFittedCands_[fitName] = inputDir.make<TProfile>(("NumFittedCands_"+(fitName)).c_str(), "; class; # of fitted tracks", 11, 0.5, 11.5, -0.5, 9.9e6);
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(11, "Num tracks exc dups passing cut");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(10, "Num tracks exc dups");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(9, "Num rejected fake tracks");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(8, "Num rejected tracks");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(7, "Num TP tracks");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(6, "Num TP track killed by cut");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(5, "TP fitted passing cut");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(4, "TP fitted");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(3, "Accepted track");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(2, "Number of Stubs");
    profNumFittedCands_[fitName]->GetXaxis()->SetBinLabel(1, "Fitted tracks including fakes");
    profNumFittedCands_[fitName]->LabelsOption("d");
 
    hisFitQinvPtMatched_[fitName] = inputDir.make<TH1F>(("FitQinvPtMatched_"+(fitName)).c_str(),"Fitted q/p_{T} for matched tracks", 100, -0.5, 0.5 );
    hisFitPhi0Matched_[fitName]   = inputDir.make<TH1F>(("FitPhi0Matched_"+(fitName)).c_str(), "Fitted #phi_{0} for matched tracks", 70, -3.5, 3.5 );
    hisFitD0Matched_[fitName]     = inputDir.make<TH1F>(("FitD0Matched_"+(fitName)).c_str(), "Fitted d_{0} for matched tracks", 100, -1., 1. );
    hisFitZ0Matched_[fitName]     = inputDir.make<TH1F>(("FitZ0Matched_"+(fitName)).c_str(), "Fitted z_{0} for matched tracks", 100, -25., 25. );
    hisFitEtaMatched_[fitName]    = inputDir.make<TH1F>(("FitEtaMatched_"+(fitName)).c_str(), "Fitted #eta for matched tracks", 70, -3.5, 3.5 );
 
    hisFitChi2Matched_[fitName]    = inputDir.make<TH1F>(("Chi2Matched_"+(fitName)).c_str(), "#Chi^{2}", 100, 0., 100. );
    hisFitChi2DofMatched_[fitName] = inputDir.make<TH1F>(("Chi2DofMatched_"+(fitName)).c_str(), "#Chi^{2}/DOF", 100, 0., 10. );
 
    hisFitQinvPtUnmatched_[fitName] = inputDir.make<TH1F>(("FitQinvPtUnmatched_"+(fitName)).c_str(), "Fitted q/p_{T} for unmatched tracks", 100, -0.5, 0.5 );
    hisFitPhi0Unmatched_[fitName]   = inputDir.make<TH1F>(("FitPhi0Unmatched_"+(fitName)).c_str(), "Fitted #phi_{0} for unmatched tracks", 70, -3.5, 3.5 );
    hisFitD0Unmatched_[fitName]     = inputDir.make<TH1F>(("FitD0Unmatched_"+(fitName)).c_str(), "Fitted d_{0} for unmatched tracks", 100, -1., 1. );
    hisFitZ0Unmatched_[fitName]     = inputDir.make<TH1F>(("FitZ0Unmatched_"+(fitName)).c_str(), "Fitted z_{0} for unmatched tracks", 100, -25., 25. );
    hisFitEtaUnmatched_[fitName]    = inputDir.make<TH1F>(("FitEtaUnmatched_"+(fitName)).c_str(), "Fitted #eta for unmatched tracks", 70, -3.5, 3.5 );
 
    hisFitChi2Unmatched_[fitName]    = inputDir.make<TH1F>(("FitChi2Unmatched_"+(fitName)).c_str(), "#Chi^{2}", 100, 0., 100. );
    hisFitChi2DofUnmatched_[fitName] = inputDir.make<TH1F>(("FitChi2DofUnmatched_"+(fitName)).c_str(), "#Chi^{2}/DOF", 100, 0., 10. );
 
    hisFitVsTrueQinvPtGoodChi2_[fitName] = inputDir.make<TH2F>(("FitVsTrueQinvPtGoodChi2_"+(fitName)).c_str(), "; TP q/p_{T}; Fitted q/p_{T} (good #chi^2)", 100, -0.5, 0.5, 100, -0.5, 0.5 );
    hisFitVsTruePhi0GoodChi2_[fitName]   = inputDir.make<TH2F>(("FitVsTruePhi0GoodChi2_"+(fitName)).c_str(), "; TP #phi_{0}; Fitted #phi_{0} (good #chi^2)", 70, -3.5, 3.5, 70, -3.5, 3.5 );
    hisFitVsTrueD0GoodChi2_[fitName]     = inputDir.make<TH2F>(("FitVsTrueD0GoodChi2_"+(fitName)).c_str(), "; TP d_{0}; Fitted d_{0} (good #chi^2)", 100, -1., 1., 100, -1., 1. );
    hisFitVsTrueZ0GoodChi2_[fitName]     = inputDir.make<TH2F>(("FitVsTrueZ0GoodChi2_"+(fitName)).c_str(), "; TP z_{0}; Fitted z_{0} (good #chi^2)" , 100, -25., 25., 100, -25., 25. );
    hisFitVsTrueEtaGoodChi2_[fitName]    = inputDir.make<TH2F>(("FitVsTrueEtaGoodChi2_"+(fitName)).c_str(), "; TP #eta; Fitted #eta (good #chi^2)", 70, -3.5, 3.5, 70, -3.5, 3.5 );
 
    hisFitVsSeedQinvPtGenCand_[fitName] = inputDir.make<TH2F>(("FitVsSeedQinvPtGenCand_"+(fitName)).c_str(), "; Seed q/p_{T} (Genuine Cand); Fitted q/p_{T}", 100, -0.5, 0.5, 100, -0.5, 0.5 );
    hisFitVsSeedPhi0GenCand_[fitName]   = inputDir.make<TH2F>(("FitVsSeedPhi0GenCand_"+(fitName)).c_str(), "; Seed #phi_{0} (Genuine Cand); Fitted #phi_{0}", 70, -3.5, 3.5, 70, -3.5, 3.5 );
    hisFitVsSeedD0GenCand_[fitName]     = inputDir.make<TH2F>(("FitVsSeedD0GenCand_"+(fitName)).c_str(), "; Seed d_{0} (Genuine Cand); Fitted d_{0}", 100, -1., 1., 100, -1., 1. );
    hisFitVsSeedZ0GenCand_[fitName]     = inputDir.make<TH2F>(("FitVsSeedZ0GenCand_"+(fitName)).c_str(), "; Seed z_{0} (Genuine Cand); Fitted z_{0}", 100, -25., 25., 100, -25., 25. );
    hisFitVsSeedEtaGenCand_[fitName]    = inputDir.make<TH2F>(("FitVsSeedEtaGenCand_"+(fitName)).c_str(), "; Seed #eta (Genuine Cand); Fitted #eta", 70, -3.5, 3.5, 70, -3.5, 3.5 ); 
 
    hisFitQinvPtResGoodChi2_[fitName] = inputDir.make<TH1F>(("FitQinvPtResGoodChi2_"+(fitName)).c_str(), "Fitted minus true q/p_{T} (good #chi^2)", 100, -0.5, 0.5 );
    hisFitPhi0ResGoodChi2_[fitName]   = inputDir.make<TH1F>(("FitPhi0ResGoodChi2_"+(fitName)).c_str(), "Fitted minus true #phi_{0} (good #chi^2)", 100, -0.1, 0.1 );
    hisFitD0ResGoodChi2_[fitName]     = inputDir.make<TH1F>(("FitD0ResGoodChi2_"+(fitName)).c_str(), "Fitted minus true d_{0} (good #chi^2)", 100,  -1., 1. );
    hisFitZ0ResGoodChi2_[fitName]     = inputDir.make<TH1F>(("FitZ0ResGoodChi2_"+(fitName)).c_str(), "Fitted minus true z_{0} (good #chi^2)", 100, -10., 10. );
    hisFitEtaResGoodChi2_[fitName]    = inputDir.make<TH1F>(("FitEtaResGoodChi2_"+(fitName)).c_str(), "Fitted minus true #eta (good #chi^2)", 100, -0.1, 0.1 );

    hisSeedQinvPtResGoodChi2_[fitName] = inputDir.make<TH1F>(("SeedQinvPtResGoodChi2_"+(fitName)).c_str(), "True minus seed q/p_{T} (good #chi^2)", 100, -0.5, 0.5 );
    hisSeedPhi0ResGoodChi2_[fitName]   = inputDir.make<TH1F>(("SeedPhi0ResGoodChi2_"+(fitName)).c_str(), "True minus seed #phi_{0} (good #chi^2)", 100, -0.1, 0.1 );
    hisSeedD0ResGoodChi2_[fitName]     = inputDir.make<TH1F>(("SeedD0ResGoodChi2_"+(fitName)).c_str(), "True minus seed d_{0} (good #chi^2)", 100,  -1., 1. );
    hisSeedZ0ResGoodChi2_[fitName]     = inputDir.make<TH1F>(("SeedZ0ResGoodChi2_"+(fitName)).c_str(), "True minus seed z_{0} (good #chi^2)", 100, -10., 10. );
    hisSeedEtaResGoodChi2_[fitName]    = inputDir.make<TH1F>(("SeedEtaResGoodChi2_"+(fitName)).c_str(), "True minus seed #eta (good #chi^2)", 100, -0.1, 0.1 );
 
    hisFitVsSeedQinvPtFakeCand_[fitName] = inputDir.make<TH2F>(("FitVsSeedQinvPtFakeCand_"+(fitName)).c_str(), "; Seed q/p_{T} (Fake Cand); Fitted q/p_{T}", 100, -0.5, 0.5, 100, -0.5, 0.5 );
    hisFitVsSeedPhi0FakeCand_[fitName]   = inputDir.make<TH2F>(("FitVsSeedPhi0FakeCand_"+(fitName)).c_str(), "; Seed #phi_{0} (Fake Cand); Fitted #phi_{0}", 70, -3.5, 3.5, 70, -3.5, 3.5 );
    hisFitVsSeedD0FakeCand_[fitName]     = inputDir.make<TH2F>(("FitVsSeedD0FakeCand_"+(fitName)).c_str(), "; Seed d_{0} (Fake Cand); Fitted d_{0}", 100, -1., 1., 100, -1., 1. );
    hisFitVsSeedZ0FakeCand_[fitName]     = inputDir.make<TH2F>(("FitVsSeedZ0FakeCand_"+(fitName)).c_str(), "; Seed z_{0} (Fake Cand); Fitted z_{0}", 100, -25., 25., 100, -25., 25. );
    hisFitVsSeedEtaFakeCand_[fitName]    = inputDir.make<TH2F>(("FitVsSeedEtaFakeCand_"+(fitName)).c_str(), "; Seed #eta (Fake Cand); Fitted #eta", 70, -3.5, 3.5, 70, -3.5, 3.5 );

    float maxEta = settings_->maxStubEta();
    hisQoverPtResVsTrueEta_[fitName]      = inputDir.make<TProfile>(("QoverPtResVsTrueEta_"+(fitName)).c_str(), "q/p_{T} resolution; |#eta|; q/p_{T} resolution", 24, 0.0, maxEta);
    hisPhi0ResVsTrueEta_[fitName]    = inputDir.make<TProfile>(("PhiResVsTrueEta_"+(fitName)).c_str(), "#phi_{0} resolution; |#eta|; #phi_{0} resolution", 24, 0.0, maxEta);
    hisEtaResVsTrueEta_[fitName]     = inputDir.make<TProfile>(("EtaResVsTrueEta_"+(fitName)).c_str(), "#eta resolution; |#eta|; #eta resolution", 24, 0.0, maxEta);
    hisZ0ResVsTrueEta_[fitName]      = inputDir.make<TProfile>(("Z0ResVsTrueEta_"+(fitName)).c_str(), "z_{0} resolution; |#eta|; z_{0} resolution", 24, 0.0, maxEta);
    hisD0ResVsTrueEta_[fitName]      = inputDir.make<TProfile>(("D0ResVsTrueEta_"+(fitName)).c_str(), "d_{0} resolution; |#eta|; d_{0} resolution", 24, 0.0, maxEta);

    hisQoverPtResVsTrueZ0_[fitName]      = inputDir.make<TProfile>(("QoverPtResVsTrueZ0_"+(fitName)).c_str(), "q/p_{T} resolution; z_{0}; q/ p_{T} resolution", 50, -25.,25.);
    hisPhi0ResVsTrueZ0_[fitName]    = inputDir.make<TProfile>(("PhiResVsTrueZ0_"+(fitName)).c_str(), "#phi_{0} resolution; z_{0}; #phi_{0} resolution", 50, -25.,25.);
    hisEtaResVsTrueZ0_[fitName]     = inputDir.make<TProfile>(("EtaResVsTrueZ0_"+(fitName)).c_str(), "#eta resolution; z_{0}; #eta resolution", 50, -25.,25.);
    hisZ0ResVsTrueZ0_[fitName]      = inputDir.make<TProfile>(("Z0ResVsTrueZ0_"+(fitName)).c_str(), "z_{0} resolution; z_{0}; z_{0} resolution", 50, -25.,25.);
    hisD0ResVsTrueZ0_[fitName]      = inputDir.make<TProfile>(("D0ResVsTrueZ0_"+(fitName)).c_str(), "d_{0} resolution; z_{0}; d_{0} resolution", 50, -25.,25.);

    float maxAbsQoverPt = 1./settings_->houghMinPt(); // Max. |q/Pt| covered by  HT array.
    hisQoverPtResVsTrueInvPt_[fitName]      = inputDir.make<TProfile>(("QoverPtResVsTrueInvPt_"+(fitName)).c_str(), "q/p_{T} resolution; 1/p_{T}; q/p_{T} resolution", 25, 0.0, maxAbsQoverPt);
    hisPhi0ResVsTrueInvPt_[fitName]    = inputDir.make<TProfile>(("PhiResVsTrueInvPt_"+(fitName)).c_str(), "#phi_{0} resolution; 1/p_{T}; #phi_{0} resolution", 25, 0.0, maxAbsQoverPt);
    hisEtaResVsTrueInvPt_[fitName]     = inputDir.make<TProfile>(("EtaResVsTrueInvPt_"+(fitName)).c_str(), "#eta resolution; 1/p_{T}; #eta resolution", 25, 0.0, maxAbsQoverPt);
    hisZ0ResVsTrueInvPt_[fitName]      = inputDir.make<TProfile>(("Z0ResVsTrueInvPt_"+(fitName)).c_str(), "z_{0} resolution; 1/p_{T}; z_{0} resolution", 25, 0.0, maxAbsQoverPt);
    hisD0ResVsTrueInvPt_[fitName]      = inputDir.make<TProfile>(("D0ResVsTrueInvPt_"+(fitName)).c_str(), "d_{0} resolution; 1/p_{T}; d_{0} resolution", 25, 0.0, maxAbsQoverPt);

    hisTrueFittedChiSquaredVsTrueEta_[fitName]      = inputDir.make<TH2F>(("TrueFittedChiSqauredVsEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} values and TP's #eta; #Chi^{2}; #eta", 50000, -0.5, 50000.5, 250, -2.50, 2.50 );
    hisTrueFittedChiSquaredDofVsTrueEta_[fitName]   = inputDir.make<TH2F>(("TrueFittedChiSqauredDofVsEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} Degrees of Freedom and TP's #eta; #Chi^{2}; #eta", 5001, -0.5, 5000.5, 250, -2.50, 2.50 );
    hisTrueFittedChiSquaredVsFittedEta_[fitName]    = inputDir.make<TH2F>(("TrueFittedChiSqauredVsFittedEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} values and fitted #eta; #Chi^{2}; #eta", 50000, -0.5, 50000.5, 250, -2.50, 2.50 );
    hisTrueFittedChiSquaredDofVsFittedEta_[fitName] = inputDir.make<TH2F>(("TrueFittedChiSqauredDofVsFittedEta_"+(fitName)).c_str(), "Comparison of #Chi^{2} Degrees of Freedom and fitted  #eta; #Chi^{2}; #eta", 5001, -0.5, 5000.5, 250, -2.50, 2.50 );
    hisFittedChiSquaredFunctionOfStubs_[fitName]    = inputDir.make<TH2F>(("FittedChiSquaredFunctionOfStubs_"+(fitName)).c_str(), "Comparison of #Chi^{2} values as a function of stubs produced by TP; # of stubs; #Chi^{2}", 21, -0.5, 20.5, 500000, -0.5, 500000.5 );
    hisFittedChiSquaredDofFunctionOfStubs_[fitName] = inputDir.make<TH2F>(("FittedChiSquaredDofFunctionOfStubs_"+(fitName)).c_str(), "Comparison of #Chi^{2} Degrees of Freedom as a function of stubs produced by TP; # of stubs; #Chi^{2} Degrees of Freedom", 21, -0.5, 20.5, 5001, -0.5, 5000.5 );

    hisTrueEtaMatchedGoodChi2_[fitName]    = inputDir.make<TH1F>(("TrueEtaMatchedGoodChi2_"+(fitName)).c_str(), "True #eta for matched tracks (good #chi^2)", 70, -3.5, 3.5 );
    hisTrueEtaMatchedBadChi2_[fitName]     = inputDir.make<TH1F>(("TrueEtaMatchedBadChi2_"+(fitName)).c_str(), "True #eta for matched tracks (bad #chi^2)", 70, -3.5, 3.5 );
    hisStubPurityMatchedGoodChi2_[fitName] = inputDir.make<TH1F>(("FracMatchedStubsMatchedGoodChi2_"+(fitName)).c_str(), "Purity of stubs on matched tracks (good #chi^2)", 102, -0.01, 1.01 );
    hisStubPurityMatchedBadChi2_[fitName]  = inputDir.make<TH1F>(("FracMatchedStubsMatchedBadChi2_"+(fitName)).c_str(), "Purity of stubs on matched tracks (bad #chi^2)", 102, -0.01, 1.01 );

    profChi2DofVsInvPtPERF_[fitName] = inputDir.make<TProfile>(("Chi2DofVsInvPtPERF_"+(fitName)).c_str(), "; q/Pt; Mean sqrt(#Chi^{2}/DOF) of perfect tracks", 16, 0.,  1./settings_->houghMinPt(), 0., 25.);
    profBigChi2DofVsInvPtPERF_[fitName] = inputDir.make<TProfile>(("BigChi2DofVsInvPtPERF_"+(fitName)).c_str(), "; q/Pt; Frac perfect tracks with #Chi^{2}/DOF > 10", 16, 0.,  1./settings_->houghMinPt());
    hisD0TPBigChi2DofPERF_[fitName]   = inputDir.make<TH1F>(("D0TPBigChi2DofPERF"+(fitName)).c_str(),"; True d0 for high Pt tracks with bad #chi2/DOF", 100,0.0,0.5);
    hisD0TPSmallChi2DofPERF_[fitName] = inputDir.make<TH1F>(("D0TPSmallChi2DofPERF"+(fitName)).c_str(),"; True d0 for high Pt tracks with good #chi2/DOF",100,0.0,0.5);

    hisNumMatchedStubsKilledVsKilled_[fitName] = inputDir.make<TH2F>(("NumMatchedStubsKilledVsKilled_"+(fitName)).c_str(), "; All stubs killed by fit; Good stubs killed by fit", 10,-0.5,9.5,10,-0.5,9.5);
    profTrksKilledByFit_[fitName] = inputDir.make<TProfile>(("TrksKilledByFit_"+(fitName)).c_str(), "Track category; Fraction of tracks killed by fit", 2,0.5,2.5,-999.,999.);
    profTrksKilledByFit_[fitName]->GetXaxis()->SetBinLabel(1, "matched");
    profTrksKilledByFit_[fitName]->GetXaxis()->SetBinLabel(2, "unmatched");
    hisNumStubsVsPurity_[fitName] = inputDir.make<TH2F>(("NumStubsVsPurity_"+(fitName)).c_str(), "; Number of stubs; Purity", 30, 0.0, 30.0, 100, 0.0, 1.0);
    
    // Histograms specific to Linear Regression track fitter.
    if (fitName.find("LinearRegression") != string::npos) {
      hisNumIterations_[fitName] = inputDir.make<TH1F>(("NumIterations_"+(fitName)).c_str(), "; Number of Iterations", 30, -0.5, 29.5);
      hisFailingState_[fitName] = inputDir.make<TH1F>(("FailingState_"+(fitName)).c_str(), "; State which lost matching", 10, 0.0, 10);
      hisTotalStateCalls_[fitName] = inputDir.make<TProfile>(("TotalStateCalls_"+(fitName)).c_str(), "; state; Total number of State calls", 10, 0.0, 10);
      hisRelativeStateCalls_[fitName] = inputDir.make<TProfile>(("RelativeStateCalls_"+(fitName)).c_str(), "; state; Number of State calls relative to total State calls", 10, 0.0, 10);
    }
    hisNumStubsFitKills_[fitName] = inputDir.make<TH1F>(("NumStubsFitKills_"+(fitName)).c_str(), "; Stubs per track killed by fit", 30, -0.5, 29.5);
    hisNumStubsFitKillsVsPurity_[fitName] = inputDir.make<TH2F>(("NumStubsFitKillsVsPurity_"+(fitName)).c_str(), "; Stubs per track killed by fit; Purity", 30, 0.0, 30, 100, -0.01, 1.01 );
    hisNumStubsFitKillsVsPurityMatched_[fitName] = inputDir.make<TH2F>(("NumStubsFitKillsVsPurityMatched_"+(fitName)).c_str(), "; Stubs per track killed bu fit; Purity", 30, 0.0, 30.0, 100, -0.01, 1.01 );
    hisNumStubsFitKillsVsPurityUnmatched_[fitName] = inputDir.make<TH2F>(("NumStubsFitKillsVsPurityUnmatched_"+(fitName)).c_str(), "; Stubs per track killed by fit; Purity", 30, 0.0, 30.0, 100, -0.01, 1.01 );

    // Duplicate track histos.
    profDupFitTrksVsEta_[fitName]   = inputDir.make<TProfile>(("DupFitTrksVsEta_"+(fitName)).c_str() ,"; #eta; Fraction of duplicate tracks",12,0.,3.);
    profDupFitTrksVsInvPt_[fitName] = inputDir.make<TProfile>(("DupFitTrksVsInvPt_"+(fitName)).c_str() ,"; 1/Pt; Fraction of duplicate tracks",16,0.,maxAbsQoverPt);
    // Fake track histos.
    profFakeFitTrksVsEta_[fitName]  = inputDir.make<TProfile>(("FakeFitTrksVsEta_"+(fitName)).c_str() ,"; #eta; Fraction of fake tracks",12,0.,3.);
    profFakeFitTrksVsZ0_[fitName]  = inputDir.make<TProfile>(("FakeFitTrksVsZ0_"+(fitName)).c_str() ,"; z_{0}; Fraction of fake tracks",50,-25.,25.);
    profFitFracTrueStubsVsLayer_[fitName] = inputDir.make<TProfile>(("FitFracTrueStubsVsLayer_"+(fitName)).c_str() ,";Layer ID; fraction of true stubs",30,0.5,30.5);
    profFitFracTrueStubsVsEta_[fitName] = inputDir.make<TProfile>(("FitFracTrueStubsVsEta_"+(fitName)).c_str() ,";#eta; fraction of true stubs",24,0.,3.);

    // Histos for tracking efficiency vs. TP kinematics. (Binning must match similar histos in bookTrackCands()).
    hisFitTPinvptForEff_[fitName] = inputDir.make<TH1F>(("FitTPinvptForEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for effi. measurement);",16,0.,maxAbsQoverPt);
    hisFitTPetaForEff_[fitName]   = inputDir.make<TH1F>(("FitTPetaForEff_"+(fitName)).c_str(),"; #eta of TP (used for effi. measurement);",20,-3.,3.);
    hisFitTPphiForEff_[fitName]   = inputDir.make<TH1F>(("FitTPphiForEff_"+(fitName)).c_str(),"; #phi of TP (used for effi. measurement);",20,-M_PI,M_PI);

    // Histo for efficiency to reconstruct track perfectly (no incorrect hits). (Binning must match similar histos in bookTrackCands()).
    hisPerfFitTPinvptForEff_[fitName]   = inputDir.make<TH1F>(("PerfFitTPinvptForEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for perf. effi. measurement);",16,0.,maxAbsQoverPt);
    hisPerfFitTPetaForEff_[fitName]   = inputDir.make<TH1F>(("PerfFitTPetaForEff_"+(fitName)).c_str(),"; #eta of TP (used for perfect effi. measurement);",20,-3.,3.);

    // Histos for algorithmic tracking efficiency vs. TP kinematics. (Binning must match similar histos in bookTrackCands()).
    hisFitTPinvptForAlgEff_[fitName] = inputDir.make<TH1F>(("FitTPinvptForAlgEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for alg. effi. measurement);",16,0.,maxAbsQoverPt);
    hisFitTPetaForAlgEff_[fitName]   = inputDir.make<TH1F>(("FitTPetaForAlgEff_"+(fitName)).c_str(),"; #eta of TP (used for alg. effi. measurement);",20,-3.,3.);
    hisFitTPphiForAlgEff_[fitName]   = inputDir.make<TH1F>(("FitTPphiForAlgEff_"+(fitName)).c_str(),"; #phi of TP (used for alg. effi. measurement);",20,-M_PI,M_PI);

    // Histo for efficiency to reconstruct track perfectly (no incorrect hits). (Binning must match similar histos in bookTrackCands()).
    hisPerfFitTPinvptForAlgEff_[fitName] = inputDir.make<TH1F>(("PerfFitTPinvptForAlgEff_"+(fitName)).c_str() ,"; 1/Pt of TP (used for perf. alg. effi. measurement);",16,0.,maxAbsQoverPt);
    hisPerfFitTPetaForAlgEff_[fitName]   = inputDir.make<TH1F>(("PerfFitTPetaForAlgEff_"+(fitName)).c_str(),"; #eta of TP (used for perf. alg. effi. measurement);",20,-3.,3.);

    // Histos for algorithmic tracking efficiency vs. TP production point. (Binning must match similar histos in bookTrackCands()).
    hisFitTPd0ForAlgEff_[fitName]  = inputDir.make<TH1F>(("FitTPd0ForAlgEff_"+(fitName)).c_str() ,"; Production point x of TP (used for alg. effi. measurement);",50,0.,1.);
    hisFitTPz0ForAlgEff_[fitName]  = inputDir.make<TH1F>(("FitTPz0ForAlgEff_"+(fitName)).c_str() ,"; Production point y of TP (used for alg. effi. measurement);",50,0.,25.);

    // Histo for algorithmic tracking efficiency vs sector number (to check if looser cuts are needed in certain regions)
    unsigned int nPhi = settings_->numPhiSectors();
    unsigned int nEta = settings_->numEtaRegions();
    hisFitTPphisecForAlgEff_[fitName]  = inputDir.make<TH1F>(("FitTPphisecForAlgEff_"+(fitName)).c_str() ,"; #phi sector of TP (used for alg. effi. measurement);",nPhi,-0.5,nPhi-0.5);
    hisFitTPetasecForAlgEff_[fitName]  = inputDir.make<TH1F>(("FitTPetasecForAlgEff_"+(fitName)).c_str() ,"; #eta sector of TP (used for alg. effi. measurement);",nEta,-0.5,nEta-0.5);
    hisPerfFitTPphisecForAlgEff_[fitName]  = inputDir.make<TH1F>(("PerfFitTPphisecForAlgEff_"+(fitName)).c_str() ,"; #phi sector of TP (used for perf. alg. effi. measurement);",nPhi,-0.5,nPhi-0.5);
    hisPerfFitTPetasecForAlgEff_[fitName]  = inputDir.make<TH1F>(("PerfFitTPetasecForAlgEff_"+(fitName)).c_str() ,"; #eta sector of TP (used for perf. alg. effi. measurement);",nEta,-0.5,nEta-0.5);
  }
}

//=== Fill histograms for studying track fitting.

void Histos::fillTrackFitting( const InputData& inputData, const std::vector<std::pair<std::string,L1fittedTrack>>& mFittedTracks, float chi2dofCutPlots ) {  
 
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

  const vector<TP>&  vTPs = inputData.getTPs();

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

    for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : mFittedTracks ){

      std::string algoName(fittedTrackPair.first); // Get fitting algo name
      const L1fittedTrack& fitTrk = fittedTrackPair.second; // Get fitted track
      //IRT
      //      const TP*   assocTP =  fitTrk.getL1track3D().getMatchedTP(); // Get the TP the fitted track matches to, if any.
      const TP*   assocTP =  fitTrk.getMatchedTP(); // Get the TP the fitted track matches to, if any.

      if ( !fitTrk.accepted() ) continue;

      // Does this match desired TP?
      if (assocTP == &tp)  tpRecoed[algoName] = true;
      if (assocTP == &tp && fitTrk.chi2dof() <= chi2dofCutPlots) tpRecoedPass[algoName] = true;

      // Does this match desired TP & useForAlgEff
      if (assocTP == &tp && tp.useForAlgEff()) {
      tpRecoedForAlgEff[algoName] = true; // If it does, set tpRecoed to true
      // Also note if TP was reconstructed perfectly (no incorrect hits on reco track).
      if (fitTrk.getPurity() == 1.) tpRecoedForAlgEffPerfect[algoName] = true;
      if (fitTrk.chi2dof() <= chi2dofCutPlots) tpRecoedForAlgEffPassCut[algoName] = true; // Does track with assoc. TP pass the cut 
      if (fitTrk.getPurity() == 1. && fitTrk.chi2dof() <= chi2dofCutPlots) tpRecoedForAlgEffPerfectPassCut[algoName] = true; // Is the track passing the cut perfectly recoed?
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

  for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : mFittedTracks ){
   
    std::string j (fittedTrackPair.first);
    const L1fittedTrack& fitTrk = fittedTrackPair.second;

    // Get original HT track candidate prior to fit for comparison.
    const L1track3D& htTrk = fitTrk.getL1track3D();

    // Get matched truth particle, if any.
    const TP* tp = fitTrk.getMatchedTP();

    // Increment nFittedTrack and nStubsOnTrack counters
    nFittedTracks[j] += 1;
    if( fitTrk.chi2dof() <= chi2dofCutPlots && fitTrk.accepted() ) nStubsOnTrack[j] += fitTrk.getNumStubs();
 
    if ( fitTrk.accepted() ) {

      nTracksGenuine[j] += 1;
      if (tp != nullptr) {
	  nTracksGenuineTP[j] += 1;
	  hisFitVsSeedQinvPtGenCand_[j]->Fill( htTrk.qOverPt(), fitTrk.qOverPt() );
	  hisFitVsSeedPhi0GenCand_[j]->Fill( htTrk.phi0(), fitTrk.phi0() );
	  hisFitVsSeedD0GenCand_[j]->Fill( htTrk.d0(), fitTrk.d0() );
	  hisFitVsSeedZ0GenCand_[j]->Fill( htTrk.z0(), fitTrk.z0() );
	  hisFitVsSeedEtaGenCand_[j]->Fill( htTrk.eta(), fitTrk.eta() );
	  
	  if ( fitTrk.chi2dof() <= chi2dofCutPlots ){
	    nTracksGenuinePass[j] += 1;
	  }

	  // Check if chi2/NDF is well behaved for perfectly reconstructed tracks.
          if (fitTrk.getPurity() == 1.) {
            profChi2DofVsInvPtPERF_[j]->Fill(fabs(tp->qOverPt()), sqrt(fitTrk.chi2dof()));
	    profBigChi2DofVsInvPtPERF_[j]->Fill(fabs(tp->qOverPt()), (fitTrk.chi2dof() > 10));
	    // Are high Pt tracks sensitive to d0 impact parameter?
	    if (tp->pt() > 10.) {
              if (fitTrk.chi2dof() > 10.) {
    	        hisD0TPBigChi2DofPERF_[j]->Fill(fabs(tp->d0()));
 	      } else {
	        hisD0TPSmallChi2DofPERF_[j]->Fill(fabs(tp->d0()));
	      }
	    }
	  }
	}
    
	if (tp == nullptr){
	  nTracksFake[j] += 1;
	  if ( fitTrk.chi2dof() > chi2dofCutPlots ){
	    nTracksFakeCut[j] += 1;
	  }
	  hisFitVsSeedQinvPtFakeCand_[j]->Fill( htTrk.qOverPt(), fitTrk.qOverPt() );
	  hisFitVsSeedPhi0FakeCand_[j]->Fill( htTrk.phi0(), fitTrk.phi0() );
	  hisFitVsSeedD0FakeCand_[j]->Fill( htTrk.d0(), fitTrk.d0() );
	  hisFitVsSeedZ0FakeCand_[j]->Fill( htTrk.z0(), fitTrk.z0() );
	  hisFitVsSeedEtaFakeCand_[j]->Fill( htTrk.eta(), fitTrk.eta() );
	}
    }
    else if ( !fitTrk.accepted() ) {
      nRejectedTracks[j] += 1;
      if ( tp == nullptr ) nRejectedFake[j] += 1;
    }

    if ( fitTrk.accepted() && fitTrk.chi2dof() < chi2dofCutPlots ){
      // Study fake tracks.
      bool fakeTrk = (tp == nullptr);
      profFakeFitTrksVsEta_[j]->Fill(fabs(fitTrk.eta()), fakeTrk);
      // Study incorrect hits on good tracks.
      if ( ! fakeTrk) {
	const vector<const Stub*> stubs = fitTrk.getStubs();
	for (const Stub* s : stubs) {
	  // Was this stub produced by correct truth particle?
	  const set<const TP*> stubTPs = s->assocTPs();
	  bool trueStub = (stubTPs.find(tp) != stubTPs.end());
	  profFitFracTrueStubsVsLayer_[j]->Fill(s->layerId(), trueStub);
	  profFitFracTrueStubsVsEta_[j]->Fill(fabs(s->eta()), trueStub);
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
    hisNumMatchedStubsKilledVsKilled_[j]->Fill( fitTrk.getNumKilledStubs(), fitTrk.getNumKilledMatchedStubs() );  

    // --- Study purity against number of stubs to augment track fitter killing stubs due to large residual study
    hisNumStubsVsPurity_[j]->Fill( fitTrk.getNumStubs(), fitTrk.getPurity() );

    // Old way before accepted Kalman came along. 
    profTrksKilledByFit_[j]->Fill(ibin, !fitTrk.accepted());

    // Histograms specific to Linear Regression track fitter.
    if (j.find("LinearRegression") != string::npos) {

      int numIterations;
      std::string lostMatchingState;
      std::unordered_map< std::string, int > stateCalls;
      fitTrk.getInfoLR( numIterations, lostMatchingState, stateCalls );
      hisNumIterations_[j]->Fill( numIterations );
      if ( tp != tpHT )
	hisFailingState_[j]->Fill( lostMatchingState.c_str(), 1. );
      int totalCalls = 0;
      for ( auto state : stateCalls ) {
	totalCalls += state.second;
	hisTotalStateCalls_[ j ]->Fill( state.first.c_str(), state.second );
      }
      for ( auto state : stateCalls )
	hisRelativeStateCalls_[ j ]->Fill( state.first.c_str(), state.second / float(totalCalls) );
    }

    hisNumStubsFitKills_[j]->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs() ); 
    hisNumStubsFitKillsVsPurity_[j]->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs(), fitTrk.getPurity() ); 
    if (fitTrk.accepted()) hisNumStubsFitKillsVsPurityMatched_[j]->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs(), fitTrk.getPurity() );
    else hisNumStubsFitKillsVsPurityUnmatched_[j]->Fill( fitTrk.getL1track3D().getNumStubs() - fitTrk.getNumStubs(), fitTrk.getPurity() );

    if ( !fitTrk.accepted() ) continue; // If a rejected track, do not make plots for these.

    // Fill fitted parameter histograms and histograms of seed parameters against fitted parameters

    // Seed track parameter distributions
    hisSeedQinvPt_[j]->Fill( htTrk.qOverPt() );
    hisSeedPhi0_[j]->Fill( htTrk.phi0() );
    hisSeedD0_[j]->Fill( htTrk.d0() );
    hisSeedZ0_[j]->Fill( htTrk.z0() );
    hisSeedEta_[j]->Fill( htTrk.eta() );
    // Fitted track parameter distributions & chi2, separately for tracks that do/do not match a truth particle
    if ( tp != nullptr){
 
      hisFitQinvPtMatched_[j]->Fill( fitTrk.qOverPt() );
      hisFitPhi0Matched_[j]->Fill( fitTrk.phi0() );
      hisFitD0Matched_[j]->Fill( fitTrk.d0() );
      hisFitZ0Matched_[j]->Fill( fitTrk.z0() );
      hisFitEtaMatched_[j]->Fill( fitTrk.eta() );
 
      hisFitChi2Matched_[j]->Fill( fitTrk.chi2() );
      hisFitChi2DofMatched_[j]->Fill( fitTrk.chi2dof() );
 
    } else {
 
      hisFitQinvPtUnmatched_[j]->Fill( fitTrk.qOverPt() );
      hisFitPhi0Unmatched_[j]->Fill( fitTrk.phi0() );
      hisFitD0Unmatched_[j]->Fill( fitTrk.d0() );
      hisFitZ0Unmatched_[j]->Fill( fitTrk.z0() );
      hisFitEtaUnmatched_[j]->Fill( fitTrk.eta() );
 
      hisFitChi2Unmatched_[j]->Fill( fitTrk.chi2() );
      hisFitChi2DofUnmatched_[j]->Fill( fitTrk.chi2dof() );
 
    }
       
    // If there is an associated tracking particle, fill up histograms of TP's parameters against fitted parameters
    if ( tp != nullptr ){

      // Do seperately for those with good/poor chi2.
      if ( fitTrk.chi2dof() <= chi2dofCutPlots ){
	// Fitted vs True parameter distribution 2D plots
	hisFitVsTrueQinvPtGoodChi2_[j]->Fill( tp->qOverPt(), fitTrk.qOverPt() );
	hisFitVsTruePhi0GoodChi2_[j]->Fill( tp->phi0(), fitTrk.phi0( ));
	hisFitVsTrueD0GoodChi2_[j]->Fill( tp->d0(), fitTrk.d0() );
	hisFitVsTrueZ0GoodChi2_[j]->Fill( tp->z0(), fitTrk.z0() );
	hisFitVsTrueEtaGoodChi2_[j]->Fill( tp->eta(), fitTrk.eta() );
	// Residuals between fitted and true helix params as 1D plot.
	hisFitQinvPtResGoodChi2_[j]->Fill( fitTrk.qOverPt() - tp->qOverPt());
	hisFitPhi0ResGoodChi2_[j]->Fill( reco::deltaPhi(fitTrk.phi0(), tp->phi0()) );
	hisFitD0ResGoodChi2_[j]->Fill( fitTrk.d0() - tp->d0() );
	hisFitZ0ResGoodChi2_[j]->Fill( fitTrk.z0() - tp->z0() );
	hisFitEtaResGoodChi2_[j]->Fill( fitTrk.eta() - tp->eta() );  
	// Residuals between true and seed helix params as 1D plot.
	hisSeedQinvPtResGoodChi2_[j]->Fill( tp->qOverPt() - htTrk.qOverPt());
	hisSeedPhi0ResGoodChi2_[j]->Fill( reco::deltaPhi(tp->phi0(), htTrk.phi0()) );
	hisSeedD0ResGoodChi2_[j]->Fill( tp->d0() - htTrk.d0() );
	hisSeedZ0ResGoodChi2_[j]->Fill( tp->z0() - htTrk.z0() );
	hisSeedEtaResGoodChi2_[j]->Fill( tp->eta() - htTrk.eta() );  

	// Understand which matched tracks have good/bad chi2.
	hisTrueEtaMatchedGoodChi2_[j]->Fill( tp->eta() );
	hisStubPurityMatchedGoodChi2_[j]->Fill( fitTrk.getPurity() );
      } else {
 
	// Plot rapidity of matched tracks with bad chi2.
	hisTrueEtaMatchedBadChi2_[j]->Fill( tp->eta() );
	hisStubPurityMatchedBadChi2_[j]->Fill( fitTrk.getPurity() );
	if (fitTrk.getPurity() > 0.99) {    
	  // These tracks have no bad hits. So why do they have bad chi2?
	  /*
	    cout<<"WIERD: chi2/dof "<<fitTrk.chi2dof()<<" "<<fitTrk.getNumStubs()<<" "<<fitTrk.getNumLayers()<<endl;
	    cout<<"WIERD: ht "<<htTrk.pt()<<" "<<htTrk.eta()<<" "<<htTrk.phi0()<<" "<<htTrk.d0()<<" "<<htTrk.z0()<<endl;
	    cout<<"WIERD: fit "<<fitTrk.pt()<<" "<<fitTrk.eta()<<" "<<fitTrk.phi0()<<" "<<fitTrk.d0()<<" "<<fitTrk.z0()<<endl;
	    cout<<"WIERD: tp "<<tp->pt()<<" "<<tp->eta()<<" "<<tp->phi0()<<" "<<tp->d0()<<" "<<tp->z0()<<endl;
	    cout<<endl;
	  */
	}
      }

      // Plot helix parameter resolution against eta or Pt.

      if (fitTrk.chi2dof() <= chi2dofCutPlots) {
        hisQoverPtResVsTrueEta_[j]->Fill( std::abs(tp->eta()), std::abs( fitTrk.qOverPt() - tp->qOverPt() ) );
        hisPhi0ResVsTrueEta_[j]->Fill( std::abs(tp->eta()), std::abs(reco::deltaPhi(fitTrk.phi0(), tp->phi0()) ) );
        hisEtaResVsTrueEta_[j]->Fill( std::abs(tp->eta()), std::abs( fitTrk.eta() - tp->eta() ) );
        hisZ0ResVsTrueEta_[j]->Fill( std::abs(tp->eta()), std::abs( fitTrk.z0() - tp->z0() ) );
        hisD0ResVsTrueEta_[j]->Fill( std::abs(tp->eta()), std::abs( fitTrk.d0() - tp->d0() ) );

        hisQoverPtResVsTrueInvPt_[j]->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.qOverPt() - tp->qOverPt() ) );
        hisPhi0ResVsTrueInvPt_[j]->Fill( std::abs(tp->qOverPt()), std::abs(reco::deltaPhi(fitTrk.phi0(), tp->phi0()) ) );
        hisEtaResVsTrueInvPt_[j]->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.eta() - tp->eta() ) );
        hisZ0ResVsTrueInvPt_[j]->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.z0() - tp->z0() ) );
        hisD0ResVsTrueInvPt_[j]->Fill( std::abs(tp->qOverPt()), std::abs( fitTrk.d0() - tp->d0() ) );

        hisQoverPtResVsTrueZ0_[j]->Fill( std::abs(tp->z0()), std::abs( fitTrk.qOverPt() 
          - tp->z0() ) );
        hisPhi0ResVsTrueZ0_[j]->Fill( std::abs(tp->z0()), std::abs(reco::deltaPhi(
          fitTrk.phi0(), tp->phi0()) ) );
        hisEtaResVsTrueZ0_[j]->Fill( std::abs(tp->z0()), std::abs( fitTrk.eta() - tp->
          eta() ) );
        hisZ0ResVsTrueZ0_[j]->Fill( std::abs(tp->z0()), std::abs( fitTrk.z0() - tp->
          z0() ) );
        hisD0ResVsTrueZ0_[j]->Fill( std::abs(tp->z0()), std::abs( fitTrk.d0() - tp->
          d0() ) );

      }

      // Plot chi^2 vs eta, and # stubs vs eta.
      hisTrueFittedChiSquaredVsTrueEta_[j]->Fill( fitTrk.chi2(), fitTrk.getL1track3D().eta() );
      hisTrueFittedChiSquaredDofVsTrueEta_[j]->Fill( fitTrk.chi2dof(), fitTrk.getL1track3D().eta() );
      hisTrueFittedChiSquaredVsFittedEta_[j]->Fill( fitTrk.chi2(), fitTrk.eta() );
      hisTrueFittedChiSquaredDofVsFittedEta_[j]->Fill( fitTrk.chi2dof(), fitTrk.eta() );
 
      hisFittedChiSquaredFunctionOfStubs_[j]->Fill( fitTrk.getStubs().size(), fitTrk.chi2() );
      hisFittedChiSquaredDofFunctionOfStubs_[j]->Fill( fitTrk.getStubs().size(), fitTrk.chi2dof() );


 
    }
  }

  for ( const string& fitterName : settings_->trackFitters() ){

    profNumFittedCands_[fitterName]->Fill(1.0, nFittedTracks[fitterName]);
    profNumFittedCands_[fitterName]->Fill(2.0, nStubsOnTrack[fitterName]);
    profNumFittedCands_[fitterName]->Fill(3.0, nTracksGenuine[fitterName]);
    profNumFittedCands_[fitterName]->Fill(4.0, nTracksGenuineTP[fitterName]);
    profNumFittedCands_[fitterName]->Fill(5.0, nTracksGenuinePass[fitterName]);
    profNumFittedCands_[fitterName]->Fill(6.0, nTracksFakeCut[fitterName]);
    profNumFittedCands_[fitterName]->Fill(7.0, nTracksFake[fitterName]);
    profNumFittedCands_[fitterName]->Fill(8.0, nRejectedTracks[fitterName]);
    profNumFittedCands_[fitterName]->Fill(9.0, nRejectedFake[fitterName]);
    profNumFittedCands_[fitterName]->Fill(10.0, nTracksExcDups[fitterName]);
    profNumFittedCands_[fitterName]->Fill(11.0, nTracksExcDupsPass[fitterName]);
  }

  //=== Study duplicate tracks. IRT

  for (const TP& tp: vTPs) {

    for (auto &fitName : settings_->trackFitters()) {

      unsigned int nMatch = 0;
      vector<const L1fittedTrack*> fitTrkStore;

      for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : mFittedTracks ){
   
	std::string j (fittedTrackPair.first);
	const L1fittedTrack& fitTrk = fittedTrackPair.second;

	if (j == fitName) {

	  if (fitTrk.accepted() && fitTrk.chi2dof() < chi2dofCutPlots) {

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
        profDupFitTrksVsEta_[fitName]->Fill(fabs(tp.eta()), dup); 
        profDupFitTrksVsInvPt_[fitName]->Fill(fabs(tp.qOverPt()), dup); 
/*
	if (dup) {
	  cout<<"DUPFIT "<<fitName<<" : TP pt = "<<fabs(tp.qOverPt())<<" eta = "<<fabs(tp.eta())<<" :";
	  unsigned int nT = 0;
	  unsigned int nTok = 0;
	  for (const L1fittedTrack* trkTmp : fitTrkStore) {
	    bool consistentCell = trkTmp->consistentHTcell();
	    if (consistentCell) nTok++;
	    cout<<" T"<<++nT<<" consistent = "<<consistentCell;   
	  }
	  bool noDupIn1stPass = (nTok <= 1); // check if this would have been a duplicate if only 1st pass of algo50 duplicate removal had been run.
	  cout<<" noDupIn1stPass = "<<noDupIn1stPass<<endl;
	}
*/
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

	for ( const std::pair<std::string,L1fittedTrack>& fittedTrackPair : mFittedTracks ){
   
	  const std::string j (fittedTrackPair.first);
	  const L1fittedTrack& fitTrk = fittedTrackPair.second;
    
	  if (j == fitName) {

	    if (fitTrk.accepted() && fitTrk.chi2dof() < chi2dofCutPlots) {

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
	  hisFitTPinvptForEff_[fitName]->Fill(1./tp.pt());
	  hisFitTPetaForEff_[fitName]->Fill(tp.eta());
	  hisFitTPphiForEff_[fitName]->Fill(tp.phi0());
	  // Also plot efficiency to perfectly reconstruct the track (no fake hits)
	  if (tpRecoedPerfect) {
	    hisPerfFitTPinvptForEff_[fitName]->Fill(1./tp.pt());
	    hisPerfFitTPetaForEff_[fitName]->Fill(tp.eta());
	  }
	  if (tp.useForAlgEff()) { // Check TP is good for algorithmic efficiency measurement.
	    hisFitTPinvptForAlgEff_[fitName]->Fill(1./tp.pt());
	    hisFitTPetaForAlgEff_[fitName]->Fill(tp.eta());
	    hisFitTPphiForAlgEff_[fitName]->Fill(tp.phi0());
	    // Plot also production point of all good reconstructed TP.
	    hisFitTPd0ForAlgEff_[fitName]->Fill(fabs(tp.d0()));
	    hisFitTPz0ForAlgEff_[fitName]->Fill(fabs(tp.z0()));
	    // Plot sector number to understand if looser cuts are needed in certain regions.
	    hisFitTPphisecForAlgEff_[fitName]->Fill(iPhiSec_TP);
	    hisFitTPetasecForAlgEff_[fitName]->Fill(iEtaReg_TP);
	    // Also plot efficiency to perfectly reconstruct the track (no fake hits)
	    if (tpRecoedPerfect) {
	      hisPerfFitTPinvptForAlgEff_[fitName]->Fill(1./tp.pt());
	      hisPerfFitTPetaForAlgEff_[fitName]->Fill(tp.eta());
  	      hisPerfFitTPphisecForAlgEff_[fitName]->Fill(iPhiSec_TP);
	      hisPerfFitTPetasecForAlgEff_[fitName]->Fill(iEtaReg_TP);
	    }
	  }
	}
      }
    }
  }
}

//=== Produce plots of tracking efficiency prior to track fit (run at end of job).

void Histos::plotTrackEfficiency() {

  TFileDirectory inputDir = fs_->mkdir("Effi");
  // Plot tracking efficiency
  makeEfficiencyPlot(inputDir, graphEffVsInvPt_, hisRecoTPinvptForEff_, hisTPinvptForEff_, "EffVsInvPt", "; 1/Pt; Tracking efficiency" );
  makeEfficiencyPlot(inputDir, graphEffVsInvPt_, hisRecoTPetaForEff_, hisTPetaForEff_, "EffVsEta","; #eta; Tracking efficiency" );

  // std::cout << "Made first graph" << std::endl;
  // graphEffVsEta_ = inputDir.make<TGraphAsymmErrors>(hisRecoTPetaForEff_, hisTPetaForEff_);
  // graphEffVsEta_->SetNameTitle("EffVsEta","; #eta; Tracking efficiency");
  makeEfficiencyPlot(inputDir, graphEffVsInvPt_, hisRecoTPphiForEff_, hisTPphiForEff_,"EffVsPhi","; #phi; Tracking efficiency");

  // Also plot efficiency to reconstruct track perfectly.
  makeEfficiencyPlot( inputDir, graphPerfEffVsInvPt_, hisPerfRecoTPinvptForEff_, hisTPinvptForEff_, "PerfEffVsInvPt","; 1/Pt; Tracking perfect efficiency");
  makeEfficiencyPlot( inputDir, graphPerfEffVsEta_, hisPerfRecoTPetaForEff_, hisTPetaForEff_, "PerfEffVsEta","; #eta; Tracking perfect efficiency");

  // Plot algorithmic tracking efficiency
  makeEfficiencyPlot( inputDir, graphAlgEffVsInvPt_, hisRecoTPinvptForAlgEff_, hisTPinvptForAlgEff_, "AlgEffVsInvPt","; 1/Pt; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffVsEta_, hisRecoTPetaForAlgEff_, hisTPetaForAlgEff_, "AlgEffVsEta","; #eta; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffVsPhi_, hisRecoTPphiForAlgEff_, hisTPphiForAlgEff_, "AlgEffVsPhi","; #phi; Tracking efficiency");

  makeEfficiencyPlot( inputDir, graphAlgEffVsD0_, hisRecoTPd0ForAlgEff_, hisTPd0ForAlgEff_, "AlgEffVsD0","; d0; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffVsZ0_, hisRecoTPz0ForAlgEff_, hisTPz0ForAlgEff_, "AlgEffVsZ0","; z0; Tracking efficiency");

  makeEfficiencyPlot( inputDir, graphAlgEffVsPhiSec_, hisRecoTPphisecForAlgEff_, hisTPphisecForAlgEff_, "AlgEffVsPhiSec","; #phi sector; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffVsEtaSec_, hisRecoTPetasecForAlgEff_, hisTPetasecForAlgEff_, "AlgEffVsEtaSec","; #eta sector; Tracking efficiency");

  // Also plot algorithmic efficiency to reconstruct track perfectly.
  makeEfficiencyPlot( inputDir, graphPerfAlgEffVsInvPt_, hisPerfRecoTPinvptForAlgEff_, hisTPinvptForAlgEff_, "PerfAlgEffVsInvPt","; 1/Pt; Tracking perfect efficiency");
  makeEfficiencyPlot( inputDir, graphPerfAlgEffVsEta_, hisPerfRecoTPetaForAlgEff_, hisTPetaForAlgEff_, "PerfAlgEffVsEta","; #eta; Tracking perfect efficiency");

  // cout << "graphPerfAlgEffVsPhiSec_" << endl;
  // graphPerfAlgEffVsPhiSec_ = inputDir.make<TGraphAsymmErrors>(hisPerfRecoTPphisecForAlgEff_, hisTPphisecForAlgEff_);
  // graphPerfAlgEffVsPhiSec_->SetNameTitle("PerfAlgEffVsPhiSec","; #phi sector; Tracking perfect efficiency");

  // cout << "graphPerfAlgEffVsEtaSec_"<< endl;
  // graphPerfAlgEffVsEtaSec_ = inputDir.make<TGraphAsymmErrors>(hisPerfRecoTPetasecForAlgEff_, hisTPetasecForAlgEff_);
  // graphPerfAlgEffVsEtaSec_->SetNameTitle("PerfAlgEffVsEtaSec","; #eta sector; Tracking perfect efficiency");
}

//=== Produce plots of tracking efficiency after track fit (run at end of job).

void Histos::plotTrackEffAfterFit(string fitName) {
  TFileDirectory inputDir = fs_->mkdir(("Effi_"+(fitName)).c_str());
  std::cout << "Tracking eff after fit" << std::endl;
  // Plot tracking efficiency
  makeEfficiencyPlot( inputDir, graphEffFitVsInvPt_[fitName], hisFitTPinvptForEff_[fitName], hisTPinvptForEff_, ("EffFitVsInvPt_"+(fitName)).c_str(),"; 1/Pt; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphEffFitVsEta_[fitName], hisFitTPetaForEff_[fitName], hisTPetaForEff_, ("EffFitVsEta_"+(fitName)).c_str(),"; #eta; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphEffFitVsPhi_[fitName], hisFitTPphiForEff_[fitName], hisTPphiForEff_, ("EffFitVsPhi_"+(fitName)).c_str(),"; #phi; Tracking efficiency");

  // Also plot efficiency to reconstruct track perfectly.
  makeEfficiencyPlot( inputDir, graphPerfEffFitVsInvPt_[fitName], hisPerfFitTPinvptForEff_[fitName], hisTPinvptForEff_, ("PerfEffFitVsInvPt_"+(fitName)).c_str(),"; 1/Pt; Tracking perfect efficiency");
  makeEfficiencyPlot( inputDir, graphPerfEffFitVsEta_[fitName], hisPerfFitTPetaForEff_[fitName], hisTPetaForEff_, ("PerfEffFitVsEta_"+(fitName)).c_str(),"; #eta; Tracking perfect efficiency");

  // Plot algorithmic tracking efficiency
  makeEfficiencyPlot( inputDir, graphAlgEffFitVsInvPt_[fitName], hisFitTPinvptForAlgEff_[fitName], hisTPinvptForAlgEff_, ("AlgEffFitVsInvPt_"+(fitName)).c_str(),"; 1/Pt; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffFitVsEta_[fitName], hisFitTPetaForAlgEff_[fitName], hisTPetaForAlgEff_, ("AlgEffFitVsEta_"+(fitName)).c_str(),"; #eta; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffFitVsPhi_[fitName], hisFitTPphiForAlgEff_[fitName], hisTPphiForAlgEff_, ("AlgEffFitVsPhi_"+(fitName)).c_str(),"; #phi; Tracking efficiency");

  makeEfficiencyPlot( inputDir, graphAlgEffFitVsD0_[fitName], hisFitTPd0ForAlgEff_[fitName], hisTPd0ForAlgEff_, ("AlgEffFitVsD0_"+(fitName)).c_str(),"; d0; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffFitVsZ0_[fitName], hisFitTPz0ForAlgEff_[fitName], hisTPz0ForAlgEff_, ("AlgEffFitVsZ0_"+(fitName)).c_str(),"; z0; Tracking efficiency");

  makeEfficiencyPlot( inputDir, graphAlgEffFitVsPhiSec_[fitName], hisFitTPphisecForAlgEff_[fitName], hisTPphisecForAlgEff_, ("AlgEffFitVsPhiSec_"+(fitName)).c_str(),"; #phi sector; Tracking efficiency");
  makeEfficiencyPlot( inputDir, graphAlgEffFitVsEtaSec_[fitName], hisFitTPetasecForAlgEff_[fitName], hisTPetasecForAlgEff_, ("AlgEffFitVsEtaSec_"+(fitName)).c_str(),"; #eta sector; Tracking efficiency");

  // Also plot algorithmic efficiency to reconstruct track perfectly.
  makeEfficiencyPlot( inputDir, graphPerfAlgEffFitVsInvPt_[fitName], hisPerfFitTPinvptForAlgEff_[fitName], hisTPinvptForAlgEff_, ("PerfAlgEffFitVsInvPt_"+(fitName)).c_str(),"; 1/Pt; Tracking perfect efficiency");
  makeEfficiencyPlot( inputDir, graphPerfAlgEffFitVsEta_[fitName], hisPerfFitTPetaForAlgEff_[fitName], hisTPetaForAlgEff_, ("Perf AlgEffFitVsEta_"+(fitName)).c_str(),"; #eta; Tracking perfect efficiency");

  // cout << "graphPerfAlgEffFitVsPhiSec_" << endl;
  // graphPerfAlgEffFitVsPhiSec_[fitName] = inputDir.make<TGraphAsymmErrors>(hisPerfFitTPphisecForAlgEff_[fitName], hisTPphisecForAlgEff_);
  // graphPerfAlgEffFitVsPhiSec_[fitName]->SetNameTitle(("PerfAlgEffFitVsPhiSec_"+(fitName)).c_str(),"; #phi sector; Tracking perfect efficiency");
  
  // cout << "graphPerfAlgEffFitVsEtaSec_ "<< endl;
  // graphPerfAlgEffFitVsEtaSec_[fitName] = inputDir.make<TGraphAsymmErrors>(hisPerfFitTPetasecForAlgEff_[fitName], hisTPetasecForAlgEff_);

  // graphPerfAlgEffFitVsEtaSec_[fitName]->SetNameTitle(("PerfAlgEffFitVsEtaSec_"+(fitName)).c_str(),"; #eta sector; Tracking perfect efficiency");
}

void Histos::makeEfficiencyPlot( TFileDirectory &inputDir, TGraphAsymmErrors* outputEfficiency, TH1F* pass, TH1F* all, TString name, TString title ) {
  TEfficiency tempEfficiency( *pass, *all );
  tempEfficiency.Draw();
  gPad->Update();
  cout << name << endl;
  outputEfficiency = inputDir.make<TGraphAsymmErrors>( *tempEfficiency.GetPaintedGraph() );
  outputEfficiency->SetNameTitle(name,title);
}


//=== Book efficiency & fake rate histograms and print track-finding performance summary.

void Histos::endJobAnalysis() {

  // Produce plots of tracking efficiency using track candidates found prior to track fit.
  this->plotTrackEfficiency();

  // // Produce more plots of tracking efficiency using track candidates after track fit.
  for (auto &fitName : settings_->trackFitters()) {
    this->plotTrackEffAfterFit(fitName);
  }

  //--- Print summary of track-finding performance (prior to helix fit).

  float numTrackCands = profNumTrackCands_->GetBinContent(1); // No. of track cands
  float numTrackCandsErr = profNumTrackCands_->GetBinError(1); // No. of track cands uncertainty
  float numTPTrackCandsIncDups = profNumTrackCands_->GetBinContent(2); // Ditto, counting only those matched to TP
  float numTPTrackCandsExcDups = profNumTrackCands_->GetBinContent(6); // Ditto, but excluding duplicates
  float numTrackCandsSF = profNumTrackCands_->GetBinContent(8);
  float numTrackCandsSFErr = profNumTrackCands_->GetBinError(8);
  float numFakeTracks = numTrackCands - numTPTrackCandsIncDups;
  float numExtraDupTracks = numTPTrackCandsIncDups - numTPTrackCandsExcDups;
  float fracFake = numFakeTracks/(numTrackCands + 1.0e-6);
  float fracDup = numExtraDupTracks/(numTrackCands + 1.0e-6);

  float numStubsOnTracks = profStubsOnTracks_->GetBinContent(1);
  float meanStubsPerTrack = numStubsOnTracks/(numTrackCands + 1.0e-6); //protection against demoninator equals zero.
  unsigned int numRecoTPforAlg = hisRecoTPinvptForAlgEff_->GetEntries();
  unsigned int numTPforAlg     = hisTPinvptForAlgEff_->GetEntries();
  float algEff = float(numRecoTPforAlg)/(numTPforAlg + 1.0e-6); //protection against demoninator equals zero.
  float algEffErr = sqrt(algEff*(1-algEff)/(numTPforAlg + 1.0e-6)); // uncertainty
  float algPerfEff = float(numPerfRecoTPforAlg_)/(numTPforAlg + 1.0e-6); //protection against demoninator equals zero.
  float algPerfEffErr = sqrt(algPerfEff*(1-algPerfEff)/(numTPforAlg + 1.0e-6)); // uncertainty
  cout<<"========================================================================="<<endl;
  cout<<"               TRACK-FINDING SUMMARY (before track fit)                  "<<endl;
  cout<<"Number of track candidates before rz filter = "<< numTrackCandsSF<<" +/- "<<numTrackCandsSFErr << endl;
  cout<<"Number of track candidates found per event = "<<numTrackCands<<" +- "<<numTrackCandsErr<<endl;
  cout<<"                     with mean stubs/track = "<<meanStubsPerTrack<<endl; 
  cout<<"Fraction of track cands that are fake = "<<fracFake<<endl;
  cout<<"Fraction of track cands that are genuine, but extra duplicates = "<<fracDup<<endl;
  cout<<"Algorithmic tracking efficiency = "<<numRecoTPforAlg<<"/"<<numTPforAlg<<" = "<<algEff<<" +- "<<algEffErr<<endl;
  cout<<"Perfect algorithmic tracking efficiency = "<<numPerfRecoTPforAlg_<<"/"<<numTPforAlg<<" = "<<algPerfEff<<" +- "<<algPerfEffErr<<" (no incorrect hits)"<<endl;
  cout<<"========================================================================="<<endl;

  //--- Print summary of track-finding performance after helix fit, for each track fitting algorithm used.
   
  for (auto &fitName : settings_->trackFitters()) {

    float numFittedTracks = profNumFittedCands_[fitName]->GetBinContent(1); // mean fitted tracks/event, no chi2/ndf cut.
    float numStubsOnTrack = profNumFittedCands_[fitName]->GetBinContent(2); // Number of stubs passing chi2/ndf cut.
    float numFittedTracksAccepted = profNumFittedCands_[fitName]->GetBinContent(3); // # of accepted fitted tracks
    float numGenuineFittedTracksTP = profNumFittedCands_[fitName]->GetBinContent(4); // # of accepted fitted tracks with assoc. TP
    float numGenuineFittedTracksPass = profNumFittedCands_[fitName]->GetBinContent(5); // # of accepted fitted tracks with assoc. TP that passed chi2nof cut
    float numGenuineFittedTracksPassErr = profNumFittedCands_[fitName]->GetBinError(5); // Error on # of accepted fitted tracks with assoc. TP that passed chi2nof cut
    float numFakeFittedTracksCut = profNumFittedCands_[fitName]->GetBinContent(6); // # of accepted fitted tracks with no assoc TP killed by chi2nof cut
    float numFakeFittedTracksCutErr = profNumFittedCands_[fitName]->GetBinError(6); // Error on # of accepted fitted tracks with no assoc TP killed by chi2nof cut
    float numFakeFittedTracks = profNumFittedCands_[fitName]->GetBinContent(7); // Number of tracks/accepted tracks  (same thing hopefully) - tracks with min 5 layers with assoc TP = number of fakes!
    float numFakeFittedTracksErr = profNumFittedCands_[fitName]->GetBinError(7); // Error on number of tracks/accepted tracks  (same thing hopefully) - tracks with min 5 layers with assoc TP = number of fakes!
    float numRejectedTracks = profNumFittedCands_[fitName]->GetBinContent(8); // Number of tracks which were rejected due to stub killing
    float numRejectedFakeTracks = profNumFittedCands_[fitName]->GetBinContent(9); // Number of fake tracks which were rejected due to stub killing
    float numTracksExcDups = profNumFittedCands_[fitName]->GetBinContent(10); // Number of track candidates excluding duplicates.
    float numTracksExcDupsPass = profNumFittedCands_[fitName]->GetBinContent(11); // Number of track candidates excluding duplicates passing chi2/ndf cut.

    float numFakeTracksAcceptedPass = numFakeFittedTracks-numFakeFittedTracksCut; // Number of fake (accepted) tracks not killed by the chi2/ndf cut.
    float numFittedTracksPass = numGenuineFittedTracksPass+numFakeTracksAcceptedPass;
    float numFittedTracksPassErr = sqrt(numGenuineFittedTracksPassErr*numGenuineFittedTracksPassErr + numFakeFittedTracksErr*numFakeFittedTracksErr + numFakeFittedTracksCutErr*numFakeFittedTracksCutErr);

    float numStubsPerFittedTrack = numStubsOnTrack/(numFittedTracksPass + 1.0e-6);

    float fracDupFitted = (numGenuineFittedTracksTP-numTracksExcDups)/(numFittedTracks + 1.0e-6);
    float fracDupFittedPass = (numGenuineFittedTracksPass-numTracksExcDupsPass)/(numFittedTracksPass + 1.0e-6);

    float fracFittedNoTP = numFakeFittedTracks/(numFittedTracksAccepted + 1.0e-6);
    float fracGenuineFittedPass = numGenuineFittedTracksPass/(numGenuineFittedTracksTP + 1.0e-6);
    float fracFakeFittedCut = numFakeFittedTracksCut/(numFakeFittedTracks + 1.0e-6);
    float fracRejectedFake = numRejectedFakeTracks/(numRejectedTracks + 1.0e-6);
    float fracTracksTooFewStubs = numRejectedTracks/(numFittedTracks + 1.0e-6);

    // Pre-cut algo and perfect algo efficiencies and errors
    float fitEff = float(numFitAlgEff_[fitName])/(numTPforAlg + 1.0e-6); //protection against demoninator equals zero.
    float fitEffErr = sqrt(fitEff*(1-fitEff)/(numTPforAlg + 1.0e-6)); //protection against demoninator equals zero.
    float fitPureEff = float(numFitPerfAlgEff_[fitName])/(numTPforAlg + 1.0e-6); //protection against demoninator equals zero.
    float fitPureEffErr = sqrt(fitPureEff*(1-fitPureEff)/(numTPforAlg + 1.0e-6)); //protection against demoninator equals zero.

    // Post-cut algo and perfect algo efficiencies and errors
    float fitEffPass = float(numFitAlgEffPass_[fitName])/(numTPforAlg + 1.0e-6); //protection against demoninator equals zero.
    float fitEffPassErr = sqrt(fitEffPass*(1-fitEffPass)/(numTPforAlg + 1.0e-6)); //protection against demoninator equals zero.
    float fitPureEffPass = float(numFitPerfAlgEffPass_[fitName])/(numTPforAlg + 1.0e-6); //protection against demoninator equals zero.
    float fitPureEffPassErr = sqrt(fitPureEffPass*(1-fitPureEffPass)/(numTPforAlg + 1.0e-6)); //protection against demoninator equals zero.

    cout << "                    GENERAL FITTING SUMMARY FOR: " << fitName << endl;
    cout << "Number of fitted track candidates found per event = " << numFittedTracksPass << " +- " << numFittedTracksPassErr << endl;
    cout << "                            with mean stubs/track = " << numStubsPerFittedTrack << endl;
    cout << "Fraction of fitted tracks that are fake = " << float(numFakeTracksAcceptedPass)/(numFittedTracksPass + 1.0e-6) << endl;
    cout << "Fraction of fitted tracks that are genuine, but extra duplicates (post-cut) = " << fracDupFittedPass << endl;
    cout << "Algorithmic fitting efficiency (post-cut) = " << numFitAlgEffPass_[fitName] << "/" << numTPforAlg << " = " << fitEffPass << " +- " << fitEffPassErr << endl;
    cout << "Perfect algorithmic fitting efficiency(post-cut) = " << numFitPerfAlgEffPass_[fitName] << "/" <<  numTPforAlg << " = " << fitPureEffPass << " +- " << fitPureEffPassErr << " (no incorrect hits)" << endl;

    if ( settings_->detailedFitOutput() ){
      cout << endl<< "More detailed information about helix fit:"<<endl<<endl;
      cout << "Fraction of fitted tracks that are genuine, but extra duplicates (pre-cut)= " << numGenuineFittedTracksTP-numTracksExcDups << "/" << numFittedTracks << " = " << fracDupFitted << endl;
      cout << "Algorithmic fitting efficiency (pre-cut) = " << numFitAlgEff_[fitName] << "/" << numTPforAlg << " = " << fitEff << " +- " << fitEffErr << endl;
      cout << "Perfect algorithmic fitting efficiency(pre-cut) = " << numFitPerfAlgEff_[fitName] << "/" <<  numTPforAlg << " = " << fitPureEff << " +- " << fitPureEffErr << " (no incorrect hits)" << endl;
      cout << "----------------Candidates with enough stubs after fit-------------------" << endl;
      cout << "Fraction of fitted tracks which are fake (before chi2/dof cut) = " << numFakeFittedTracks << "/" << numFittedTracksAccepted << " = " << fracFittedNoTP << endl;
      cout << "Fraction of genuine fitted tracks that pass chi2/dof cut = " << numGenuineFittedTracksPass << "/" << numGenuineFittedTracksTP << " = " << fracGenuineFittedPass << endl;
      cout << "Fraction of fake fitted tracks killed by chi2/dof cut = " << numFakeFittedTracksCut << "/" << numFakeFittedTracks << " = " << fracFakeFittedCut << endl;
      cout << "----------------Candidates with too few stubs after fit------------------" << endl;
      cout << "Fraction of fitted tracks with too few stubs = " << numRejectedTracks << "/" << numFittedTracks << " = " << fracTracksTooFewStubs << endl;;
      cout << "Fraction of fitted tracks with too few stubs which are fakes = " << numRejectedFakeTracks << "/" << numRejectedTracks << " = " << fracRejectedFake << endl;
    }
    cout << "=========================================================================" << endl;
  }

  // Check that stub filling was consistent with known limitations of firmware design.

  cout<<endl<<"Max. |gradients| of stub lines in HT arrays are: r-phi = "<<HTrphi::maxLineGrad()<<", r-z = "<<HTrz::maxLineGrad()<<endl;

  if (HTrphi::maxLineGrad() > 1. || HTrz::maxLineGrad() > 1.) {

    cout<<"WARNING: Line |gradient| exceeds 1, which firmware will not be able to cope with! Please adjust HT array size to avoid this."<<endl;

  } else if (HTrphi::fracErrorsTypeA() > 0. || HTrz::fracErrorsTypeA() > 0.) {

    cout<<"WARNING: Despite line gradients being less than one, some fraction of HT columns have filled cells with no filled neighbours in W, SW or NW direction. Firmware will object to this! ";
    cout<<"This fraction = "<<HTrphi::fracErrorsTypeA()<<" for r-phi HT & "<<HTrz::fracErrorsTypeA()<<" for r-z HT"<<endl; 

  } else if (HTrphi::fracErrorsTypeB() > 0. || HTrz::fracErrorsTypeB() > 0.) {

    cout<<"WARNING: Despite line gradients being less than one, some fraction of HT columns recorded individual stubs being added to more than two cells! Thomas firmware will object to this! "; 
    cout<<"This fraction = "<<HTrphi::fracErrorsTypeB()<<" for r-phi HT & "<<HTrz::fracErrorsTypeB()<<" for r-z HT"<<endl;   
  }

  // Check for presence of common MC bug.

  float meanShared = hisFracStubsSharingClus0_->GetMean();
  if (meanShared > 0.01) cout<<endl<<"WARNING: You are using buggy MC. A fraction "<<meanShared<<" of stubs share clusters in the module seed sensor, which front-end electronics forbids."<<endl;

 
  // Vertex Efficiency
  TFileDirectory inputDir= fs_->mkdir("VertexEfficiency");

  PVefficiencyVsTrueZ0_ = inputDir.make<TEfficiency>(*hisRecoPrimaryVertexVsTrueZ0_,*hisPrimaryVertexTrueZ0_);
  PVefficiencyVsTrueZ0_->SetNameTitle("PVefficiencyVsTrueZ0_", "Primary Vertex Finding Efficiency; true z_{0}; Efficiency");

  PVefficiencyVsNoTracks_ = inputDir.make<TEfficiency>(*hisNoRecoPVTracks_,*hisNoRealPVTracks_);
  PVefficiencyVsNoTracks_->SetNameTitle("PVefficiencyVsNoTracks_", "Primary Vertex Finding Efficiency; True No. Tracks; Efficiency");
    
  PVefficiencyVspT_ = inputDir.make<TEfficiency>(*hisRecoPVPt_,*hisRealPVPt_);
  PVefficiencyVspT_->SetNameTitle("PVefficiencyVspT_", "Primary Vertex Finding Efficiency; True pT; Efficiency");

  PVefficiencyVsGoodEtaTracks_ = inputDir.make<TEfficiency>(*hisNoRecoTracksGoodEta_,*hisNoRealTracksGoodEta_);
  PVefficiencyVsGoodEtaTracks_->SetNameTitle("PVefficiencyVsGoodEtaTracks_","Primary Vertex Finding Efficiency; Tracks good eta; Efficiency");
    
  PVefficiencyVsBadEtaTracks_ = inputDir.make<TEfficiency>(*hisNoRecoTracksBadEta_,*hisNoRealTracksBadEta_);
  PVefficiencyVsBadEtaTracks_->SetNameTitle("PVefficiencyVsBadEtaTracks_","Primary Vertex Finding Efficiency; Tracks bad eta; Efficiency");

  tdrPVefficiencyVsTrueZ0_ = inputDir.make<TEfficiency>(*hisTDRPrimaryVertexVsTrueZ0_,*hisPrimaryVertexTrueZ0_);
  tdrPVefficiencyVsTrueZ0_->SetNameTitle("tdrPVefficiencyVsTrueZ0_", "Primary Vertex Finding Efficiency (Technical Proposal Algo); true z_{0}; Efficiency");
    
  PVefficiencyVsNoTracksGoodPercentage50_ = inputDir.make<TEfficiency> (*hisNoRecoTracksPerc50_,*hisNoRealTracksPerc50_);
  PVefficiencyVsNoTracksGoodPercentage50_->SetNameTitle("PVefficiencyVsNoTracksGoodPercentage50_", "PV Efficiency; No. tracks; Efficiency");
    
  PVefficiencyVsNoTracksGoodPercentage90_ = inputDir.make<TEfficiency> (*hisNoRecoTracksPerc90_,*hisNoRealTracksPerc90_);
  PVefficiencyVsNoTracksGoodPercentage90_->SetNameTitle("PVefficiencyVsNoTracksGoodPercentage90_", "PV Efficiency; No. tracks; Efficiency");
    

  PVefficiencyz0Perc50 = inputDir.make<TEfficiency> (*hisz0RecoPerc50_, *hisz0RealPerc50_);
  PVefficiencyz0Perc50->SetNameTitle("PVefficiencyz0Perc50","PV Efficiency; z0; Efficiency");
    
  PVefficiencyz0Perc90 = inputDir.make<TEfficiency> (*hisz0RecoPerc90_, *hisz0RealPerc90_);
  PVefficiencyz0Perc90->SetNameTitle("PVefficiencyz0Perc90","PV Efficiency; z0; Efficiency");
                                       
  PVefficiencypTPerc50 = inputDir.make<TEfficiency> (*hisRecopTPerc50_, *hisRealpTPerc50_);
  PVefficiencypTPerc50->SetNameTitle("PVefficiencypTPerc50","PV Efficiency; pT; Efficiency");
                                       
  PVefficiencypTPerc90 = inputDir.make<TEfficiency> (*hisRecopTPerc90_, *hisRealpTPerc90_);
  PVefficiencypTPerc90->SetNameTitle("PVefficiencypTPerc90","PV Efficiency; pT; Efficiency");
    
  PVefficiencyNoLayers = inputDir.make<TEfficiency> (*numLayersRecoTracks_, *numLayersRealTracks_);
  PVefficiencyNoLayers->SetNameTitle("PVefficiencyNoLayers","PV Efficiency; No. Layers; Efficiency");
  
  PVefficiencyNoPSLayers = inputDir.make<TEfficiency> (*numPSLayersRecoTracks_, *numPSLayersRealTracks_);
  PVefficiencyNoPSLayers->SetNameTitle("PVefficiencyNoPSLayers","PV Efficiency; No. PS Layers; Efficiency");
    
  PVefficiencyLayerID = inputDir.make<TEfficiency> (*recoTracksLayersId_, *realTracksLayersId_);
  PVefficiencyLayerID->SetNameTitle("PVefficiencyLayerID","PV Efficiency; First Layer Hit; Efficiency");
    
    
    
    

  cout << "==================== VERTEX RECONSTRUCTION ======================" << endl;
  
  if(settings_->vx_algoId()==0) cout << "GAP ALGORITHM with vertex gap = "<< settings_->vx_distance() << " cm "<< endl;
  if(settings_->vx_algoId()==1) cout << "SIMPLE MERGE CLUSTERING ALGORITHM with mininal vertex distance = "<< settings_->vx_distance() << " cm "<< endl;
  if(settings_->vx_algoId()==2) cout << "DBSCAN ALGORITHM with vertex gap = "<< settings_->vx_distance() << " cm "<< endl;
  if(settings_->vx_algoId()==3) cout << "PVR ALGORITHM with vertex gap = "<< settings_->vx_distance() << " cm "<< endl;
  if(settings_->vx_algoId()==4) cout << "AVR ALGORITHM with vertex gap = "<< settings_->vx_distance() << " cm "<< endl;

  cout << "Average no. Reconstructed Vertices: "<< hisNoRecoVertices_->GetMean() << "(" << hisNoRecoVertices_->GetMean()*100./(hisNoPileUpVertices_->GetMean()+1.)<<"%)"<< endl;
  cout << "Average ratio of matched tracks in primary vertex "<< hisRatioMatchedTracksInPV_->GetMean()*100 <<"%, Technical Proposal Algo : "<< hisRatioMatchedTracksInTDRPV_->GetMean()*100 << " % "<< endl;
  cout << "Averate ratio of fake tracks in primary vertex "<< hisFakeTracksRateInPV_->GetMean()*100 << "%, Technical Proposal Algo : "<< hisFakeTracksRateInTDRPV_->GetMean()*100 << " % "<< endl;
  cout << "Average PV z0 resolution "<< hisRecoVertexZ0Resolution_->GetMean() << " cm , Technical Proposal Algo : "<< hisTDRVertexZ0Resolution_->GetMean() << " cm "<< endl;

  float recoPVeff = double(hisRecoPrimaryVertexVsTrueZ0_->GetEntries())/double(hisPrimaryVertexTrueZ0_->GetEntries());
  float numRecoPV = double(hisRecoPrimaryVertexVsTrueZ0_->GetEntries());
  float numPVs = double(hisPrimaryVertexTrueZ0_->GetEntries());
  float numTDRPV = double(hisTDRPrimaryVertexVsTrueZ0_->GetEntries());

  float recoPVeff_err = sqrt((numRecoPV+1)*(numRecoPV+2)/((numPVs+2)*(numPVs+3)) - (numRecoPV+1)*(numRecoPV+1)/((numPVs+2)*(numPVs+2)) );

  float tdrPVeff = double(hisTDRPrimaryVertexVsTrueZ0_->GetEntries())/double(hisPrimaryVertexTrueZ0_->GetEntries());
  float tdrPVeff_err = sqrt((numTDRPV+1)*(numTDRPV+2)/((numPVs+2)*(numPVs+3)) - (numTDRPV+1)*(numTDRPV+1)/((numPVs+2)*(numPVs+2)) );

  cout << "PrimaryVertex Finding Efficiency = "<< recoPVeff << " +/- "<<recoPVeff_err << " Technical Proposal Algo "<< tdrPVeff << " +/- "<<tdrPVeff_err << endl;

  
  cout << "================= L1TrkMET Trigger =============================" << endl;
  for (unsigned int i = 0; i < 4; ++i){
    float genmet = 25. + i*25.;
    float met_steps = (genmet-10.)/10;

    cout << "********** GenMET > " << genmet << " GeV << ************" << endl;
    ostringstream name;
    name << "GenMET" <<genmet<<".txt";
    float sigEvents = double(noSignalEvents[i]);
    float bkgEvents = double(noBackgroundEvents[i]);

    hisMETevents_[i]->Fill(0.5, sigEvents);
    hisMETevents_[i]->Fill(21.5, bkgEvents);
    hisMETevents_[i]->GetXaxis()->SetBinLabel(1,"Signal Events");
    hisMETevents_[i]->GetXaxis()->SetBinLabel(22,"Bkg Events");

    for (unsigned int j = 0; j < 10; ++j){
      float cutmet = 10. + j*met_steps;
      cout << "... L1TrkMet > "<<cutmet<<" GeV ..."<< endl;
      float efficiency = 0.;
      float rejection = 0.;
      float efficiency_tdr = 0.;
      float rejection_tdr =  0.;
      float sigma_eff = 0., sigma_eff_tdr = 0., sigma_rej = 0., sigma_rej_tdr = 0.; 

      float recoEvents = double(noRecoSignalEvents[i][j]);
      float tdrEvents = double(noTDRSignalEvents[i][j]);

      if(noSignalEvents[i] > 0){
        efficiency = double(noRecoSignalEvents[i][j])/double(noSignalEvents[i]);
        efficiency_tdr = double(noTDRSignalEvents[i][j])/double(noSignalEventsTDR[i]);
     
        sigma_eff = sqrt( (recoEvents+1)*(recoEvents+2)/((sigEvents+2)*(sigEvents+3)) - (recoEvents+1)*(recoEvents+1)/((sigEvents+2)*(sigEvents+2)) );
        sigma_eff_tdr = sqrt( (tdrEvents+1)*(tdrEvents+2)/((sigEvents+2)*(sigEvents+3)) - (tdrEvents+1)*(tdrEvents+1)/((sigEvents+2)*(sigEvents+2)) );
      }

      float recoBkgEvents = double(noRecoBackgroundEvents[i][j]);
      float tdrBkgEvents = double(noTDRBackgroundEvents[i][j]);

      hisMETevents_[i]->Fill(1.5+j, recoEvents);
      hisMETevents_[i]->Fill(11.5+j, tdrEvents);
      hisMETevents_[i]->Fill(22.5+j, recoBkgEvents);
      hisMETevents_[i]->Fill(32.5+j, tdrBkgEvents);

      ostringstream label ;
      label << "RecoSignals"<<cutmet;
      hisMETevents_[i]->GetXaxis()->SetBinLabel(j+2,label.str().c_str());
      label.clear();
      label.str("");
      label << "TPSignals"<<cutmet;
      hisMETevents_[i]->GetXaxis()->SetBinLabel(j+12,label.str().c_str());
      label.clear();
      label.str("");
      label << "RecoBkg"<<cutmet;
      hisMETevents_[i]->GetXaxis()->SetBinLabel(j+23,label.str().c_str());
      label.clear();
      label.str("");
      label << "TPBkg"<<cutmet;
      hisMETevents_[i]->GetXaxis()->SetBinLabel(j+33,label.str().c_str());


      if(noBackgroundEvents[i] > 0){
        rejection = double(noRecoBackgroundEvents[i][j])/double(noBackgroundEvents[i]);
        rejection_tdr = double(noTDRBackgroundEvents[i][j])/double(noBackgroundEventsTDR[i]);

        sigma_rej = sqrt( (recoBkgEvents+1)*(recoBkgEvents+2)/((bkgEvents+2)*(bkgEvents+3)) - (recoBkgEvents+1)*(recoBkgEvents+1)/((bkgEvents+2)*(bkgEvents+2)) );
        sigma_rej_tdr = sqrt( (tdrBkgEvents+1)*(tdrBkgEvents+2)/((bkgEvents+2)*(bkgEvents+3)) - (tdrBkgEvents+1)*(tdrBkgEvents+1)/((bkgEvents+2)*(bkgEvents+2)) );
               
      }

      

      cout << "Signal Efficiency : "<< efficiency << " +/- "<< sigma_eff << " ("<< noRecoSignalEvents[i][j]<<"/"<<noSignalEvents[i]<<") Bkg Rejection Efficiency : "<< rejection << " +/- " << sigma_rej << "("<< noRecoBackgroundEvents[i][j]<<"/"<<noBackgroundEvents[i]<< ")"<<endl;
      cout << "TP Algo: Signal Efficiency : "<< efficiency_tdr << " +/- "<< sigma_eff_tdr <<" ("<< noTDRSignalEvents[i][j]<<"/"<<noSignalEvents[i]<<") Bkg Rejection Efficiency : "<< rejection_tdr <<" +/- " << sigma_rej_tdr << " ("<< noTDRBackgroundEvents[i][j]<<"/"<<noBackgroundEvents[i]<< ")"<<endl;

      grMET_[i]->SetPoint(j, efficiency, rejection);
      grMET_[i]->SetPointError(j, sigma_eff, sigma_rej);
    
      grMET_tdr_[i]->SetPoint(j, efficiency_tdr, rejection_tdr);
      grMET_tdr_[i]->SetPointError(j, sigma_eff_tdr, sigma_rej_tdr);

    }
  }
}
