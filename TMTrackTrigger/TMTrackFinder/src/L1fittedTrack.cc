#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"

// Function for merging two tracks into a single track, used by by KillDupTracks.h for duplicate track removal.

L1fittedTrack L1fittedTrack::mergeTracks(const L1fittedTrack B) const {

  vector<const Stub*> aStubs=this->getStubs(), bStubs=B.getStubs();

  // Since a "set" contains no duplicates, add stubs to a set to eliminate the duplicates.
  set<const Stub*> mStubs;
  mStubs.insert(aStubs.begin(), aStubs.end());
  mStubs.insert(bStubs.begin(), bStubs.end());

  // Now copy the set back to the required vector for output.
  vector<const Stub*> mergedStubs;
  for (const Stub* s: mStubs) {
    mergedStubs.push_back(s);
  }

  // N.B. This defines the HT cell location as that of the first track, meaning that the merged tracks depends
  // on which track is first and which is second. This will make it hard to get identical results from hardware 
  // & software.
  return L1fittedTrack(settings_, l1track3D_, mergedStubs, 
           qOverPt_, d0_, phi0_, z0_, tanLambda_, chi2_, nHelixParam_,
           iPhiSec_, iEtaReg_, accepted_);
}

// Digitize track and degrade helix parameter resolution according to effect of digitisation.

void L1fittedTrack::digitizeTrack(const string& fitterName){
  if (settings_->enableDigitize()) {
    if (! digitizedTrack_ ) {
      digitizedTrack_ = true;
     
      bool consistent = this->consistentHTcell();
      int  mbinhelix  = int(this->getCellLocationRphi().first) - floor(settings_->houghNbinsPt()/2);
      int  cbinhelix  = int(this->getCellLocationRphi().second) - floor(settings_->houghNbinsPhi()/2);
      int  mBinHT     = int(this->getCellLocationRphiHT().first) - floor(settings_->houghNbinsPt()/2);
      int  cBinHT     = int(this->getCellLocationRphiHT().second) - floor(settings_->houghNbinsPhi()/2);

      if(matchedTP_ != nullptr){
        digitalTrack_.init(fitterName,
         iPhiSec_, iEtaReg_, mBinHT, cBinHT, mbinhelix, cbinhelix, 
         qOverPt_, phi0_,tanLambda_, z0_, chi2_, nLayers_, consistent, this->accepted(),
         matchedTP_->pt(), matchedTP_->eta(), 
         matchedTP_->index(), matchedTP_->useForAlgEff(), matchedTP_->useForEff(), matchedTP_->pdgId(), 
         matchedTP_->phi0(), matchedTP_->tanLambda(), matchedTP_->z0(), matchedTP_->qOverPt());
      } else{
        digitalTrack_.init(fitterName,
         iPhiSec_, iEtaReg_, mBinHT, cBinHT, mbinhelix, cbinhelix,       
         qOverPt_, phi0_, tanLambda_, z0_, chi2_, nLayers_, consistent, this->accepted(), 
         0, 0, -1, 0, 0, 0, 0, 0, 0, 0);
      }

      // Digitize track
      digitalTrack_.makeDigitalTrack();

      // Convert digitized track params back to floating point with degraded resolution.
      qOverPt_   = digitalTrack_.qOverPt();
      phi0_      = digitalTrack_.phi0();
      z0_        = digitalTrack_.z0();
      tanLambda_ = digitalTrack_.tanLambda();
      chi2_      = digitalTrack_.chisquared();

    }
  }
}

