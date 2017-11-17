#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

//#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TRandom.h"

#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/DataCorrection.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/DeadModuleDB.h"

#include <iostream>

using namespace std;

//=== Store useful info about this stub.

Stub::Stub(TTStubRef ttStubRef, unsigned int index_in_vStubs, const Settings* settings, 
           const TrackerGeometry*  trackerGeometry, const TrackerTopology*  trackerTopology) : 
  TTStubRef(ttStubRef), 
  settings_(settings), 
  index_in_vStubs_(index_in_vStubs), 
  digitalStub_(settings),
  digitizedForGPinput_(false), // notes that stub has not yet been digitized for GP input.
  digitizedForHTinput_(false), // notes that stub has not yet been digitized for HT input.
  digitizedForSFinput_(false)  // notes that stub has not yet been digitized for seed filter input.
{
  // Get coordinates of stub.
  const TTStub<Ref_Phase2TrackerDigi_> *ttStubP = ttStubRef.get(); 

  DetId geoDetId = ttStubRef->getDetId();
  for (auto gd=trackerGeometry->dets().begin(); gd != trackerGeometry->dets().end(); gd++) 
  {
    DetId detid = (*gd)->geographicalId();
    if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue; // only run on OT
    if(!trackerTopology->isLower(detid) ) continue; // loop on the stacks: choose the lower arbitrarily
    DetId stackDetid = trackerTopology->stack(detid); // Stub module detid

    if ( ttStubRef->getDetId() != stackDetid ) continue;
    geoDetId = detid;
    break;
  }
  const GeomDetUnit* det0 = trackerGeometry->idToDetUnit( geoDetId );
  // To get other module, can do this
  // const GeomDetUnit* det1 = trackerGeometry->idToDetUnit( trackerTopology->partnerDetId( geoDetId ) );

  const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
  const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );
  MeasurementPoint measurementPoint = ttStubRef->getClusterRef(0)->findAverageLocalCoordinatesCentered();
  LocalPoint clustlp   = topol->localPosition(measurementPoint);
  GlobalPoint pos  =  theGeomDet->surface().toGlobal(clustlp);

  phi_ = pos.phi();
  r_   = pos.perp();
  z_   = pos.z();

  if (r_ < settings_->trackerInnerRadius() || r_ > settings_->trackerOuterRadius() || fabs(z_) > settings_->trackerHalfLength()) {
    throw cms::Exception("Stub: Stub found outside assumed tracker volume. Please update tracker dimensions specified in Settingsm.h!")<<" r="<<r_<<" z="<<z_<<" "<<ttStubRef->getDetId().subdetId()<<endl;
  }

  // Set info about the module this stub is in
  this->setModuleInfo(trackerGeometry, trackerTopology, geoDetId);
  // Uncertainty in stub coordinates due to strip or pixel length in r-z.
  // float sensorSpacing = sqrt( (moduleMaxR_ - moduleMinR_) * (moduleMaxR_ - moduleMinR_) + (
    // moduleMaxZ_ - moduleMinZ_) * (moduleMaxZ_ - moduleMinZ_) );
  if (barrel_) {
    rErr_ = 0.;
    zErr_ = 0.5*stripLength_; 
    
    if((fabs(z_) >= 11.8542 and layerId_ == 1) or (fabs(z_) >= 26.7164 and layerId_ == 2) or (fabs(z_) >= 36.2523 and layerId_ == 3)){
      rErr_ = 0.5*(moduleMaxR_ - moduleMinR_);
      zErr_ = 0.5*(moduleMaxZ_ - moduleMinZ_);
    }
    // cout << "barrel ";
  } else {
    rErr_ = 0.5*stripLength_; 
    zErr_ = 0.;
    // cout << "endcap ";
  }
  // cout << "r "<< r_ << " rErr_ "<< rErr_ << " zErr_ "<< zErr_ << endl;

  // Get the coordinates of the two clusters that make up this stub, measured in units of strip pitch, and measured
  // in the local frame of the sensor. They have a granularity  of 0.5*pitch.
  for (unsigned int iClus = 0; iClus <= 1; iClus++) { // Loop over two clusters in stub.  
    localU_cluster_[iClus] = ttStubP->getClusterRef(iClus)->findAverageLocalCoordinatesCentered().x();
    localV_cluster_[iClus] = ttStubP->getClusterRef(iClus)->findAverageLocalCoordinatesCentered().y();
  }

  // Get location of stub in module in units of strip number (or pixel number along finest granularity axis).
  // Range from 0 to (nStrips - 1) inclusive.
  // N.B. Since iphi is integer, this degrades the granularity by a factor 2. This seems silly, but track fit wants it.
  iphi_ = localU_cluster_[0]; // granularity 1*strip (unclear why we want to degrade it ...)

  // Get stub bend that is available in front-end electronics, where bend is displacement between 
  // two hits in stubs in units of strip pitch.
  bendInFrontend_ = ttStubRef->getTriggerBend();
  bool isEndcap = GeomDetEnumerators::isEndcap( trackerGeometry->idToDetUnit( geoDetId )->subDetector() );
  if (isEndcap && pos.z() > 0) bendInFrontend_ *= -1;

  // EJC Bend in barrel seems to be flipped
  // Below can be reduced to one line, but suprised that this has to change over entire
  // barrel, and not just for e.g. plus/minus z for positive/negative bend
  // so left as two lines to leave clues for any future bugs...
  if ( barrel_ && pos.z() >= 0 ) bendInFrontend_ = (-1)*bendInFrontend_;
  if ( barrel_ && pos.z() < 0 ) bendInFrontend_ = (-1)*bendInFrontend_;    


  // Get stub bend that is available in off-detector electronics, allowing for degredation of 
  // bend resolution due to bit encoding by FE chip if required.
  bool rejectStub = false;          // indicates if bend is outside window assumed in DataCorrection.h
  numMergedBend_ = 1;               // Number of bend values merged into single degraded one.
  if (settings->bendResReduced()) {
    float degradedBend;       // degraded bend
    this->degradeResolution(bendInFrontend_, geoDetId,
          degradedBend, rejectStub, numMergedBend_); // sets value of last 3 arguments.
    bend_ = degradedBend;
  } else {
    bend_ = bendInFrontend_;
  }

  // Estimate track Pt and phi0 based on stub bend info, and angle in r-phi projection of stub direction to sensor plane.
  // float sensorSpacing = barrel_ ? (moduleMaxR_ - moduleMinR_) : (moduleMaxZ_ - moduleMinZ_);
  // EJC Above not true for tilted modules
  float sensorSpacing = sqrt( (moduleMaxR_ - moduleMinR_) * (moduleMaxR_ - moduleMinR_) + (moduleMaxZ_ - moduleMinZ_) * (moduleMaxZ_ - moduleMinZ_) );
  
  pitchOverSep_ = stripPitch_/sensorSpacing;
  // IRT - use stub (r,z) instead of module (r,z). Logically correct but has negligable effect on results.
  // dphiOverBend_ = barrel_  ?  pitchOverSep_  :  pitchOverSep_*fabs(z_)/r_;
  // EJC Above not true for tilted modules
  double alpha = this->moduleTilt();
  dphiOverBendCorrection_ = fabs( cos( fabs(this->theta()) - alpha ) / sin( this->theta() ) );
  dphiOverBend_ = pitchOverSep_ * dphiOverBendCorrection_;    

  dphi_ = bend_ * dphiOverBend();

  // Fill frontendPass_ flag, indicating if frontend readout electronics will output this stub.
  this->setFrontend(rejectStub); 

  // Calculate bin range along q/Pt axis of r-phi Hough transform array consistent with bend of this stub.
  this->calcQoverPtrange();

  // Initialize class used to produce digital version of stub, with original stub parameters pre-digitization.
  digitalStub_.init(phi_, r_, z_, dphi(), this->rhoParameter(), min_qOverPt_bin_, max_qOverPt_bin_, layerId_, this->layerIdReduced(), bend_, stripPitch_, sensorSpacing, rErr_, zErr_);
}

//=== Calculate bin range along q/Pt axis of r-phi Hough transform array consistent with bend of this stub.

void Stub::calcQoverPtrange() {
  // First determine bin range along q/Pt axis of HT array 
  const int nbinsPt = (int) settings_->houghNbinsPt(); // Use "int" as nasty things happen if multiply "int" and "unsigned int".
  const int min_array_bin = 0;
  const int max_array_bin = nbinsPt - 1;  
  // Now calculate range of q/Pt bins allowed by bend filter.
  float qOverPtMin = this->qOverPtOverBend() * (this->bend() - this->bendRes());
  float qOverPtMax = this->qOverPtOverBend() * (this->bend() + this->bendRes());
  const float houghMaxInvPt = 1./settings_->houghMinPt();
  const float qOverPtBinSize = (2. * houghMaxInvPt)/settings_->houghNbinsPt();
  // Convert to bin number along q/Pt axis of HT array.
  // N.B. The terms involving "0.5" here have the effect that the cell is accepted if the q/Pt at its centre is
  // consistent with the stub bend. This gives the same behaviour for the "daisy chain" firmware, which uses
  // this bin range, and for the systolic/2-c-bin firmwares which instead use the calculation in HTcell::bendFilter().
  // If you choose to remove the "0.5" terms here, which loosens the bend filter cut, then I recommend that you 
  // tighten up the value of the "BendResolution" config parameter by about 0.05 to compensate.
  // Decision to remove them taken in softare & GP firmware on 9th August 2016.
  //int min_bin = std::floor(  0.5 + (qOverPtMin + houghMaxInvPt)/qOverPtBinSize);
  //int max_bin = std::floor( -0.5 + (qOverPtMax + houghMaxInvPt)/qOverPtBinSize);
  int min_bin = std::floor((qOverPtMin + houghMaxInvPt)/qOverPtBinSize);
  int max_bin = std::floor((qOverPtMax + houghMaxInvPt)/qOverPtBinSize);
  // Limit it to range of HT array.
  min_bin = max(min_bin, min_array_bin);
  max_bin = min(max_bin, max_array_bin);
  // If min_bin > max_bin at this stage, it means that the Pt estimated from the bend is below the range we wish to find tracks in.
  // Keep min_bin > max_bin, so such stubs can be identified, but set both variables to values inside the allowed range.
  if (min_bin > max_bin) {
    min_bin = max_array_bin;
    max_bin = min_array_bin;
    //if (frontendPass_) throw cms::Exception("Stub: m bin calculation found low Pt stub not killed by FE electronics cuts")<<qOverPtMin<<" "<<qOverPtMax<<endl;
  }
  min_qOverPt_bin_ = (unsigned int) min_bin;
  max_qOverPt_bin_ = (unsigned int) max_bin;
}

//=== Digitize stub for input to Geographic Processor, with digitized phi coord. measured relative to closest phi sector.
//=== (This approximation is valid if their are an integer number of digitisation bins inside each phi octant).
//=== However, you should also call digitizeForHTinput() before accessing digitized stub data, even if you only care about that going into GP! Otherwise, you will not identify stubs assigned to more than one octant.

void Stub::digitizeForGPinput(unsigned int iPhiSec) {
  if (settings_->enableDigitize()) {

    // Save CPU by not redoing digitization if stub was already digitized for this phi sector.
    if ( ! (digitizedForGPinput_ && digitalStub_.iGetOctant(iPhiSec) == digitalStub_.iDigi_Octant()) ) {

      // Digitize
      digitalStub_.makeGPinput(iPhiSec);

      // Replace stub coordinates with those degraded by digitization process.
      phi_  = digitalStub_.phi();
      r_    = digitalStub_.r();
      z_    = digitalStub_.z();
      bend_ = digitalStub_.bend();

      // If the Stub class contains any data members that are not input to the GP, but are derived from variables that
      // are, then be sure to update these here too, unless Stub.h uses the check*() functions to declare them invalid.

      dphiOverBend_ = barrel_  ?  pitchOverSep_  :  pitchOverSep_*fabs(z_)/r_;

      // Note that stub has been digitized for GP input
      digitizedForGPinput_ = true;
    }
    digitizedForHTinput_ = false;
  }
}

//=== Digitize stub for input to Hough transform, with digitized phi coord. measured relative to specified phi sector.

void Stub::digitizeForHTinput(unsigned int iPhiSec) {

  if (settings_->enableDigitize()) {

    // Save CPU by not redoing digitization if stub was already digitized for this phi sector.
    if ( ! (digitizedForHTinput_ && iPhiSec == digitalStub_.iDigi_PhiSec()) ) {

      // Call digitization for GP in case not already done. (Needed for variables that are common to GP & HT).
      this->digitizeForGPinput(iPhiSec);

      // Digitize
      digitalStub_.makeHTinput(iPhiSec);

      // Since GP and HT use same digitisation in r and z, don't bother updating their values.
      // (Actually, the phi digitisation boundaries also match, except for systolic array, so could skip updating phi too).

      // Replace stub coordinates and bend with those degraded by digitization process. (Don't bother with r & z, as already done by GP digitisation).
      phi_  = digitalStub_.phi();

      if (settings_->firmwareType() >= 1 && settings_->firmwareType() < 99) {
  // If using daisy chain firmware, then recalculate bin range along q/Pt axis of r-phi Hough transform array 
  // consistent with bend of this stub, since it depends on r & z which have now been digitized.
  // (This recalculation should really be done in DigitalStub::makeHTinput(), but too lazy to move it there ...).
  this->calcQoverPtrange();

      } else {

  // If using Thomas/Systolic firmware, then variables dphi & rho are digitized.
  dphi_ = digitalStub_.dphi();
  float rho  = digitalStub_.rho();
  this->setRhoParameter(rho);
      }

      // If the Stub class contains any data members that are not input to the HT, but are derived from variables that
      // are, then be sure to update these here too, unless Stub.h uses the check*() functions to declare them invalid. 
      // - currently none.

      // Note that stub has been digitized.
      digitizedForHTinput_ = true;
    }
  }
}

//=== Digitize stub for input to r-z Seed Filter.

void Stub::digitizeForSFinput() {
  if (settings_->enableDigitize()) {

    if ( ! digitizedForSFinput_) {
      // Digitize variables specific to seed filter if not already done.
      digitalStub_.makeSFinput();

      // Replace stub (r,z) uncertainties, estimated from half-pixel/strip-length, by those degraded by the digitization process. 
      rErr_ = digitalStub_.rErr();
      zErr_ = digitalStub_.zErr();
      // Must also replace stub r coordinate, as seed filter works with digitized r instead of digitized rT.
      r_    = digitalStub_.r();

      digitizedForSFinput_ = true;
    }
  }
}

//=== Digitize stub for input to r-z Seed Filter.

void Stub::digitizeForDRinput(unsigned int stubId) {
  if (settings_->enableDigitize()) {
    
      // Digitize variables specific to seed filter if not already done.
      digitalStub_.makeDRinput(stubId);
      // digitizedForDRinput_ = true;
    
  }
}


//===  Restore stub to pre-digitized state. i.e. Undo what function digitize() did.

void Stub::reset_digitize() {
  if (settings_->enableDigitize()) {
    // Save CPU by not undoing digitization if stub was not already digitized.
    if (digitizedForGPinput_ || digitizedForHTinput_) {

      // Replace stub coordinates and bend with original coordinates stored prior to any digitization.
      phi_  = digitalStub_.orig_phi();
      r_    = digitalStub_.orig_r();
      z_    = digitalStub_.orig_z();
      bend_ = digitalStub_.orig_bend();

      // Also restore original uncertainties in stub coordinates (estimated from strip or pixel half-length).
      rErr_ = digitalStub_.orig_rErr();
      zErr_ = digitalStub_.orig_zErr();

      // Variables dphi & rho are not used with daisy-chain firmware.
      if (settings_->firmwareType() == 0 || settings_->firmwareType() == 99)  {
        dphi_ = digitalStub_.orig_dphi();
        float rho  = digitalStub_.orig_rho();
        this->setRhoParameter(rho);
      }

      // Note that stub is (no longer) digitized.
      digitizedForGPinput_ = false;
      digitizedForHTinput_ = false;
      digitizedForSFinput_ = false;

      // If the Stub class contains any data members that are not input to the GP or HT, but are derived from 
      // variables that are, then be sure to update these here too.

      dphiOverBend_ = barrel_  ?  pitchOverSep_  :  pitchOverSep_*fabs(z_)/r_;

      if (settings_->firmwareType() == 0 || settings_->firmwareType() == 99) {
  // Recalculate bin range along q/Pt axis of r-phi Hough transform array consistent with bend of this stub,
  // since it depends on dphi which is no longer digitized. Not needed with daisy-chain firmware, since this range
  // is transmitted to HT hardware along optical link.
  this->calcQoverPtrange();
      }
    }
  }
}

//=== Degrade assumed stub bend resolution.
//=== Also return boolean indicating if stub bend was outside assumed window, so stub should be rejected
//=== and return an integer indicating how many values of bend are merged into this single one.

void Stub::degradeResolution(float bend, const DetId& detId,
           float& degradedBend, bool& reject, unsigned int& num) const {
  if (barrel_) {
    // unsigned int layer = detId.iLayer();
    DataCorrection::ConvertBarrelBend( bend, layerId_,
               degradedBend, reject, num);
  } else {
    // unsigned int ring = detId.iRing();
    DataCorrection::ConvertEndcapBend( bend, endcapRing_,
               degradedBend, reject, num);
  }
}


//=== Set flag indicating if stub will be output by front-end readout electronics 
//=== (where we can reconfigure the stub window size and rapidity cut).
//=== Argument indicates if stub bend was outside window size encoded in DataCorrection.h
//=== Note that this should run on quantities as available inside front-end chip, which are not
//=== degraded by loss of bits or digitisation.

void Stub::setFrontend(bool rejectStub) {
  frontendPass_ = true; // Did stub pass cuts applied in front-end chip
  stubFailedDataCorrWindow_ = false; // Did it only fail cuts corresponding to windows encoded in DataCorrection.h?
  // Don't use stubs at large eta, since it is impossible to form L1 tracks from them, so they only contribute to combinatorics.
  if ( fabs(this->eta()) > settings_->maxStubEta() ) frontendPass_ = false;
  // Don't use stubs whose Pt is significantly below the Pt cut used in the L1 tracking, allowing for uncertainty in q/Pt due to stub bend resolution.
  if (settings_->killLowPtStubs()) {
    const float qOverPtCut = 1./settings_->houghMinPt();
    // Apply this cut in the front-end electronics.
    if (fabs(this->bendInFrontend()) - this->bendResInFrontend() > qOverPtCut/this->qOverPtOverBend()) frontendPass_ = false;
    // Reapply the same cut using the degraded bend information available in the off-detector electronics.
    // The reason is  that the bend degredation can move the Pt below the Pt cut, making the stub useless to the off-detector electronics.
    if (fabs(this->bend())           - this->bendRes()           > qOverPtCut/this->qOverPtOverBend()) frontendPass_ = false;
  } 
  // Don't use stubs whose bend is outside the window encoded into DataCorrection.h
  if (rejectStub) {
    if (frontendPass_) stubFailedDataCorrWindow_ = true;
    frontendPass_ = false;
  }

  // Kill stubs in tracker regions declared dead by the user.
  if (settings_->deadSimulateFrac() > 0.) { // Is option to emulate dead modules enabled?
    const DeadModuleDB dead;
    if (dead.killStub(this)) {
      static TRandom randomGenerator;
      if (randomGenerator.Rndm() < settings_->deadSimulateFrac()) frontendPass_ = false;
    }
  }
}

//=== Note which tracking particle(s), if any, produced this stub.
//=== The 1st argument is a map relating TrackingParticles to TP.

void Stub::fillTruth(const map<edm::Ptr< TrackingParticle >, const TP* >& translateTP, edm::Handle<TTStubAssMap> mcTruthTTStubHandle, edm::Handle<TTClusterAssMap> mcTruthTTClusterHandle){

  TTStubRef ttStubRef(*this); // Cast to base class

  //--- Fill assocTP_ info. If both clusters in this stub were produced by the same single tracking particle, find out which one it was.

  bool genuine =  mcTruthTTStubHandle->isGenuine(ttStubRef); // Same TP contributed to both clusters?
  assocTP_ = nullptr;

  // Require same TP contributed to both clusters.
  if ( genuine ) {
    edm::Ptr< TrackingParticle > tpPtr = mcTruthTTStubHandle->findTrackingParticlePtr(ttStubRef);
    if (translateTP.find(tpPtr) != translateTP.end()) {
      assocTP_ = translateTP.at(tpPtr);
      // N.B. Since not all tracking particles are stored in InputData::vTPs_, sometimes no match will be found.
    }
  }

  // Fill assocTPs_ info.

  if (settings_->stubMatchStrict()) {

    // We consider only stubs in which this TP contributed to both clusters.
    if (assocTP_ != nullptr) assocTPs_.insert(assocTP_);

  } else {

    // We consider stubs in which this TP contributed to either cluster.

    for (unsigned int iClus = 0; iClus <= 1; iClus++) { // Loop over both clusters that make up stub.
       const TTClusterRef& ttClusterRef = ttStubRef->getClusterRef(iClus);

      // Now identify all TP's contributing to either cluster in stub.
      vector< edm::Ptr< TrackingParticle > > vecTpPtr = mcTruthTTClusterHandle->findTrackingParticlePtrs(ttClusterRef);

      for (edm::Ptr< TrackingParticle> tpPtr : vecTpPtr) {
  if (translateTP.find(tpPtr) != translateTP.end()) {
    assocTPs_.insert( translateTP.at(tpPtr) );
    // N.B. Since not all tracking particles are stored in InputData::vTPs_, sometimes no match will be found.
  }
      }
    }
  }

  //--- Also note which tracking particles produced the two clusters that make up the stub

  for (unsigned int iClus = 0; iClus <= 1; iClus++) { // Loop over both clusters that make up stub.
    const TTClusterRef& ttClusterRef = ttStubRef->getClusterRef(iClus);

    bool genuineCluster =  mcTruthTTClusterHandle->isGenuine(ttClusterRef); // Only 1 TP made cluster?
    assocTPofCluster_[iClus] = nullptr;

    // Only consider clusters produced by just one TP.
    if ( genuineCluster ) {
      edm::Ptr< TrackingParticle > tpPtr = mcTruthTTClusterHandle->findTrackingParticlePtr(ttClusterRef);

      if (translateTP.find(tpPtr) != translateTP.end()) {
  assocTPofCluster_[iClus] = translateTP.at(tpPtr);
  // N.B. Since not all tracking particles are stored in InputData::vTPs_, sometimes no match will be found.
      }
    }
  }

  // Sanity check - is truth info of stub consistent with that of its clusters?
  // Commented this out, as it throws errors for unknown reason with iErr=1. Apparently, "genuine" stubs can be composed of two clusters that are
  // not "genuine", providing that one of the TP that contributed to each cluster was the same.
  /*
  unsigned int iErr = 0;
  if (this->genuine()) { // Stub matches truth particle
    if ( ! ( this->genuineCluster()[0] && (this->assocTPofCluster()[0] == this->assocTPofCluster()[1]) ) ) iErr = 1;
  } else {
    if ( ! ( ! this->genuineCluster()[0] || (this->assocTPofCluster()[0] != this->assocTPofCluster()[1]) )  ) iErr = 2;
  }
  if (iErr > 0) {
    cout<<" DEBUGA "<<(this->assocTP() == nullptr)<<endl;
    cout<<" DEBUGB "<<(this->assocTPofCluster()[0] == nullptr)<<" "<<(this->assocTPofCluster()[1] == nullptr)<<endl;
    cout<<" DEBUGC "<<this->genuineCluster()[0]<<" "<<this->genuineCluster()[1]<<endl;
    if (this->assocTPofCluster()[0] != nullptr) cout<<" DEBUGD "<<this->assocTPofCluster()[0]->index()<<endl;
    if (this->assocTPofCluster()[1] != nullptr) cout<<" DEBUGE "<<this->assocTPofCluster()[1]->index()<<endl;
    //    throw cms::Exception("Stub: Truth info of stub & its clusters inconsistent!")<<iErr<<endl;
  }
  */
}

//=== Estimated phi angle at which track crosses a given radius rad, based on stub bend info. Also estimate uncertainty on this angle due to endcap 2S module strip length.
//=== N.B. This is identical to Stub::beta() if rad=0.

pair <float, float> Stub::trkPhiAtR(float rad) const { 
  float rStubMax = r_ + rErr_; // Uncertainty in radial stub coordinate due to strip length.
  float rStubMin = r_ - rErr_;
  float trkPhi1 = (phi_ + dphi()*(1. - rad/rStubMin));
  float trkPhi2 = (phi_ + dphi()*(1. - rad/rStubMax));
  float trkPhi    = 0.5*    (trkPhi1 + trkPhi2);
  float errTrkPhi = 0.5*fabs(trkPhi1 - trkPhi2); 
  return pair<float, float>(trkPhi, errTrkPhi);
}


//=== Note if stub is a crazy distance from the tracking particle trajectory that produced it.
//=== If so, it was probably produced by a delta ray.

bool Stub::crazyStub() const {

  bool crazy;
  if (assocTP_ == nullptr) {
    crazy = false; // Stub is fake, but this is not crazy. It happens ...
  } else {
    // Stub was produced by TP. Check it lies not too far from TP trajectory.
    crazy = fabs( reco::deltaPhi(phi_, assocTP_->trkPhiAtStub( this )) )  >  settings_->crazyStubCut();
  } 
  return crazy;
}

//=== Get reduced layer ID (in range 1-7), which can be packed into 3 bits so simplifying the firmware).

unsigned int Stub::layerIdReduced() const {
  // Don't bother distinguishing two endcaps, as no track can have stubs in both.
  unsigned int lay = (layerId_ < 20) ? layerId_ : layerId_ - 10; 

  // No genuine track can have stubs in both barrel layer 6 and endcap disk 11 etc., so merge their layer IDs.
  // WARNING: This is tracker geometry dependent, so may need changing in future ...
  if (lay == 6) lay = 11; 
  if (lay == 5) lay = 12; 
  if (lay == 4) lay = 13; 
  if (lay == 3) lay = 15; 
  // At this point, the reduced layer ID can have values of 1, 2, 11, 12, 13, 14, 15. So correct to put in range 1-7.
  if (lay > 10) lay -= 8;

  if (lay < 1 || lay > 7) throw cms::Exception("Stub: Reduced layer ID out of expected range");

  return lay;
}


//=== Set info about the module that this stub is in.

void Stub::setModuleInfo(const TrackerGeometry* trackerGeometry, const TrackerTopology* trackerTopology, const DetId& detId) {

  idDet_ = detId();

  // Get min & max (r,phi,z) coordinates of the centre of the two sensors containing this stub.
  const GeomDetUnit* det0 = trackerGeometry->idToDetUnit( detId );
  const GeomDetUnit* det1 = trackerGeometry->idToDetUnit( trackerTopology->partnerDetId( detId ) );

  float R0 = det0->position().perp();
  float R1 = det1->position().perp();
  float PHI0 = det0->position().phi();
  float PHI1 = det1->position().phi();
  float Z0 = det0->position().z();
  float Z1 = det1->position().z();
  moduleMinR_   = std::min(R0,R1);
  moduleMaxR_   = std::max(R0,R1);
  moduleMinPhi_ = std::min(PHI0,PHI1);
  moduleMaxPhi_ = std::max(PHI0,PHI1);
  moduleMinZ_   = std::min(Z0,Z1);
  moduleMaxZ_   = std::max(Z0,Z1);

  // Note if module is PS or 2S, and whether in barrel or endcap.
  psModule_ = trackerGeometry->getDetectorType( detId ) == TrackerGeometry::ModuleType::Ph2PSP; // From https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/Geometry/TrackerGeometryBuilder/README.md
  barrel_ = detId.subdetId()==StripSubdetector::TOB || detId.subdetId()==StripSubdetector::TIB;

  //  cout<<"DEBUG STUB "<<barrel_<<" "<<psModule_<<"  sep(r,z)=( "<<moduleMaxR_ - moduleMinR_<<" , "<<moduleMaxZ_ - moduleMinZ_<<" )    stub(r,z)=( "<<0.5*(moduleMaxR_ + moduleMinR_) - r_<<" , "<<0.5*(moduleMaxZ_ + moduleMinZ_) - z_<<" )"<<endl;

  // Encode layer ID.
  if (barrel_) {
    layerId_ = trackerTopology->layer( detId ); // barrel layer 1-6 encoded as 1-6
  } else {
    // layerId_ = 10*detId.iSide() + detId.iDisk(); // endcap layer 1-5 encoded as 11-15 (endcap A) or 21-25 (endcapB)
    // EJC This seems to give the same encoding
    layerId_ = 10*trackerTopology->side( detId ) + trackerTopology->tidWheel( detId );
  }

  // Note module ring in endcap
  // endcapRing_ = barrel_  ?  0  :  detId.iRing();
  endcapRing_ = barrel_  ?  0  :  trackerTopology->tidRing( detId );

  // Get sensor strip or pixel pitch using innermost sensor of pair.

  const PixelGeomDetUnit* unit = reinterpret_cast<const PixelGeomDetUnit*>( det0 );
  const PixelTopology& topo = unit->specificTopology();
  const Bounds& bounds = det0->surface().bounds();

  std::pair<float, float> pitch = topo.pitch();
  stripPitch_ = pitch.first; // Strip pitch (or pixel pitch along shortest axis)
  stripLength_ = pitch.second;  //  Strip length (or pixel pitch along longest axis)
  nStrips_ = topo.nrows(); // No. of strips in sensor
  sensorWidth_ = bounds.width(); // Width of sensitive region of sensor (= stripPitch * nStrips).

  outerModuleAtSmallerR_ = false;
  if ( barrel_ && det0->position().perp() > det1->position().perp() ) {
    outerModuleAtSmallerR_ = true;
  }


  sigmaPerp_ = stripPitch_/sqrt(12.); // resolution perpendicular to strip (or to longest pixel axis)
  sigmaPar_  = stripLength_/sqrt(12.); // resolution parallel to strip (or to longest pixel axis)
}
