#include "TMTrackTrigger/TMTrackFinder/interface/Sector.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"

#include "DataFormats/Math/interface/deltaPhi.h"

using namespace std;

//=== Initialise
 
void Sector::init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg) {
  settings_ = settings;

  // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when assigning stubs to sectors?
  handleStripsPhiSec_ = settings->handleStripsPhiSec();
  handleStripsEtaSec_ = settings->handleStripsEtaSec();

  //===  Characteristics of this eta region.

  iPhiSec_ = iPhiSec;
  iEtaReg_ = iEtaReg;
  // Using lines of specified rapidity drawn from centre of CMS, determine the z coords at which
  // they cross the radius rChosenRofZ_.
  etaMin_ = settings->etaRegions()[iEtaReg];
  etaMax_ = settings->etaRegions()[iEtaReg + 1];
  chosenRofZ_ = settings->chosenRofZ();
  // Get range in z of tracks covered by this sector at chosen radius from beam-line
  zOuterMin_ = chosenRofZ_ / tan( 2. * atan(exp(-etaMin_)) );
  zOuterMax_ = chosenRofZ_ / tan( 2. * atan(exp(-etaMax_)) );
  beamWindowZ_ = settings->beamWindowZ(); // Assumed half-length of beam-spot

  // IRT
  // If rapidity line leaves tracker endcap before reaching r = rChosenOfZ_, try something different.
  /*
  trackerOuterRadius_ = settings->trackerOuterRadius();
  trackerInnerRadius_ = settings->trackerInnerRadius();
  trackerHalfLength_  = settings->trackerHalfLength();
  if (fabs(zOuterMin_) > trackerHalfLength_) {
    float scale = trackerHalfLength_/fabs(zOuterMin_);
    zOuterMin_ *= scale;
    rOuterMin_ *= scale;
  }
  if (fabs(zOuterMax_) > trackerHalfLength_) {
    float scale = trackerHalfLength_/fabs(zOuterMax_);
    zOuterMax_ *= scale;
    rOuterMax_ *= scale;
  }
  */

  //=== Characteristics of this phi region.

  phiCentre_ = 2.*M_PI * (0.5 + float(iPhiSec)) / float(settings->numPhiSectors()) - M_PI; // Centre of sector in phi
  sectorHalfWidth_ = M_PI / float(settings->numPhiSectors()); // Sector half width excluding overlaps.
  chosenRofPhi_     = settings->chosenRofPhi();
  useStubPhi_       = settings->useStubPhi();
  minPt_            = settings->houghMinPt(); // Min Pt covered by  HT array.
  useStubPhiTrk_    = settings->useStubPhiTrk();
  assumedPhiTrkRes_ = settings->assumedPhiTrkRes();
  calcPhiTrkRes_    = settings->calcPhiTrkRes();

  //=== Check if subsectors in eta are being used within each sector.
  numSubSecsEta_    = settings->numSubSecsEta();
  /*
  // If subsectors have equal width in z50, do this.
  float subSecWidth = (zOuterMax_ - zOuterMin_)/float(numSubSecsEta_); 
  for (unsigned int i = 0; i < numSubSecsEta_; i++) {
    zOuterMinSub_.push_back( zOuterMin_ +  i     *subSecWidth);
    zOuterMaxSub_.push_back( zOuterMin_ + (i + 1)*subSecWidth);
  }
  */
  // If subsectors have equal width in rapidity, do this.
  float subSecWidth = (etaMax_ - etaMin_)/float(numSubSecsEta_); 
  for (unsigned int i = 0; i < numSubSecsEta_; i++) {
    float subSecEtaMin = etaMin_ + i * subSecWidth;
    float subSecEtaMax = subSecEtaMin + subSecWidth;
    float subSecZmin = chosenRofZ_ / tan( 2. * atan(exp(-subSecEtaMin)) );
    float subSecZmax = chosenRofZ_ / tan( 2. * atan(exp(-subSecEtaMax)) );
    zOuterMinSub_.push_back( subSecZmin );
    zOuterMaxSub_.push_back( subSecZmax );
  }
}

//=== Check if stub is inside this eta region.

bool Sector::insideEta( const Stub* stub ) const {
  // Lower edge of this eta region defined by line from (r,z) = (0,-beamWindowZ) to (chosenRofZ_, zOuterMin_).
  // Upper edge of this eta region defined by line from (r,z) = (0, beamWindowZ) to (chosenRofZ_, zOuterMax_).

  bool inside = this->insideEtaRange(stub, zOuterMin_, zOuterMax_);
  return inside;
}


//=== Check if stub is within subsectors in eta that sector may be divided into.

vector<bool> Sector::insideEtaSubSecs( const Stub* stub) const {

  if (settings_->enableDigitize() && numSubSecsEta_ == 2) {
    // Use (complicated) digitized firmware emulation
    return subEtaFwCalc(stub->digitalStub().iDigi_Rt(), stub->digitalStub().iDigi_Z());

  } else {
    // Use (simpler) floating point calculation.

    vector<bool> insideVec;

    // Loop over subsectors.
    for (unsigned int i = 0; i < numSubSecsEta_; i++) {
      bool inside = this->insideEtaRange(stub, zOuterMinSub_[i], zOuterMaxSub_[i]);
      insideVec.push_back(inside);
    }

    return insideVec;
  }
}

//=== Check if stub is within eta sector or subsector that is delimated by specified zTrk range.

bool Sector::insideEtaRange( const Stub* stub, float zRangeMin, float zRangeMax) const {
  // Lower edge of this eta region defined by line from (r,z) = (0,-beamWindowZ) to (chosenRofZ_, zRangeMin).
  // Upper edge of this eta region defined by line from (r,z) = (0, beamWindowZ) to (chosenRofZ_, zRangeMax).

  float zMin, zMax;

  bool inside;

  if ( ! handleStripsEtaSec_) {
    //--- Don't modify algorithm to allow for uncertainty in stub (r,z) coordinates caused by 2S module strip length?

    // Calculate z coordinate of lower edge of this eta region, evaluated at radius of stub.
    zMin = ( zRangeMin * stub->r() - beamWindowZ_ * fabs(stub->r() - chosenRofZ_) ) / chosenRofZ_;
    // Calculate z coordinate of upper edge of this eta region, evaluated at radius of stub.
    zMax = ( zRangeMax * stub->r() + beamWindowZ_ * fabs(stub->r() - chosenRofZ_) ) / chosenRofZ_;

    // zMin = ( zRangeMin * stub->r() - beamWindowZ_ * fabs(stub->r() - rOuterMin_) ) / rOuterMin_;
    // zMax = ( zRangeMax * stub->r() + beamWindowZ_ * fabs(stub->r() - rOuterMax_) ) / rOuterMax_;

    inside = (stub->z() > zMin && stub->z() < zMax);

  } else {
    //--- Do modify algorithm to allow for uncertainty in stub (r,z) coordinates caused by 2S module strip length?

    float stubMinR = stub->r() - stub->rErr(); 
    float stubMaxR = stub->r() + stub->rErr(); 
    float stubMinZ = stub->z() - stub->zErr(); 
    float stubMaxZ = stub->z() + stub->zErr(); 

    // Calculate z coordinate of lower edge of this eta region, evaluated at radius of stub.
    float rStubA = (zRangeMin + beamWindowZ_) >= 0 ? stubMinR : stubMaxR; // stub r coordinate uncertain (especially in endcap), so use one which gives most -ve zMin.
    zMin = -beamWindowZ_ + (rStubA / chosenRofZ_) * (zRangeMin + beamWindowZ_);

    // Calculate z coordinate of upper edge of this eta region, evaluated at radius of stub.
    float rStubB = (zRangeMax - beamWindowZ_) >= 0 ? stubMaxR : stubMinR; // stub r coordinate uncertain (especially in endcap), so use one which gives most +ve zMax.
    zMax =  beamWindowZ_ + (rStubB / chosenRofZ_) * (zRangeMax - beamWindowZ_);

    inside = (stubMaxZ > zMin && stubMinZ < zMax);
  }

  return inside;
}

//=== Check if stub is inside this phi region.

bool Sector::insidePhi( const Stub* stub ) const {

  // N.B. The logic here for preventing a stub being assigned to > 2 sectors seems overly agressive.
  // But attempts at improving it have failed ...

  bool okPhi    = true;
  bool okPhiTrk = true;

  if (useStubPhi_) {
    float delPhi = reco::deltaPhi(stub->phi(), phiCentre_); // Phi difference between stub & sector in range -PI to +PI.
    float tolerancePhi = stub->phiDiff(chosenRofPhi_, minPt_); // How much stub phi might differ from track phi because of track curvature.
    float outsidePhi = fabs(delPhi) - sectorHalfWidth_ - tolerancePhi; // If > 0, then stub is not compatible with being inside this sector. 
    if (outsidePhi > 0) okPhi = false;
  }

  if (useStubPhiTrk_) {
    // Estimate either phi0 of track from stub info, or phi of the track at radius chosenRofPhi_.
    float phiTrk = stub->trkPhiAtR( chosenRofPhi_ ).first; // N.B. This equals stub->beta() if chosenRofPhi_ = 0.
    float delPhiTrk = reco::deltaPhi(phiTrk, phiCentre_); // Phi difference between stub & sector in range -PI to +PI.
    float tolerancePhiTrk = assumedPhiTrkRes_ * (2*sectorHalfWidth_); // Set tolerance equal to nominal resolution assumed in phiTrk
    if (calcPhiTrkRes_) {
      // Calculate uncertainty in phiTrk due to poor resolution in stub bend
      float phiTrkRes = stub->trkPhiAtRres( chosenRofPhi_ );
      // Reduce tolerance if this is smaller than the nominal assumed resolution.
      tolerancePhiTrk = min(tolerancePhiTrk, phiTrkRes);
    }
    float outsidePhiTrk = fabs(delPhiTrk) - sectorHalfWidth_ - tolerancePhiTrk; // If > 0, then stub is not compatible with being inside this sector.

    // Modify algorithm to allow for uncertainty due to 2S module strip length, if requested.
    if (handleStripsPhiSec_) {
      float chosenStubPhiErr = stub->trkPhiAtR( chosenRofPhi_ ).second; // The "Err" here is uncertainty due to 2S strip length.
      outsidePhiTrk -= chosenStubPhiErr;
    }

    if (outsidePhiTrk > 0) okPhiTrk = false;
  }

  return (okPhi && okPhiTrk);
}

//=== For performance studies, note which stubs on given tracking particle are inside the sector.
//=== Returns two booleans for each stub, indicating if they are in phi & eta sectors respectively.
//=== AND them together to get (eta,phi) sector decision.

unordered_map<const Stub*, pair<bool, bool> > Sector::stubsInside ( const TP& tp) const {
  unordered_map<const Stub*, pair<bool, bool> > inside;
  // Loop over stubs produced by tracking particle
  const vector<const Stub*>& assStubs= tp.assocStubs();
  for (const Stub* stub: assStubs) {
    // Check if this stub is inside sector
    inside[stub] = pair<bool, bool>(this->insidePhi(stub), this->insideEta(stub));
  }
  return inside;
}

//=== Count number of stubs in given tracking particle which are inside this (phi,eta) sector;
//=== or inside it if only the eta cuts are applied; or inside it if only the phi cuts are applied.
//=== The results are returned as the 3 last arguments of the function.

void Sector::numStubsInside( const TP& tp, 
                             unsigned int& nStubsInsideEtaPhi, unsigned int& nStubsInsideEta, 
			     unsigned int& nStubsInsidePhi) const 
{
  nStubsInsideEtaPhi = 0;
  nStubsInsideEta    = 0;
  nStubsInsidePhi    = 0;
  for (const auto& iter: this->stubsInside(tp) ) {
    bool insidePhi = iter.second.first;
    bool insideEta = iter.second.second;
    if (insidePhi && insideEta) nStubsInsideEtaPhi++;
    if (insideEta)              nStubsInsideEta++;
    if (insidePhi)              nStubsInsidePhi++;
  }
}

// Digitize a floating point number to 2s complement integer. (Kristian Harder)

Long64_t Sector::forceBitWidth( const float value, const UInt_t nBits) const {

  // slightly hand-waving treatment of 2s complement
  Long64_t sign = 1;
  if (value<0) sign = -1;
  Long64_t iValue = Long64_t(fabs(value));
  Long64_t mask =  (Long64_t(1)<<nBits)-Long64_t(1);
  Long64_t result = sign*(iValue & mask);
  return result;
}

//=== Check if stub is within subsectors in eta that sector may be divided into. Uses digitized calculation corresponding to GP firmware. (Kristian Harder)

vector<bool> Sector::subEtaFwCalc(const int rT, const int z) const {

  // these should be depending on parameters!
  // HACK BEGIN
  float numPhiSecPerOct = 2.;
  float numPhiSec = 8* numPhiSecPerOct;
  float cBins = 64.;
  float cBinBase = 2.0 * M_PI / numPhiSec / cBins;
  float Bfield = 3.8112;
  float Bconv = Bfield*0.00015;
  float ptCut = 3.0;
  float mBins = 32.;
  float mBinBase = 2.0*Bconv /ptCut/mBins;
  float rTBase = cBinBase/mBinBase / (1<<9);
  float zBase = 1./0.64;
  // HACK END
  
  // Ian - perhaps better to set above using cfg params?
  //Bfield = settings_->getBfield();
  //Bconv  = settings_->invPtToDphi() / 10.;  // extra factor of 10 as hardware uses mm.
  //ptCut  = settings_->houghMinPt();
  //cBins  = settings_->houghNbinsPhi();
  //mBins  = settings_->houghNbinsPt();
  //numPhiSec = settings_->numPhiSectors();

  // unit transformations: firmware uses mm, software uses cm
  float BeamWindow=10.*beamWindowZ_;
  float T_rphi=10.*chosenRofPhi_;
  float T_rz=10.*chosenRofZ_;

  // actual algorithm as used in firmware, mostly using same variable names
  float Beam_over_T = BeamWindow/T_rz;
  float Beam_over_T_base = 1./(1<<24);
  Long64_t bot = forceBitWidth(Beam_over_T*rTBase/zBase/Beam_over_T_base,25);
  Long64_t bw = forceBitWidth(BeamWindow/zBase/Beam_over_T_base,48);
  float tanlSecBase=1./(1<<16);
  float etaSecMid = (settings_->etaRegions()[iEtaReg_] + settings_->etaRegions()[iEtaReg_ + 1]) / 2.0;
  float tanlSecMid = 1.0 / tan(2.0*atan(exp(-etaSecMid))); 
  Long64_t int_tanlSecMid = forceBitWidth(int(tanlSecMid * rTBase / zBase / tanlSecBase),31);
  Long64_t tanlSec_Mid = forceBitWidth(int_tanlSecMid,25);
  Long64_t r=forceBitWidth(rT+T_rphi/rTBase,12);
  Long64_t g = forceBitWidth(bot*r-bw,48);
  Long64_t absg = abs(g);
  Long64_t shift_g = forceBitWidth(absg>>24,24);
  Long64_t tlsr=forceBitWidth(tanlSec_Mid*r,37);
  Long64_t shift_tlsr = forceBitWidth(tlsr>>16,21);

  vector<bool> insideVec;
  insideVec.push_back(z<=(shift_tlsr+shift_g));
  insideVec.push_back(z>=(shift_tlsr-shift_g));
  return insideVec;
}
