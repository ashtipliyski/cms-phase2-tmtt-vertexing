#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/DigitalTrack.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include <map>

//=== Note configuration parameters.

DigitalTrack::DigitalTrack(const Settings* settings) :

  // Check DigitalTrack correctly initialized.
  ranInit_ (false),
  ranMake_ (false),

  // Digitization configuration parameters
  settings_(settings),

  // Number of phi sectors and phi octants.
  numPhiSectors_ (settings->numPhiSectors()),
  numPhiOctants_ (8),
  // Phi sector and phi octant width (radians)
  phiSectorWidth_(2.*M_PI / float(numPhiSectors_)), 
  phiOctantWidth_(2.*M_PI / float(numPhiOctants_)), 
  // Radius from beamline with respect to which stub r coord. is measured.
  chosenRofPhi_  (settings->chosenRofPhi()),

  // Number of q/Pt bins in Hough  transform array.
  nbinsPt_       ((int) settings->houghNbinsPt()),
  invPtToDPhi_   (settings->invPtToDphi())
{
}

//=== Get digitisation configuration parameters for the specific track fitter being used here.

void DigitalTrack::getDigiCfg(const string& fitterName) {
  if (fitterName == "SimpleLR") {
    // SimpleLR track fitter
    skipTrackDigi_  = settings_->slr_skipTrackDigi();
    oneOver2rBits_  = settings_->slr_oneOver2rBits();
    oneOver2rRange_ = settings_->slr_oneOver2rRange();
    phi0Bits_       = settings_->slr_phi0Bits();
    phi0Range_      = settings_->slr_phi0Range();
    z0Bits_         = settings_->slr_z0Bits();
    z0Range_        = settings_->slr_z0Range();
    tanLambdaBits_  = settings_->slr_tanlambdaBits();
    tanLambdaRange_ = settings_->slr_tanlambdaRange();
    chisquaredBits_ = settings_->slr_chisquaredBits();
    chisquaredRange_= settings_->slr_chisquaredRange();
  } else {
    // KF track fitter
    // Also used for all other fitters, though unlikely to be correct them them ...
    if (fitterName == "KF4ParamsComb") {
      skipTrackDigi_  = settings_->kf_skipTrackDigi();
    } else {
      skipTrackDigi_  = settings_->other_skipTrackDigi(); // Allows to skip digitisation for other fitters
    }
    oneOver2rBits_  = settings_->kf_oneOver2rBits();
    oneOver2rRange_ = settings_->kf_oneOver2rRange();
    phi0Bits_       = settings_->kf_phi0Bits();
    phi0Range_      = settings_->kf_phi0Range();
    z0Bits_         = settings_->kf_z0Bits();
    z0Range_        = settings_->kf_z0Range();
    tanLambdaBits_  = settings_->kf_tanlambdaBits();
    tanLambdaRange_ = settings_->kf_tanlambdaRange();
    chisquaredBits_ = settings_->kf_chisquaredBits();
    chisquaredRange_= settings_->kf_chisquaredRange();
  }

  // Calculate multipliers to digitize the floating point numbers.
  oneOver2rMult_ = pow(2.,oneOver2rBits_)/oneOver2rRange_;
  phi0Mult_      = pow(2.,phi0Bits_)/phi0Range_;
  z0Mult_        = pow(2.,z0Bits_)/z0Range_;
  tanLambdaMult_ = pow(2.,tanLambdaBits_)/tanLambdaRange_;
  chisquaredMult_= pow(2.,chisquaredBits_)/chisquaredRange_;
}

//=== Initialize track with original, floating point track params

void DigitalTrack::init(const string& fitterName, unsigned int iPhiSec, unsigned int iEtaReg, int mbin, int cbin, int mBinhelix, int cBinhelix, float qOverPt_orig, float phi0_orig, float tanLambda_orig, float z0_orig, float chisquared_orig, unsigned int nLayers, bool consistent, bool accepted, float tp_pt, float tp_eta, int tp_index, bool tp_useForAlgEff, bool tp_useForEff, int tp_pdgId, float tp_phi0, float tp_tanLambda, float tp_z0, float tp_qoverpt) {

  ranInit_ = true; // Note we ran init().

  fitterName_     = fitterName;

  // Get digitisation parameters for this particular track fitter.
  this->getDigiCfg(fitterName);

  phiSectorCentre_ = phiSectorWidth_ * (0.5 + double(iPhiSec)) - M_PI; 

  qOverPt_orig_   = qOverPt_orig;
  oneOver2r_orig_ = qOverPt_orig*invPtToDPhi_;
  phi0_orig_      = phi0_orig;
  phi0rel_orig_   = reco::deltaPhi(phi0_orig_, phiSectorCentre_);
  tanLambda_orig_ = tanLambda_orig;
  z0_orig_        = z0_orig;
  chisquared_orig_= chisquared_orig;
  nlayers_        = nLayers;
  iPhiSec_        = iPhiSec;
  iEtaReg_        = iEtaReg;
  mBin_           = mbin;
  cBin_           = cbin;
  mBinhelix_      = mBinhelix;
  cBinhelix_      = cBinhelix;

  consistent_     = consistent;
  accepted_       = accepted;
  tp_tanLambda_   = tp_tanLambda;
  tp_pt_          = tp_pt;
  tp_eta_         = tp_eta;
  tp_index_       = tp_index;
  tp_useForAlgEff_= tp_useForAlgEff;
  tp_useForEff_   = tp_useForEff;
  tp_pdgId_       = tp_pdgId;
  tp_phi0_        = tp_phi0;
  tp_z0_          = tp_z0;
  tp_qoverpt_     = tp_qoverpt;
}

//=== Digitize track

void DigitalTrack::makeDigitalTrack() {

  if (! ranInit_) throw cms::Exception("DigitalTrack: You forgot to call init() before makeDigitalTrack()!");

  ranMake_ = true; // Note we ran makeDigitalTrack()

  if (skipTrackDigi_) {
    // Optionally skip track digitisaton if done internally inside track fitting code, so
    // retain original helix params.
    iDigi_oneOver2r_  = 0;
    iDigi_phi0rel_    = 0;
    iDigi_tanLambda_  = 0;
    iDigi_z0_         = 0;
    iDigi_chisquared_ = 0;

    oneOver2r_    = oneOver2r_orig_; 
    qOverPt_      = qOverPt_orig_;
    phi0rel_      = phi0rel_orig_;
    phi0_         = phi0_orig_;
    tanLambda_    = tanLambda_orig_;
    z0_           = z0_orig_;
    chisquared_   = chisquared_orig_;

  } else {

    //--- Digitize variables

    iDigi_oneOver2r_  = floor(oneOver2r_orig_*oneOver2rMult_)+1;
    iDigi_phi0rel_    = floor(phi0rel_orig_*phi0Mult_);
    iDigi_tanLambda_  = floor(tanLambda_orig_*tanLambdaMult_);
    iDigi_z0_         = floor(z0_orig_*z0Mult_);
    iDigi_chisquared_ = floor(chisquared_orig_*chisquaredMult_);

    // If fitted declared track invalid, it will have set its chi2 to very large number. 
    // So truncate it at maximum allowed by digitisation range.
    if ( ! accepted_ ) iDigi_chisquared_ = pow(2.,chisquaredBits_) - 1;

    // if(settings_->digitizeSLR()){
    //   mBinhelix_ = floor(iDigi_1over2r_/pow(2,5));
    //   cBinhelix_ = floor(iDigi_phiT_/pow(2,7));
    // }

    //--- Determine floating point track params from digitized numbers (so with degraded resolution).

    oneOver2r_    = (iDigi_oneOver2r_)/oneOver2rMult_;
    qOverPt_      = oneOver2r_/invPtToDPhi_;
    phi0rel_      = (iDigi_phi0rel_)/phi0Mult_;
    phi0_         = reco::deltaPhi(phi0rel_, -phiSectorCentre_);
    tanLambda_    = (iDigi_tanLambda_)/tanLambdaMult_;
    z0_           = (iDigi_z0_)/z0Mult_;
    chisquared_   = (iDigi_chisquared_)/chisquaredMult_;

    // Check that track coords. are within assumed digitization range.
    // (Don't bother with this check if using Linear Regression or Chi2 fitters without Seed Filter,
    //  since they are not expected to work well in this scenario)
    // (Also, only do this check for accepted tracks, since SimpleLR currently has bug affected rejected tracks)
    if (fitterName_ == "KF4ParamsComb" || settings_->useSeedFilter()) {
      if (accepted_) {
        this->checkInRange();
      }
    }

    // Check that digitization followed by undigitization doesn't change results too much.
    // (Don't bother with this check if using Linear Regression or Chi2 fitters without Seed Filter,
    //  since they are not expected to work well in this scenario)
    if (fitterName_ == "KF4ParamsComb" || settings_->useSeedFilter()) {
      this->checkAccuracy();
    }
  }
}

//=== Check that stub coords. are within assumed digitization range.

void DigitalTrack::checkInRange() const {
  if (fabs(oneOver2r_orig_) >= 0.5*oneOver2rRange_)   throw cms::Exception("DigitalTrack: Track oneOver2r is out of assumed digitization range.")<<" |oneOver2r| = " <<fabs(oneOver2r_orig_) <<" > "<<0.5*oneOver2rRange_<<"; Fitter="<<fitterName_<<"; track accepted = "<<accepted_<<endl;  
  if (fabs(phi0rel_orig_) >= 0.5*phi0Range_)   throw cms::Exception("DigitalTrack: Track phi0rel is out of assumed digitization range.")<<" |phi0rel| = " <<fabs(phi0rel_orig_) <<" > "<<0.5*phi0Range_<<"; Fitter="<<fitterName_<<"; track accepted = "<<accepted_<<endl;  
  if (fabs(z0_orig_) >= 0.5*z0Range_)   throw cms::Exception("DigitalTrack:  Track z0 is out of assumed digitization range.")<<" |z0| = " <<fabs(z0_orig_) <<" > "<<0.5*z0Range_<<"; Fitter="<<fitterName_<<"; track accepted = "<<accepted_<<endl;  
  if (fabs(tanLambda_orig_) >= 0.5*tanLambdaRange_)   throw cms::Exception("DigitalTrack: Track tanLambda is out of assumed digitization range.")<<" |tanLambda| = " <<fabs(tanLambda_orig_) <<" > "<<0.5*tanLambdaRange_<<"; Fitter="<<fitterName_<<"; track accepted = "<<accepted_<<endl;  
  if (accepted_) { // Tracks declared invalid by fitter can have very large original chi2.
    if (chisquared_orig_ >= chisquaredRange_ or chisquared_orig_ < 0.)   throw cms::Exception("DigitalTrack: Track chisquared is out of assumed digitization range.")<<" chisquared = " <<chisquared_orig_ <<" > "<<chisquaredRange_<<" or < 0"<<"; Fitter="<<fitterName_<<"; track accepted = "<<accepted_<<endl;  
  }
}

//=== Check that digitisation followed by undigitisation doesn't change significantly the stub coordinates.

void DigitalTrack::checkAccuracy() const {
  float TA = qOverPt_- qOverPt_orig_;
  float TB = reco::deltaPhi(phi0_, phi0_orig_);
  float TC = z0_ - z0_orig_;
  float TD = tanLambda_ - tanLambda_orig_;
  float TE = chisquared_  - chisquared_orig_;
  if ( ! accepted_ ) TE = 0.; // don't apply this check to tracks declared invalid by the fitter.

  static map<string, unsigned int> nErr; // Count precision errors from each fitter.
  if (nErr.find(fitterName_) == nErr.end()) nErr[fitterName_] = 0; // Initialize error count.
  const  unsigned int maxErr = 20;  // Print error message only this number of times.
  if (nErr[fitterName_] < maxErr) {
    if (fabs(TA) > 0.01 || fabs(TB) > 0.001 || fabs(TC) > 0.05 || fabs(TD) > 0.001 || fabs(TE) > 0.5) {
      nErr[fitterName_]++;
      cout<<"WARNING: DigitalTrack lost precision: "<<fitterName_<<" accepted="<<accepted_<<" "<<TA<<" "<<TB<<" "<<TC<<" "<<TD<<" "<<TE<<endl;
    }
  }
}
