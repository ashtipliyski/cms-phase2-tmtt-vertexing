#include "TMTrackTrigger/TMTrackFinder/interface/HTrphi.h"
#include "TMTrackTrigger/TMTrackFinder/interface/InputData.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <vector>

//=== The r-phi Hough Transform array for a single (eta,phi) sector.
//===
//=== Its axes are (q/Pt, phiTrk), where phiTrk is the phi at which the track crosses a 
//=== user-configurable radius from the beam-line.

using namespace std;

// Maximum |gradient| of line corresponding to any stub. Should be less than the value of 1.0 assumed by the firmware.
float        HTrphi::maxLineGradient_ = 0.;
// Error count when stub added to cell which does not lie NE, E or SE of stub added to previous HT column.
unsigned int HTrphi::numErrorsTypeA_ = 0;
// Error count when stub added to more than 2 cells in one HT column (problem only for Thomas' firmware).
unsigned int HTrphi::numErrorsTypeB_ = 0;
// Error count normalisation
unsigned int HTrphi::numErrorsNormalisation_ = 0;

//=== Initialise
 
void HTrphi::init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg, 
		  float etaMinSector, float etaMaxSector, float phiCentreSector) {

  HTbase::init(settings, iPhiSec, iEtaReg);

  invPtToDphi_         = settings->invPtToDphi();

  int nCellsHT = settings->houghNcellsRphi(); // Total number of required cells in HT array (if > 0)

  //--- Specification of HT q/Pt axis.

  maxAbsQoverPtAxis_  = 1./settings->houghMinPt(); // Max. |q/Pt| covered by  HT array.
  nBinsQoverPtAxis_   = settings->houghNbinsPt();  // No. of bins in HT array in q/Pt.
  if (nCellsHT > 0) nBinsQoverPtAxis_ = 1; // Will calculate number of bins automatically. Initialize it to non-zero value.
  binSizeQoverPtAxis_ = 2*maxAbsQoverPtAxis_ / nBinsQoverPtAxis_;

  //--- Specification of HT phiTrk axis

  // N.B. phiTrk corresponds to phi where track crosses radius = chosenRofPhi_.
  chosenRofPhi_       = settings->chosenRofPhi();
  phiCentreSector_    = phiCentreSector; // Centre of phiTrk sector.
  maxAbsPhiTrkAxis_   = M_PI / float(settings->numPhiSectors()); // Half-width of phiTrk axis in HT array.
  nBinsPhiTrkAxis_    = settings->houghNbinsPhi(); // No. of bins in HT array phiTrk
  if (nCellsHT > 0) nBinsPhiTrkAxis_ = 1; // Will calculate number of bins automatically. Initialize it to non-zero value.
  binSizePhiTrkAxis_  = 2*maxAbsPhiTrkAxis_ / nBinsPhiTrkAxis_;

  // Did user specify number of cells required in HT array? If so, determine number of bins along
  // array axes such that their product equals required number of cells, and that their ratio gives
  // a maximum line |gradient| of stubs crossing the array of 1.0.
  if (nCellsHT > 0) { 
    // Get line gradient with current array axes.
    float currentLineGrad = this->calcMaxLineGradArray();
    // Calculate new number of bins on each axis to meet constraint.
    float fact = nBinsQoverPtAxis_ * currentLineGrad / nBinsPhiTrkAxis_;
    nBinsQoverPtAxis_ = ceil( sqrt(nCellsHT * fact) );
    nBinsPhiTrkAxis_  = int ( sqrt(nCellsHT / fact) );
    // And recalculate bin size accordingly.
    binSizeQoverPtAxis_ = 2*maxAbsQoverPtAxis_ / nBinsQoverPtAxis_;
    binSizePhiTrkAxis_  = 2*maxAbsPhiTrkAxis_  / nBinsPhiTrkAxis_;
  }

  // Note max. |gradient| that the line corresponding to any stub in any of the r-phi HT arrays could have.
  // Firmware assumes this should not exceed 1.0;
  HTrphi::maxLineGradient_ = max( HTrphi::maxLineGradient_, this->calcMaxLineGradArray());

  // Optionally merge 2x2 neighbouring cells into a single cell at low Pt, to reduce efficiency loss due to 
  // scattering.
  enableMerge2x2_  = settings->enableMerge2x2();
  minInvPtToMerge2x2_ = 1./(settings->maxPtToMerge2x2());
  if (minInvPtToMerge2x2_ > maxAbsQoverPtAxis_) enableMerge2x2_ = false;

  // Merging cells only allowed if HT array dimensions are even.
  if (enableMerge2x2_ && (nBinsQoverPtAxis_%2 != 0 || nBinsPhiTrkAxis_%2 != 0)) throw cms::Exception("HTrphi: You are not allowed to set EnableMerge2x2=True if you have an odd number of bins in r-phi HT array ")<<nBinsQoverPtAxis_<<" "<<nBinsPhiTrkAxis_<<endl;

  //--- Other options used when filling the HT.

  // Don't fill all HT cells nominally crossed by line corresponding to stub.
  killSomeHTCellsRphi_ = settings->killSomeHTCellsRphi();
  // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when fill stubs in r-phi HT?
  handleStripsRphiHT_  = settings->handleStripsRphiHT();

  // Used to kill excess stubs or tracks that can't be transmitted within time-multiplexed period.
  nReceivedStubs_ = 0;
  busyInputSectorKill_     = settings_->busyInputSectorKill();   // Kill excess stubs going fron GP to HT?
  busyInputSectorNumStubs_ = settings_->busyInputSectorNumStubs(); // Max. num. of stubs that can be sent from GP to HT within TM period
  busySectorKill_          = settings_->busySectorKill();        // Kill excess tracks flowing out of HT?
  busySectorNumStubs_      = settings_->busySectorNumStubs();    // Max. num. of stubs that can be sent out of HT within TM period
  busySectorEachCharge_    = settings_->busySectorEachCharge();  // +ve & -ve tracks output on separate optical link?
  busySectorMbinRanges_    = settings_->busySectorMbinRanges();  // or individual m bin (=q/Pt) ranges to be output to optical links. 
  busySectorUseMbinRanges_ = ( ! busySectorMbinRanges_.empty() ); // m bin ranges option disabled if vector empty. 
  busySectorMbinOrder_     = settings_->busySectorMbinOrder();   // Specifies which m bins should be grouped together by BusySectorMbinRanges. If empty, then they are grouped in order 0,1,2,3,4,5 ...
  busySectorUseMbinOrder_  = ( ! busySectorMbinOrder_.empty() );

  bool rescaleMbins = false;
  if (busySectorUseMbinRanges_) {
    // Check if the total number of bins specified in cfg option BusySectorMbinRanges corresponds
    // to the number of m bins (q/Pt) in the HT. If not, determine how much the ranges must be scaled
    // to make this true.
    unsigned int nTotalBins = 0;
    for (unsigned int j = 0; j < busySectorMbinRanges_.size(); j++) {
      nTotalBins += busySectorMbinRanges_[j];
    }
    rescaleMbins = (nTotalBins != nBinsQoverPtAxis_);
    // No rescaling allowed with MBinOrder option.
    if (rescaleMbins && busySectorUseMbinOrder_) throw cms::Exception("HTrphi: BusySectorUserMbin error");
    float rescaleFactor = rescaleMbins  ?  float(nBinsQoverPtAxis_)/float(nTotalBins)  :  1. ;
    // Find lower and upper inclusive limits of each m bin range to be sent to a separate optical link.
    busySectorMbinLow_.resize( busySectorMbinRanges_.size() );
    busySectorMbinHigh_.resize( busySectorMbinRanges_.size() );
    float mBinSum = 0.;
    for (unsigned int i = 0; i < busySectorMbinRanges_.size(); i++) {
      busySectorMbinLow_[i]  = std::round(mBinSum);
      busySectorMbinHigh_[i] = std::round(mBinSum + rescaleFactor * busySectorMbinRanges_[i]) - 1;
      mBinSum += rescaleFactor * busySectorMbinRanges_[i];
    }
  }

  //--- Options for duplicate track removal after running HT.
  unsigned int dupTrkAlgRphi = settings->dupTrkAlgRphi();
  HTbase::killDupTrks_.init(settings, dupTrkAlgRphi);

  // Resize HT array to suit these specifications, and initialise each cell with configuration parameters.
  HTbase::htArray_.resize(nBinsQoverPtAxis_, nBinsPhiTrkAxis_, false);

  const bool isRphiHT = true;
  for (unsigned int i = 0; i < nBinsQoverPtAxis_; i++) {
    for (unsigned int j = 0; j < nBinsPhiTrkAxis_; j++) {
      pair<float, float> helix = this->helix2Dconventional(i, j); // Get track params at centre of cell.
      float qOverPt = helix.first;
      HTbase::htArray_(i,j).init( settings, isRphiHT, iPhiSec, iEtaReg, 
				  etaMinSector, etaMaxSector, qOverPt, i); // Calls HTcell::init()
    }
  }

  static bool first = true;
  if (first) {
    first = false;
    cout<<"=== R-PHI HOUGH TRANSFORM AXES RANGES: abs(q/Pt) < "<<maxAbsQoverPtAxis_<<" & abs(track-phi) < "<<maxAbsPhiTrkAxis_<<" ==="<<endl<<endl;
    cout<<"=== R-PHI HOUGH TRANSFORM ARRAY SIZE: q/Pt bins = "<<nBinsQoverPtAxis_<<" & track-phi bins = "<<nBinsPhiTrkAxis_<<endl; 
    if (busySectorKill_ && busySectorUseMbinRanges_ && rescaleMbins) {
      cout<<"=== R-PHI HOUGH TRANSFORM WARNING: Rescaled m bin ranges specified by cfg parameter BusySectorMbinRanges, as they were inconsistent with total number of m bins in HT."<<endl;
      cout<<"=== Rescaled values for BusySectorMbinRanges =";
      for (unsigned int i = 0; i < busySectorMbinRanges_.size(); i++) {
        cout<<" "<<(busySectorMbinHigh_[i] - busySectorMbinLow_[i] + 1);
      }
      cout<<endl;
    }
  }
}

//=== Add stub to HT array.
//=== If eta subsectors are being used within each sector, specify which ones the stub is compatible with.

void HTrphi::store(const Stub* stub, const vector<bool>& inEtaSubSecs) {
  // Optionally, only store stubs that can be sent from GP to HT within TM period.
  if ( ( ! busyInputSectorKill_) || (nReceivedStubs_ <  busyInputSectorNumStubs_) ) {

    nReceivedStubs_++;

    // Loop over q/Pt related bins in HT array.
    for (unsigned int i = 0; i < nBinsQoverPtAxis_; i++) {

      // In this q/Pt bin, find the range of phi bins that this stub is consistent with.
      pair<unsigned int, unsigned int> iRange = this->iPhiRange( stub, i);
      unsigned int iPhiTrkBinMin = iRange.first;
      unsigned int iPhiTrkBinMax = iRange.second;

      // Store stubs in these cells.
      for (unsigned int j = iPhiTrkBinMin; j <= iPhiTrkBinMax; j++) {  

	bool canStoreStub = true;
	unsigned int iStore = i;
	unsigned int jStore = j;

	// Optionally merge 2x2 neighbouring cells into a single cell at low Pt, to reduce efficiency loss
	// due to scattering.
	if (enableMerge2x2_) {
	  // Check if this cell is merged with its neighbours (as in low Pt region).
	  if (this->mergedCell(i, j)) {
	    // Get location of cell that this cell is merged into (iStore, jStore).
	    // Calculation assumes HT array has even number of bins in both dimensions.
	    if (i%2 == 1) iStore = i - 1;
	    if (j%2 == 1) jStore = j - 1;
	    // If this stub was already stored in this merged 2x2 cell, then don't store it again.
	    if (HTbase::htArray_(iStore, jStore).stubStoredInCell( stub )) canStoreStub = false;
	  }
	}

	if (canStoreStub) 
    HTbase::htArray_(iStore, jStore).store( stub, inEtaSubSecs ); // Calls HTcell::store()
      }

      // Check that limitations of firmware would not prevent stub being stored correctly in this HT column.
      this->countFirmwareErrors(i, iPhiTrkBinMin, iPhiTrkBinMax);
    }
  }
}

//=== Determine the m-bin (q/pt) range the specified track is in. (Used if outputting each m bin range on a different opto-link).

unsigned int HTrphi::getMbinRange(const L1track2D& trk) const {
  if (busySectorUseMbinRanges_) {
    unsigned int mBin = trk.getCellLocation().first;
    unsigned int mBinOrder;
    if (busySectorUseMbinOrder_) {
      // User wants to group bins in a wierd order.
      mBinOrder = 99999;
      for (unsigned int k = 0; k < busySectorMbinOrder_.size(); k++) {
	if (mBin == busySectorMbinOrder_[k]) mBinOrder = k;
      }
      if (mBinOrder == 99999) throw cms::Exception("HTrphi::getMbinRange() mBinOrder calculation wrong.");
    } else {
      // User grouping bins in numerical order 0,1,2,3,4,5...
      mBinOrder = mBin;
    }
    for (unsigned int i = 0; i < busySectorMbinRanges_.size(); i++) {
      if (mBinOrder >= busySectorMbinLow_[i] && mBinOrder <= busySectorMbinHigh_[i]) return i;
    }
    throw cms::Exception("HTrphi::getMbinRange() messed up");
  } else {
    return 0;
  }
}

//=== For a given Q/Pt bin, find the range of phi bins that a given stub is consistent with.
//=== Return as a pair (min bin, max bin)
//=== If it range lies outside the HT array, then the min bin will be set larger than the max bin.

pair<unsigned int, unsigned int> HTrphi::iPhiRange( const Stub* stub, unsigned int iQoverPtBin, bool debug) const {

  // Note q/Pt value corresponding to centre of this bin.
  float qOverPtBin    = -maxAbsQoverPtAxis_ + (iQoverPtBin + 0.5) * binSizeQoverPtAxis_;
  // Note change in this q/Pt value needed to reach either edge of the bin. 
  float qOverPtBinVar = 0.5*binSizeQoverPtAxis_;

  // Reducing effective bin width can reduce fake rate.
  //qOverPtVar = 0.4*binSizeQoverPtAxis_;

  // Calculate range of track-phi that would allow a track in this q/Pt range to pass through the stub.
  float phiTrk    = stub->phi() + invPtToDphi_ * qOverPtBin    *     (stub->r() - chosenRofPhi_);
  // The next line does the phiTrk calculation without the usual approximation, but it doesn't 
  // improve performance.
  //float phiTrk    = stub->phi() + asin(invPtToDphi_ * qOverPtBin * stub->r()) - asin(invPtToDphi_ * qOverPtBin * chosenRofPhi_);
  float phiTrkVar =               invPtToDphi_ * qOverPtBinVar * fabs(stub->r() - chosenRofPhi_);
  float phiTrkMin = phiTrk - phiTrkVar;
  float phiTrkMax = phiTrk + phiTrkVar;

  // Allow for uncertainty due to strip length if requested.
  if (handleStripsRphiHT_) {
    // Estimate uncertainty due to strip length, using first order derivative of phiTrk w.r.t. stub coords.
    // Note that barrel modules only care about zErr and endcap ones about rErr.
    float phiTrkVarStub;
    if (stub->barrel()) {
      phiTrkVarStub = 0.;
    } else {
      phiTrkVarStub = invPtToDphi_ * fabs(qOverPtBin) * stub->rErr();
    }
    phiTrkMin -= phiTrkVarStub; 
    phiTrkMax += phiTrkVarStub; 
  }

  // Allow for multiple scattering/resolution
  // phiTrkMin -= 0.005;
  // phiTrkMax += 0.005;

  float deltaPhiMin = reco::deltaPhi(phiTrkMin, phiCentreSector_); // Offset to centre of sector.
  float deltaPhiMax = reco::deltaPhi(phiTrkMax, phiCentreSector_);
  pair<float, float> phiTrkRange( deltaPhiMin, deltaPhiMax );

  // Determine which HT array cell range in track-phi this range "phiTrkRange" corresponds to.
  pair<unsigned int, unsigned int> iPhiTrkBinRange = this->HTbase::convertCoordRangeToBinRange(phiTrkRange, nBinsPhiTrkAxis_, (-maxAbsPhiTrkAxis_), binSizePhiTrkAxis_, killSomeHTCellsRphi_);

  return iPhiTrkBinRange;
}

//=== Check that limitations of firmware would not prevent stub being stored correctly in this HT column.

void HTrphi::countFirmwareErrors(unsigned int iQoverPtBin, unsigned int iPhiTrkBinMin, unsigned int iPhiTrkBinMax) {
  static unsigned int iPhiTrkBinMinLast = 0;
  static unsigned int iPhiTrkBinMaxLast = 99999;
  // Reinitialize if this is left-most column in HT array.
  if (iQoverPtBin == 0) {
    iPhiTrkBinMinLast = 0;
    iPhiTrkBinMaxLast = 99999;
  }

  // Only do check if stub is being stored somewhere in this HT column.
  if (iPhiTrkBinMax >= iPhiTrkBinMin) {
    //--- Remaining code below checks that firmware could successfully store this stub in this column.
    //   (a) Does cell lie NE, E or SE of cell filled in previous column?
    bool OK_a = (iPhiTrkBinMin + 1 >= iPhiTrkBinMinLast) && (iPhiTrkBinMax <= iPhiTrkBinMaxLast + 1);
    //   (b) Are no more than 2 cells filled in this column (problem only for Thomas' firmware)
    bool OK_b = (iPhiTrkBinMax - iPhiTrkBinMin + 1 <= 2);

    if ( ! OK_a ) numErrorsTypeA_++;
    if ( ! OK_b ) numErrorsTypeB_++;
    numErrorsNormalisation_++; // No. of times a stub is added to an HT column.

    iPhiTrkBinMinLast = iPhiTrkBinMin;
    iPhiTrkBinMaxLast = iPhiTrkBinMax;
  }
}

//=== Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
//=== The helix parameters returned will be those corresponding to the two axes of the HT array.
//=== So they might be (q/pt, phi0) or (q/pt, phi65) etc. depending on the configuration.

pair<float, float> HTrphi::helix2Dhough(unsigned int i, unsigned int j) const {

  unsigned int iStore = i;
  unsigned int jStore = j;

  // If using merged 2x2 cells in low Pt parts of array, must correct for this.
  bool merged = false;
  if (enableMerge2x2_) {
    // Check if this cell is merged with its neighbours (as in low Pt region).
    if (this->mergedCell(i, j)) {
      merged = true;
      // Get location of cell that this cell is merged into (iStore, jStore).
      // Calculation assumes HT array has even number of bins in both dimensions.
      if (i%2 == 1) iStore = i - 1;
      if (j%2 == 1) jStore = j - 1;
    }
  }

  float qOverPt = -maxAbsQoverPtAxis_ + (iStore + 0.5)*binSizeQoverPtAxis_;
  float phiTrk  = -maxAbsPhiTrkAxis_  + (jStore + 0.5)*binSizePhiTrkAxis_;

  if (merged) {
    qOverPt += 0.5*binSizeQoverPtAxis_;
    phiTrk  += 0.5*binSizePhiTrkAxis_;
  }
  
  phiTrk = reco::deltaPhi(phiTrk + phiCentreSector_, 0.); // Correct phiTrk to centre of sector, taking care of 2*pi wrapping
  return pair<float, float>(qOverPt, phiTrk); 
}

//=== Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
//=== The helix parameters returned will be always be (q/pt, phi0), irrespective of how the axes
//=== of the HT array are defined.

pair<float, float> HTrphi::helix2Dconventional(unsigned int i, unsigned int j) const {
  // Get the helix parameters corresponding to the axes definitions of the HT.
  pair<float, float> helix2Dht = this->helix2Dhough(i,j);
  // Convert to the conventionally agreed pair of helix parameters, (q/pt, phi0).
  float qOverPt = helix2Dht.first; // easy
  float phi0    = reco::deltaPhi(helix2Dht.second + invPtToDphi_*chosenRofPhi_*qOverPt, 0.); // If HT defined track phi other than at r=0, must correct to get phi0. Allow for 2*pi wrapping of phi.
  return pair<float, float>(qOverPt, phi0);
}

//=== Which cell in HT array should this TP be in, based on its true trajectory?
//=== (If TP was outside HT array, then cell index is set to either -1 or nBins in array in the corresponding dimension).

pair<int, int> HTrphi::trueCell( const TP* tp ) const {

  // Get HT axis variables corresponding to this TP.
  float qOverPt = tp->qOverPt();
  float phiTrk  = tp->trkPhiAtR( chosenRofPhi_ );
  // Measure phi relative to centre of sector.
  float deltaPhi = reco::deltaPhi(phiTrk, phiCentreSector_);
  // Convert to bin numbers inside HT array.
  int iQoverPt = floor( ( qOverPt  - ( -maxAbsQoverPtAxis_) ) / binSizeQoverPtAxis_ );
  int iPhiTrk  = floor( ( deltaPhi - ( -maxAbsPhiTrkAxis_ ) ) / binSizePhiTrkAxis_  );
  // Check if this cell was within the HT array.
  if (iQoverPt >= 0 && iQoverPt < int(nBinsQoverPtAxis_) && iPhiTrk >= 0 && iPhiTrk < int(nBinsPhiTrkAxis_)) {
    // Check if this cell is merged with its neighbours (as in low Pt region), and if so return merged cell location.
    if (this->mergedCell((unsigned int) iQoverPt, (unsigned int) iPhiTrk)) {
      if (iQoverPt%2 == 1) iQoverPt -= 1;
      if (iPhiTrk%2  == 1) iPhiTrk  -= 1;
    }
  } else {
    // TP is not in this HT array at all, so flag this by setting "outside" bin index to either -1 or nBins in array.
    iQoverPt = max(iQoverPt, -1);
    iQoverPt = min(iQoverPt, (int) nBinsQoverPtAxis_);
    iPhiTrk  = max(iPhiTrk , -1);
    iPhiTrk  = min(iPhiTrk , (int) nBinsPhiTrkAxis_);
  }
  return pair<int, int>(iQoverPt, iPhiTrk);
}

//=== Which cell in HT array should this fitted track be in, based on its fitted trajectory?
//=== (If fitted track was outside HT array, then cell index is set to either -1 or nBins in array in the corresponding dimension).

pair<int, int> HTrphi::getCell( const L1fittedTrack* fitTrk ) const {

  // Get HT axis variables corresponding to this fitted track.
  float qOverPt = fitTrk->qOverPt();
  // Convert phi0 to phi at chosen radius used in HT.
  float phiTrk  = fitTrk->phiAtChosenR();
  // Measure phi relative to centre of sector.
  float deltaPhi = reco::deltaPhi(phiTrk, phiCentreSector_);
  // Convert to bin numbers inside HT array.
  int iQoverPt = floor( ( qOverPt  - ( -maxAbsQoverPtAxis_) ) / binSizeQoverPtAxis_ );
  int iPhiTrk  = floor( ( deltaPhi - ( -maxAbsPhiTrkAxis_ ) ) / binSizePhiTrkAxis_  );
  // Check if this cell was within the HT array.
  if (iQoverPt >= 0 && iQoverPt < int(nBinsQoverPtAxis_) && iPhiTrk >= 0 && iPhiTrk < int(nBinsPhiTrkAxis_)) {
    // Check if this cell is merged with its neighbours (as in low Pt region), and if so return merged cell location.
    if (this->mergedCell((unsigned int) iQoverPt, (unsigned int) iPhiTrk)) {
      if (iQoverPt%2 == 1) iQoverPt -= 1;
      if (iPhiTrk%2  == 1) iPhiTrk  -= 1;
    }
  } else {
    // Fitted track is not in this HT array at all, so flag this by setting "outside" bin index to either -1 or nBins in array.
    iQoverPt = max(iQoverPt, -1);
    iQoverPt = min(iQoverPt, (int) nBinsQoverPtAxis_);
    iPhiTrk  = max(iPhiTrk , -1);
    iPhiTrk  = min(iPhiTrk , (int) nBinsPhiTrkAxis_);
  }
  return pair<int, int>(iQoverPt, iPhiTrk);
}

//=== Check if specified cell is merged with its 2x2 neighbours into a single cell,
//=== as it is in low Pt region.

bool HTrphi::mergedCell(unsigned int iQoverPtBin, unsigned int jPhiTrkBin) const {

  bool merge = false;

  if (enableMerge2x2_) {
    unsigned int i = iQoverPtBin;
    unsigned int j = jPhiTrkBin;

    // Number of bins to merge on each q/Pt side of array (+ve & -ve charge) must be even number.
    float fMergeBins = (maxAbsQoverPtAxis_ - minInvPtToMerge2x2_)/(2.*binSizeQoverPtAxis_);
    unsigned int numQoverPtBinsToMerge = 2 * min( (unsigned int)(floor(fMergeBins)), (nBinsQoverPtAxis_/4) );
    unsigned int iB = (nBinsQoverPtAxis_ - 1) - i; // Count backwards across array.
    if (min(i, iB) < numQoverPtBinsToMerge) merge = true;
  }

  return merge;
}

//=== Calculate maximum |gradient| that any stub's line across this HT array could have, so can check it doesn't exceed 1.

float HTrphi::calcMaxLineGradArray() const {
  // Get max. |gradient| possible in this HT array.
  float gradOuter = fabs(invPtToDphi_ * (settings_->trackerOuterRadius() - chosenRofPhi_));
  float gradInner = fabs(invPtToDphi_ * (settings_->trackerInnerRadius() - chosenRofPhi_));
  float maxGrad = max(gradOuter, gradInner);
  // Convert it to units of bin width.
  maxGrad *= binSizeQoverPtAxis_/binSizePhiTrkAxis_;
  return maxGrad;
}

//=== If requested, kill those tracks in this sector that can't be read out during the time-multiplexed period, because
//=== the HT has associated too many stubs to tracks.

vector<L1track2D> HTrphi::killTracksBusySec(const vector<L1track2D>& tracks) const {

  vector<L1track2D> outTracks;

  if (busySectorKill_) {

    unsigned int nStubsOut          = 0;     // #stubs assigned to tracks in this sector.
    unsigned int nStubsOutPosCharge = 0;     // #stubs assigned to +ve charged tracks in this sector.
    unsigned int nStubsOutNegCharge = 0;     // #stubs assigned to -ve charged tracks in this sector.
    vector<unsigned int> nStubsOutInRange(busySectorMbinRanges_.size(), 0); // #stubs assigned to each m bin range in this sector. 

    for (const L1track2D& trk : tracks) {
      bool keep = true;
      unsigned int nStubs = trk.getNumStubs(); // #stubs on this track.
      if (busySectorUseMbinRanges_) { // Are tracks from different m bin ranges output seperately to increase bandwidth?
	unsigned int mBinRange = this->getMbinRange(trk); // Which m bin range is this track in?
	nStubsOutInRange[ mBinRange ] += nStubs;
	if (nStubsOutInRange[ mBinRange ] > busySectorNumStubs_) keep = false;
      } else {
	if (busySectorEachCharge_) { // Are +ve and -ve charged tracks output on separate optical links to increase bandwidth?
	  if (trk.qOverPt() > 0) {
	    nStubsOutPosCharge += nStubs;
	    if (nStubsOutPosCharge > busySectorNumStubs_) keep = false;
	  } else {
	    nStubsOutNegCharge += nStubs;
	    if (nStubsOutNegCharge > busySectorNumStubs_) keep = false;
	  }
	} else {
	  nStubsOut += nStubs;
	  if (nStubsOut > busySectorNumStubs_) keep = false;
	}
      }

      if (keep) outTracks.push_back(trk);
    }

  } else {

    outTracks = tracks;

  }

  return outTracks;
}

//=== Define the order in which the hardware processes rows of the HT array when it outputs track candidates.
//=== Currently corresponds to highest Pt tracks first. 
//=== If two tracks have the same Pt, the -ve charge one is output before the +ve charge one.
  
vector<unsigned int> HTrphi::rowOrder(unsigned int numRows) const {
  vector<unsigned int> iOrder;

  // Logic slightly different depending on whether HT array has even or odd number of rows.
  const bool oddNumRows = (numRows%2 == 1);

  // This selects middle rows first before moving to exterior ones.
  if (oddNumRows) {
    unsigned int middleRow = (numRows - 1)/2;
    iOrder.push_back(middleRow);
    for (unsigned int i = 1; i < (numRows - 1)/2; i++) {
      iOrder.push_back(middleRow - i); // -ve charge
      iOrder.push_back(middleRow + i); // +ve charge
    }
  } else {
    unsigned int startRowPos = numRows/2;
    unsigned int startRowNeg = startRowPos - 1;
    for (unsigned int i = 0; i < numRows/2; i++) {
      iOrder.push_back(startRowNeg - i); // -ve charge
      iOrder.push_back(startRowPos + i); // +ve charge
    }
  }

  return iOrder;
}
