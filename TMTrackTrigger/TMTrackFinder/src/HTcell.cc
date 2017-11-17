#include <TMTrackTrigger/TMTrackFinder/interface/HTcell.h>
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

//=== Initialization with cfg params, boolean indicating if this is r-phi or r-z HT, 
//=== rapidity range of current sector, and estimated q/Pt of cell,
//=== and (if called from r-phi HT) the bin number of the cell along the q/Pt axis of the r-phi HT array.

void HTcell::init(const Settings* settings, bool isRphiHT, unsigned int iPhiSec, unsigned int iEtaReg,
		  float etaMinSector, float etaMaxSector, float qOverPt, unsigned int ibin_qOverPt) {
  settings_ = settings;

  // Note if if this r-phi or r-z HT.
  isRphiHT_ = isRphiHT;

  // Sector number
  iPhiSec_ = iPhiSec;
  iEtaReg_ = iEtaReg;

  // Note track q/Pt. 
  // In this case of an r-phi HT, each cell corresponds to a unique q/Pt.
  // In the case of an r-z HT, it is assumed that we know q/Pt from previously run r-phi HT.
  qOverPtCell_ = qOverPt;
  // Note bin number of cell along q/Pt axis of r-phi HT array. (Not used if r-z HT).
  ibin_qOverPt_ = ibin_qOverPt;
  // Rapidity range of sector.
  etaMinSector_ = etaMinSector;
  etaMaxSector_ = etaMaxSector;

  // The following only relevant to r-phi Hough transform.
  if (isRphiHT_) {
    invPtToDphi_   = settings->invPtToDphi();  // B*c/2E11

    // Use filter in each HT cell using only stubs which have consistent bend?
    useBendFilter_ = settings->useBendFilter();
  }

  // A filter is used each HT cell, which prevents more than the specified number of stubs being stored in the cell. (Reflecting memory limit of hardware).
  maxStubsInCell_ = settings->maxStubsInCell();

  // Note if daisy-chain firmware is in use, together with digitized stubs.
  daisyChainFirmware_ = (settings->firmwareType() >= 1 && settings->firmwareType() < 99) && settings->enableDigitize();

  // Check if subsectors are being used within each sector. These are only ever used for r-phi HT.
  numSubSecs_ = isRphiHT_   ?   settings->numSubSecsEta()  :  1;
}

//=== Termination. Search for track in this HT cell etc.

void HTcell::end(){
  // Produce list of filtered stubs by applying all requested filters (e.g. on stub bend).
  // (If no filters are requested, then filtered & unfiltered stub collections will be identical).

  // N.B. Other filters,  such as the r-z filters, which the firmware runs after the HT because they are too slow within it,
  // are not defined here, but instead inside class TrkFilterAfterRphiHT.
  // if ( vStubs_.size() > 2 ) {
  //   std::cout << "In HTcell::end()" << std::endl;
  //   std::cout << "Number of stubs coming in : " << vStubs_.size() << std::endl;    
  // }

  vFilteredStubs_ = vStubs_;
  // The bend filter is only relevant to r-phi Hough transform.
  if (isRphiHT_) {
    if (useBendFilter_) vFilteredStubs_ = this->bendFilter(vFilteredStubs_);
  }

  // Prevent too many stubs being stored in a single HT cell if requested (to reflect hardware memory limits).
  // N.B. This MUST be the last filter applied.
  if (maxStubsInCell_ <= 99) vFilteredStubs_ = this->maxStubCountFilter(vFilteredStubs_);

  // Calculate the number of layers the filtered stubs in this cell are in.
  numFilteredLayersInCell_ = this->calcNumFilteredLayers();

  if (numSubSecs_ > 1) { 
    // If using subsectors within each sector, calculate the number of layers the filters stubs in this cell are in,
    // when one considers only the subset of the stubs within each subsector.
    // Look for the "best" subsector.
    numFilteredLayersInCellBestSubSec_ = 0;
    for (unsigned int i = 0; i < numSubSecs_; i++) {
      unsigned int numLaySubSec = this->calcNumFilteredLayers(i);
      numFilteredLayersInCellBestSubSec_ = max(numFilteredLayersInCellBestSubSec_, numLaySubSec);
    }
  } else {
    // If only 1 sub-sector, then subsector and sector are identical.
    numFilteredLayersInCellBestSubSec_ = numFilteredLayersInCell_;
  }
  // if ( vStubs_.size() > 2 ) {
  //    std::cout << "Number of filtered stubs going out : " << vFilteredStubs_.size() << " " << numFilteredLayersInCell_  << std::endl;    
  // }
}

// Calculate how many tracker layers the filter stubs in this cell are in, when only the subset of those stubs
// that are in the specified subsector are counted.

unsigned int HTcell::calcNumFilteredLayers(unsigned int iSubSec) const {
  vector<const Stub*> stubsInSubSec;
  for (const Stub* s : vFilteredStubs_) {
    const vector<bool>& inSubSec = subSectors_.at(s); // Find out which subsectors this stub is in.
    if (inSubSec[iSubSec]) stubsInSubSec.push_back(s);
  }
  return Utility::countLayers( settings_, stubsInSubSec );
}


//=== Produce a filtered collection of stubs in this cell that all have consistent bend.
//=== Only called for r-phi Hough transform.

vector<const Stub*> HTcell::bendFilter( const vector<const Stub*>& stubs ) const {

  // Create bend-filtered stub collection.
  vector<const Stub*> filteredStubs;
  for (const Stub* s : stubs) {

    // Require stub bend to be consistent with q/Pt of this cell.

    if (daisyChainFirmware_) {
      // Daisy chain firmware doesn't have access to variables needed to calculate dphi of stub,
      // but instead knows integer range of q/Pt bins that stub bend is compatible with, so use these.
      if (s->min_qOverPt_bin() <= ibin_qOverPt_ && ibin_qOverPt_ <= s->max_qOverPt_bin() )  filteredStubs.push_back(s);
    } else {
      // Systolic array & 2-c-bin firmware do hace access to stub dphi, so can use it.
      // Predict track bend angle based on q/Pt of this HT cell and radius of stub.
      float predictedDphi = this->dphi( s->r() );
      // Require reconstructed and predicted values of this quantity to be consistent within estimated resolution. 
      if (fabs(s->dphi() - predictedDphi) < s->dphiRes()) filteredStubs.push_back(s);
    }
  }

  return filteredStubs;
}

//=== Filter stubs so as to prevent more than specified number of stubs being stored in one cell.
//=== This reflects finite memory of hardware.

vector<const Stub*> HTcell::maxStubCountFilter( const vector<const Stub*>& stubs ) const {
  vector<const Stub*> filteredStubs;
  // If there are too many stubs in a cell, the hardware keeps (maxStubsInCell - 1) of the first stubs in the list
  // plus the last stub.  
  if (stubs.size() > maxStubsInCell_) {
    for (unsigned int i = 0; i < maxStubsInCell_ - 1; i++) { // first stubs
      filteredStubs.push_back(stubs[i]);
    }
    filteredStubs.push_back(stubs[ stubs.size() - 1] ); // plus last stub
  } else {
    filteredStubs = stubs;
  }
  return filteredStubs;
}
