#include "TMTrackTrigger/TMTrackFinder/interface/HTpair.h"
#include "TMTrackTrigger/TMTrackFinder/interface/HTrz.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track2D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"

#include <iostream>
#include <unordered_set>

using namespace std;

//=== Initialization

void HTpair::init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg,
		  float etaMinSector, float etaMaxSector, float phiCentreSector) {

  // Store config params.
  settings_        = settings;
  enableRzHT_      = settings->enableRzHT();   // Use the r-z Hough Transform?
  iPhiSec_         = iPhiSec;                  // Sector number
  iEtaReg_         = iEtaReg;
  etaMinSector_    = etaMinSector;             // Range of eta sector
  etaMaxSector_    = etaMaxSector;             // Range of eta sector
  phiCentreSector_ = phiCentreSector;          // Centre of phi sector

  // Initialize r-phi Hough transform array.
  htArrayRphi_.init(settings_, iPhiSec_, iEtaReg_, etaMinSector_, etaMaxSector_, phiCentreSector_);

  // Initialize any track filters (e.g. r-z) run after the r-phi Hough transform.
  rzFilters_.init(settings_, iPhiSec_, iEtaReg_, etaMinSector_, etaMaxSector_, phiCentreSector_);  

  //--- Option for duplicate track removal on collection of L1track3D produced after running all track-finding steps.
  unsigned int dupTrkAlgRzSeg = settings->dupTrkAlgRzSeg();
  killDupTrks_.init(settings, dupTrkAlgRzSeg);
}

//=== Add stub to r-phi HT array.
//== If eta subsectors are being used within each sector, specify which ones the stub is compatible with.

void HTpair::store( const Stub* stub, const vector<bool>& inEtaSubSecs) {
  htArrayRphi_.store(stub, inEtaSubSecs);
}

//=== Causes r-phi HT to search for tracks. 

void HTpair::end() {
  htArrayRphi_.end();
}

//=== Create list of all 3D track candidates in "vecTracks3D_", obtained either by combining r-phi and r-z HTs, 
//=== or by just using r-phi HT and guessing r-z track parameters from centre of sector.
//=== Optionally runs track filters (e.g. r-z filter) after r-phi HT if requested, which may improve r-z track parameter estimate. 
//=== Optionally also runs duplicate track removal on the final 3D track collection.

void HTpair::make3Dtracks() {

  vector<L1track3D> vecTracks3D;

  // Get tracks found by r-phi HT.
  const vector<L1track2D>& vecTracksRphi = htArrayRphi_.trackCands2D();
  // if ( htArrayRphi_.trackCands2D().size() > 0 ) {
  //   std::cout << "In HTpair make3Dtracks" << std::endl;
  //   std::cout << "Number of 2d tracks : " << vecTracksRphi.size() << std::endl;    
  // }


  // Run requested track filters (e.g. r-z filters) to clean up the tracks by removing inconsistent stubs, and killing
  // some tracks altogether if this procedure leaves them with too few stubs.
  // The r-z filters may also add a good estimate of the r-z helix parameters to each track.
  const vector<L1track2D>& vecTracksRphiFilt = rzFilters_.filterTracks(vecTracksRphi);
  // if ( htArrayRphi_.trackCands2D().size() > 0 ) {
  //   std::cout << "After rz filter : " << vecTracksRphiFilt.size() << std::endl;
  // }

  // Loop over track candidates found by r-phi HT.

  for (const L1track2D& trkRphi : vecTracksRphiFilt) {
    const vector<const Stub*>& stubsOnTrkRphi = trkRphi.getStubs(); // stubs assigned to track 

    if (enableRzHT_) {

      // --- Run r-z HT on stubs assigned to each track by r-phi HT.

      float qOverPt = trkRphi.getHelix2D().first; // Estimated q/Pt of this track from r-phi HT.
      HTrz  htArrayRz; 
      htArrayRz.init(settings_, iPhiSec_, iEtaReg_, etaMinSector_, etaMaxSector_, qOverPt);
      // Loop over stubs on each track and pass them to r-z HT.
      for (const Stub* s : stubsOnTrkRphi) {
	htArrayRz.store( s );
      }
      htArrayRz.end();

      // Loop over tracks found by r-z HT obtained using stubs on tracks found by r-phi HT..
      const vector<L1track2D>& trackCandsRz = htArrayRz.trackCands2D();
      for (const L1track2D& trkRz : trackCandsRz) {

	// Create 3D track (N.B. Set stubs equal to those on r-z track, which are filtered with respect to those on the r-phi track by the r-z HT).
	// The L1track3D class automatically finds the associated truth Tracking Particle (if any).
	L1track3D trk3D(settings_, trkRz.getStubs(), 
			trkRphi.getCellLocation(), trkRphi.getHelix2D(),
			trkRz.getCellLocation()  , trkRz.getHelix2D(),
			iPhiSec_                 , iEtaReg_);
        // Add to list of stored 3D tracks.
	// N.B. The stubs on each r-phi track can give track cands in multiple cells of r-z HT,
	// in which case, several 3D tracks may be created for each r-phi track.
        vecTracks3D.push_back( trk3D );
      }

    } else {

      // --- Don't run  r-z HT. 

      // Estimate r-z track parameters.

      // If an r-z filter was run after the r-phi HT, it might provide an estimate of track's (z0,eta).
      float z0, tan_lambda;
      bool validRZ = trkRphi.getTrkEstZ0andTanLam(z0, tan_lambda);

      // If not, then estimate r-z track parameters from centre of rapidity sector.
      if ( ! validRZ ) {
        z0 = 0.;
        // float etaCentreSector = 0.5*(etaMinSector_ + etaMaxSector_);
        // float theta = 2. * atan(exp(-etaCentreSector));
        // tan_lambda = 1./tan(theta);
        tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector_))) + 1/tan(2*atan(exp(-etaMaxSector_))));
      }

      pair<float, float> helixRz(z0, tan_lambda);

      // Set cell location in r-z HT to crazy large number.
      pair<unsigned int, unsigned int> cellRz(99999, 99999);

      // Create 3D track (N.B. Set stubs equal to those on r-phi track, since no r-z HT was run to filter them).
      // The L1track3D class automatically finds the associated truth Tracking Particle (if any).
      L1track3D trk3D(settings_, stubsOnTrkRphi, 
		      trkRphi.getCellLocation(), trkRphi.getHelix2D(),
		      cellRz                   , helixRz,
		      iPhiSec_                 , iEtaReg_);
      // Add to list of stored 3D tracks.
      vecTracks3D.push_back( trk3D );
    }
  }
  // if ( htArrayRphi_.trackCands2D().size() > 0 ) {
  //   std::cout << "3D tracks after rz : " << vecTracks3D.size() << std::endl;
  // }

  // Make copy of tracks before duplicate removal for debug purposes.
  vecTracksRZfilter_ = vecTracks3D;

  // Run duplicate track removal on all the 3D tracks found in this sector & store final 3D track collection.
  vecTracks3D_ = killDupTrks_.filter( vecTracks3D );
}


//=== Get number of stubs assigned to 3D track candidates.

unsigned int HTpair::numStubsOnTrackCands3D() const {

  unsigned int nStubs = 0;

  // Loop over track candidates
  for (const L1track3D& trk : vecTracks3D_) {
    nStubs += trk.getStubs().size();
  }

  return nStubs;
}

//=== Get all 3D track candidates that were associated to the given tracking particle.
//=== (If the vector is empty, then the tracking particle was not reconstructed in this sector).

vector<const L1track3D*> HTpair::assocTrackCands3D(const TP& tp) const {

  vector<const L1track3D*> assocRecoTrk;

  // Loop over track candidates, looking for those associated to given TP.
  for (const L1track3D& trk : vecTracks3D_) {
    if (trk.getMatchedTP() != nullptr) {
      if (trk.getMatchedTP()->index() == tp.index()) assocRecoTrk.push_back(&trk); 
    }
  }

  return assocRecoTrk;
}

//=== Get the number of r-phi cells that a given set of 3D track candidates came from.

unsigned int HTpair::numRphiCells(const vector<const L1track3D*>& vectorTrk3D) const {

  // Unique set will only store one copy of something, even if multiple copies are added to it.
  unordered_set<unsigned int> rphiCellAddresses;
  for (const L1track3D* trk : vectorTrk3D) {
    // Store location in r-phi HT of cell that provided stubs for this 3D track.
    pair<unsigned int, unsigned int> cellRphi = trk->getCellLocationRphi();
    // Unfortunately, can't store "std::pair" inside "unordered_set", so must encode cell address as a single integer.
    rphiCellAddresses.insert( 10000*cellRphi.first + cellRphi.second); // This is safe providing HT array has dimensions < 10000 x 10000.
  }
  return rphiCellAddresses.size(); // This is the number of r-phi cells with a unique address.
}
