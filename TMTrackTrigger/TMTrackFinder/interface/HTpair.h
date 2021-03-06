#ifndef __HTpair_H__
#define __HTpair_H__

#include "TMTrackTrigger/TMTrackFinder/interface/HTrphi.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TrkRZfilter.h"
#include "TMTrackTrigger/TMTrackFinder/interface/KillDupTrks.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"

#include "boost/numeric/ublas/matrix.hpp"
#include <vector>
#include <utility>

using  boost::numeric::ublas::matrix;

class Settings;
class Stub;
class TP;

using namespace std;

//=== This is a linked pair of r-phi and r-z Hough Transform arrays for a single (eta,phi) sector.
//=== It allows track reconstruction in 3D, and gives access to these 3D tracks.
//===
//=== The use of the r-z HT is optional. If it is not used, 3D tracks are still produced but with
//=== only a rough estimate of their r-z track parameters based on the eta sectors they are assigned to.
//===
//=== If requested, it will also run some track filters (e.g. r-z filters) on the track candidates found
//=== by the r-phi HT. These not only clean up the tracks, but may also provide an estimate of the r-z track parameters.
//=== 
//=== To create 3D tracks, call the sequence init(), store(), end(), make3Dtracks().

class HTpair {

public:
  
  HTpair() {}
  ~HTpair(){}

  //=== Main routines to do tracking.

  // Initialization
  void init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg, 
      float etaMinSector, float etaMaxSector, float phiCentreSector);

  // Add stub to r-phi HT array.
  // If eta subsectors are being used within each sector, specify which ones the stub is compatible with.
  void store( const Stub* stub, const vector<bool>& inEtaSubSecs);

  // Causes r-phi HT to search for tracks. 
  void end();

  // Create list of all 3D track candidates in "vecTracks3D_", obtained either by combining r-phi and r-z HTs, 
  // or by just using r-phi HT and guessing r-z track parameters from centre of sector.
  // Optionally runs track filters (e.g. r-z filter) after r-phi HT if requested, which may improve r-z track parameter estimate. 
  // Optionally also runs duplicate track removal on the final 3D track collection.
   void make3Dtracks();

  //=== Access to intermediate tracking.

  // Access to filled r-phi Hough transform array. (The r-z one is not permanently stored).
  const HTrphi& getRphiHT()    const {return htArrayRphi_;}
  // Ditto, but allowing the r-phi Hough transform to be changed. (Used by class MuxHToutputs).
  HTrphi& getRphiHT_nonConst()       {return htArrayRphi_;}

  // Access to track r-z filters.
  const TrkRZfilter& getRZfilters() const {return rzFilters_;}

  //=== Info about 3D track candidates found.

  // Get list of all 3D track candidates found, obtained either by combining r-phi and r-z HTs, 
  // or by just using r-phi HT and guessing r-z track parameters from centre of sector.
  // Optionally runs track filters (e.g. r-z filter) after r-phi HT if requested, which may improve r-z track parameter estimate.
  // Each L1track3D object gives access to stubs on each track and helix parameters, 
  // and also to the associated truth tracking particle.
  const vector<L1track3D>& trackCands3D() const {return vecTracks3D_;}

  // For debugging, return 3D tracks after r-z filter/HT, but before any duplicate removal run on 3D tracks.
  const vector<L1track3D>& trackCandsRZfilter() const {return vecTracksRZfilter_;}

  // Number of 3D track candidates found.
  unsigned int numTrackCands3D() const {return vecTracks3D_.size();}

  // Get number of stubs assigned to 3D track candidates.
   unsigned int   numStubsOnTrackCands3D() const;

  // Get all 3D track candidates associated to the given tracking particle.
  // (If the vector is empty, then the tracking particle was not reconstructed in this sector).
   vector<const L1track3D*> assocTrackCands3D(const TP& tp) const;

  // Get number of 3D track candidates associated to the given tracking particle.
  unsigned int numAssocTrackCands3D(const TP& tp) const {return this->assocTrackCands3D( tp ).size();}

  // Get the number of r-phi HT cells that a given set of 3D track candidates came from.
   unsigned int   numRphiCells(const vector<const L1track3D*>& trk3D) const;

private:

  // Configuration parameters
  const Settings* settings_;
  bool  enableRzHT_;       // Use the r-z Hough Transform?
  unsigned int iPhiSec_;   // Sector number.
  unsigned int iEtaReg_;
  float etaMinSector_;     // Range of eta sector
  float etaMaxSector_;     // Range of eta sector
  float phiCentreSector_;  // Phi angle of centre of this (eta,phi) sector.

  // r-phi Hough transform
  HTrphi htArrayRphi_; 
  // Don't bother storing r-z HT array, as would need one for each r-phi track cand, taking too much memory.

  // Track filter(s), such as r-z filters, run after the r-phi Hough transform.
  TrkRZfilter rzFilters_;

  // Contains algorithm used for duplicate track removal.
  KillDupTrks<L1track3D> killDupTrks_;

  // List of all found 3D track candidates and their associated properties.
  vector<L1track3D> vecTracks3D_;

  vector<L1track3D> vecTracksRZfilter_;

  // Since r-z HT array is not stored, store instead this variable used to debug it.
  vector<float> fracCellsWithNoNeighboursRz_;
};
#endif

