#ifndef __TrkRZfilter_H__
#define __TrkRZfilter_H__

#include "TMTrackTrigger/TMTrackFinder/interface/L1track2D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/HTrphi.h"

#include <vector>

class Settings;
class Stub;

using namespace std;

//=== This class runs filters in track candidates previously found by the r-phi Hough transform
//=== e.g. Filters requiring the stubs to be consistent with a straight line in the r-z plane.
//===
//=== The filtering removes inconsistent stubs from the track candidates, & also kills some track candidates
//=== altogether if the filter leaves them with too few stubs.
//===
//=== Some r-z filters may also add an estimate of the r-z helix parameters to the selected track candidates.
//===
//=== It does NOT contain filters such as the bend filter, which are so simple that the firmware can run them 
//=== INSIDE the r-phi HT. Simple filters of this kind are in class HTcell.

class TrkRZfilter {

public:

  TrkRZfilter() {}
  ~TrkRZfilter() {}
  
  struct SortStubsInLayer{
    inline bool operator() (const Stub *stub1, const Stub *stub2){
      return(fabs(stub1->layerId()) < fabs(stub2->layerId()));
    }
  };
  // Initialize configuration parameters, and note sector number, eta range covered by sector and phi coordinate of its centre.
  void init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg, 
	    float etaMinSector, float etaMaxSector, float phiCentreSector);

  // Filters track candidates (found by the r-phi Hough transform), removing inconsistent stubs from the tracks, 
  // also killing some of the tracks altogether if they are left with too few stubs.
  // Also adds an estimate of r-z helix parameters to the selected track objects, if the filters used provide this.
  vector<L1track2D> filterTracks(const vector<L1track2D>& tracks);

  //=== Extra information about each track input to filter. (Only use after you have first called filterTracks).

  // Number of seed combinations considered by the ZTrk Filter for each input track.
  vector<unsigned int> numZtrkSeedCombsPerTrk() const {return numZtrkSeedCombsPerTrk_; }
  // Number of seed combinations considered by the Seed Filter for each input track.
  vector<unsigned int> numSeedCombsPerTrk() const {return numSeedCombsPerTrk_; }
  vector<unsigned int> numGoodSeedCombsPerTrk() const {return numGoodSeedCombsPerTrk_; } // Only counts seeds compatible with beam-spot.


private:

  //--- Filters returning filtered stubs based on input ones.
 
  // Produce a filtered collection of stubs from the input ones (on original track) that all have consistent rapidity
  vector<const Stub*> etaFilter ( const vector<const Stub*>& stubs, float trkQoverPt ) const;
  // Produce a filtered collection of stubs from the input ones (on original track) that all have consistent zR.
  vector<const Stub*> zTrkFilter (const vector<const Stub*>& stubs, float trkQoverPt );
  // Produce a filtered collection of stubs from the input ones (on original track)that are consistent with a straight line in r-z using tracklet algo.
  vector<const Stub*> seedFilter (const vector<const Stub*>& stubs, float trkQoverPt, bool print );

private:

  //=== Configuration parameters

  const Settings* settings_;

  unsigned int iPhiSec_; // Sector number.
  unsigned int iEtaReg_;
  float etaMinSector_; // rapidity range of this sector.
  float etaMaxSector_;
  float chosenRofZ_;   // Radius used to defined zTrkMinSector and zTrkMaxSector.
  float zTrkMinSector_; // corresponding range of this sector specified as z coordinate of track at given radius.
  float zTrkMaxSector_; 
  float phiCentreSector_; // phi coordinate of its centre.

  // Use filter in each HT cell using only stubs which have consistent rapidity?
  bool   useEtaFilter_;
  // Use filter in each HT cell using only stubs which have consistent zR
  bool   useZTrkFilter_;
  // Filter stubs in cell using a tracklet-like algorithm
  bool   useSeedFilter_;

  // Options for Ztrk filter
  float chosenRofZFilter_;

  // Options for Seed filter.
  float seedResolution_;
  bool  keepAllSeed_;

  // Number of seed combinations considered by the ZTrk Filter, for each input track.
  vector<unsigned int>  numZtrkSeedCombsPerTrk_;

  // Number of seed combinations considered by the Seed Filter, for each input track.
  vector<unsigned int>  numSeedCombsPerTrk_;
  vector<unsigned int>  numGoodSeedCombsPerTrk_;
  unsigned int maxSeedCombinations_;
  unsigned int maxGoodSeedCombinations_;
  unsigned int maxSeedsPerStub_;
  bool         zTrkSectorCheck_;

  float beamWindowZ_; // Assumed length of beam spot in z.

  // Track (z0, tan_lambda) estimate from r-z filter if available.
  bool   estValid_;
  float  estZ0_;
  float  estTanLambda_;

  // For debugging
  unsigned int minNumMatchLayers_;
};
#endif

