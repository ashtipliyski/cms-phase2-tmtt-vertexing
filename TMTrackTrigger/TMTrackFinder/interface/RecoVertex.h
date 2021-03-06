#ifndef __RecoVertex_H__
#define __RecoVertex_H__

#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"

#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"

#include <vector>

class RecoVertex {

public:
  // Fill useful info about tracking particle.
  RecoVertex(){z0_ = -999.; pT_ = -9999.; met_ = -999.;}
  ~RecoVertex(){}

  /// Tracking Particles in vertex    
  vector<const L1fittedTrack*> tracks()    const { return    tracks_;    }
  /// Tracking Particles in vertex    
  std::set< const TP* > trueTracks()    const { return    trueTracks_;    }
  /// Number of tracks originating from this vertex
  unsigned int      numTracks() const { return  tracks_.size();}
  /// Number of true particles assigned to this vertex
  unsigned int      numTrueTracks() const {return trueTracks_.size();}
  /// Assign fitted track to this vertex
  void              insert(const L1fittedTrack* fitTrack)     { tracks_.push_back(fitTrack); if(fitTrack->getMatchedTP()!= nullptr and fitTrack->getMatchedTP()->physicsCollision()) trueTracks_.insert(fitTrack->getMatchedTP());}
  /// Compute vertex parameters
  void              computeParameters();
  /// Set z0 position
  void              setZ(double z)    {z0_ = z;}
  /// Sum ot fitted tracks transverse momentum
  double            pT()        const {return pT_;}
  /// Vertex z0 position
  double            z0()        const {return z0_;}
  /// Vertex z0 width
  double            z0width()   const {return z0width_;}
  /// Clear track vector
  void              clear()     {tracks_.clear(); trueTracks_.clear();}
  /// True if primary vertex
  bool              PrimaryVertex() const { return pv_;}
  /// Set primary vertex tag
  void              isPrimary(bool is) { pv_ = is;}
  /// Contain high-pT track?
  bool              hasHighPt() const { return highPt_;}
  /// Number of high-pT tracks (pT > 10 GeV)
  unsigned int      numHighPtTracks() const { return numHighPtTracks_;}
  /// Vertex MET
  double            met()       const {return met_;}
  /// chi2 of the tracks
  vector<float>     tracksChi2dof() const;

private:

  double            z0_;
  double            z0width_;
  double            pT_;
  double            met_;
  double            metX_;
  double            metY_;

  std::vector<const L1fittedTrack*>   tracks_;
  std::set< const TP* >   trueTracks_;
  bool              pv_;
  bool              highPt_;
  unsigned int      numHighPtTracks_;
  
};

#endif
