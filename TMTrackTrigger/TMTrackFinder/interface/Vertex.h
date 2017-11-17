#ifndef __Vertex_H__
#define __Vertex_H__

#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"

#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"

#include <vector>

class Stub;

class Vertex {

public:
  // Fill useful info about tracking particle.
  Vertex(){}

  Vertex(double vz):vz_(vz){}

  ~Vertex(){}

  /// Tracking Particles in vertex    
  const vector<TP>& tracks()    const { return    tracks_;    }
  /// Number of tracks originating from this vertex
  unsigned int      numTracks() const { return  tracks_.size();}
  /// Assign TP to this vertex
  void              insert(TP tp)     { tracks_.push_back(tp);}
  /// Compute vertex parameters
  void              computeParameters();
  /// Sum ot fitted tracks transverse momentum
  double            pT()        const {return pT_;}
  /// Vertex z0 position
  double            z0()        const {return z0_;}
  /// Vertex z0 width
  double            z0width()   const {return z0width_;}
  /// Vertex z position
  double            vz()        const {return vz_;}
  /// Vertec MET
  double            met()       const {return met_;}
  /// Number tracks with eta lower than 0.9
  int               noTracksBelow15eta() const;
  /// Number tracks with eta bigger than 0.9
  int               noTracksAbove15eta() const;
  /// Percentage of tracks with eta lower than 0.9
  float            percentageTracksGoodEta() const;
  /// No. layers
  vector<unsigned int> numLayersTracks() const;
  /// No PS layers
  vector<unsigned int> numPSLayersTracks() const;
  // Layers ID
  vector<int> layerIDTracks() const;

    
    
private:

  double            vz_;
  double            z0_;
  double            z0width_;
  double            pT_;
  double            met_;
  double            metX_;
  double            metY_;


  std::vector<TP>   tracks_;

  
};

#endif
