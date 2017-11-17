#ifndef __INPUTDATA_H__
#define __INPUTDATA_H__

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Vertex.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <vector>

class Settings;

using namespace std;

//=== Unpacks stub & tracking particle (truth) data into user-friendlier format in Stub & TP classes.
//=== Also makes B-field available to Settings class.

class InputData {

public:
  
  InputData(const edm::Event& iEvent, const edm::EventSetup& iSetup, Settings* settings,
  const edm::EDGetTokenT<TrackingParticleCollection> tpInputTag,
  const edm::EDGetTokenT<DetSetVec> stubInputTag,
  const edm::EDGetTokenT<TTStubAssMap> stubTruthInputTag,
  const edm::EDGetTokenT<TTClusterAssMap> clusterTruthInputTag
   );

  // Sort Tracking Particles by vertex z position
  struct SortVertexByZ0{
    inline bool operator() (const Vertex vertex0, const Vertex vertex1){
      return(vertex0.z0() < vertex1.z0());
    }
  };

  // Get tracking particles
  const vector<TP>&          getTPs()      const {return vTPs_;}
  // Get stubs that would be output by the front-end readout electronics 
  const vector<const Stub*>& getStubs()    const {return vStubs_;}
  /// Get Primary vertex information
  const Vertex&              getPrimaryVertex()       const {return vertex_;}
  /// Get PileUp Vertices information
  const vector<Vertex>&      getPileUpVertices()      const {return vertices_;}
  /// Get reconstructable pile-up vertices information
  const vector<Vertex>&      getRecoPileUpVertices()      const {return recoVertices_;}
  /// Generated MET
  const float                GenMET()          const { return genMET_;}
  const float                GenMET_PU()          const { return genMET_PU_;}

  //--- of minor importance ...

  // Get number of stubs prior to applying tighted front-end readout electronics cuts specified in section StubCuts of Analyze_Defaults_cfi.py. (Only used to measure the efficiency of these cuts).
  const vector<Stub>&        getAllStubs() const {return vAllStubs_;}

private:
  // const edm::EDGetTokenT<TrackingParticleCollection> inputTag;

  // Can optionally be used to sort stubs by bend.
  struct SortStubsInBend {
     inline bool operator() (const Stub* stub1, const Stub* stub2) {
        return(fabs(stub1->bend()) < fabs(stub2->bend()));
     }
  };

private:

  vector<TP> vTPs_; // tracking particles
  Vertex vertex_;
  vector<Vertex> vertices_;
  vector<Vertex> recoVertices_;
  vector<const Stub*> vStubs_; // stubs that would be output by the front-end readout electronics.

  //--- of minor importance ...

  vector<Stub> vAllStubs_; // all stubs, even those that would fail any tightened front-end readout electronic cuts specified in section StubCuts of Analyze_Defaults_cfi.py. (Only used to measure the efficiency of these cuts).
  float genMET_;
  float genMET_PU_;
};
#endif

