#ifndef __TMTRACKPRODUCER_H__
#define __TMTRACKPRODUCER_H__

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"


#include <vector>
#include <map>
#include <string>

using namespace std;

class Settings;
class Histos;
class TrackFitGeneric;

class TMTrackProducer : public edm::EDProducer {

public:
  explicit TMTrackProducer(const edm::ParameterSet&);	
  ~TMTrackProducer(){}

private:

  typedef std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > TTTrackCollection;
  typedef std::vector< const L1fittedTrack* > FitTrackCollection;
  typedef std::vector< L1fittedTrack > TrackCollection;
  typedef std::vector< TP > TrackParticleCollection;

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

private:
  const edm::EDGetTokenT<TrackingParticleCollection> tpInputTag;
  const edm::EDGetTokenT<DetSetVec> stubInputTag;
  const edm::EDGetTokenT<TTStubAssMap> stubTruthInputTag;
  const edm::EDGetTokenT<TTClusterAssMap> clusterTruthInputTag;

  Settings *settings_;
  Histos   *hists_;
  map<string, TrackFitGeneric*> fitterWorkerMap_;
};
#endif

