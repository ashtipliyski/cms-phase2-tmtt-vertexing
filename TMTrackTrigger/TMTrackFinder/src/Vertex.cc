#include "TMTrackTrigger/TMTrackFinder/interface/Vertex.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"


void Vertex::computeParameters(){
	pT_ = 0.;
	z0_ = 0.;
	met_ = 0.;
	metX_ = 0.;
	metY_ = 0.;
	float z0square = 0.;
	for(TP track : tracks_){
		pT_ += track.pt();
		z0_ += track.z0();
		z0square += track.z0()*track.z0();
		metX_ += track.pt()*cos(track.phi0());
		metY_ += track.pt()*sin(track.phi0());
	}
	met_ = sqrt(metX_*metX_ + metY_*metY_);
	z0_ /= tracks_.size();
	z0square /= tracks_.size();
	z0width_ = sqrt(fabs(z0_*z0_ - z0square));
}

int Vertex::noTracksBelow15eta() const{
  int noTracks15 = 0;
  for(TP track : tracks_){
    if((track.eta()<0.9) && (track.eta()>-0.9))
      noTracks15 = noTracks15 + 1;
  }
   return noTracks15;
}

int Vertex::noTracksAbove15eta() const{
  int nTracks15 = 0;
  for(TP track : tracks_){
    if((track.eta()>0.9) || (track.eta()<-0.9))
      nTracks15 = nTracks15 + 1;
  }
   return nTracks15;
}

float Vertex::percentageTracksGoodEta() const{
    float percent = 0.0;
    unsigned int noTrk = 0;

    for(TP track : tracks_){
        if((track.eta()<0.9) && (track.eta()>-0.9))
            noTrk = noTrk + 1;
    }

    if (tracks_.size() != 0)
      percent = (float)noTrk/tracks_.size();

    return percent;
}


vector<unsigned int> Vertex::numLayersTracks() const{
    vector<unsigned int> noLayers;
    for(TP tracks : tracks_){
        noLayers.push_back(tracks.numLayers());
    }
    
    return noLayers;
}

vector<unsigned int> Vertex::numPSLayersTracks() const{
  vector<unsigned int> noPSLayers;
  for(TP tracks : tracks_){
    noPSLayers.push_back(tracks.numPSLayers());
  }

  return noPSLayers;
}

vector<int> Vertex::layerIDTracks() const{
    int ok=0;
    vector<int> firstLayer;
    for(TP tracks : tracks_){
        for(unsigned int k=0;k<tracks.assocStubs().size();k++){
            if (tracks.assocStubs()[k]->layerId()==1){
                ok=1;
            }
        }
        firstLayer.push_back(ok);
    }
    
    return firstLayer;
}

