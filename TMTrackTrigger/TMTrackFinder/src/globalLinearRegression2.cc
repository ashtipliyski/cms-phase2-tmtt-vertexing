///=== This is the global Linear Regression for 4 helix parameters track fit algorithm.

///=== Written by: Thomas Schuh

#include "TMTrackTrigger/TMTrackFinder/interface/globalLinearRegression2.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <vector>
#include <algorithm>
#include <limits>

void globalLinearRegression2::bookHistos() {
  TFileDirectory inputDir = fs_->mkdir("globalLinearRegression2Internal");

  hisNumIter_      = inputDir.make<TH2F>("NumIter"     ,"; State; No. of iterations for Step",10,0.5,10.5,20,-0.5,19.5);

  hisGoodStubsKilled_ = inputDir.make<TH2F>("GoodStubsKilled","; State; Layer; No. good stubs in killed",10,0.5,10.5,30,0.5,30.5);
  hisBadStubsKilled_  = inputDir.make<TH2F>("BadStubsKilled" ,"; State; Layer; No. good stubs in killed",10,0.5,10.5,30,0.5,30.5);

  hisGoodStubsResA_ = inputDir.make<TH1F>("GoodStubsResA","; Good stub residual for Step A",200,0.,10.);
  hisBadStubsResA_  = inputDir.make<TH1F>("BadStubsResA" ,"; Bad stub residual for Step A",200,0.,10.);
  hisGoodStubsResB_ = inputDir.make<TH1F>("GoodStubsResB","; Good stub residual for Step B",200,0.,10.);
  hisBadStubsResB_  = inputDir.make<TH1F>("BadStubsResB" ,"; Bad stub residual for Step B",200,0.,10.);
  hisGoodStubsResC_ = inputDir.make<TH1F>("GoodStubsResC","; Good stub residual for Step C",200,0.,10.);
  hisBadStubsResC_  = inputDir.make<TH1F>("BadStubsResC" ,"; Bad stub residual for Step C",200,0.,10.);
  hisGoodStubsResD_ = inputDir.make<TH1F>("GoodStubsResD","; Good stub residual for Step D",200,0.,10.);
  hisBadStubsResD_  = inputDir.make<TH1F>("BadStubsResD" ,"; Bad stub residual for Step D",200,0.,10.);
  hisGoodStubsResE_ = inputDir.make<TH1F>("GoodStubsResE","; Good stub residual for Step E",200,0.,10.);
  hisBadStubsResE_  = inputDir.make<TH1F>("BadStubsResE" ,"; Bad stub residual for Step E",200,0.,10.);

  hisGoodStubsPhiRes_ = inputDir.make<TH1F>("GoodStubsPhiRes","; #phi residual of good stubs",200,0.,10.);
  hisBadStubsPhiRes_  = inputDir.make<TH1F>("BadStubsPhiRes" ,"; #phi residual of bad stubs",200,0.,10.);
  hisGoodStubsZRes_   = inputDir.make<TH1F>("GoodStubsZRes"  ,"; z residual of good stubs",200,0.,10.);
  hisBadStubsZRes_    = inputDir.make<TH1F>("BadStubsZRes"   ,"; z residual of bad stubs",200,0.,10.);

  hisPerfTrackPhiRes_ = inputDir.make<TH1F>("PerfTrackPhiRes","; #phi residual of good stubs",200,0.,10.);
  hisBadTrackPhiRes_  = inputDir.make<TH1F>("BadTrackPhiRes" ,"; #phi residual of bad stubs",200,0.,10.);
  hisPerfTrackZRes_   = inputDir.make<TH1F>("PerfTrackZRes"  ,"; z residual of good stubs",200,0.,10.);
  hisBadTrackZRes_    = inputDir.make<TH1F>("BadTrackZRes"   ,"; z residual of bad stubs",200,0.,10.);

  hisPerfTrackBigPhiRes_ = inputDir.make<TH1F>("PerfTrackBigPhiRes","; biggest #phi residual of good stubs",200,0.,10.);
  hisBadTrackBigPhiRes_  = inputDir.make<TH1F>("BadTrackBigPhiRes" ,"; biggest #phi residual of bad stubs",200,0.,10.);
  hisPerfTrackBigZRes_   = inputDir.make<TH1F>("PerfTrackBigZRes"  ,"; biggest z residual of good stubs",200,0.,10.);
  hisBadTrackBigZRes_    = inputDir.make<TH1F>("BadTrackBigZRes"   ,"; biggest z residual of bad stubs",200,0.,10.);  
  hisPerfTrackBigRes_    = inputDir.make<TH1F>("PerfTrackBigRes"   ,"; biggest residual of good stubs",200,0.,10.);
  hisBadTrackBigRes_     = inputDir.make<TH1F>("BadTrackBigRes"    ,"; biggest residual of bad stubs",200,0.,10.);
}

void globalLinearRegression2::fillHistos(const Stub* s) {
  bool goodStubLost = false;
  if (goodHTCandidate_) {
    const set<const TP*> stubTPs = s->assocTPs();
    goodStubLost = (stubTPs.find(tpHT_) != stubTPs.end());
  }

  if (goodStubLost) {
    hisGoodStubsKilled_->Fill(t_state_names[ state_ ], s->layerId(), 1.);
  } else {
    hisBadStubsKilled_->Fill(t_state_names[ state_ ], s->layerId(), 1.);
  }
  if (state_ == t_state::StepA) {
    if (goodStubLost) {
      hisGoodStubsResA_->Fill(largestResid_.second.second);
    } else {
      hisBadStubsResA_->Fill(largestResid_.second.second);
    }
  } else if (state_ == t_state::StepB) {
    if (goodStubLost) {
      hisGoodStubsResB_->Fill(largestResid_.second.second);
    } else {
      hisBadStubsResB_->Fill(largestResid_.second.second);
    }
  } else if (state_ == t_state::StepC) {
    if (goodStubLost) {
      hisGoodStubsResC_->Fill(largestResid_.second.second);
    } else {
      hisBadStubsResC_->Fill(largestResid_.second.second);
    }
  } else if (state_ == t_state::StepD) {
    if (goodStubLost) {
      hisGoodStubsResD_->Fill(largestResid_.second.second);
    } else {
      hisBadStubsResD_->Fill(largestResid_.second.second);
    }
  }	
}

void globalLinearRegression2::fillEndHistos() {
  bool perfTrack = false;
  if (goodCandidate_) {
    tp_ = Utility::matchingTP(settings_, stubs_, nMatchedLayersBest_, matchedStubsBest_);
    nMatchedPSLayers_ = Utility::countLayers(settings_, matchedStubsBest_, true, true);  
    unsigned int wrongStubs = stubs_.size() - matchedStubsBest_.size();
    perfTrack = (tp_ != nullptr && nMatchedPSLayers_ >= 2 && wrongStubs == 0);
  }

  if (tp_ == nullptr) {
    for (const Stub* s : stubs_) {
      hisBadTrackPhiRes_->Fill(residRPhi_[s]);
      hisBadTrackZRes_->Fill(residRZ_[s]);
    }
  } else {
    if (perfTrack) {
      for (const Stub* s : stubs_) {
	hisPerfTrackPhiRes_->Fill(residRPhi_[s]);
	hisPerfTrackZRes_->Fill(residRZ_[s]);
      }
    } else {
      for (const Stub* s : stubs_) {
	if (s->assocTP() == tp_) {
	  hisGoodStubsPhiRes_->Fill(residRPhi_[s]);
	  hisGoodStubsZRes_->Fill(residRZ_[s]);
	} else { 
	  hisBadStubsPhiRes_->Fill(residRPhi_[s]);
	  hisBadStubsZRes_->Fill(residRZ_[s]);
	}
      }
    }
  }

  float maxPhiRes = 0.;
  float maxZRes = 0;
  float maxRes = 0.;
  for (const Stub* s : stubs_) {
    maxPhiRes = max(maxPhiRes, residRPhi_[s]);
    maxZRes   = max(maxZRes, residRZ_[s]);
    maxRes    = max(maxRes, float(0.5*(2*residRPhi_[s] + residRZ_[s])));
  }
  if (perfTrack) {
    hisPerfTrackBigPhiRes_->Fill(maxPhiRes);
    hisPerfTrackBigZRes_->Fill(maxZRes);
    hisPerfTrackBigRes_->Fill(maxRes);
  } else {
    hisBadTrackBigPhiRes_->Fill(maxPhiRes);
    hisBadTrackBigZRes_->Fill(maxZRes);
    hisBadTrackBigRes_->Fill(maxRes);
  }
}

void globalLinearRegression2::initRun() {

  // IRT
  makeHistos_ = settings_->LRFillInternalHists();
  maxIterations_ = settings_->maxIterationsLR();

  invPtToDPhi_ = - settings_->invPtToDphi(); // turns HT +ive m into +ive RPhi slope
  numPhiSectors_ = settings_->numPhiSectors();
  etaRegions_ = settings_->etaRegions();
  minStubLayers_ = settings_->minStubLayers();
  minPSLayers_ = settings_->minPSLayers();
  houghMinPt_ = settings_->houghMinPt();
  beamWindowZ_ = settings_->beamWindowZ();
  debug_ = settings_->debug() == 7;
  chosenRofPhi_ = settings_->chosenRofPhi();
  chosenRofZ_ = settings_->chosenRofZ();
  houghNbinsPt_ = settings_->houghNbinsPt();
  houghNbinsPhi_ = settings_->houghNbinsPhi();
  binSizeQoverPtAxis_ = std::fabs( 2. / houghMinPt_ / (double)houghNbinsPt_ * invPtToDPhi_ );
  binSizePhiTrkAxis_ = 2. * M_PI / (double)numPhiSectors_ / (double)houghNbinsPhi_;
  phiSectorSize_ = 2. * M_PI / (double)numPhiSectors_;
  qOverPtSectorSize_ = std::fabs( 2. * houghMinPt_ / invPtToDPhi_ );
  zSectorSize_ = -1.;

  // IRT - this was buggy. fixed and moved to initFit()
  //for ( unsigned int iEtaReg = 0; iEtaReg < etaRegions_.size() - 1; iEtaReg++ )
  //  zSectorSize_ = std::max( zSectorSize_, (double)std::abs( std::sinh( etaRegions_[ iEtaReg ] ) - std::sinh( etaRegions_[ iEtaReg + 1 ] ) ) * chosenRofZ_ );
  //tanLambdaSectorSize_ = ( zSectorSize_ / 2. + beamWindowZ_ ) / chosenRofZ_;

  digitizeLR_ = settings_->digitizeLR();
  PhiPrecision_ = settings_->PhiPrecision();
  RPrecision_ = settings_->RPrecision();
  ZPrecision_ = settings_->ZPrecision();
  ZSlopeWidth_ = settings_->ZSlopeWidth();
  ZInterceptWidth_ = settings_->ZInterceptWidth();
  HTMBinBase_ = settings_->invPtToInvR() / houghMinPt_ / houghNbinsPt_;
  HTCBinBase_ = 2. * M_PI / numPhiSectors_ / houghNbinsPhi_;
  PhiBase_ = HTCBinBase_;
  PhiBase_ *= std::pow( 2., std::floor( std::log2( PhiPrecision_ / PhiBase_ ) ) );
  RPhiBase_ = HTCBinBase_ / HTMBinBase_;
  RPhiBase_ *= std::pow( 2., std::floor( std::log2( RPrecision_ / RPhiBase_ ) ) );
  ZSlopeBase_ = ( zSectorSize_ + 2. * beamWindowZ_ ) / chosenRofZ_ * std::pow( 2., - (int)ZSlopeWidth_ );
  ZInterceptBase_ = zSectorSize_ * std::pow( 2., - (int)ZInterceptWidth_ );
  ZBase_ = ZInterceptBase_;
  ZBase_ *=  std::pow( 2, std::floor( std::log2( ZPrecision_ / ZBase_ ) ) );
  RZBase_ = ZInterceptBase_ / ZSlopeBase_;
  RZBase_ *= std::pow( 2, std::floor( std::log2( RPrecision_ / RZBase_ ) ) );
  //for ( auto m : { "globalZ", "PhiTrack", "ZTrack", "PhiResid", "ZResid", "RPhi", "RZ", "Barrel", "Phi", "Z", "InterceptZ", "SlopeZ", "InterceptPhi", "SlopePhi", "Pixel", "RPhiSum", "RZSum", "PhiSum", "ZSum", "RPhiTimesPhiSum", "RZTimesZSum", "RPhiSquareSum", "RZSquareSum", "NumeratorSlopePhi", "NumeratorInterceptPhi", "DenominatorPhi", "NumeratorSlopeZ", "NumeratorInterceptZ", "DenominatorZ" } )
    //minMax_[ m ] = std::make_pair( std::numeric_limits< float >::infinity(), - std::numeric_limits< float >::infinity() );

  if (makeHistos_) this->bookHistos();
};

L1fittedTrack globalLinearRegression2::fit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg ) {

  // IRT - easily configurable cut (ignore name).
  //const float myCut = settings_->dupMaxZ0Scan();
  const float myCut = 1.4;

  initFit( l1track3D, iPhiSec, iEtaReg );
  onlySeedForRZHelix_ = false;
  onlySeedForRPhiHelix_ = false;
  onlySeedKilling_ = false;
  onlyNoSeedKilling_ = false;

  // Note if number of iterations must be truncated.
  const int numIterDesired = NumStubs_ + 1 - minStubLayersRed_;
  bool mustTruncate = (numIterDesired > int(maxIterations_)); 
  // Calculate number of iterations required for states B-D.
  const int numIterB = 0;
  // If state A is truncated, then state B runs killing stubs, so meaning that state C can no longer use last helix fit from state A. 
  const int numIterC = mustTruncate ? int(nLayers_ - minStubLayersRed_) : max(int(nLayers_ - minStubLayersRed_) - 1, 0);
  const int numIterD = 1;
  // Calculate max. number of iterations for state A.
  const int maxIterA = max(int(maxIterations_) - numIterB - numIterC - numIterD, 0);
  //cout<<"START "<<NumStubs_<<" "<<nLayers_<<" "<<maxIterA<<" "<<numIterB<<" "<<numIterC<<" "<<numIterD<<endl;

  unsigned int nIterationsLast = 0;

  while ( true ) {
    if ( not checkValidity() )
      return createTrack( l1track3D );
    switch ( state_ ) {
      // Reduce tracks to one stub per layer, refitting helix at each iteration.
      case t_state::StepA :
	onlyPSForRZ_ = false;
	onlyRPhiResiduals_ = false;
	onlyRZResiduals_ = false;
	onlyPSKilling_ = false;
	onlySSKilling_ = false;
	useRelativResiduals_ = true;
	calcHelix_ = true;
        residualCut_ = 0.;
        killLargestResid_ = (nIterations_ < maxIterA);
        break;
      // Reduce tracks to one stub per layer, saving time by not refitting helix.
      case t_state::StepB :
	onlyPSForRZ_ = false;
	onlyRPhiResiduals_ = false;
	onlyRZResiduals_ = false;
	onlyPSKilling_ = false;
	onlySSKilling_ = false;
	useRelativResiduals_ = true;
	calcHelix_ = false;
	// This allows for case where no iterations were performed in StepA.
	if (nIterations_ == 0) calcHelix_ = true;
        residualCut_ = 0.;
        killLargestResid_ = true;
        break;
      // Reduce tracks to 4 layers
      case t_state::StepC :
	onlyPSForRZ_ = true;
	onlyRPhiResiduals_ = false;
	onlyRZResiduals_ = false;
	onlyPSKilling_ = false;
	onlySSKilling_ = (nPSLayers_ == 2);
	useRelativResiduals_ = false;
	calcHelix_ = true;
        residualCut_ = 0.;
        killLargestResid_ = (nLayers_ > minStubLayersRed_);
        break;
      // Do final fit to track and kill if it is poor quality (done in 1 iteration).
      case t_state::StepD :
	onlyPSForRZ_ = true;
	onlyRPhiResiduals_ = false;
	onlyRZResiduals_ = false;
	onlyPSKilling_ = false;
	onlySSKilling_ = false;
	useRelativResiduals_ = false;
	calcHelix_ = true;
        residualCut_ = myCut;
        killLargestResid_ = true;
        break;
      default:
	cout<<"MISTAKE"<<endl;
	exit(1);
    }

    // IRT
    bool killedStub = false;
    if (killLargestResid_) {
       stateCalls_[ state_ ]++;
       if ( calcHelix_ ) calcHelix();
       //cout<<"ITER "<<nIterations_<<" "<<state_<<" "<<nLayers_<<" "<<NumStubs_<<endl;
       calcResidual();
       findLargestResidual();
       if ( largestResid_.second.second > residualCut_ ) {
         // Staying within this state.
	 killedStub = true;
         killLargestResidual();
       }
    }

    if (not killedStub) {
      // We have finished with this state.
      if (makeHistos_) {
	unsigned int nIter = nIterations_ - nIterationsLast;
	nIterationsLast = nIterations_;
	hisNumIter_->Fill( t_state_names[ state_ ], nIter, 1.);
      }
      if ( state_ < lastState_ ) {
	// Move on to next state.
        state_ = static_cast< t_state >( state_ + 1 );
      } else {
	// Finished
        break;
      }
    }
  }

  checkValidity();
  calcChiSq();

  if (makeHistos_) this->fillEndHistos();

  // Don't bother with cuts on z0 and Pt.
  bool insidePhiSec = this->sectorCheck( interceptRPhi_, phiSectorSize_ );
  bool insideEtaSec = this->sectorCheck( interceptRZ_, zSectorSize_ );
  if ( not (insidePhiSec && insideEtaSec) ) valid_ = false;

  return createTrack( l1track3D );
}

L1fittedTrack globalLinearRegression2::createTrack( const L1track3D& l1track3D ) {

  trackParams_["qOverPt"] = helixRPhi_.first / invPtToDPhi_;
  trackParams_["phi0"] = helixRPhi_.second - helixRPhi_.first * chosenRofPhi_;
  trackParams_["t"] = helixRZ_.first;
  trackParams_["z0"] =  helixRZ_.second - helixRZ_.first * chosenRofZ_;
  L1fittedTrack fitTrk( settings_, l1track3D, stubs_, trackParams_["qOverPt"], 0, trackParams_["phi0"], trackParams_["z0"], trackParams_["t"], chiSq_, nPar_, iPhiSec_, iEtaReg_, valid_);
  std::string lostMatchingState = t_state_names[ lostMatchingState_ ];
  std::unordered_map< std::string, int > stateCalls;
  for ( auto state : stateCalls_ )
    stateCalls[ t_state_names[ state.first ] ] = state.second;
  fitTrk.setInfoLR( nIterations_, lostMatchingState, stateCalls ); // Set extra info for debugging/histogramming.
  return fitTrk;

}

void globalLinearRegression2::initFit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg ) {

  iPhiSec_ = iPhiSec;
  iEtaReg_ = iEtaReg;

  // IRT - make z sector size dependent on sector number.
  zSectorSize_ = std::fabs( std::sinh( etaRegions_[ iEtaReg ] ) - std::sinh( etaRegions_[ iEtaReg + 1 ] ) ) * chosenRofZ_;
  tanLambdaSectorSize_ = ( zSectorSize_ + 2.*beamWindowZ_ ) / chosenRofZ_;

  phiCentre_ = phiSectorSize_ * ( (double)iPhiSec_ +.5 ) - M_PI;
  zCentre_ = ( std::sinh( etaRegions_[ iEtaReg_ ] ) + std::sinh( etaRegions_[ iEtaReg_ + 1 ] ) ) / 2. * chosenRofZ_;
  HTHelixRPhi_ = l1track3D.getHelixRphi();
  HTHelixRPhi_.first *= invPtToDPhi_;
  HTHelixRPhi_.second = reco::deltaPhi( HTHelixRPhi_.second - phiCentre_ + HTHelixRPhi_.first * chosenRofPhi_, 0. );
  HTHelixRZ_ = std::make_pair( zCentre_ / chosenRofZ_, 0. );
  stubs_ = l1track3D.getStubs();
  NumStubs_ = stubs_.size();
  minStubLayersRed_ = Utility::numLayerCut("FIT", settings_, iPhiSec, iEtaReg, fabs(l1track3D.qOverPt()), l1track3D.eta());
  minStubLayersRed_ = 4;
  nLayers_ = Utility::countLayers( settings_, stubs_ );
  nPSLayers_ = Utility::countLayers( settings_, stubs_, false, true );
  if ( not checkValidity() )
    return;
  stubMap_.clear();
  for ( auto stub : stubs_ ) {
    layerID_ = stub->layerIdReduced();
    stubMap_[ layerID_ ].push_back( stub );
    /*RPhi_ = rPhi( stub );
    Phi_ = phi( stub );
    RZ_ = rZ( stub );
    Z_ = z( stub );
    minMax_["RPhi"] = std::make_pair( std::min( minMax_["RPhi"].first, RPhi_ ), std::max( minMax_["RPhi"].second, RPhi_ ) );
    minMax_["RZ"] = std::make_pair( std::min( minMax_["RZ"].first, RZ_ ), std::max( minMax_["RZ"].second, RZ_ ) );
    minMax_["Phi"] = std::make_pair( std::min( minMax_["Phi"].first, Phi_ ), std::max( minMax_["Phi"].second, Phi_ ) );
    minMax_["Z"] = std::make_pair( std::min( minMax_["Z"].first, Z_ ), std::max( minMax_["Z"].second, Z_ ) );
    minMax_["globalZ"] = std::make_pair( std::min( minMax_["globalZ"].first, stub->z() ), std::max( minMax_["globalZ"].second, stub->z() ) );
    if ( stub->psModule() )
      minMax_["Pixel"].second = max( minMax_["Pixel"].second, stub->r() );
    else
      minMax_["Pixel"].first = min( minMax_["Pixel"].first, stub->r() );
    if ( stub->barrel() )
      minMax_["Barrel"].second = max( minMax_["Barrel"].second, std::fabs( stub->z() ) );
    else
      minMax_["Barrel"].first = min( minMax_["Barrel"].first, std::fabs( stub->z() ) );*/
  }
  nIterations_ = 0;
  // IRT
  //state_ = generalCleaning;
  //state_ = seedFinding;
  state_ = firstState_;
  stateCalls_.clear();
  tpHT_ = Utility::matchingTP( settings_, stubs_, nMatchedLayersBest_, matchedStubsBest_ );
  nMatchedPSLayers_ = Utility::countLayers(settings_, matchedStubsBest_, true, true);  
  goodCandidate_ = (l1track3D.getMatchedTP() != nullptr && nMatchedPSLayers_ >= 2);
  // IRT
  goodHTCandidate_ = goodCandidate_;
  lostMatchingState_ = firstState_;

  oldNumStubs_ = 0;
  if ( debug_ )
    initDebug( l1track3D );

}

double globalLinearRegression2::digi( double val, double base ) {
  if ( digitizeLR_ )
    return std::floor( val / base ) * base;
  return val;
}

double globalLinearRegression2::rPhi( const Stub* stub ) {
  return digi( stub->r() - chosenRofPhi_, RPhiBase_ );
}

double globalLinearRegression2::rZ( const Stub* stub ) {
  return digi( rPhi( stub ) + chosenRofPhi_ - chosenRofZ_, RZBase_ );
}

double globalLinearRegression2::phi( const Stub* stub ) {
  return digi( reco::deltaPhi( stub->phi() - phiCentre_ - HTHelixRPhi_.second - HTHelixRPhi_.first * rPhi( stub ), 0. ), PhiBase_ );
}

double globalLinearRegression2::z( const Stub* stub ) {
  return digi( stub->z() - zCentre_ - HTHelixRZ_.second - HTHelixRZ_.first * rZ( stub ), RZBase_ );
}

bool globalLinearRegression2::checkValidity() {

  valid_ = true;
  if ( nLayers_ < minStubLayersRed_ )
    valid_ = false;
  if ( nPSLayers_ < minPSLayers_ )
    valid_ = false;
  return valid_;

}

void globalLinearRegression2::calcHelix( ) {

  // IRT
  if ( oldNumStubs_ != NumStubs_ ) nIterations_++;
  oldNumStubs_ = NumStubs_;

  resetSums();
  for ( auto layer : stubMap_ ) if ( layer.second.size() > 0 ) {
    phiMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    rPhiMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    rZMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    zMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    sumOnlyRZ_ = true;
    sumOnlyRPhi_ = true;
    for ( auto stub : layer.second ) {
      RPhi_ = rPhi( stub );
      Phi_ = phi( stub );
      RZ_ = rZ( stub );
      Z_ = z( stub );
      seedingStub_ = isSeedingStub( stub );
      rPhiMinMax_ = std::make_pair( std::min( rPhiMinMax_.first, RPhi_ ), std::max( rPhiMinMax_.second, RPhi_ ) );
      rZMinMax_ = std::make_pair( std::min( rZMinMax_.first, RZ_ ), std::max( rZMinMax_.second, RZ_ ) );
      if ( seedingStub_ or not onlySeedForRPhiHelix_ ) {
        sumOnlyRZ_ = false;
        phiMinMax_ = std::make_pair( std::min( phiMinMax_.first, Phi_ ), std::max( phiMinMax_.second, Phi_ ) );
      }

      // IRT
      //if ( seedingStub_ or not onlySeedForRZHelix_ ) {
      //  sumOnlyRPhi_ = false;
      //  zMinMax_ = std::make_pair( std::min( zMinMax_.first, Z_ ), std::max( zMinMax_.second, Z_ ) );
      //}
      if ( stub->psModule() or not onlyPSForRZ_ ) {
        sumOnlyRPhi_ = false;
        zMinMax_ = std::make_pair( std::min( zMinMax_.first, Z_ ), std::max( zMinMax_.second, Z_ ) );
      }

    }
    RPhi_ = ( rPhiMinMax_.first + rPhiMinMax_.second ) / 2.;
    RZ_ = ( rZMinMax_.first + rZMinMax_.second ) / 2.;
    Phi_ = ( phiMinMax_.first + phiMinMax_.second ) / 2.;
    Z_ = ( zMinMax_.first + zMinMax_.second ) / 2.;
    updateSums();
  }
  calcLinearParameter();
  helixRPhi_ = std::make_pair( HTHelixRPhi_.first + slopeRPhi_, reco::deltaPhi( HTHelixRPhi_.second + interceptRPhi_ + phiCentre_, 0. ) );
  helixRZ_ = std::make_pair( HTHelixRZ_.first + slopeRZ_, HTHelixRZ_.second + interceptRZ_ + zCentre_ );

}

bool globalLinearRegression2::isSeedingStub( const Stub* stub ) {

  if ( not stub->psModule() )
    return false;
  // Tracker layer ID number (1-6 = barrel layer; 11-15 = endcap A disk; 21-25 = endcap B disk)
  if ( not stub->barrel() ){
    if ( stub->layerId() < 11 or stub->layerId() > 23 )
      return false;
    if ( stub->layerId() > 13 and stub->layerId() < 20 )
      return false;
  }
  return true;

}

void globalLinearRegression2::calcResidual() {

  RPhiResiduals_.clear();
  RZResiduals_.clear();
  residuals_.clear();
  relativResiduals_.clear();

  // IRT
  residRPhi_.clear();
  residRZ_.clear();
  relResidRPhi_.clear();
  relResidRZ_.clear();

  for ( auto layer : stubMap_ ) if ( layer.second.size() > 0 ) {
    minRPhiResid_ = std::numeric_limits< double >::infinity();
    minRZResid_ = std::numeric_limits< double >::infinity();
    minResid_ = std::numeric_limits< double >::infinity();
    for ( auto stub : layer.second ) {
      Phi_ = interceptRPhi_ + slopeRPhi_ * rPhi( stub );
      //minMax_["PhiTrack"] = std::make_pair( min( minMax_["PhiTrack"].first, Phi_ ), max( minMax_["PhiTrack"].second, Phi_ ) );
      RPhiResid_ = phi( stub ) - Phi_;
      //minMax_["PhiResid"] = std::make_pair( min( minMax_["PhiResid"].first, RPhiResid_ ), max( minMax_["PhiResid"].second, RPhiResid_ ) );
      RPhiResid_ = std::fabs( RPhiResid_ );
      RPhiSigma_ = 0.001;
      if ( not stub->barrel() )
        RPhiSigma_ = 0.002;
      RPhiResid_ /= RPhiSigma_;
      RPhiResiduals_[ layer.first ].push_back( RPhiResid_ );
      Z_ = interceptRZ_ + slopeRZ_ * rZ( stub );
      //minMax_["ZTrack"] = std::make_pair( min( minMax_["ZTrack"].first, Z_ ), max( minMax_["ZTrack"].second, Z_ ) );
      RZResid_ = z( stub ) - Z_;
      //minMax_["ZResid"] = std::make_pair( min( minMax_["ZResid"].first, RZResid_ ), max( minMax_["ZResid"].second, RZResid_ ) );
      RZResid_ = std::fabs( RZResid_ );

      // IRT
      // New Thomas code
      //RZSigma_ = stub->stripLength();
      //if ( stub->psModule() and ( state_ == generalCleaning or state_ == seedFinding ) ) {
      //  RZSigma_ = 2.;
      //}
      // Old Thomas code
      RZSigma_ = 5. / std::sqrt( 12 );
      if ( onlyPSForRZ_ and stub->psModule() )
        RZSigma_ = stub->sigmaPar() * 2.;

      if ( not stub->barrel() )
        RZSigma_ *= std::fabs( HTHelixRZ_.first );
      RZResid_ /= RZSigma_;
      RZResiduals_[ layer.first ].push_back( RZResid_ );
      resid_ = ( RPhiResid_ + RZResid_ ) / 2.;
      if ( onlyRPhiResiduals_ )
        resid_ = RPhiResid_;
      if ( onlyRZResiduals_ )
        resid_ = RZResid_;
      residuals_[ layer.first ].push_back( resid_ );
      minRPhiResid_ = std::min( minRPhiResid_, RPhiResid_ );
      minRZResid_ = std::min( minRZResid_, RZResid_ );
      minResid_ = std::min( minResid_, resid_ );
      // IRT for histograms
      residRPhi_[stub] = RPhiResid_;
      residRZ_[stub]   = RZResid_;
      minRPhiResid_ = std::min( minRPhiResid_, RPhiResid_ );
      minRZResid_   = std::min( minRZResid_, RZResid_ );
    }
    for ( unsigned int i = 0; i < layer.second.size(); i++ ){
      relativResiduals_[ layer.first ].push_back( residuals_[ layer.first ][ i ] - minResid_ );
    }
    // IRT for histograms.
    for ( const Stub* s : layer.second) {
      relResidRPhi_[s] = residRPhi_[s] - minRPhiResid_;
      relResidRZ_[s]   = residRZ_[s]   - minRZResid_;
    }
  }

}

void globalLinearRegression2::findLargestResidual() {

  std::map< unsigned int, std::vector< double > > residuals = residuals_;
  if ( useRelativResiduals_ )
    residuals = relativResiduals_;
  largestResid_.second.second = -1.;
  for ( auto layer : residuals ) for ( unsigned int i = 0; i < layer.second.size(); i++ ) {
    const Stub* stub( stubMap_[ layer.first ][ i ] );

    // IRT
    if (onlyPSKilling_ && onlySSKilling_) {
      cout<<"Stupid!"<<endl;
      exit(1);
    }
    if ( onlyPSKilling_ and not stub->psModule() )
      continue;
    if ( onlySSKilling_ and stub->psModule() )
      continue;

    seedingStub_ = isSeedingStub( stub );
    if ( ( onlySeedKilling_ and not seedingStub_ ) or ( onlyNoSeedKilling_ and seedingStub_ ) )
      continue;
    resid_ = residuals[ layer.first ][ i ];
    if ( resid_ > largestResid_.second.second )
      largestResid_ = std::make_pair( layer.first, std::make_pair( stub, resid_ ) );
  }

}

void globalLinearRegression2::killLargestResidual() {

  if ( goodCandidate_ )
    iterDebug();
  layerID_ = largestResid_.first;
  const Stub* stub( largestResid_.second.first );

  // IRT
  if (makeHistos_) this->fillHistos(stub);

  stubMap_[ layerID_ ].erase( std::remove( stubMap_[ layerID_ ].begin(), stubMap_[ layerID_ ].end(), stub ), stubMap_[ layerID_ ].end() );
  stubs_.erase( std::remove( stubs_.begin(), stubs_.end(), stub ), stubs_.end() );
  nLayers_ = Utility::countLayers( settings_, stubs_ );
  nPSLayers_ = Utility::countLayers( settings_, stubs_, false, true );
  // IRT
  //oldNumStubs_ = NumStubs_;
  NumStubs_--;
  tp_ = Utility::matchingTP( settings_, stubs_, nMatchedLayersBest_, matchedStubsBest_ );
  nMatchedPSLayers_ = Utility::countLayers(settings_, matchedStubsBest_, true, true);  
  if ( goodCandidate_ and (tpHT_ != tp_ or nMatchedPSLayers_ < 2) ) {
    goodCandidate_ = false;
    lostMatchingState_ = state_; 
  }
}

void globalLinearRegression2::calcChiSq() {

  chiSq_ = 0.;
  for( auto layer : residuals_ ) for( auto resid : layer.second )
    chiSq_ += resid;
  chiSq_ /= NumStubs_;

}

void globalLinearRegression2::resetSums() {

  NRPhi_ = 0.;
  NRZ_ = 0.;
  RPhiTimesPhiSum_ = 0.;
  RZTimesZSum_ = 0.;
  RPhiSum_ = 0.;
  RZSum_ = 0.;
  PhiSum_ = 0.;
  ZSum_ = 0.;
  RPhiSquareSum_ = 0.;
  RZSquareSum_ = 0.;

}

void globalLinearRegression2::updateSums() {

  if ( not sumOnlyRZ_ ) {
    NRPhi_ ++;
    RPhiTimesPhiSum_ += RPhi_ * Phi_;
    RPhiSum_ += RPhi_;
    PhiSum_ += Phi_;
    RPhiSquareSum_ += RPhi_ * RPhi_;
    /*minMax_["RPhiSum"] = std::make_pair( min( minMax_["RPhiSum"].first, RPhiSum_ ), max( minMax_["RPhiSum"].second, RPhiSum_ ) );
    minMax_["PhiSum"] = std::make_pair( min( minMax_["PhiSum"].first, PhiSum_ ), max( minMax_["PhiSum"].second, PhiSum_ ) );
    minMax_["RPhiTimesPhiSum"] = std::make_pair( min( minMax_["RPhiTimesPhiSum"].first, RPhiTimesPhiSum_ ), max( minMax_["RPhiTimesPhiSum"].second, RPhiTimesPhiSum_ ) );
    minMax_["RPhiSquareSum"] = std::make_pair( min( minMax_["RPhiSquareSum"].first, RPhiSquareSum_ ), max( minMax_["RPhiSquareSum"].second, RPhiSquareSum_ ) );*/
  }
  if ( not sumOnlyRPhi_ ) {
    NRZ_ ++;
    RZTimesZSum_ += RZ_ * Z_;
    RZSum_ += RZ_;
    ZSum_ += Z_;
    RZSquareSum_ += RZ_ * RZ_;
    /*minMax_["RZSum"] = std::make_pair( min( minMax_["RZSum"].first, RZSum_ ), max( minMax_["RZSum"].second, RZSum_ ) );
    minMax_["ZSum"] = std::make_pair( min( minMax_["ZSum"].first, ZSum_ ), max( minMax_["ZSum"].second, ZSum_ ) );
    minMax_["RZTimesZSum"] = std::make_pair( min( minMax_["RZTimesZSum"].first, RZTimesZSum_ ), max( minMax_["RZTimesZSum"].second, RZTimesZSum_ ) );
    minMax_["RZSquareSum"] = std::make_pair( min( minMax_["RZSquareSum"].first, RZSquareSum_ ), max( minMax_["RZSquareSum"].second, RZSquareSum_ ) );*/
  }

}

void globalLinearRegression2::calcLinearParameter() {

  denominatorRPhi_ = NRPhi_ * RPhiSquareSum_ - RPhiSum_ * RPhiSum_;
  denominatorRZ_ = NRZ_ * RZSquareSum_ - RZSum_ * RZSum_;
  slopeRPhi_ = ( NRPhi_ * RPhiTimesPhiSum_ - RPhiSum_ * PhiSum_ );
  slopeRZ_ = ( NRZ_ * RZTimesZSum_ - RZSum_ * ZSum_ );
  interceptRPhi_ = digi( RPhiSquareSum_ * PhiSum_ - RPhiSum_ * RPhiTimesPhiSum_, RPhiBase_ * RPhiBase_ * PhiBase_ * std::pow( 2., 5 ) );
  interceptRZ_ = digi( RZSquareSum_ * ZSum_ - RZSum_ * RZTimesZSum_, RPhiBase_ * RPhiBase_ * PhiBase_ * std::pow( 2., 4 ) );
  /*minMax_["NumeratorSlopePhi"] = std::make_pair( min( minMax_["NumeratorSlopePhi"].first, slopeRPhi_ ), max( minMax_["NumeratorSlopePhi"].second, slopeRPhi_ ) );
  minMax_["NumeratorInterceptPhi"] = std::make_pair( min( minMax_["NumeratorInterceptPhi"].first, interceptRPhi_ ), max( minMax_["NumeratorInterceptPhi"].second, interceptRPhi_ ) );
  minMax_["DenominatorPhi"] = std::make_pair( min( minMax_["DenominatorPhi"].first, denominatorRPhi_ ), max( minMax_["DenominatorPhi"].second, denominatorRPhi_ ) );
  minMax_["NumeratorSlopeZ"] = std::make_pair( min( minMax_["NumeratorSlopeZ"].first, slopeRZ_ ), max( minMax_["NumeratorSlopeZ"].second, slopeRZ_ ) );
  minMax_["NumeratorInterceptZ"] = std::make_pair( min( minMax_["NumeratorInterceptZ"].first, interceptRZ_ ), max( minMax_["NumeratorInterceptZ"].second, interceptRZ_ ) );
  minMax_["DenominatorZ"] = std::make_pair( min( minMax_["DenominatorZ"].first, denominatorRZ_ ), max( minMax_["DenominatorZ"].second, denominatorRZ_ ) );
  */
  if ( digitizeLR_ ) {
    int shift = (int)std::ceil( std::log2( denominatorRPhi_ / RPhiBase_ / RPhiBase_ ) ) - 10;
    if ( shift > 0 ) {
      denominatorRPhi_ /= std::pow( 2., shift );
      slopeRPhi_ /= std::pow( 2., shift );
      interceptRPhi_ /= std::pow( 2., shift );
    }
  }
  reciprocalDenominatorRPhi_ = digi( 1. / denominatorRPhi_, std::pow( 2., -16 ) / RPhiBase_ / RPhiBase_ ) ;
  slopeRPhi_ *= reciprocalDenominatorRPhi_;
  interceptRPhi_ *= reciprocalDenominatorRPhi_;
  // IRT - removed these constraints as they degrade resolution & lose efficiency.
  //slopeRPhi_ = sectorConstrain( slopeRPhi_, qOverPtSectorSize_ );
  //interceptRPhi_ = sectorConstrain( interceptRPhi_, phiSectorSize_ );
  if ( digitizeLR_ ) {
    int shift = (int)std::ceil( std::log2( denominatorRZ_ / RZBase_ / RZBase_ ) ) - 10;
    if ( shift > 0 ) {
      denominatorRZ_ /= std::pow( 2., shift );
      slopeRZ_ /= std::pow( 2., shift );
      interceptRZ_ /= std::pow( 2., shift );
    }
  }
  reciprocalDenominatorRZ_ = digi( 1. / denominatorRZ_, std::pow( 2., -16 ) / RZBase_ / RZBase_ );
  slopeRZ_ *= reciprocalDenominatorRZ_;
  interceptRZ_ *= reciprocalDenominatorRZ_;
  // IRT - removed these constraints as they degrade resolution & lose efficiency.
  //slopeRZ_ = sectorConstrain( slopeRZ_, tanLambdaSectorSize_ );
  //interceptRZ_ = sectorConstrain( interceptRZ_, zSectorSize_ );
  /*minMax_["SlopePhi"] = std::make_pair( min( minMax_["SlopePhi"].first, slopeRPhi_ ), max( minMax_["SlopePhi"].second, slopeRPhi_ ) );
  minMax_["InterceptPhi"] = std::make_pair( min( minMax_["InterceptPhi"].first, interceptRPhi_ ), max( minMax_["InterceptPhi"].second, interceptRPhi_ ) );
  minMax_["SlopeZ"] = std::make_pair( min( minMax_["SlopeZ"].first, slopeRZ_ ), max( minMax_["SlopeZ"].second, slopeRZ_ ) );
  minMax_["InterceptZ"] = std::make_pair( min( minMax_["InterceptZ"].first, interceptRZ_ ), max( minMax_["InterceptZ"].second, interceptRZ_ ) );
*/
}

double globalLinearRegression2::sectorConstrain( double val, double limit ) {

  if ( val > limit / 2. )
    return limit / 2.;
  else if ( val < - limit / 2. )
    return - limit / 2.;
  else
    return val;

}


bool globalLinearRegression2::sectorCheck( double val, double limit ) {

  bool insideSector = (fabs(val) < limit/2.);
  return insideSector;
}

void globalLinearRegression2::initDebug( const L1track3D& l1track3D ) {

  goodCandidate_ = tpHT_ != nullptr;
  if ( goodCandidate_ ) {
    assocStubs_ = tpHT_->assocStubs();
    tpHelixRPhi_ = std::make_pair( tpHT_->qOverPt() * invPtToDPhi_, reco::deltaPhi( tpHT_->phi0() - phiCentre_ + tpHT_->qOverPt() * invPtToDPhi_ * chosenRofPhi_, 0. ) );
    tpHelixRZ_ = std::make_pair( tpHT_->tanLambda(), tpHT_->z0() - zCentre_ + tpHT_->tanLambda() * chosenRofZ_ );
    debugStream_.str("");
    debugStream_.clear();
    debugStream_ << std::setprecision( 6 );
    //debugStream_ << "HTHelixRPhi " << HTHelixRPhi_.first << " " << HTHelixRPhi_.second << " HTHelixRZ " << HTHelixRZ_.first << " " << HTHelixRZ_.second << std::endl;
    debugStream_ << "tpHelixRPhi " << tpHelixRPhi_.first - HTHelixRPhi_.first << " " << tpHelixRPhi_.second - HTHelixRPhi_.second;
    debugStream_ << " tpHelixRZ " << tpHelixRZ_.first - HTHelixRZ_.first << " " << tpHelixRZ_.second - HTHelixRZ_.second << std::endl;
    debugStream_ << HTHelixRPhi_.first << " " << reco::deltaPhi( HTHelixRPhi_.second + phiCentre_ - helixRPhi_.first * chosenRofPhi_, 0. ) << " ";
    debugStream_ << HTHelixRZ_.first << " " << HTHelixRZ_.second + zCentre_ - helixRZ_.first * chosenRofZ_ << std::endl;
    debugStream_ << tpHT_->qOverPt() * invPtToDPhi_ << " " << reco::deltaPhi( tpHT_->phi0(), 0. ) << " ";
    debugStream_ << tpHT_->tanLambda() << " " << tpHT_->z0() << std::endl;
    debugStream_ << "tpStubs: ( layerID, ps, barrel, rPhi, phi, rZ, z )" << std::endl;
    for ( auto s : assocStubs_ ) {
      debugStream_ << s->layerIdReduced() << " " << s->psModule() << " " << s->barrel() << " " << std::setfill(' ') << std::setw( 8 ) << std::fixed << rPhi( s ) << " " << phi( s ) << " " << rZ( s ) << " " << std::setw( 8 ) << z( s );
      debugStream_ << " " << s->r() << " " << s->phi() << " " << s->z() << " ";
      for ( auto as : stubs_ ) if ( s == as )
        debugStream_ << "#this one is present in the candidate";
      debugStream_ << std::endl;
    }
    debugStream_ << "htStubs:" << std::endl;
    for ( auto s : stubs_ ) {
      debugStream_ << s->layerIdReduced() << " " << s->psModule() << " " << s->barrel() << " " << std::setfill(' ') << std::setw( 8 ) << std::fixed << rPhi( s ) << " " << phi( s ) << " " << rZ( s ) << " " << std::setw( 8 ) << z( s );
      debugStream_ << " " << s->r() << " " << s->phi() << " " << s->z() << " ";
      for ( auto as : assocStubs_ ) if ( s == as )
        debugStream_<< "#that is a matched stub.";
      debugStream_<< std::endl;
    }
  }

}

void globalLinearRegression2::iterDebug() {

  debugStream_ << "state: " << state_ << std::endl;
  debugStream_ << std::setprecision( 6 );
  debugStream_ << "helix fit " << slopeRPhi_ << " " << interceptRPhi_ << " " << slopeRZ_ << " " << interceptRZ_ << std::endl;
  debugStream_ << std::setprecision( 3 );
  debugStream_ << "residuals:" << std::endl;
  for ( auto layer : stubMap_ ) for ( unsigned int is = 0; is < layer.second.size(); is++ ) {
    const Stub* stub( stubMap_[ layer.first ][ is ] );
    debugStream_ << layer.first << " " << stub->psModule() << " " << stub->barrel() << " " << std::fixed  << std::setfill(' ') << std::setprecision( 3 ) << std::setw( 8 ) << rPhi( stub ) << " " << std::setw( 8 ) << phi( stub ) << " "  << std::setw( 8 ) << rZ( stub ) << " " << std::setw( 8 ) << z( stub ) << " | " << std::setw( 8 ) << RPhiResiduals_[ layer.first ][ is ] << " " << std::setw( 8 ) <<  RZResiduals_[ layer.first ][ is ] << " | " << std::setw( 8 ) << relativResiduals_[ layer.first ][ is ] << " | " << std::setw( 8 ) << residuals_[ layer.first ][ is ] << std::endl;
  }
  const Stub* stub( largestResid_.second.first );
  debugStream_ << "killed stub " << stub->layerIdReduced() << " " << stub->psModule() << " " << stub->barrel() << " " << std::setfill(' ') << std::setw( 6 ) << std::setprecision( 3 ) << std::fixed << rPhi( stub ) << " " << phi( stub ) << " " << rZ( stub ) << " " << std::setw( 8 ) << z( stub ) << " ";
  for ( auto as : assocStubs_ ) if ( as == stub )
    debugStream_ << "that was a matche one ...";
  debugStream_<< std::endl;

}
