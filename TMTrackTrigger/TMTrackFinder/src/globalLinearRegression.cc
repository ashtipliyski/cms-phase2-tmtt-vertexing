///=== This is the global Linear Regression for 4 helix parameters track fit algorithm.

///=== Written by: Thomas Schuh

#include "TMTrackTrigger/TMTrackFinder/interface/globalLinearRegression.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include <vector>
#include <algorithm>
#include <limits>


void globalLinearRegression::initRun() {

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

};

L1fittedTrack globalLinearRegression::fit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg ) {

  initFit( l1track3D, iPhiSec, iEtaReg );
  onlyPSForRZHelix_ = false;
  calcHelix();
  onlyPSForRZHelix_ = true;
  while ( true ) {
    if ( not checkValidity() )
      return createTrack( l1track3D );
    switch ( state_ ) {
      case generalCleaning :
        useRelativResiduals_ = false;
        residualCut_ = 1.;
        killLargestResid_ = true;
        break;
      case layerPurging :
        useRelativResiduals_ = true;
        residualCut_ = 0.;
        break;
      case layerKilling :
        onlySSKilling_ = nPSLayers_ == minPSLayers_;
        useRelativResiduals_ = false;
        residualCut_ = 0.;
        killLargestResid_ = nLayers_ > minStubLayersRed_;
        break;
      case candidateKilling :
        onlySSKilling_ = false;
        residualCut_ = 1.4;
        killLargestResid_ = true;
        break;
      default :
        std::cout << "MISTAKE" << std::endl; std::exit(1);
    }
    stateCalls_[ state_ ]++;
    calcResidual();
    findLargestResidual();
    if ( killLargestResid_ and largestResid_.second.second > residualCut_ )
      killLargestResidual();
    else if ( state_ < candidateKilling ) {
      calcHelix();
      state_ = static_cast< t_state >( state_ + 1 );
    } else
      break;
  }
  checkValidity();
  calcChiSq();

  // Declare track to be invalid if its direction is not within sector.
  // Don't bother with cuts on z0 and Pt.
  bool insidePhiSec = this->sectorCheck( interceptRPhi_, phiSectorSize_ );
  bool insideEtaSec = this->sectorCheck( interceptRZ_, zSectorSize_ );
  if ( not (insidePhiSec && insideEtaSec) ) valid_ = false;

  return createTrack( l1track3D );

}

L1fittedTrack globalLinearRegression::createTrack( const L1track3D& l1track3D ) {

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

void globalLinearRegression::initFit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg ) {

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
  state_ = static_cast< t_state >( 0 );
  stateCalls_.clear();
  tpHT_ = Utility::matchingTP( settings_, stubs_, nMatchedLayersBest_, matchedStubsBest_ );
  goodCandidate_ = l1track3D.getMatchedTP() != nullptr;
  oldNumStubs_ = 0;
  if ( debug_ )
    initDebug( l1track3D );

}

double globalLinearRegression::digi( double val, double base ) {
  if ( digitizeLR_ )
    return std::floor( val / base ) * base;
  return val;
}

double globalLinearRegression::rPhi( const Stub* stub ) {
  return digi( stub->r() - chosenRofPhi_, RPhiBase_ );
}

double globalLinearRegression::rZ( const Stub* stub ) {
  return digi( rPhi( stub ) + chosenRofPhi_ - chosenRofZ_, RZBase_ );
}

double globalLinearRegression::phi( const Stub* stub ) {
  return digi( reco::deltaPhi( stub->phi() - phiCentre_ - HTHelixRPhi_.second - HTHelixRPhi_.first * rPhi( stub ), 0. ), PhiBase_ );
}

double globalLinearRegression::z( const Stub* stub ) {
  return digi( stub->z() - zCentre_ - HTHelixRZ_.second - HTHelixRZ_.first * rZ( stub ), RZBase_ );
}

bool globalLinearRegression::checkValidity() {

  valid_ = true;
  if ( nLayers_ < minStubLayersRed_ )
    valid_ = false;
  if ( nPSLayers_ < minPSLayers_ )
    valid_ = false;
  return valid_;

}

void globalLinearRegression::calcHelix( ) {

  if ( oldNumStubs_ != NumStubs_ )
    nIterations_++;
  resetSums();
  for ( auto layer : stubMap_ ) if ( layer.second.size() > 0 ) {
    phiMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    rPhiMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    rZMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    zMinMax_ = std::make_pair( std::numeric_limits< double >::infinity(), - std::numeric_limits< double >::infinity() );
    sumOnlyRPhi_ = true;
    for ( auto stub : layer.second ) {
      RPhi_ = rPhi( stub );
      Phi_ = phi( stub );
      RZ_ = rZ( stub );
      Z_ = z( stub );
      rPhiMinMax_ = std::make_pair( std::min( rPhiMinMax_.first, RPhi_ ), std::max( rPhiMinMax_.second, RPhi_ ) );
      rZMinMax_ = std::make_pair( std::min( rZMinMax_.first, RZ_ ), std::max( rZMinMax_.second, RZ_ ) );
      phiMinMax_ = std::make_pair( std::min( phiMinMax_.first, Phi_ ), std::max( phiMinMax_.second, Phi_ ) );
      if ( stub->psModule() or not onlyPSForRZHelix_ ) {
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

void globalLinearRegression::calcResidual() {

  RPhiResiduals_.clear();
  RZResiduals_.clear();
  residuals_.clear();
  relativResiduals_.clear();
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
      RZSigma_ = stub->stripLength();
      if ( stub->psModule() and state_ == generalCleaning ) {
        RZSigma_ = 1.;
      }
      if ( not stub->barrel() )
        RZSigma_ *= std::fabs( HTHelixRZ_.first );
      RZResid_ /= RZSigma_;
      RZResiduals_[ layer.first ].push_back( RZResid_ );
      resid_ = ( RPhiResid_ + RZResid_ ) / 2.;
      residuals_[ layer.first ].push_back( resid_ );
      minRPhiResid_ = std::min( minRPhiResid_, RPhiResid_ );
      minRZResid_ = std::min( minRZResid_, RZResid_ );
      minResid_ = std::min( minResid_, resid_ );
    }
    for ( unsigned int i = 0; i < layer.second.size(); i++ )
      relativResiduals_[ layer.first ].push_back( residuals_[ layer.first ][ i ] - minResid_ );
  }

}

void globalLinearRegression::findLargestResidual() {

  std::map< unsigned int, std::vector< double > > residuals = residuals_;
  if ( useRelativResiduals_ )
    residuals = relativResiduals_;
  largestResid_.second.second = -1.;
  for ( auto layer : residuals ) for ( unsigned int i = 0; i < layer.second.size(); i++ ) {
    const Stub* stub( stubMap_[ layer.first ][ i ] );
    if ( stub->psModule() and onlySSKilling_ )
      continue;
    resid_ = residuals[ layer.first ][ i ];
    if ( resid_ > largestResid_.second.second )
      largestResid_ = std::make_pair( layer.first, std::make_pair( stub, resid_ ) );
  }

}

void globalLinearRegression::killLargestResidual() {

  if ( goodCandidate_ )
    iterDebug();
  layerID_ = largestResid_.first;
  const Stub* stub( largestResid_.second.first );
  stubMap_[ layerID_ ].erase( std::remove( stubMap_[ layerID_ ].begin(), stubMap_[ layerID_ ].end(), stub ), stubMap_[ layerID_ ].end() );
  stubs_.erase( std::remove( stubs_.begin(), stubs_.end(), stub ), stubs_.end() );
  nLayers_ = Utility::countLayers( settings_, stubs_ );
  nPSLayers_ = Utility::countLayers( settings_, stubs_, false, true );
  NumStubs_--;
  tp_ = Utility::matchingTP( settings_, stubs_, nMatchedLayersBest_, matchedStubsBest_ );
  if ( goodCandidate_ and tpHT_ != tp_ ) {
    goodCandidate_ = false;
    lostMatchingState_ = state_; 
 }

}

void globalLinearRegression::calcChiSq() {

  chiSq_ = 0.;
  for( auto layer : residuals_ ) for( auto resid : layer.second )
    chiSq_ += resid;
  chiSq_ /= NumStubs_;

}

void globalLinearRegression::resetSums() {

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

void globalLinearRegression::updateSums() {

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

void globalLinearRegression::calcLinearParameter() {

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

double globalLinearRegression::sectorConstrain( double val, double limit ) {

  if ( val > limit / 2. )
    return limit / 2.;
  else if ( val < - limit / 2. )
    return - limit / 2.;
  else
    return val;

}

bool globalLinearRegression::sectorCheck( double val, double limit ) {

  bool insideSector = (fabs(val) < limit/2.);
  return insideSector;
}

void globalLinearRegression::initDebug( const L1track3D& l1track3D ) {

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

void globalLinearRegression::iterDebug() {

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
