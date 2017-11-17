///=== This is the global Linear Regression for 4 helix parameters track fit algorithm.

///=== Written by: Thomas Schuh

#ifndef __GLOBALLINEARREGRESSION__
#define __GLOBALLINEARREGRESSION__

#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"

#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include <vector>
#include <sstream>
#include <string>


class globalLinearRegression : public TrackFitGeneric {

public:

    globalLinearRegression(const Settings* settings, const uint nPar) : TrackFitGeneric( settings ), settings_( settings ), nPar_( nPar ) {};

    virtual ~globalLinearRegression() {};

    virtual void initRun();

    L1fittedTrack fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg);

    std::string getParams() { return "globalLinearRegression"; };

protected:

  void initFit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg );
  L1fittedTrack createTrack( const L1track3D& l1track3D );
  void initDebug( const L1track3D& l1track3D );
  void iterDebug();
  bool checkValidity();
  void calcHelix();
  void calcResidual();
  void findLargestResidual();
  void killLargestResidual();
  void calcChiSq();
  void resetSums();
  void updateSums();
  void calcLinearParameter();
  double rPhi( const Stub* stub );
  double rZ( const Stub* stub );
  double phi( const Stub* stub );
  double z( const Stub* stub );
  double digi( double val, double base );
  double sectorConstrain( double val, double base );
  bool sectorCheck( double val, double base );

  // LinearRegression(const Settings* settings, const uint nPar)
  const Settings* settings_;
  uint nPar_;
  int nIterations_;
  enum t_state { generalCleaning, layerPurging, layerKilling, candidateKilling };
  const char* t_state_names[ candidateKilling + 1 ] = { "generalCleaning", "layerPurging", "layerKilling", "candidateKilling" };
  t_state state_;
  t_state lostMatchingState_;
  std::map< t_state, int > stateCalls_;

  // fit()
  bool onlyPSForRZHelix_;
  bool onlySSKilling_;
  bool useRelativResiduals_;
  double residualCut_;
  bool killLargestResid_;
  
  // initRun()
  double invPtToDPhi_;
  unsigned int numPhiSectors_;
  std::vector< double > etaRegions_;
  double zSectorSize_;
  double phiSectorSize_;
  double qOverPtSectorSize_;
  double tanLambdaSectorSize_;
  double chosenRofPhi_;
  double chosenRofZ_;
  unsigned int minPSLayers_;
  unsigned int minStubLayers_;
  unsigned int minStubLayersRed_;
  double houghMinPt_;
  double beamWindowZ_;
  bool digitizeLR_;
  double PhiPrecision_;
  double RPrecision_;
  double ZPrecision_;
  unsigned int ZSlopeWidth_;
  unsigned int ZInterceptWidth_;
  double HTMBinBase_;
  double HTCBinBase_;
  double PhiBase_;
  double RPhiBase_;
  double ZSlopeBase_;
  double ZInterceptBase_;
  double ZBase_;
  double RZBase_;
  bool debug_;

  // initFit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg )
  unsigned int iPhiSec_;
  unsigned int iEtaReg_;
  double phiCentre_;
  double zCentre_;
  std::vector< const Stub* > stubs_, assocStubs_;
  std::map< unsigned int, std::vector< const Stub* > > stubMap_;
  std::pair< double, double > HTHelixRPhi_;
  std::pair< double, double > HTHelixRZ_;
  unsigned int houghNbinsPt_;
  unsigned int houghNbinsPhi_;
  double binSizeQoverPtAxis_;
  double binSizePhiTrkAxis_;
  std::stringstream debugStream_;
  bool goodCandidate_;
  std::pair< double, double > tpHelixRPhi_;
  std::pair< double, double > tpHelixRZ_;
  unsigned int NumStubs_;
  unsigned int oldNumStubs_;
  unsigned int layerID_;
  unsigned int nLayers_;
  unsigned int nPSLayers_;

  // checkValidity( const L1track3D& l1track3D )
  bool valid_;

  // calcHelixRPhi(), calcHelixRZ()
  std::pair< double, double > rPhiMinMax_;
  std::pair< double, double > rZMinMax_;
  std::pair< double, double > phiMinMax_;
  std::pair< double, double > zMinMax_;
  double RPhi_;
  double RZ_;
  double Phi_;
  double Z_;
  bool seedingStub_;
  bool sumOnlyRZ_;
  bool sumOnlyRPhi_;
  double slopeRPhi_;
  double slopeRZ_;
  double interceptRPhi_;
  double interceptRZ_;
  std::pair< double, double > helixRPhi_;
  std::pair< double, double > helixRZ_;

  // calcResidualRPhi(), calcResidualRZ()
  double RPhiSigma_;
  double RZSigma_;
  double resid_;
  double RPhiResid_;
  double RZResid_;
  double minRPhiResid_;
  double minRZResid_;
  double minResid_;
  std::map< unsigned int, std::vector< double > > residuals_;
  std::map< unsigned int, std::vector< double > >  relativResiduals_;
  std::map< unsigned int, std::vector< double > > RPhiResiduals_;
  std::map< unsigned int, std::vector< double > > RZResiduals_;
  std::pair< unsigned int, std::pair< const Stub*, double > > largestResid_;

  // calcChiSq()
  float chiSq_;

  // createTrack( const L1track3D& l1track3D )
  std::map< std::string, double > trackParams_;

  // resetSums(), updateSums()
  int NRPhi_;
  int NRZ_;
  double RPhiTimesPhiSum_, RZTimesZSum_, RPhiSum_, RZSum_, PhiSum_, ZSum_, RPhiSquareSum_, RZSquareSum_;

  // calcLinearParameter()
  double denominatorRPhi_;
  double denominatorRZ_;
  double reciprocalDenominatorRPhi_;
  double reciprocalDenominatorRZ_;

  // initDebug()
  unsigned int nMatchedLayersBest_;
  std::vector< const Stub* > matchedStubsBest_;
  const TP* tp_;
  const TP* tpHT_;

};
#endif
