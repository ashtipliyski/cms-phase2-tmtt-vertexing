///=== This is the simple linear regression with 4 helix parameters (qOverPt, phiT, z0, tanLambda) track fit algorithm.

///=== Written by: Davide Cieri (davide.cieri@stfc.ac.uk)

#ifndef __SIMPLELR__
#define __SIMPLELR__

#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"

#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include <vector>
#include <sstream>
#include <string>


class SimpleLR : public TrackFitGeneric {

public:
	SimpleLR(const Settings* settings) : TrackFitGeneric(settings), settings_(settings) {};

	virtual ~SimpleLR() {};

	virtual void initRun();

    L1fittedTrack fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg);

    std::string getParams() { return "SimpleLR"; }

protected:
	const Settings* settings_;

	double phiSectorWidth_;
	double phiSectorCentre;

	double phiMult_;
	double rTMult_;
	double zMult_;
	double qOverPtMult_;
	double phiTMult_;
	double z0Mult_;
	double tanLambdaMult_;
	double numeratorPtMult_;
	double numeratorZ0Mult_;
	double numeratorLambdaMult_;
	double numeratorPhiMult_;
	double denominatorMult_; 
	double chi2Mult_;
	double chi2cut_;
	double invPtToDPhi_;
	double chosenRofPhi_;

	bool                 digitize_;
	unsigned int         dividerBitsHelix_;
	unsigned int         dividerBitsHelixZ_;
	unsigned int         dividerBitsChi2_;
	unsigned int         ShiftingBits_;
	unsigned int 	     shiftingBitsPt_;
	unsigned int 	     shiftingBitsz0_;
	unsigned int         shiftingBitsLambda_;

	
};

#endif