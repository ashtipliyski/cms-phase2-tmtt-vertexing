///=== This is the Kalman Combinatorial Filter for 4 helix parameters track fit algorithm.
 
 
#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsComb.h"
#include "TMTrackTrigger/TMTrackFinder/interface/kalmanState.h"
#include "TMTrackTrigger/TMTrackFinder/interface/StubCluster.h"
#define CKF_DEBUG
 
static unsigned nlayer_eta[25] = 
  { 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 7, 7, 7,
    7, 7, 7, 7, 6, 6, 6, 6, 6, 6};

static double matx_outer[25] = {
  0.16, 0.17, 0.18, 0.19, 0.20, 
  0.21, 0.26, 0.22, 0.26, 0.38,
  0.41, 0.40, 0.44, 0.50, 0.54,
  0.60, 0.44, 0.48, 0.60, 0.68,
  0.50, 0.48, 0.64, 0.39, 0.20
};

static double matx_inner[25] = {
  0.14, 0.1, 0.1, 0.1, 0.1, 
  0.1, 0.1, 0.1, 0.1, 0.1, 
  0.12, 0.1, 0.1, 0.1, 0.15,
  0.20, 0.25, 0.25, 0.3, 0.3,
  0.35, 0.40, 0.40, 0.6, 0.6
};

static double wrapRadian( double t ){

  if( t > 0 ){
    while( t > M_PI ) t-= 2*M_PI; 
  }
  else{
    while( t < - M_PI ) t+= 2*M_PI; 
  }
  return t;
}



KF4ParamsComb::KF4ParamsComb(const Settings* settings, const uint nPar, const string &fitterName ) : L1KalmanComb(settings, nPar, fitterName ){

  hdxmin[0] = -1.1e-4;
  hdxmax[0] = +1.1e-4;
  hdxmin[1] = -6.e-3;
  hdxmax[1] = +6.e-3;
  hdxmin[2] = -4.1;
  hdxmax[2] = +4.1;
  hdxmin[3] = -6.;
  hdxmax[3] = +6.;

  hxmin[0] = -0.3 * 0.0057;
  hxmax[0] = +0.3 * 0.0057;
  hxmin[1] = -0.3;
  hxmax[1] = +0.3;
  hxmin[2] = -120;
  hxmax[2] = +120;
  hxmin[3] = -6.;
  hxmax[3] = +6.;

  hddMeasmin[1] = -1.e1;
  hddMeasmax[1] = +1.e1;

  hresmin[1] = -0.5;
  hresmax[1] = +0.5;

  hresmin[1] = -10.;
  hresmax[1] = +10.;


  hxaxtmin[0] = -1.e-3;
  hxaxtmax[0] = +1.e-3;
  hxaxtmin[1] = -1.e-1;
  hxaxtmax[1] = +1.e-1;
  hxaxtmin[2] = -10.;
  hxaxtmax[2] = +10.;
  hxaxtmin[3] = -1.e-0;
  hxaxtmax[3] = +1.e-0;
}


std::map<std::string, double> KF4ParamsComb::getTrackParams(const kalmanState *state )const{

  std::vector<double> x = state->xa();
  std::map<std::string, double> y;
  y["qOverPt"] = x.at(INV2R) / getSettings()->invPtToInvR() * 2.; 
  y["phi0"] = wrapRadian( x.at(PHI0) + sectorPhi() );
  y["z0"] = x.at(Z0);
  y["t"] = x.at(T);
  y["d0"] = 0;
  return y;
}
 
/* The Kalman measurement matrix
 * Here I always measure phi(r), and z(r) */
TMatrixD KF4ParamsComb::H(const StubCluster* stubCluster)const{
  TMatrixD h(2, 4);
  h(PHI,INV2R) = -stubCluster->r();
  h(PHI,PHI0) = 1;
  h(Z,Z0) = 1;
  h(Z,T) = stubCluster->r();
  return h;
}


TMatrixD KF4ParamsComb::dH(const StubCluster* stubCluster)const{

  double dr(0);
  if(stubCluster->layerId() > 10){
    dr = stubCluster->sigmaZ();
  }

  TMatrixD h(2, 4);
  h(PHI,INV2R) = -dr;
  h(Z,T) = dr;

  return h;
}
 
/* Seed the state vector */
std::vector<double> KF4ParamsComb::seedx(const L1track3D& l1track3D)const{

  std::vector<double> x(nPar_);
  x[INV2R] = getSettings()->invPtToInvR() * l1track3D.qOverPt()/2;
  x[PHI0]  = wrapRadian( l1track3D.phi0() - sectorPhi() );
  x[Z0]    = l1track3D.z0();
  x[T]     = l1track3D.tanLambda();
    
  return x;
}

/* Seed the covariance matrix */
TMatrixD KF4ParamsComb::seedP(const L1track3D& l1track3D)const{
  TMatrixD p(4,4);

  double c = getSettings()->invPtToInvR() / 2; 

  if ( getSettings()->numEtaRegions() == 18 ) { 
      
    // optimised for 18x2 with additional error factor in pt/phi to avoid pulling towards wrong HT params
    p(INV2R,INV2R) = 0.0157 * 0.0157 * c * c * 4; 
    p(PHI0,PHI0) = 0.0051 * 0.0051 * 4; 
    p(Z0,Z0) = 5.0 * 5.0; 
    p(T,T) = 0.25 * 0.25;
      
  } else {
      
    // choose large errors
    p(INV2R,INV2R) = 0.0157 * 0.0157 * c * c * 10; 
    p(PHI0,PHI0) = 0.0051 * 0.0051 * 10; 
    p(Z0,Z0) = 5.0 * 5.0; 
    p(T,T) = 0.25 * 0.25 * 10;
      
  }

  return p;
}

/* The forecast matrix
 * (here equals identity matrix) */
TMatrixD KF4ParamsComb::F(const StubCluster* stubCluster, const kalmanState *state )const{
  TMatrixD F(4,4);
  for(int n = 0; n < 4; n++)
    F(n, n) = 1;
  return F;
}

/* the vector of measurements */
std::vector<double> KF4ParamsComb::d(const StubCluster* stubCluster )const{
  std::vector<double> meas;
  meas.resize(2);
  meas[0] = wrapRadian( stubCluster->phi() - sectorPhi() );
  meas[1] = stubCluster->z();
  return meas;
}


TMatrixD KF4ParamsComb::PddMeas(const StubCluster* stubCluster, const kalmanState *state )const{

  double inv2R = 1.05 * getSettings()->invPtToInvR() / ( 2 * state->candidate().pt() ) ; // factor 1.05 improves eff, alternatively use state->xa().at(INV2R)
  double inv2R2 = inv2R * inv2R;

  double tanl = 0.9 * state->xa().at(T);  // factor of 0.9 improves rejection
  double tanl2 = tanl * tanl; 

  TMatrixD p(2,2);

  double vphi(0);
  double vz(0);

  // consider error due to integerisation only for z (r in encap) coord when enabled
  double err_digi2(0);
  if (getSettings()->enableDigitize()) err_digi2 = 0.15625 * 0.15625 / 12.0;

  double a = stubCluster->sigmaX() * stubCluster->sigmaX();
  double b = stubCluster->sigmaZ() * stubCluster->sigmaZ() + err_digi2;
  double invr2 = (1.0 / stubCluster->r()) * (1.0 / stubCluster->r());

  if ( stubCluster->barrel() ) {

    vphi = (a * invr2);
    vz = b;

  } else {

    vphi = (a * invr2) + (b * inv2R2);
    vz = (b * tanl2);

  }

  p(0,0) = vphi;
  p(1,1) = vz;

  return p;

}

/* State uncertainty */
TMatrixD KF4ParamsComb::PxxModel( const kalmanState *state, const StubCluster *stubCluster )const
{

  TMatrixD p(4,4);

  if( getSettings()->kalmanMultiScattFactor() ){

    unsigned i_eta = abs( stubCluster->eta() / 0.1 );
    if( i_eta > 24 ) i_eta = 24;
    double dl = matx_outer[i_eta] / nlayer_eta[i_eta];

    unsigned stub_itr = state->nextLayer();

    const kalmanState * last_update_state = state->last_update_state();
    unsigned last_itr(1);
    if( last_update_state ) last_itr = last_update_state->nextLayer();
    dl = ( stub_itr - last_itr ) * dl; 

    if( dl ){
      std::map<std::string, double> y = getTrackParams( state );
      double dtheta0 = 1./sqrt(3) * 0.0136 * fabs(y["qOverPt"]) * sqrt(dl)*( 1+0.038*log(dl) ); 
      dtheta0 *= getSettings()->kalmanMultiScattFactor();
      p(1,1) = dtheta0 * dtheta0; 
    }
  }

  return p;
}


std::string KF4ParamsComb::getParams(){
  return "KF4ParamsComb";
}


bool KF4ParamsComb::isGoodState( const kalmanState &state )const
{

  unsigned nStubLayers = state.nStubLayers();
  bool goodState( true );

  // todo : make configurable

  // state parameter selections
  if( nStubLayers >= 2 ){
      
    double z0=fabs( state.xa()[Z0] ); 
    if( z0 > 15. ) goodState = false;
      
    double pt=fabs( getSettings()->invPtToInvR() / (2*state.xa()[INV2R]) ); 
    if( pt < 2.9 ) goodState = false;
      
  }

  // chi2 selections
  if( nStubLayers == 2 ) {
    if (state.chi2() > 15.0) goodState=false;
  }
  else if ( nStubLayers == 3 ) {
    if (state.chi2() > 100.0) goodState=false;
  }
  else if( nStubLayers == 4 ) {
    if (state.chi2() > 320.0) goodState=false;
  }

  return goodState;
}

