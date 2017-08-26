#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  n_x_ = 5;//state vector dimension
  x_ = VectorXd(n_x_);
  //first neasurement
  x_ << 1,1,1,1,1;



  // initial covariance matrix
  P_ = MatrixXd(n_x_,n_x_);
  //initial covariance matrix
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .2;//originally 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .4;//originally 30

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.0175;//0.075 in lecture

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.1;//0.1 in lecture

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_= false;
  //P_ was originally 0.1 for the identity
 // P_<<0.1,0,0,0,0,
  //	  0,0.1,0,0,0,
//	  0,0,1,0,0,
//	  0,0,0,1,0,
//	  0,0,0,0,1;
  // below values came from Forum discussion
  P_<< 0.0043, -0.0013, 0.0030, -0.0022,-0.0020,
	 -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
	 0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
	 -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
	 -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

  H_laser_ = MatrixXd(2,n_x_);
  H_laser_.fill(0);
  H_laser_ << 1,0,0,0,0,
		      0,1,0,0,0;

  //measurement covariance matrix for the laser
  R_laser_ = MatrixXd(2,2);
  R_laser_.fill(0.0);
  R_laser_ << std_laspx_*std_laspx_,0,
		      0,std_laspy_*std_laspy_;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ =3 - n_x_;

  time_us_ = 0;
  ///* Weights of sigma points
  weights_= VectorXd(2*n_aug_+1);
  weights_.fill(0.0);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_[0] = weight_0;
  for(int i=1;i<2*n_aug_+1;i++){
	  double weight = 0.5/(n_aug_+lambda_);
	  weights_[i] = weight;
  }
  //NIS variables
  nis_laser_=0;
  nis_radar_= 0;
  //predicted sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_,2 * n_aug_+1);
  Xsig_pred_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if(!is_initialized_) {
	        // Initialize x_, P_, previous_time, anything else needed.
	        if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
	            // Initialize laser
	            x_(0) = meas_package.raw_measurements_[0];
	            x_(1) = meas_package.raw_measurements_[1];
	            x_(2) = 0;
	            x_(3) = 0;
	            x_(4) = 0;
	            std::cout<<"x_(0) = "<<x_(0)<<"\n";
	            std::cout<<"x_(1) = "<<x_(1)<<"\n";
	        } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	            // Initialize radar
	            float rho = meas_package.raw_measurements_[0];
	            float phi = meas_package.raw_measurements_[1];
	            float rho_dot = meas_package.raw_measurements_[2];
	            float px = rho * cos(phi);
	            float py = rho * sin(phi);
	            x_<<px,py,0,0,0;
	            std::cout<<"px = "<<px<<"\n";
	            std::cout<<"py = "<<py<<"\n";
	        }

	        time_us_ = meas_package.timestamp_;
	        is_initialized_ = true;
	        return;
	    }

	    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	    time_us_ = meas_package.timestamp_;
        std::cout<<"A x_ = "<<x_<<"\n";
	    Prediction(dt);

	    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
	        UpdateLidar(meas_package);
	        std::cout<<"LASER\n";
	    } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	    	UpdateRadar(meas_package);
	    	std::cout<<"RADAR\n";
	    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  std::cout<<"Prediction\n";
  MatrixXd Xsig_aug = MatrixXd();

  //get augmented sigma points
  AugmentedSigmaPoints(&Xsig_aug);
  // get sigma points
  Xsig_pred_ = MatrixXd();

  SigmaPointPrediction(delta_t, Xsig_aug, &Xsig_pred_);

  //get mean and covariance
  PredictMeanAndCovariance(&x_, &P_);
  std::cout<<"end of Prediction\n";
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the lidar NIS.
  */
  std::cout<<"UpdateLidar\n";
  VectorXd z = VectorXd(2);
  z<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1];

  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt*Si;
  std::cout<<"K = "<<K<<"\n";
  //new estimate
  std::cout<<"y = "<<y<<"\n";
  std::cout<<"x_ = "<<x_<<"\n";
  x_ = x_+(K*y);
  std::cout<<"x_ = "<<x_<<"\n";
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);
  P_ = (I-K*H_laser_)*P_;
  //update NIS for Lidar
  nis_laser_ = calculateNIS(z,z_pred,S);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */
  std::cout<<"UpdateRadar\n";
  MatrixXd ZSig_out = MatrixXd();
  VectorXd z_out = VectorXd();
  MatrixXd S_out = MatrixXd();

  ZSig_out.fill(0.0);
  z_out.fill(0.0);
  S_out.fill(0.0);

  PredictRadarMeasurement(&ZSig_out, &z_out,&S_out);

  float rho = meas_package.raw_measurements_[0];
  float phi = meas_package.raw_measurements_[1];
  float rho_dot = meas_package.raw_measurements_[2];

  VectorXd z = VectorXd(3);


  z << rho, phi, rho_dot;

  UpdateRadarState(ZSig_out, z_out, S_out, z);

}

void UKF::PredictRadarMeasurement(MatrixXd* ZSig_out, VectorXd* z_out, MatrixXd* S_out){
	std::cout<<"PredictRadar\n";
	int n_z = 3;//sets the measurement dimension
	const int nu_sigma_points_ = 2*n_aug_+1;
	MatrixXd Zsig = MatrixXd(n_z, nu_sigma_points_);

	//transform sigma points into measurement space
	for (int i = 0; i<nu_sigma_points_;i++){
		//extrace values for better readibility
		double p_x = Xsig_pred_(0,i);
		double p_y = Xsig_pred_(1,i);
		double v = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;
		//measurement model
		Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);//r
		Zsig(1,i) = atan2(p_y,p_x);//phi
		Zsig(2,i) = (p_x*v1 + p_y*v2)/sqrt(p_x*p_x+p_y*p_y);
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	for(int i=0;i<nu_sigma_points_;i++){
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);
	for(int i=0;i<nu_sigma_points_;i++){
		//residual
		VectorXd z_diff = Zsig.col(i)-z_pred;
		//normalize angle
		z_diff(1) = normalize_angle(z_diff(1));
		S=S+weights_(i)*z_diff*z_diff.transpose();

	}
	//add in measurement noise covariance
	MatrixXd R = MatrixXd(n_z,n_z);

	R<<     std_radr_*std_radr_,0,0,
			0,std_radphi_*std_radphi_,0,
			0,0,std_radrd_*std_radrd_;
	S=S+R;

	*ZSig_out = Zsig;
	*z_out = z_pred;
	*S_out = S;


}

void UKF::UpdateRadarState(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S, VectorXd& z){
    std::cout<<"UpdateRadarState\n";
	//set measurement dimension r,phi,r_dot
	int n_z=3;
	const int nu_sigma_points_ = 2 *n_aug_ +1;
    std::cout<<"x_ = "<<x_<<"\n";
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_,n_z);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i=0; i<nu_sigma_points_;i++){

		//residual
		VectorXd z_diff = Zsig.col(i)-z_pred;
		//angle normalize
		z_diff(1)=normalize_angle(z_diff(1));
		//state difference
		VectorXd x_diff = Xsig_pred_.col(i)-x_;
		//angle normalize again
		x_diff(3)=normalize_angle(x_diff(3));

		Tc = Tc + weights_[i]*x_diff*z_diff.transpose();
		//std::cout<<"z_diff(1) = "<<z_diff(1)<<"\n";
	}

	//Kalman gain K;
	MatrixXd K = Tc *S.inverse();
	//residual
	VectorXd z_diff = z-z_pred;
	//angle normalization
	z_diff(1) = normalize_angle(z_diff(1));
	//update state mean and covariance matrix

	x_ = x_ + K*z_diff;

	P_=P_-K*S*K.transpose();
	//calculate NIS for radar
	nis_radar_ = calculateNIS(z,z_pred, S);

}

double UKF::normalize_angle(double theta){
	//std::cout<<"Normalize Angle\n";
	while (theta>M_PI)  theta-=2.0*M_PI;
	while (theta<-M_PI) theta+=2.0*M_PI;
	return theta;
}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out){
	std::cout<<"PredictMean\n";

	//create vector for predicted state
	VectorXd x = VectorXd(n_x_);
	x.fill(0.0);

	//create covariant matrix for prediction
	MatrixXd P = MatrixXd(n_x_,n_x_);
	P.fill(0.0);
	std::cout<<"x = "<<x<<"\n";
	for(int i= 0;i<2*n_aug_ + 1;i++){
		x = x + weights_(i) * Xsig_pred_.col(i);
		std::cout<<"weights_("<<i<<") = "<<weights_(i)<<"\n";
		std::cout<<"Xsig_pred_.col("<<i<<") = "<<Xsig_pred_.col(i)<<"\n";
	}
	std::cout<<"x = "<<x<<"\n";
	//predicted state covariance matrix
	for(int i = 0;i<2*n_aug_+1;i++){
		//state difference
		VectorXd x_diff = Xsig_pred_.col(i)-x;
		//angle normalization
		x_diff(3) = normalize_angle(x_diff(3));

		P = P + weights_(i) * x_diff * x_diff.transpose();
	}
	*x_out = x;
	*P_out = P;

}
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){
  //create augmented mean state
  std::cout<<"Augmented\n";

  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //make square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_,2*n_aug_+1);

  Xsig_aug.col(0) = x_aug;

  for (int i = 0;i<n_aug_;i++){
	  Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_)*L.col(i);
	  Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_)*L.col(i);
  }

  *Xsig_out = Xsig_aug;



}


void UKF::SigmaPointPrediction(double delta_t,MatrixXd& Xsig_aug, MatrixXd* Xsig_out){

  //make matrix with predicted sigma points as columns
  std::cout<<"Sigma\n";
  MatrixXd Xsig_pred = MatrixXd(n_x_,2*n_aug_+1);


  //predicted sigma points
  for(int i = 0;i<2*n_aug_+1;i++){
	  //extract values for better readability
	  double p_x      = Xsig_aug(0,i);
	  double p_y      = Xsig_aug(1,i);
	  double v        = Xsig_aug(2,i);
	  double yaw      = Xsig_aug(3,i);
	  double yawd     = Xsig_aug(4,i);
	  double nu_a     = Xsig_aug(5,i);
	  double nu_yawdd = Xsig_aug(6,i);

	  //prediced state values
	  double px_p,py_p;

	  //avoiding division by zero
	  if (fabs(yawd) >0.001){//fabs computes absolute value
		  px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t)-sin(yaw));
		  py_p = p_y * v/yawd * (cos(yaw)- cos(yaw+yawd*delta_t));
	  }
	  else {
		  px_p = p_x + v*delta_t*cos(yaw);
		  py_p = p_y + v*delta_t*sin(yaw);
	  }

	  double v_p = v;
	  double yaw_p = yaw + yawd*delta_t;
	  double yawd_p = yawd;

	  //add noise
	  px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
	  px_p = px_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
	  v_p = v_p + nu_a*delta_t;

	  yaw_p = yaw_p + 05.*nu_yawdd*delta_t*delta_t;
	  yawd_p= yawd_p + nu_yawdd*delta_t;

	  //write predicted sigma points into right column
	  Xsig_pred(0,i) = px_p;
	  Xsig_pred(1,i) = py_p;
	  Xsig_pred(2,i) = v_p;
	  Xsig_pred(3,i) = yaw_p;
	  Xsig_pred(4,i) = yawd_p;

	  //write out result
	  *Xsig_out = Xsig_pred;

  }
}

double UKF::calculateNIS(VectorXd z, VectorXd z_pred, MatrixXd S){
	double nis = 0;
	VectorXd z_diff = z-z_pred;
	nis = z_diff.transpose()*S*z_diff;
	std::cout<<"nis = "<<nis<<"\n";
	return nis;
}
