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
  P_<< 1.0,0,0,0,0,
	   0,1.0,0,0,0,
	   0,0,1.0,0,0,
	   0,0,0,1.0,0,
	   0,0,0,0,1.0;


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;//0.075 in lecture

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;//0.1 in lecture

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_= true;
  use_laser_= true;
  use_radar_= true;

  ///* Weights of sigma points
  weights_= VectorXd::Zero(2*n_aug_+1);


  ///* Augmented state dimension
  int n_aug_ = 7;

  ///* Sigma point spreading parameter
  double lambda_ =3 - n_x_;

  //predicted sigma point matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_,2 * n_aug_+1)

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
	            x_(0) = meas_package.raw_measurements_(0);
	            x_(1) = meas_package.raw_measurements_(1);
	        } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	            // Initialize radar
	            float ro = meas_package.raw_measurements_(0);
	            float phi = meas_package.raw_measurements_(1);
	            float rho_dot = meas_package.raw_measurements_(2);
	            float px = ro * cos(phi);
	            float py = = ro * sin(phi);
	            x_<<px,py,0,0,0;
	        }

	        time_us_ = meas_package.timestamp_;
	        is_initialized_ = true;
	        return;
	    }

	    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	    time_us_ = meas_package.timestamp_;

	    Prediction(dt);

	    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
	        UpdateLidar(meas_package);
	    } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	    	UpdateRadar(meas_package);
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

  MatrixXd Xsig_aug = MatrixXd();

  //get augmented sigma points
  AugmentedSigmaPoints(&Xsig_aug);
  // get sigma points
  Xsig_pred = MatrixXd();
  SigmaPointPrediction(delta_t, Xsig_aug, &Xsig_pred_);

  //get mean and covariance
  PredictMeanandCovariance(&x_, &P_);


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
}
void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out){
	//create vector fro predicted state
	VectorXd x = VectorXd(n_x_);
	//create covariant matrix for prediction
	MatrixXd P = MatrixXd(n_x_,n_x_);

	//predicted state mean
	x.fill(0.0);
	for(int i = 0;i<2*n_aug_ + 1;i++){
		x = x+ weights_(i)*Xsig_pred_.col(i);
	}
	//predicted state covariance matrix
	P.fill(0.0);
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
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //make square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_,2*n_aug_+1);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0;i<n_aug_;i++){
	  Xsig_aug.col(i+1)    = x_aug + sqrt(lambda_+n_aug_)*L.col(i);
	  Xsig_aug.col(i+1+n_aug_)    = x_aug - sqrt(lambda_+n_aug_)*L.col(i);
  }

  *Xsig_out = Xsig_aug;



}


void UKF::SigmaPointPrediction(double delta_t,MatrixXd& Xsig_aug, MatrixXd* Xsig_out){

  //make matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_,2*n_aug_+1);

  //predicted sigma points
  for(int i = 0;i<2*n_aug_+1;i++){
	  //extrace values for better readability
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
	  if (fabs(yawd) >0.001){
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
