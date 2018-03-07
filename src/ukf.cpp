#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  // cout << "Initializing" << endl;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.9;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  is_initialized_ = false;

  ///* time when the state is true, in us
  time_us_ = 0.0;

  ///* State dimension
  n_x_ = 5 ;

  ///* Augmented state dimension
  n_aug_ = 7 ;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;
  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ +1);

  //  Lidar Measurement Matrix
  H_ = MatrixXd(2, n_x_);
  H_.fill(0.0);
  H_(0, 0) = 1;
  H_(1, 1) = 1;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ *std_laspx_, 0,
  0, std_laspy_ *std_laspy_;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);;

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if( (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
      (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    ) {
      // cout << " Start - ProcessMeasurement " << endl;
      if (!is_initialized_) {

        time_us_ = meas_package.timestamp_;
        x_.fill(0.0);

        P_ = MatrixXd::Identity(5, 5);

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
            float ro     = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
            x_(0) = ro     * cos(phi);
            x_(1) = ro     * sin(phi);
        }else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }
        is_initialized_ = true;
        return;
      }
      float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
      time_us_ = meas_package.timestamp_;

      Prediction(dt);

      if (meas_package.sensor_type_ == MeasurementPackage::LASER  && use_laser_) {
        UpdateLidar(meas_package);
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        UpdateRadar(meas_package);
      }

      // cout << " End Initialization" << endl;
  }
  return;
}

/** Generating Sigma Points */
void UKF::GenerateSigmaPoints(MatrixXd *Xsig_out) {
  // cout << "Generating Sigma Points" << endl;
  MatrixXd Xsig( n_x_, 2 * n_x_ + 1 );

  MatrixXd A = P_.llt().matrixL();
  Xsig.col(0) = x_;

  for (int i = 0; i < n_x_; i++) {
    Xsig.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

  *Xsig_out = Xsig;
}


void UKF::AugmentSigmaPoints(MatrixXd *Xsig_out) {
  // cout << "Augment Sigma Points" << endl;
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  lambda_ = 3 - n_aug_;

  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_ )= P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for(int i=0; i < n_aug_; i++){
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  *Xsig_out = Xsig_aug;
}


void UKF::PredictSigmaPoints(MatrixXd *Xsig_aug, double delta_t) {
  // cout << "PredictSigmaPoints "<< endl;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  for(int i = 0; i< 2*n_aug_+1; i++) {
    //extract values for better readability
    double p_x = (*Xsig_aug)(0,i);
    double p_y = (*Xsig_aug)(1,i);
    double v = (*Xsig_aug)(2,i);
    double yaw = (*Xsig_aug)(3,i);
    double yawd = (*Xsig_aug)(4,i);
    double nu_a = (*Xsig_aug)(5,i);
    double nu_yawdd = (*Xsig_aug)(6,i);

    double px_p, py_p;

    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::SetWeights() {
  weights_.fill(0.0);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }
}

void UKF::PredictStateMean(){
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
}

void UKF::PredictStateCovarianceMatrix(){
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig(n_x_, n_x_);
  GenerateSigmaPoints(&Xsig);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ +1);
  AugmentSigmaPoints(&Xsig_aug);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  PredictSigmaPoints(&Xsig_aug, delta_t);
  SetWeights();
  PredictStateMean();
  PredictStateCovarianceMatrix();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void UKF::NormalizeAngles(VectorXd *v, int idx){
  while ((*v)(idx)> M_PI) {
      (*v)(idx) -= 2.*M_PI;
  }
  while ((*v)(idx)<-M_PI) {
      (*v)(idx) += 2.*M_PI;
  }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  VectorXd z = meas_package.raw_measurements_;
  int n_z = 3; // for rho, phi and phi_dot

  MatrixXd Zsig(n_z, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = v * cos(yaw);
    double v2 = v * sin(yaw);

    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
  }

  // mean vec
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngles(&z_diff, 1);
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_,0,0,
      0, std_radphi_*std_radphi_,0,
      0,0, std_radrd_*std_radrd_;
  S = S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngles(&z_diff, 1);
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngles(&x_diff, 1);
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;
  NormalizeAngles(&z_diff, 1);
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
