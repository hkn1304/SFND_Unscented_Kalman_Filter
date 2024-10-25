#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

  // set state dimension
  int n_x = 5;
  // set augmented dimension
  int n_aug = 7;
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  // define spreading parameter
  double lambda = 3 - n_x;

    // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_ = 0.2;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  previous_timestamp_=0.0;

  // initial state vector
  x_ = VectorXd::Zero(n_x);

  // initial covariance matrix
  P_ = MatrixXd::Zero(n_x, n_x);

  P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

    // create augmented mean vector
  x_aug = VectorXd::Zero(n_aug);

  // create augmented state covariance
  P_aug = MatrixXd::Zero(n_aug, n_aug);

  // create sigma point matrix
  Xsig_aug = MatrixXd::Zero(n_aug, 2 * n_aug + 1);

  // create sigma point matrix
  Xsig = MatrixXd::Zero(n_x, 2 * n_x + 1);

      // create matrix for sigma points in measurement space
  Zsig = MatrixXd::Zero(n_z, 2 * n_aug + 1);

  // mean predicted measurement
  z_pred_ = VectorXd::Zero(n_z);

  z_ = VectorXd::Zero(n_z);

  // Kalman gain matrix K
  K = MatrixXd::Zero(n_x, n_z);
  
  // measurement covariance matrix S
  S = MatrixXd::Zero(n_z,n_z);

  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd::Zero(n_x, 2 * n_aug + 1);

  weights_ = VectorXd::Zero(2*n_aug+1);
    // set weights
    weights_.fill(1.0 / (2.0 * (lambda + n_aug)));
    weights_(0) = lambda/(lambda+n_aug);

  R = MatrixXd::Zero(3, 3); // 3x3 matrix for radar
  Rl = MatrixXd::Zero(2, 2); // 2x2 matrix for lidar

  // R = (Eigen::Array3d(pow(std_radr_, 2), pow(std_radphi_, 2), pow(std_radrd_, 2)).matrix()).asDiagonal();
  R <<  pow(std_radr_,2),     0,          0,
                0,      pow(std_radphi_,2), 0,
                0,      0,      pow(std_radrd_,2);

  Rl << pow(std_laspx_,2), 0,
                0, pow(std_laspy_,2);

  H = Eigen::MatrixXd::Zero(n_z-1,n_x) ; 

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 1;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    //cout << "Kalman Filter Initialization " << endl;

    // set the state with the initial location and zero velocity
    x_ << 5.7441,
          1.3800,
          2.2049,
          0.5015,
          0.3528;

    //previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  else if ((!is_initialized_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
        // set the state with the initial location and zero velocity
      double rho=   meas_package.raw_measurements_[0]; 
      double phi=   meas_package.raw_measurements_[1]; 
      double rhodot=meas_package.raw_measurements_[2]; 
              
    x_ << 5.7441,
          1.3800,
          2.2049,
          0.5015,
          0.3528;
    //previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

    // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  std::cout << "Delta_t:" << delta_t << std::endl;

  this->Prediction(delta_t/*VectorXd* z_out, MatrixXd* S_out*/);
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    z_ = VectorXd(n_z-1);
    this->UpdateLidar(meas_package);
  }


  this->PredictRadarMeasurement(/*VectorXd* z_out, MatrixXd* S_out*/);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    z_ = VectorXd(n_z);
    this->UpdateRadar(meas_package);
  }

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  this->AugmentedSigmaPoints();
  this->SigmaPointPrediction(delta_t/*MatrixXd* Xsig_out*/);
  this->PredictMeanAndCovariance(/*VectorXd* x_out, MatrixXd* P_out*/);
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

    z_ <<     meas_package.raw_measurements_[0], // px
              meas_package.raw_measurements_[1]; // py

              
    H << 1,0,0,0,0,
         0,1,0,0,0;

    VectorXd y = z_ - H * x_;
    MatrixXd Sl = H * P_ * H.transpose() + Rl;
    MatrixXd Kl =  P_ * H.transpose() * Sl.inverse();
    MatrixXd I = MatrixXd::Identity(n_x,n_x); 
    // new state
    x_ = x_ + (Kl * y);
    P_ = (I - Kl * H) * P_;

    std::cout << "x=" << std::endl <<  x_ << std::endl;
    std::cout << "P=" << std::endl <<  P_ << std::endl;

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   * 
   *
    meas_package.sensor_type_ = MeasurementPackage::RADAR;
    meas_package.raw_measurements_ = VectorXd(3);
    meas_package.raw_measurements_ << marker.rho, marker.phi, marker.rho_dot;
    meas_package.timestamp_ = timestamp;
  */

    this->PredictRadarMeasurement(/*VectorXd* z_out, MatrixXd* S_out*/);

    double rho=   meas_package.raw_measurements_[0]; 
    double phi=   meas_package.raw_measurements_[1]; 
    double rhodot=meas_package.raw_measurements_[2]; 

    z_ <<  rho,   // px
            phi,   // py
            rhodot;// v

    // z_ <<  rho * cos(phi),   // px
    //         rho * sin(phi),   // py
    //         rhodot / cos(phi);// v
                  

    this->UpdateRadarState(/*VectorXd* x_out, MatrixXd* P_out*/);


}

void UKF::AugmentedSigmaPoints(/*MatrixXd* Xsig_out*/) {

  /**
   * Student part begin
   */
  // create augmented mean state
  x_aug.head(5)= x_;
  x_aug.tail(n_aug - n_x).setZero();

  // create augmented covariance matrix
  // P_aug.topLeftCorner(n_x, n_x)=P_;
  // // Create a 2x2 diagonal matrix by scaling the identity matrix
  // P_aug.block(5, 5, 2, 2) = Eigen::Matrix2d::Identity() * Eigen::Vector2d(pow(std_a,2), pow(std_yawdd,2)).asDiagonal();
    P_aug.setZero();
    P_aug.topLeftCorner(n_x, n_x) = P_;
    P_aug(n_x, n_x) = 0.04;//std_a_ * std_a_;
    P_aug(n_x + 1, n_x + 1) = 0.04;//std_yawdd_ * std_yawdd_;
    std::cout << "P_aug= " << P_aug << std::endl;


  // create square root matrix
  MatrixXd A_aug= P_aug.llt().matrixL();
  // create augmented sigma points
    // Create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    double sqrt_lambda_n_aug = sqrt(lambda + n_aug);
    for (int i = 0; i < n_aug; ++i) {
        Xsig_aug.col(i + 1) = x_aug + sqrt_lambda_n_aug * A_aug.col(i);
        Xsig_aug.col(i + 1 + n_aug) = x_aug - sqrt_lambda_n_aug * A_aug.col(i);
    }
    std::cout << "Xsig_aug= " << Xsig_aug << std::endl;
  /**
   * Student part end
   */

  // write result
  //*Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(double delta_t/*MatrixXd* Xsig_out*/) {

//double delta_t = 0.1; // time diff in sec

  // Loop over sigma points
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
      // Extract values for better readability
      double px      = Xsig_aug(0,i);
      double py      = Xsig_aug(1,i);
      double v       = Xsig_aug(2,i);
      double psi     = Xsig_aug(3,i);
      double psidot  = Xsig_aug(4,i);
      double nu_a    = Xsig_aug(5,i);
      double nu_psi  = Xsig_aug(6,i);

      // Predicted state values
      double px_p, py_p;

      // Avoid division by zero
      if (fabs(psidot) > 0.001) {
          px_p = px + (v/psidot) * (sin(psi + psidot * delta_t) - sin(psi));
          py_p = py + (v/psidot) * (-cos(psi + psidot * delta_t) + cos(psi));
      } else {
          px_p = px + v * cos(psi) * delta_t;
          py_p = py + v * sin(psi) * delta_t;
      }

      double v_p = v;
      double psi_p = psi + psidot * delta_t;
      double psidot_p = psidot;

      // Add noise
      px_p += 0.5 * delta_t * delta_t * cos(psi) * nu_a;
      py_p += 0.5 * delta_t * delta_t * sin(psi) * nu_a;
      v_p  += delta_t * nu_a;
      psi_p += 0.5 * delta_t * delta_t * nu_psi;
      psidot_p += delta_t * nu_psi;

      // Write predicted sigma points into the correct column
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = psi_p;
      Xsig_pred_(4,i) = psidot_p;
}

  /**
   * Student part end
   */

  // print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;

  // write result
  //*Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(/*VectorXd* x_out, MatrixXd* P_out*/) {

  /**
   * Student part begin
   */


    // For predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; ++i) {
        x_ += weights_(i) * Xsig_pred_.col(i);
    }
  // predict state covariance matrix
    
    P_.fill(0.0);
    // SUMMARY:
    // First take each column in Xsig_pred and subtract from x
    // Then multiply with its transpose so that n_x by n_x element obtained
    // Last multiply with column of the weight corresponding to column of Xsig_pred
    // in this case first column with first element, second column with second element of weight
    // Iterate through all rows and columns of the matrix
    for (int col = 0; col < Xsig_pred_.cols(); ++col) {
        // Access each element using (row, col)
        VectorXd x_diff = (Xsig_pred_.col(col)- x_);
            //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ += weights_(col)*x_diff* x_diff.transpose();
    } 

  /**
   * Student part end
   */

  // print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x_ << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P_ << std::endl;

  // // write result
  // *x_out = x_;
  // *P_out = P_;
}

void UKF::PredictRadarMeasurement(/*VectorXd* z_out, MatrixXd* S_out*/) {

  /**
   * Student part begin
   */

  // transform sigma points into measurement space
  for (int i=0; i < 2 * n_aug + 1; ++i)
  {
      double px= Xsig_pred_.col(i)(0);
      double py= Xsig_pred_.col(i)(1);
      double v= Xsig_pred_.col(i)(2);
      double psi= Xsig_pred_.col(i)(3);
      double psidot= Xsig_pred_.col(i)(4);
      
      Zsig(0,i)   = sqrt( pow(px,2)+pow(py,2) );
      Zsig(1,i)   = atan(py/px);
      Zsig(2,i)= (px*cos(psi)*v+py*sin(psi)*v) / sqrt( pow(px,2)+pow(py,2) ); 

  }
  // calculate mean predicted measurement

    // Iterate through all rows and columns of the matrix
    Eigen::VectorXd zp= Eigen::VectorXd::Zero(n_z);
    for (int row = 0; row < Zsig.rows(); ++row) {

        for (int col = 0; col < Zsig.cols(); ++col) {
            // Access each element using (row, col)
            zp(row) += Zsig(row, col)*weights_(col);
        }
        
    }
    z_pred_= zp;
    
    
  // predict state covariance matrix
    MatrixXd Zp = MatrixXd::Zero(n_z, n_z);
    // SUMMARY:
    // First take each column in Xsig_pred and subtract from x
    // Then multiply with its transpose so that n_x by n_x element obtained
    // Last multiply with column of the weight corresponding to column of Xsig_pred
    // in this case first column with first element, second column with second element of weight
    // Iterate through all rows and columns of the matrix
    for (int col = 0; col < Zsig.cols(); ++col) {
        // Access each element using (row, col)
        Zp += weights_(col)*(Zsig.col(col)- z_pred_)*( (Zsig.col(col)- z_pred_) ).transpose();
    } 

    
  // calculate innovation covariance matrix S
  S= Zp+R;
  /**
   * Student part end
   */

  // // print result
  // std::cout << "z_pred: " << std::endl << z_pred_ << std::endl;
  // std::cout << "S: " << std::endl << S << std::endl;

  // // write result
  // *z_out = z_pred_;
  // *S_out = S;
}

void UKF::UpdateRadarState(/*VectorXd* x_out, MatrixXd* P_out*/) {

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x, n_z);

  /**
   * Student part begin
   */

  // calculate cross correlation matrix

    // SUMMARY:
    // First take each column in Xsig_pred and subtract from x
    // Then multiply with its transpose so that n_x by n_x element obtained
    // Last multiply with column of the weight corresponding to column of Xsig_pred
    // in this case first column with first element, second column with second element of weight
    // Iterate through all rows and columns of the matrix
    for (int col = 0; col < Xsig_pred_.cols(); ++col) {
        // Access each element using (row, col)
        
        VectorXd z_diff = Zsig.col(col) - z_pred_;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
        // state difference
        VectorXd x_diff = Xsig_pred_.col(col) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc += weights_(col)*x_diff*z_diff.transpose();
    } 

    
  // calculate Kalman gain K;

  K = Tc*S.inverse();
  
  // update state mean and covariance matrix
  VectorXd z_diff = z_-z_pred_;

  // Angle normalization for the bearing angle
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;


  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S*K.transpose();
  /**
   * Student part end
   */

  // // print result
  // std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  // std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

  // // write result
  // *x_out = x_;
  // *P_out = P_;
}