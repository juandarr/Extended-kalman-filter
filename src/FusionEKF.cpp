#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  

  // initializing matrices
  x = VectorXd(4);
  F = MatrixXd(4,4);
  P = MatrixXd(4,4);
  I = MatrixXd::Identity(2, 2);
  Q = MatrixXd(4,4);

  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  

  //object state
  x << 0 , 0 , 0, 0;

  //State transition matrix
  F << 1 , 0, 1, 0,
       0, 1, 0, 1,
      0, 0, 1 , 0,
      0 , 0, 0, 1;

  //state covariance matrix
  P << 1, 0, 0, 0,
      0 , 1 ,0, 0,
      0, 0 , 1000, 0,
      0, 0, 0, 1000;

  /**
   * Set the process and measurement noises
  */

  // process covariance matrix
  Q << 0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0, 
       0, 0, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix - laser
  H_laser_ << 1 , 0, 0, 0,
              0 , 1 , 0 , 0;

  //Initialize enhanced kalman filter with the matrices used for lidar data
  ekf_.Init(x, P, F, H_laser_, R_laser_, Q);

  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * - Initialize the state ekf_.x_ with the first measurement.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      //Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]),
                measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]),
                0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      //Initialize state.
      // set the state with the initial location and zero velocity
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1], 
              0,
              0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /**
   * Prediction
   */

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // 1. Modify the F matrix so that the time is integrated
  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   */
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // 2. Set the process covariance matrix Q
  /** Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << (pow(dt,4)/4)*noise_ax, 0, (pow(dt,3)/2)*noise_ax, 0 ,
            0, (pow(dt,4)/4)*noise_ay, 0, (pow(dt,3)/2)*noise_ay,  
            (pow(dt,3)/2)*noise_ax, 0, pow(dt,2)*noise_ax, 0 ,
            0, (pow(dt,3)/2)*noise_ay, 0, pow(dt,2)*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
     * Radar updates
     */ 
    //Define the measurement matrix Hj with the jacobian for radar data
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    //Define the measurement covariance matrix for radar data
    ekf_.R_ = R_radar_;
    //Call the enhanced kalman filter function to update the state given new radar measurements
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    /** 
     * Laser updates
     */
    //Define the measurement matrix H for laser data
    ekf_.H_ = H_laser_;
    //Define the measurement covariance matrix for laser data
    ekf_.R_ = R_laser_;
    //Call the kalman filter function to update the state given new laser measurements
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "Type of sensor: " << measurement_pack.sensor_type_ << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
