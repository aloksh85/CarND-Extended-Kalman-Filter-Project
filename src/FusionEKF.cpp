#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

VectorXd compute_hx(VectorXd &x)
{
  VectorXd h_x(3);
  h_x<<0.0,0.0,0.0;

  h_x(0) = sqrt(pow(x(0),2) + pow(x(1),2));
  h_x(1) = atan2(x(1),x(0));

  if (h_x(0) < 0.1)
    h_x(2) = 0.0;
  else
    h_x(2) = ((x(0)*x(2)) + (x(1)*x(3)))/h_x(0);

  return h_x;
}


/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  init_counter_ = 0;


  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225,0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,0,0,
        0, 0.09, 0,
        0, 0, 0.09;

  H_laser_ << 1.0,0.0,0.0,0.0,
              0.0,1.0,0.0,0.0;

  Hj_ <<0,0,0,0,
  	    0,0,0,0,
	      0,0,0,0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

	ekf_.F_ = MatrixXd(4,4);
	ekf_.F_<< 1.0,0.0,1.0,0.0,
						0.0,1.0,0.0,1.0,
						0.0,0.0,1.0,0.0,
						0.0,0.0,0.0,1.0;

  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1.0,0.0,0.0,0.0,
             0.0,1.0,0.0,0.0,
             0.0,0.0,1000.0,0.0,
             0.0,0.0,0.0,1000.0;



  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ <<  1.0,0.0,1.0,0.0,
              0.0,1.0,0.0,1.0,
              1.0,0.0,1.0,0.0,
              0.0,1.0,0.0,1.0;

  ekf_.R_ = MatrixXd(2,2);
  ekf_.R_ << 1.0,0.0,
             0.0,1.0;

  ekf_.H_ = MatrixXd(2,4);
  ekf_.H_<<1.0,0.0,0.0,0.0,
           0.0,1.0,0.0,0.0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd radar_measure = measurement_pack.raw_measurements_;

      VectorXd x = ekf_.x_;
      ekf_.x_(0) = radar_measure(0)*cos(radar_measure(1));
      ekf_.x_(1) = radar_measure(0)*sin(radar_measure(1));

      if(init_counter_ > 0) {
        double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
        ekf_.x_(2) = (ekf_.x_(0) - x(0))/dt;
        ekf_.x_(3) = (ekf_.x_(1) - x(1))/dt;
    }
    else {
        ekf_.x_(2) = 0.0;
        ekf_.x_(3) = 0.0;
    }
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
      VectorXd lidar_measure  = measurement_pack.raw_measurements_;

      VectorXd x = ekf_.x_;
      ekf_.x_(0) = lidar_measure(0);
      ekf_.x_(1) = lidar_measure(1);

      if(init_counter_ > 0) {
        double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
        ekf_.x_(2) = (ekf_.x_(0) - x(0))/dt;
        ekf_.x_(3) = (ekf_.x_(1) - x(1))/dt;
    }
    else {
      ekf_.x_(2) = 0.0;
      ekf_.x_(3) = 0.0;
    }

    }

    // done initializing, no need to predict or update
    if(init_counter_ >= 0) {
       is_initialized_ = true;
      cout<<"initialized"<<endl;
    }
    else {
        init_counter_++;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  // convert to seconds
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  double noise_ax = 15.0;
  double noise_ay = 15.0;


  ekf_.Q_(0,0) = (dt_4*noise_ax)/4.0;
  ekf_.Q_(0,2) = (dt_3*noise_ax)/2.0;
  ekf_.Q_(1,1) = (dt_4*noise_ay)/4.0;
  ekf_.Q_(1,3) = (dt_3*noise_ay)/2.0;
  ekf_.Q_(2,0) = (dt_3*noise_ax)/2.0;
  ekf_.Q_(2,2) = (dt_2*noise_ax);
  ekf_.Q_(3,1) = (dt_3*noise_ay)/2.0;
  ekf_.Q_(3,3) = (dt_2*noise_ay);

	ekf_.F_(0,2) = dt;
	ekf_.F_(1,3) = dt;

  ekf_.Predict();
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**

   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    VectorXd h_x = compute_hx(ekf_.x_);
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.RadarUpdate(measurement_pack.raw_measurements_,
                     h_x,
                     Hj_,
                     R_radar_);
    }

    else {
    // Laser updates
    auto z = measurement_pack.raw_measurements_;
    ekf_.LidarUpdate(z,R_laser_,H_laser_);
  }

  previous_timestamp_ = measurement_pack.timestamp_;
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
