#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0,0,0,0;
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // Laser measurement mapping matrix
  MatrixXd H_in = MatrixXd(2,4);
  H_in << 1, 0, 0, 0,
          0, 1, 0, 0;

  // state vector
  VectorXd x_in = VectorXd(4);
  x_in << 0, 0, 0, 0;

  // state covariance matrix
  MatrixXd P_in = MatrixXd(4,4);
  P_in << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;

  // state transistion matrix
  MatrixXd F_in = MatrixXd(4,4);
  F_in << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

  // process covariance matrix
  MatrixXd Q_in = MatrixXd(4,4);
  Q_in << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;

  //Initializes Kalman filter
  ekf_.Init(x_in, P_in, F_in, H_in, R_laser_, R_radar_, Q_in);


  // Radar measurement covariance matrix
  ekf_.R_radar_ = MatrixXd(3,3);
  ekf_.R_radar_ <<  0.09, 0, 0,
                    0, 0.0009, 0,
                    0, 0, 0.09;

  // Laser measurement covariance matrix
  ekf_.R_laser_ = MatrixXd(2,2);
  ekf_.R_laser_ <<  0.0225, 0,
                    0, 0.0225;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF()
{

}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  // printf("Measurement type %s\n", (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)?"radar":"laser");
  if (!is_initialized_)
  {
    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rate = measurement_pack.raw_measurements_[2];

      // Convert polar measuremnts to Cartesian coordinates
      float px = ro * cos(phi);
      float py = ro * sin(phi);
      float vx = rate * cos(phi);
      float vy = rate * sin(phi);

      ekf_.x_ << px, py, vx, vy;
    }
    else
    {
      if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
      {
        float px = measurement_pack.raw_measurements_[0];
        float py = measurement_pack.raw_measurements_[1];
        float vx = 0.0; //Not moving
        float vy = 0.0; //Not moving

        ekf_.x_ << px, py, vx, vy;
      }
    }

    // Initialize the covariance matrix
    float pVar = 100;
    float vVar = 10;
    // float pVar = 1000;
    // float vVar = 100;

    ekf_.P_ << pVar, 0, 0, 0,
             0, pVar, 0, 0,
             0, 0, vVar, 0,
             0, 0, 0, vVar;

    previous_timestamp_ = measurement_pack.timestamp_;

    //Initialized!
    is_initialized_ = true;
  }
  else
  {
    // Compute elapsed time from previous state estimate.
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / (1000*1000.0);
    previous_timestamp_ = measurement_pack.timestamp_;

    //Update transistion matrix
    // cout << "dt: " << dt << "\n";
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;

    float dt2 = dt * dt;
    float dt3 = dt2 * dt;
    float dt4 = dt3 * dt;

    float qx = ekf_.noise_ax;
    float qy = ekf_.noise_ay;

    //Update process covariance matrix
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt4/4*qx,  0,          dt3/2*qx,  0,
                0,          dt4/4*qy,  0,          dt3/2*qy,
                dt3/2*qx,  0,          dt2*qx,    0,
                0,          dt3/2*qy,  0,           dt2*qy;


    //Predict!
    ekf_.Predict();

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      //Validate measurement_pack
      if (measurement_pack.raw_measurements_[0] != 0.0 && measurement_pack.raw_measurements_[1] != 0.0 && measurement_pack.raw_measurements_[2] != 0.0)
      {
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      }
    }
    else
    {
      //Validate measurement_pack
      if (measurement_pack.raw_measurements_[0] != 0.0 && measurement_pack.raw_measurements_[1] != 0.0)
      {
        ekf_.Update(measurement_pack.raw_measurements_);
      }
    }
  }
}
