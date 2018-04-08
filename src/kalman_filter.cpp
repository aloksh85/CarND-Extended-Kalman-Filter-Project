#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  auto y = z - H_*x_;
  auto S = H_*P_*H_.transpose() + R_;
  auto K = P_*H_.transpose()*S.inverse();

  x_ += K*y;
  auto I =  MatrixXd::Identity(4,4);
  P_ =(I-K*H_)*P_;
}


void KalmanFilter::LidarUpdate(const VectorXd &z,
                              const MatrixXd &R,
                              const MatrixXd& H) {

  VectorXd y = z - H*x_;
  MatrixXd S = H*P_*H.transpose() + R;
  MatrixXd K = P_*H.transpose()*S.inverse();
  x_ += K*y;

  auto I =  MatrixXd::Identity(4,4);
  P_ = (I-K*H)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
    *
  */
  Predict();
  Update(z);
}
