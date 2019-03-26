#define _USE_MATH_DEFINES

#include "kalman_filter.h"
#include "math.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

const double epsilon = 0.00001;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  x_ = F_ * x_; // + u
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  MatrixXd I; // Identity matrix
  VectorXd y;
  MatrixXd S;
  MatrixXd K;

  I = MatrixXd::Identity(P_.rows(), P_.cols());

  y = z - H_ * x_;
  S = H_ * P_ * H_.transpose() + R_;
  K = P_ * H_.transpose() * S.inverse();
  
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  MatrixXd I; // Identity matrix
  VectorXd y;
  MatrixXd S;
  MatrixXd K;

  I = MatrixXd::Identity(P_.rows(), P_.cols());

  y = z - RadarMeasurement(x_);
  S = H_ * P_ * H_.transpose() + R_;
  K = P_ * H_.transpose() * S.inverse();
  
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

VectorXd KalmanFilter::RadarMeasurement(const VectorXd &x_state) {
   double px = x_state(0), py = x_state(1), vx = x_state(2), vy = x_state(3);

   double rho = sqrt(px*px + py*py);
   double phi;
   if (fabs(px) <= epsilon) {
      if (py < 0) {
         phi = - M_PI / 2;
      } else {
         phi = M_PI / 2;
      }
   } else {
     phi = atan2(py, px);
     while (phi > M_PI) {
        phi -= 2 * M_PI;
     }
     while (phi < -M_PI) {
        phi += 2 * M_PI;
     }
   }
   double rhodot = (px*vx + py*vy) / (rho + epsilon);

   VectorXd meas = VectorXd(3);
   meas << rho, phi, rhodot;
   return meas;
}