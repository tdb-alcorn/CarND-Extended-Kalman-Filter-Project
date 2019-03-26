#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residuals = estimations[i] - ground_truth[i];
    VectorXd r;

    r = residuals.array()*residuals.array();
    rmse += r;
  }

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);

  double px = x_state(0), py = x_state(1), vx = x_state(2), vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  double rho2 = px*px+py*py;
  double rho = sqrt(rho2);
  double rho3 = rho*rho2;

  // check division by zero
  if (fabs(rho2) < 0.00001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px/rho), (py/rho), 0, 0,
      -(py/rho2), (px/rho2), 0, 0,
      py*(vx*py - vy*px)/rho3, px*(px*vy - py*vx)/rho3, px/rho, py/rho;

  return Hj;
}
