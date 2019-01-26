#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /*
  * Calculates the RMSE between the estimations and the ground truth
  * */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  
  if (estimations.size()==0)
  {
      std::cout << "Vector size for estimations is zero"<< std::endl;
      return rmse;
  }
  
  if (estimations.size() != ground_truth.size())
  {
      std::cout << "Vector size for estimations and ground truth don't match"<< std::endl;
      return rmse;
  }
  
  
  // TODO: accumulate squared residuals
  for (int i=0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd diff = estimations[i] - ground_truth[i];
    
    diff = diff.array() * diff.array();
    rmse += diff;
  }
  
  // TODO: calculate the mean
  rmse = rmse/estimations.size();
  
  // TODO: calculate the squared root
  rmse = rmse.array().sqrt();
  
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // TODO: YOUR CODE HERE 

  // check division by zero
  float pxy2 = px*px + py*py;
  if (pxy2==0)
  {
      std::cout << "Divide by zero error" << std::endl;
      return Hj;
  }
  
  Hj << px/sqrt(pxy2), py/sqrt(pxy2), 0, 0, 
        -py/pxy2, px/pxy2, 0, 0, 
        py*(vx*py-vy*px)/pow(pxy2, 1.5),px*(vy*px-vx*py)/pow(pxy2, 1.5), px/sqrt(pxy2), py/sqrt(pxy2);
  // compute the Jacobian matrix

  return Hj;
}
