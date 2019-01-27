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
  
  
  //accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    //difference between estimation and ground truth values at index i
    VectorXd diff = estimations[i] - ground_truth[i];
    //square of the difference element wise
    diff = diff.array() * diff.array();
    //add square of the differentce sum to the rmse accumulator
    rmse += diff;
  }
  
  //calculate the mean
  rmse = rmse/estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculates the jacobian matrix
   */
  MatrixXd Hj(3,4);

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check division by zero
  float pxy2 = px*px + py*py;

   // check division by zero - 0.0001
  if (fabs(pxy2) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }
  
  // compute the Jacobian matrix
  Hj << px/sqrt(pxy2), py/sqrt(pxy2), 0, 0, 
        -py/pxy2, px/pxy2, 0, 0, 
        py*(vx*py-vy*px)/pow(pxy2, 1.5),px*(vy*px-vx*py)/pow(pxy2, 1.5), px/sqrt(pxy2), py/sqrt(pxy2);
  
  return Hj;
}
