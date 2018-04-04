#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size

	if(estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];
		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();
	//calculate the squared root
	rmse = rmse.array().sqrt();
	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3,4);

	Hj <<0,0,0,0,
	     0,0,0,0,
	     0,0,0,0;

	//recover state parameters

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	 if (fabs(px) < 0.001 && fabs(py) < 0.001)
	 {
	    std:: cerr<<"x and y of state are zero. Cannot compute Jacobian"<<std::endl;
	    return Hj;
	 }
	//compute the Jacobian matrix
    float px_2 = px*px;
    float py_2 = py*py;
    float denom1 =  pow(px_2+py_2,1.5);
    float denom2 = pow(px_2 + py_2, 0.5);

    Hj(0,0) = px/denom2;
    Hj(0,1) = py/denom2;
    Hj(1,0) = -py/(px_2+py_2);
    Hj(1,1) = px/(px_2+py_2);
    Hj(2,0) = (py*((vx*py) - (px*vy)))/denom1;\
    Hj(2,1) = (px*((vy*px)-(vx*py)))/denom1;
    Hj(2,2) =  px/denom2;
    Hj(2,3) =  py/denom2;

	return Hj;
}
