#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
	outfile_.open("data/obj_pose-laser-radar-ukf-output.txt");
	outfile2_.open("data/full-output.txt");
}

Tools::~Tools() {
	outfile_.close();
	outfile2_.close();
}

ofstream* Tools::getOStream(){
	return &outfile_;
}

ofstream* Tools::getOStream2(){
	return &outfile2_;
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse.fill(0);

  if(estimations.size()==0){
	 cerr<<"Nothing to estimate. Vector empty"<<endl;
	 return rmse;
  }

  if(estimations.size() != ground_truth.size()){
	  cerr << "The vector sizes differ."<<endl;
  }
  for(int i=0;i<estimations.size();++i){
	  for(int j=0;j<4;j++){
		  float d = estimations[i][j]-ground_truth[i][j];
		  rmse[j] += d*d;
	  }
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;

}
