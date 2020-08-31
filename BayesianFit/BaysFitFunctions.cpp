#include "BaysFitFunctions.h"
#include "ParseCSV.h"
#include "UBCUtil.h"
#include "EllipseImplicitFit.h"
#include "Ellipse.h"
#include "Ellipse3D.h"
#include <Eigen/Dense>

#include "ParseMathematica.h"

using namespace Eigen;

Eigen::MatrixXd changeBasis(double inX, int k){
	Eigen::MatrixXd retMat(k+1,1);

	for(size_t i = 0; i < k+1; i++) 
		retMat(i,0) = pow(inX,i);
	
	return retMat;
}

Eigen::MatrixXd getPtsXYAlongCols(MatrixXd inPtsAlongColsHom) {

	//Eigen::MatrixXd PtsXYAlongCols = EllipseImplicitFit::getPtsXYAlongCols();



	EllipseImplicitFit eif(inPtsAlongColsHom);
	//Eigen::MatrixXd PtsXYAlongCols = EllipseImplicitFit::getPtsXYAlongCols();
	//return PtsXYAlongCols;
}

// The scenario is that we are given a data set of exchangeable pairs 
// D = {(y 1 , x 1 ), ..., (y N , x N )}.
// Collect Y = [y 1 · · · y N ] and X = [x 1 · · · x N ].
Eigen::MatrixXd buildXMatIn(Eigen::MatrixXd inXvec, int k){
	
	int n = inXvec.rows();

	Eigen::MatrixXd XMat(k+1, n);

	//cout<< "start" << endl;
	for(size_t col = 0; col < n; col++) {
			XMat.block(0,col,k+1,1) = changeBasis(inXvec(col,0), k);
			//cout<< XMat.block(0,col,k+1,1) << endl;
	}
	return XMat;
}

double computeAlpha(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){
	double alpha = 0.0;
	
	//n -> sample size
	int n = inXmat.cols();
	//m -> x-features (d)
	int m = inXmat.rows();
	//d -> is one in our case
	int d = inYmat.rows();
	Eigen::MatrixXd variance = varianceSolve ((inYmat).transpose());
	//double variance = varianceSolve1D((inYmat).transpose());
	//cout << "variance: " << variance <<	endl;

	Eigen::MatrixXd buttom = ( inYmat*inXmat.transpose() 
		* (inXmat*inXmat.transpose()).inverse() 
		* inXmat*inYmat.transpose() );

	alpha = double(m*d) / ((variance.inverse() * buttom).trace()- double(m*d));
	return alpha;
}

Eigen::MatrixXd computeSy_x(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){
	
	//n -> sample size
	int n = inXmat.cols();
	double alpha = computeAlpha(inXmat, inYmat);

	Eigen::MatrixXd right = ( inYmat*inXmat.transpose() 
		* (inXmat*inXmat.transpose()).inverse() 
		* inXmat*inYmat.transpose() );

	Eigen::MatrixXd Sy_x = inYmat*inYmat.transpose() - (1.0/(alpha+1.0))*right;
	return Sy_x;
}

double computePy_xvalpha(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){
	
	//m -> x-features (d)
	int m = inXmat.rows();
	//d -> is one in our case
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();

	double alpha = computeAlpha(inXmat, inYmat);
	//cout << "alpha" << alpha << endl;
	
	//double v = varianceSolve1D((inYmat).transpose()); 
	Eigen::MatrixXd variance = varianceSolve ((inYmat).transpose());
	Eigen::MatrixXd Sy_x = computeSy_x(inXmat, inYmat);

	double Py_xvalpha = pow(alpha/(alpha+1.0),(m*d)/2) 
		* pow((2.0*M_PI*variance).determinant(),-n/2)
		* exp(-0.5*(variance.inverse()*Sy_x).trace());

	return Py_xvalpha;
}

vector<double> findAEv(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){
	
	//d -> is one in our case
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();
	//m -> polyOrder+1 (Feature size)
	int m = inXmat.rows();
	
	Eigen::MatrixXd variance = varianceSolve((inYmat).transpose());
	double stddev = pow(variance(0,0), 0.5);	
	double alpha = computeAlpha(inXmat, inYmat);
	vector<double> retAvec(m,0.0);
	
	
	ceres::Problem problem;
	//for (int i = 0; i < n; ++i) {
		vector<double*> parameter_blocks;
		MEE::MEECostFunction* Ecf =
			MEE::Create(
				inXmat, inYmat, stddev, alpha, &retAvec, &parameter_blocks);
		problem.AddResidualBlock(Ecf, new CauchyLoss(0.5), parameter_blocks);
	//}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);

	std::cout<<"DONE"<<endl;
	//retVal.first = summary.final_cost;
	
	//std::cout << summary.FullReport() << "\n";

	return retAvec;
}

vector<double> findALi(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){

	//d -> is one in our case
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();
	//m -> polyOrder+1 (Feature size)
	int m = inXmat.rows();
	
	Eigen::MatrixXd variance = varianceSolve((inYmat).transpose());
	double stddev = pow(variance(0,0), 0.5);
	vector<double> retAvec(m,0.0);
	
	//the problem construction is a matter of creating a CostFunction
	//for every observation.
	ceres::Problem problem;
	for (int i = 0; i < n; ++i) {
		vector<double*> parameter_blocks;
		MLE::MLECostFunction* cf =
			MLE::Create(
				inXmat.col(i), inYmat(0,i),
				stddev, &retAvec, &parameter_blocks);
		problem.AddResidualBlock(cf, new CauchyLoss(0.5), parameter_blocks);
	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);
	//retVal.first = summary.final_cost;
	
	//std::cout << summary.FullReport() << "\n";
	return retAvec;
}

double computePy_xav(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){
	
	//d -> is one in our case
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();
	//m -> feature size
	int m = inXmat.rows();

	Eigen::MatrixXd variance = varianceSolve ((inYmat).transpose());
	Eigen::MatrixXd A = (MatrixXd::Zero(d,m));

	vector<double> Avec = findALi(inXmat, inYmat);
	for (int i=0; i < m; i++) A(0,i) = Avec[i];
	
	Eigen::MatrixXd tr = 
		variance.inverse()
		*( A*inXmat*inXmat.transpose()*A.transpose() 
			- 2.0*inYmat*inXmat.transpose()*A.transpose() 
			+ inYmat*inYmat.transpose() );

	double Py_xva = pow((2.0*M_PI*variance).determinant(),-n/2)
		* exp(-0.5*(tr).trace());

	return Py_xva;
}

Eigen::MatrixXd predictY(Eigen::MatrixXd inXmat, Eigen::MatrixXd inA){	
	int n = inXmat.cols();
	int m = inXmat.rows();
	
	Eigen::MatrixXd retYmat = Eigen::MatrixXd::Zero(1,n);
	for (int i = 0; i < n; ++i) {
		retYmat.col(i) = inA * inXmat.col(i);
	}
	return retYmat;
}








