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


Eigen::MatrixXd Tsigmoid(MatrixXd inZmat) {
		
	MatrixXd retZmat = 
		MatrixXd::Zero(inZmat.rows(),inZmat.cols());

	for (int i = 0; i < inZmat.rows(); i++) {
		for (int j = 0; j < inZmat.cols(); j++) 
			retZmat(i,j) = 1.0/(1.0 + exp(-inZmat(i,j)));
	}
	return retZmat;
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

double computeAlpha(Eigen::MatrixXd inXmat
	, Eigen::MatrixXd inYmat){
	double alpha = 0.0;
	
	//n -> sample size
	int n = inXmat.cols();
	//m -> x-features (d)
	int m = inXmat.rows();
	//d -> is one in our case
	int d = inYmat.rows();
	Eigen::MatrixXd variance = 
		varianceSolve (inYmat);
	//double variance = varianceSolve1D((inYmat).transpose());
	//cout << "variance: " << variance <<	endl;

	Eigen::MatrixXd buttom = ( inYmat*inXmat.transpose() 
		* (inXmat*inXmat.transpose()).inverse() 
		* inXmat*inYmat.transpose() );

	alpha = double(m*d) / ((variance.inverse() * buttom)
		.trace()- double(m*d));
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
	Eigen::MatrixXd variance = varianceSolve (inYmat);
	Eigen::MatrixXd Sy_x = computeSy_x(inXmat, inYmat);
	//Eqn. 32
	double Py_xvalpha = pow(alpha/(alpha+1.0),(m*d)/2) 
		* pow((2.0*M_PI*variance).determinant(),-n/2)
		* exp(-0.5*(variance.inverse()*Sy_x).trace());

	return Py_xvalpha;
}

vector<double> findAEv(Eigen::MatrixXd inXmat
	, Eigen::MatrixXd inYmat) {
	
	//d output size
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();
	//m -> (Feature size)
	int m = inXmat.rows();
	//d.d
	Eigen::MatrixXd variance = MatrixXd::Zero(d,d);
	variance = varianceSolve(inYmat);

//	double stddev = pow(variance(0,0), 0.5);	
	double alpha = computeAlpha(inXmat, inYmat);
	vector<double> retAvec(m*d,0.0);
	
	//(m.m)
	Eigen::MatrixXd Sxx = (alpha + 1.0) 
		* (inXmat * inXmat.transpose()); 
	//(d.m)
	Eigen::MatrixXd Syx = inYmat * inXmat.transpose();		
	//(d.d)
	Eigen::MatrixXd Syy = inYmat * inYmat.transpose();	
	//Eigen::MatrixXd Sy_x = computeSy_x(inXmat, inYmat);

	// the problem construction is a matter of creating a 
	// CostFunction for every observation.
	ceres::Problem problem;
//	for (int i = 0; i < n; ++i) {
		vector<double*> parameter_blocks;
		MEE::MEECostFunction* Ecf = MEE::Create(
				inXmat, inYmat, variance, alpha
				, Sxx, Syx, Syy, &retAvec, &parameter_blocks);

		problem.AddResidualBlock(Ecf, new CauchyLoss(0.5)
			, parameter_blocks);
//	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);

	std::cout<<"DONE"<<endl;
	//retVal.first = summary.final_cost;
	
	//std::cout << summary.FullReport() << "\n";

	return retAvec;
}

vector<double> findALi(Eigen::MatrixXd inXmat
	, Eigen::MatrixXd inYmat){

	//d -> is one in our case
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();
	//m -> polyOrder+1 (Feature size)
	int m = inXmat.rows();

	/* solve for std of y */
	//d.d
	Eigen::MatrixXd variance = 
		varianceSolve(inYmat);
	double alpha = computeAlpha(inXmat, inYmat);
	//(m.m)
	Eigen::MatrixXd Sxx = (alpha + 1.0) 
		* (inXmat * inXmat.transpose()); 
	
	Eigen::MatrixXd stddev = Eigen::MatrixXd::Zero(d,d);
	for (int j = 0; j < d; j++) 
		stddev(j,j) = pow(variance(j,j), 0.5);

	
	// the problem construction is a matter of creating a 
	// CostFunction for every observation.
	vector<double> retAvec(m*d,0.0);
	ceres::Problem problem;
	for (int i = 0; i < n; ++i) {
		vector<double*> parameter_blocks;
		
		double cVar = 1.0 - inXmat.col(i).transpose() 
			* (Sxx - inXmat.col(i) * inXmat.col(i).transpose()).inverse()
			* inXmat.col(i);
		MatrixXd stddev_i = stddev * (1.0/cVar);

		MLE::MLECostFunction* cf = MLE::Create(
			inXmat.col(i), inYmat.col(i), stddev, &retAvec
			, &parameter_blocks);
		problem.AddResidualBlock(cf, new CauchyLoss(0.5)
			, parameter_blocks);
	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);
	//retVal.first = summary.final_cost;
	
	//std::cout << summary.FullReport() << "\n";
	return retAvec;
}

vector<double> findANLi(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){
	
	int k = 10; 
	//d -> is one in our case
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();
	//m -> polyOrder+1 (Feature size)
	int m = inXmat.rows();
	

	//the problem construction is a matter of creating a 
	//CostFunction for every observation.
	ceres::Problem problem;
	vector <double> retwvec(k+k*m,0.0);
	for (int i = 0; i < n; ++i) {
		vector<double*> parameter_blocks;
		NMLE::NMLECostFunction* cf =
			NMLE::Create(
				inXmat.col(i), inYmat(0,i), k
				, &retwvec, &parameter_blocks);
		problem.AddResidualBlock(cf, new CauchyLoss(0.5)
		, parameter_blocks);
	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);
	//retVal.first = summary.final_cost;
	
	//std::cout << summary.FullReport() << "\n";

	return retwvec;
}

double computePy_xav(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat){
	
	//d -> is one in our case
	int d = inYmat.rows();
	//n -> sample size
	int n = inXmat.cols();
	//m -> feature size
	int m = inXmat.rows();

	Eigen::MatrixXd variance = varianceSolve (inYmat);
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

Eigen::MatrixXd predictY
	(Eigen::MatrixXd inXmat, vector<double> inA) {	

	int n = inXmat.cols();
	int m = inXmat.rows();
	int d = inA.size()/m;

	MatrixXd inAmat = MatrixXd::Zero(d,m);
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < m; j++)
			inAmat(i,j) = inA[i*m + j];
	}
	MatrixXd retYmat = MatrixXd::Zero(d,n);
	// d.n  = d.m    * m.n
	retYmat = inAmat * inXmat;

	return retYmat;
}

Eigen::MatrixXd predictYNN
	(Eigen::MatrixXd inXmat, vector<double> inA) {

	int n = inXmat.cols();
	int m = inXmat.rows();
	// inA.size() = (m+1)*k
	int k = inA.size() / (m+1);
	cout << "m: " << m << "k: " << k << endl;


	MatrixXd Wmat = MatrixXd::Zero(k,m);
	MatrixXd wVec = MatrixXd::Zero(k,1);
	for (int i=0; i < k; i++) {

		for (int j=0; j < m; j++)
			Wmat(i,j) = inA[j + i*m];

		wVec(i,0) = inA[i + k*m];
	}
	
	cout << "Wmat: " << Wmat << endl;
	cout << "wVec: " << wVec << endl;


	Eigen::MatrixXd retYmat = Eigen::MatrixXd::Zero(1,n);
	Eigen::MatrixXd Zmat = Eigen::MatrixXd::Zero(k,n);
	Eigen::MatrixXd hZmat = Eigen::MatrixXd::Zero(k,n);

	for (int i = 0; i < n; i++) 
		// k.1	 =	k.m * m.1
		Zmat.col(i) = Wmat * inXmat.col(i);

	// add non-linearity by transforming Zmat
	hZmat = Tsigmoid(Zmat);
	
	// make predictions
	for (int i = 0; i < n; i++) 
		// 1.1	 			 = 1.k 							* k.1
		retYmat.col(i) = wVec.transpose() * hZmat.col(i);

	return retYmat;
}






