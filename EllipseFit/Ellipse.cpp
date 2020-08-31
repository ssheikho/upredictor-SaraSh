#include "Ellipse.h"
#include "EllipseConicConstraints.h"
#include "EllipseImplicitFit.h"

Ellipse::Ellipse(double a, double b, double alpha, double cX, double cY) :
	_a(a), _b(b), _rt2D(alpha, cX, cY) {}

Ellipse::Ellipse() : _rt2D(0.0, 0.0, 0.0) {}

Ellipse::Ellipse(MatrixXd inPtsAlongCols) : _rt2D(0.0, 0.0, 0.0) {
	double* coniVectArr =  solveConicVect(inPtsAlongCols);
	MatrixXd coniVect = MatrixXd::Zero(6,1);
	for (int i = 0; i < 6; i++)
		coniVect (i,0) = coniVectArr[i];

	MatrixXd conicImplicit = 
		arrangeConicVectIntoImplicit(coniVect);
	MatrixXd a33 = 
		aThirtyThreeFromConicSolution(coniVect);
	
	std::pair<MatrixXd, MatrixXd> solution = solveEigensystem(a33);
	double detImp = conicImplicit.determinant();
	double detA33 = a33.determinant();
	double k = -detImp / detA33;

	_a = sqrt(k / solution.first(1,0));
	_b = sqrt(k / solution.first(0,0));

	double ccBottomTerm = 4.0 *	coniVect(0,0) * coniVect(2,0)
		- coniVect(1,0) * coniVect(1,0);

	double ccXTerm = (coniVect(1,0) * coniVect(4,0)
		- 2.0 * coniVect(2,0) * coniVect(3,0)) / ccBottomTerm;
	double ccYTerm = (coniVect(3,0) * coniVect(1,0)
		- 2.0 * coniVect(0,0) * coniVect(4,0)) / ccBottomTerm;

	//tan(2θ)=B/A−C
	double alphaAngle = 
		atan2(coniVect(1,0),(coniVect(2,0)-coniVect(1,0)))/2.0;

//	double alphaAngle = 
//		atan2(solution.second(1,0), solution.second(0,0));

	_rt2D.setRT(alphaAngle, ccXTerm, ccYTerm);

}

Ellipse::Ellipse(EllipseImplicitFit eif) : _rt2D(0.0, 0.0, 0.0) {
	MatrixXd a33 = 
		aThirtyThreeFromConicSolution(eif.getConicVect());
	std::pair<MatrixXd, MatrixXd> solution = solveEigensystem(a33);
	double detImp = eif.getConicImplicit().determinant();
	double detA33 = a33.determinant();
	double k = -detImp / detA33;

	_a = sqrt(k / solution.first(1,0));
	_b = sqrt(k / solution.first(0,0));

	double ccBottomTerm = 4.0 *
		eif.getConicVect()(0,0) * eif.getConicVect()(2,0)
		- eif.getConicVect()(1,0) * eif.getConicVect()(1,0);
	double ccXTerm = (eif.getConicVect()(1,0) * eif.getConicVect()(4,0)
		- 2.0 * eif.getConicVect()(2,0) * eif.getConicVect()(3,0)) /
		ccBottomTerm;
	double ccYTerm = (eif.getConicVect()(3,0) * eif.getConicVect()(1,0)
		- 2.0 * eif.getConicVect()(0,0) * eif.getConicVect()(4,0)) /
		ccBottomTerm;
	
	double alphaAngle = atan2(eif.getConicVect()(1,0)
		, (eif.getConicVect()(2,0)-eif.getConicVect()(0,0)))/2.0;

/* same sln as above

	double alphaAngle = 
		atan2(solution.second(1,0), solution.second(0,0));

*/

	_rt2D.setRT(alphaAngle, ccXTerm, ccYTerm);
}

double* Ellipse::solveConicVect(MatrixXd inPtsAlongCols){
	// 3xn
	int nPts = inPtsAlongCols.cols();
	double* parameter_blocks = new double[6];
	for (int i = 0; i < 6; i++)
		parameter_blocks[i] = 1.0;

	ceres::Problem problem;

	// Configure the loss function.
  LossFunction* loss = new CauchyLoss(0.5);
  
	for (int i = 0; i < nPts; ++i) {
		
		ceres::CostFunction *cf =
			EllipseConicConstraints::Create(
				inPtsAlongCols(0,i),inPtsAlongCols(1,i));
		problem.AddResidualBlock(cf,
			new ceres::CauchyLoss(0.5), parameter_blocks);
	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	ceres::Solve(options, &problem, &summary);
	return parameter_blocks;
}

Ellipse::~Ellipse() {}

void Ellipse::setAB(double a, double b) {
	_a = a;
	_b = b;
}

double Ellipse::getA() {
	return _a;
}

double Ellipse::getB() {
	return _b;
}

RigidTrans2D Ellipse::getRT2D() {
	return _rt2D;
}

Matrix<double,Dynamic,Dynamic> Ellipse::getPointAtThetaH(double theta) {
	// Points on the 2D ellipse
	Matrix<double,Dynamic,Dynamic> pointOnEllipsePseudo2D(3,1);
	pointOnEllipsePseudo2D(0,0) = _a * cos(theta);
	pointOnEllipsePseudo2D(1,0) = _b * sin(theta);
	pointOnEllipsePseudo2D(2,0) = 1.0;

	//return with reference to the global/camera coords
	return _rt2D.getRT() * pointOnEllipsePseudo2D;

	//return pointOnEllipsePseudo2D;
	
}

double Ellipse::getPointXAtThetaH(double theta) {
	return getPointAtThetaH(theta)(0,0);
}

double Ellipse::getPointYAtThetaH(double theta) {
	return getPointAtThetaH(theta)(1,0);
}

double Ellipse::getPointWAtThetaH(double theta) {
	return getPointAtThetaH(theta)(2,0);
}


Matrix<double,Dynamic,Dynamic> Ellipse::getPointAtThetasH(
		Matrix<double, Dynamic, Dynamic> thetas) {
	size_t nTheta = thetas.cols();

	Matrix<double,Dynamic,Dynamic> retVals(3, nTheta);

	for(size_t i = 0; i < nTheta; i++)
		retVals.block(0,i,3,1) = getPointAtThetaH(thetas(0,i));

	return retVals;
}


