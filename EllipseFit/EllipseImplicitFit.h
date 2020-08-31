#ifndef ELLIPSE_IMPLICIT_FIT_H
#define ELLIPSE_IMPLICIT_FIT_H

#include "FitFunctions.h"
//#include "BaysFitFunctions.h"
#include "Plane.h"
#include "UBCUtil.h"

#include <cmath>

#include <Eigen/Dense>

using namespace std;

//template<typename T>
class EllipseImplicitFit {
public:

	EllipseImplicitFit(MatrixXd inPtsAlongColsHom);

	~EllipseImplicitFit();


	Plane getPlane();

	MatrixXd getRotToXY();

	MatrixXd getPtsXYAlongCols();

	MatrixXd getConicConstraintMat();

	MatrixXd getConicVect();

	MatrixXd getConicImplicit();

protected:
	MatrixXd _inPtsAlongCols;
	Plane _plane;
	double _dotProduct;
	MatrixXd _crossProduct;
	MatrixXd _rotToXY
			, _ptsXYAlongCols
			, _conicConstraintMat
			, _conicVect
			, _conicImplicit;

};

#endif
