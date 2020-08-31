#ifndef MEE_H
#define MEE_H

#include "BaysFitFunctions.h"
#include "UBCUtil.h"

#include <Eigen/Dense>

#include "ParseMathematica.h"
#include "ceres/ceres.h"
#include "ceres/dynamic_autodiff_cost_function.h"


using namespace Eigen;
using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::CauchyLoss;

//static const int _m = 1;

// The stride length of the dynamic_autodiff_cost_function evaluator.
static const int evStride = 100;

//template <typename T>
struct MEE {

	typedef DynamicAutoDiffCostFunction<MEE, evStride> MEECostFunction;

	MEE(MatrixXd inX, MatrixXd inY , MatrixXd variance
		, double alpha, MatrixXd Sxx, MatrixXd Syx, MatrixXd Syy) :
		_X(inX) //(m.1) m is # of featues
		, _y(inY)	 //(d.1) d is 1 
		, _m(_X.rows())
		, _d(_y.rows())
		, _variance(variance) //(d.d)
		, _alpha(alpha)
		, _Sxx(Sxx) //(m.m)
		, _Syx(Syx) //(d.m)
		, _Syy(Syy) //(d.d)
		, _meanA(_Syx * _Sxx.inverse()) {}

	template <typename T>
	bool operator()(T const* const* A, T* residuals) const {
/*
		T AminusMean[_d*_m];	//d.m
		T AminusMeanT[_m*_d];	//m.d
		T AminusMeanSxx[_d*_m];	//d.m
		for (int d = 0; d < _d; d++) {
			for (int m = 0; m < _m; m++) 
				AminusMean[d*_m + m] = A[d*_m + m][0] - T(_meanA(d,m));
		}
		for (int m = 0; m < _m; m++) {
			for (int d = 0; d < _d; d++) 
				AminusMeanT[m*_d + d] = AminusMean[d*_m + m];
		}
		//d.m x m.m = d.m
		for (int d = 0; d < _d; d++) {
			for (int m = 0; m < _m; m++) {
				AminusMeanSxx[d*_m + m] = T(0.0);
				for (int c = 0; c < _m; c++) 
					AminusMeanSxx[d*_m + m] += 
						AminusMean[d*_m + c] * T(_Sxx(c,m));
			}	
		}
		//d.m x m.d = d.d 
		T AminusMeanK_times_AminusMeanT[_d*_d]; //d.d
		for (int d = 0; d < _d; d++) {
			for (int m = 0; m < _d; m++) {
				AminusMeanK_times_AminusMeanT[d*_d + m] = T(0.0);
				for (int c = 0; c < _m; c++) 
					AminusMeanK_times_AminusMeanT[d*_d + m] += 
						AminusMeanSxx[d*_m + c] * AminusMeanT[c*_d+m];
			}	
		}
		MatrixXd Vinv = _variance.inverse(); //d.d
		T retMatA[_d*_d];
		T retMatY[_d*_d]; //−1/2 tr(V −1 S y|x )
		for (int d = 0; d < _d; d++) {
			for (int m = 0; m < _d; m++) {
				retMatA[d*_d + m]  = T(0.0);
				retMatY[d*_d + m]  = T(0.0);
				for (int c = 0; c < _d; c++) {
					retMatA[d*_d + m] += T(Vinv(m,c)) 
						* AminusMeanK_times_AminusMeanT[c*_d + m];
					retMatY[d*_d + m] += T(Vinv(m,c)) 
						* T(_Sy_x(c,m));
				}
			}			
		}

		T retValA (0);
		T retValY (0);
		for (int d = 0; d < _d; d++) {
			 retValA += T(retMatA[d*_d+d]);
			 retValY += T(retMatY[d*_d+d]);
		}
	residuals[0] = T(1.0/2.0) * T(retValA) 
			+ T(1.0/2.0) * T(retValY);

*/

    // Compute final prediction of Yi input vector
		T ASxx[_d*_m]; 
		//d.m x m.m = d.m
		for (int d = 0; d < _d; d++) {
			for (int m = 0; m < _m; m++) {
				ASxx[d*_m + m] = T(0.0);
				for (int c = 0; c < _m; c++) 
					ASxx[d*_m + m] += A[d*_m + c][0] * T(_Sxx(c,m));
			}	
		}

		T ASxxAtrans[_d*_d];		//d.m X m.d = (d.d) 
		T SyxAtrans[_d*_d];		//d.m X m.d = (d.d) 
		for (int d = 0; d < _d; d++) {
			for (int m = 0; m < _d; m++) {
				ASxxAtrans[d*_d + m] = T(0.0);
				for (int c = 0; c < _m; c++) {
					ASxxAtrans[d*_d + m] += 
						ASxx[d*_m + c] * A[m*_m+c][0];
					SyxAtrans[d*_d + m] += 
						 T(_Syx(d,c)) * A[m*_m+c][0];
				}
			}	
		}
		// d.d * d.d = d.d
		//Vinv * ( ASxxAtrans - T(2.0) * SyxAtrans + T(_Syy(0,0)) );
		T retMat[_d*_d];
		MatrixXd Vinv = _variance.inverse(); //d.d
		for (int d = 0; d < _d; d++) {
				retMat[d]  = T(0.0);
				for (int c = 0; c < _d; c++) {
					retMat[d] += T(Vinv(d,c))  
						* (ASxxAtrans[c*_d+d] - T(2.0) * SyxAtrans[c*_d+d] 
						+ T(_Syy(c,d)));
				}
		}			
		
		T retVal (0);
		for (int d = 0; d < _d; d++) {
			 retVal += T(retMat[d]);
		}

		residuals[0] = T(1.0/2.0) * ( retVal );

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static MEECostFunction* Create ( MatrixXd inX_i
																	, MatrixXd inY_i
																	, MatrixXd variance
																	, double alpha
																	, MatrixXd Sxx
																	, MatrixXd Syx
																	, MatrixXd Syy
																	, vector<double>* Amat
																	, vector<double*>* 
																		parameterBlocks) {

		MEE* mee = new MEE(inX_i, inY_i, variance
			, alpha, Sxx, Syx, Syy);
		MEECostFunction* cf = new MEECostFunction(mee);
		
		// Add all the parameter blocks that affect this constraint.
    parameterBlocks->clear();
		int inM = inX_i.rows();
		int inD = inY_i.rows();
		for (int d = 0; d < inD; d++) {
			for (int m = 0; m < inM; m++) {
				//d.m
    		parameterBlocks->push_back(&((*Amat)[d*inM + m]));
    		cf->AddParameterBlock(1);
			}
		}
    cf->SetNumResiduals(1);
		return (cf);
	}

	Eigen::MatrixXd _X, _y;
	int _m, _d;
	Eigen::MatrixXd _variance;
	double _alpha;
	Eigen::MatrixXd _Sxx, _Syx, _Syy, _meanA; //(d.m)

	};
#endif 
