#ifndef ELLIPSE_MODELER
#define ELLIPSE_MODELER

#include "BaysFitFunctions.h"
#include "UBCUtil.h"
#include "LinearAlgebra.h"

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

static const int MAPEstimatiorStride = 100;
struct  MAPEstimatior {

	typedef DynamicAutoDiffCostFunction<MAPEstimatior
		,  MAPEstimatiorStride> MAPEstimatiorCostFunction;

	 MAPEstimatior(MatrixXd inX_i, MatrixXd inY_i
		, MatrixXd variance, double alpha) :
		_inX_i(inX_i) //(m.1) m is # of featues
		, _inY_i(inY_i)	 //(d.1) d is 1 
		, _m(inX_i.rows())
		, _d(inY_i.rows())
		, _alpha(alpha)
		, _K(inX_i * inX_i.transpose() * (alpha+1)) //(m.m)
		, _variance(buildTensorProduct(_K.inverse(), variance)) //(d.d)
		,	_mean( (inY_i*inX_i.transpose()) * 
			((inX_i*inX_i.transpose()).inverse() / (alpha+1)) ) {}

	template <typename T>
	bool operator()(T const* const* A, T* residuals) const {

		// new Variance = [ K^−1 ⊗ V ](md.md)
/*
		for (int d = 0; d < _d; d++) {
			for (int m = 0; m < _m; m++) {
//(A[d*_m + m][0] - _mean(d,m))/pow(_variance(d*_m+m,d*_m+m),0.5);
				ASxx[d*_m + m] = T(0.0);
				for (int c = 0; c < _m; c++) 
					ASxx[d*_m + m] += A[d*_m + c][0] * T(_Sxx(c,m));
			}	
		}
*/
/*


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
*/
		T retVal (0);
		residuals[0] = T(1.0/2.0) * ( retVal );

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static MAPEstimatiorCostFunction* Create( MatrixXd inX_i
																		, MatrixXd inY_i
																		, MatrixXd variance
																		, double alpha
																		, vector<double>* Amat
																		, vector<double*>* 
																			parameterBlocks) {

		MAPEstimatior* map =
			new MAPEstimatior(inX_i, inY_i, variance, alpha);
		MAPEstimatiorCostFunction* MAPcf = 
			new MAPEstimatiorCostFunction(map);
		
		// Add all the parameter blocks that affect this constraint.
    parameterBlocks->clear();
		int inM = inX_i.rows();
		int inD = inY_i.rows();
		for (int d = 0; d < inD; d++) {
			for (int m = 0; m < inM; m++) {
				//d.m
    		parameterBlocks->push_back(&((*Amat)[d*inM + m]));
    		MAPcf->AddParameterBlock(1);
			}
		}
    MAPcf->SetNumResiduals(1);
		return (MAPcf);
	}

	Eigen::MatrixXd _inX_i, _inY_i;
	int _m, _d;
	double _alpha;
	Eigen::MatrixXd _K, _variance, _mean; 

	};
#endif 
