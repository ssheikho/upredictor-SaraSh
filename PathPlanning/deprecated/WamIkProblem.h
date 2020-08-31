/* 
IK CHAIN + CAMERA T
	Given inPts tracked camera markers, 
	the goal is to find cameraT (3 orient. 3 trans)
	that	minimize the re-Transformation error. 
*/
#ifndef WAM_IK_PROBLEM_H
#define WAM_IK_PROBLEM_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include "SpatialJacobian.h"
#include <UBCUtil.h>

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

#include <Eigen/SVD>
#include <string>
#include <fstream>
#include <cmath>        // std::abs
#include <math.h>       

#include <iostream>
#include <vector>

using namespace Eigen;
using namespace ceres;
using namespace std;

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::CauchyLoss;


template<typename T>
class WamIkProblem {
public:
	WamIkProblem( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
				inPtsAlongRows		
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles

		) :
		_fk(fk)
		, _JsWam(JsWam)
		, _inPtsAlongRows(inPtsAlongRows)
		, _nPts(_inPtsAlongRows.rows())
		, _inPtsRPi(_inPtsAlongRows.rightCols(3).transpose())
		, _inPtsRTh(_inPtsAlongRows.block(0,3,_nPts,3).transpose())
		, _inPtsRW(_inPtsAlongRows.leftCols(3).transpose())
		, _inPtsRLA(_inPtsAlongRows.block(0,6,_nPts,3).transpose())
		, _inPtsREl(_inPtsAlongRows.block(0,9,_nPts,3).transpose())
		, _inPtsRUA(_inPtsAlongRows.block(0,12,_nPts,3).transpose())
		, _inPtsRSh(_inPtsAlongRows.block(0,15,_nPts,3).transpose())
		,	_inPtsRCh(_inPtsAlongRows.block(0,18,_nPts,3).transpose())
		,	_inPtsMCh(_inPtsAlongRows.block(0,21,_nPts,3).transpose())
		,	_inPtsLCh(_inPtsAlongRows.block(0,24,_nPts,3).transpose())
		,	_inPtsLSh(_inPtsAlongRows.block(0,27,_nPts,3).transpose())
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		, _startConfigJoints(new T[getPtsBlockSize()])
		, _startConfigCamera(new T[getCameraBlockSize()])
		// inVectors WAM
		, _shoToUaVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToElVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _laToWrVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _elToWrVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToLaVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		,_shoToWrVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToThVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToPiVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _wrToThVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _wrToPiVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _wrToEeVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToEeVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		,_elOffsetVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		, _wOffsetVsWam(Eigen::Matrix<T,3
			, Eigen::Dynamic>::Zero(3,_nPts))
		// inPts WAM
		, _pShoWam(Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _pUaWam(Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _pElWam(Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _pLaWam(Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _pWrWam(Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _pThWam(Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _pPiWam(Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _EeOrientation(
			Eigen::Matrix<T, Eigen::Dynamic,3>::Zero(3*_nPts,3))
		
		, _inPtsRWrInChest(Eigen::Matrix<
			T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _inPtsRThInChest(Eigen::Matrix<
			T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _inPtsRPiInChest(Eigen::Matrix<
			T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _inPtsRLaInChest(Eigen::Matrix<
			T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _inPtsRElInChest(Eigen::Matrix<
			T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _inPtsRUaInChest(Eigen::Matrix<
			T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _inPtsRShInChest(Eigen::Matrix<
			T, 3, Eigen::Dynamic>::Zero(3,_nPts))	
		, _shoToUaVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToElVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _laToWrVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _elToWrVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToLaVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToWrVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToThVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToPiVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _wrToThVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _wrToPiVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _wrToEeVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _shoToEeVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _elOffsetVsWamInChest(
				Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3,_nPts))
		, _wOffsetVsWamInChest(
				Eigen::Matrix<T,3,Eigen::Dynamic>::Zero(3,_nPts))
		, _EeOrientationInChest(Eigen::Matrix<T, Eigen::Dynamic,3>
			::Zero(3*_nPts,3))	
{
	//initialize cerese parameter blocks (i.e. 
	//PtsBlockSize & CameraBlockSize
	for (int j = 0; j < getPtsBlockSize(); j++) 
		_startConfigJoints[j] = T(0.0);
	for (int j = 0; j < getCameraBlockSize(); j++) 
		_startConfigCamera[j] = T(0.0);

	mapInPtsHtoVecsR();
	mapInPtsHtoInPtsR();
	rotInPtsToChestPlane();
	//** For print in Mathematica **//
	/** inPts **/
	//(a) for human
	printEigenMathematica( _inPtsRPi.transpose(), cout
		, "inPtsRPi");	
	printEigenMathematica( _inPtsRTh.transpose(), cout
		, "inPtsRTh");	
	printEigenMathematica( _inPtsRW.transpose(), cout
		, "inPtsRW");	
	printEigenMathematica( _inPtsRLA.transpose(), cout
		, "inPtsRLA");	
	printEigenMathematica( _inPtsREl.transpose(), cout
		, "inPtsREl");	
	printEigenMathematica( _inPtsRUA.transpose(), cout
		, "inPtsRUA");	
	printEigenMathematica( _inPtsRSh.transpose(), cout
		, "inPtsRSh");	
	printEigenMathematica( _inPtsRCh.transpose(), cout
		, "inPtsRCh");	
	printEigenMathematica( _inPtsMCh.transpose(), cout
		, "inPtsMCh");	
	printEigenMathematica( _inPtsLCh.transpose(), cout
		, "inPtsLCh");	
	printEigenMathematica( _inPtsLSh.transpose(), cout
		, "inPtsLSh");	
	//(b) for wam, i.e. scaled to wam
	printEigenMathematica( _pPiWam.transpose(), cout
		, "inPtsRPiWam");	
	printEigenMathematica( _pThWam.transpose(), cout
		, "inPtsRThWam");	
	printEigenMathematica( _pWrWam.transpose(), cout
		, "inPtsRWWam");	
	printEigenMathematica( _pLaWam.transpose(), cout
		, "inPtsRLAWam");	
	printEigenMathematica( _pElWam.transpose(), cout
		, "inPtsRElWam");	
	printEigenMathematica( _pUaWam.transpose(), cout
		, "inPtsRUAWam");	
	printEigenMathematica( _pShoWam.transpose(), cout
		, "inPtsRShWam");	
	/** (b.2) inVectors **/
	printEigenMathematica( _shoToUaVsWam.transpose()
		, cout, "shoToUaVsWam");	
	printEigenMathematica( _shoToElVsWam.transpose()
		, cout, "shoToElVsWam");	
	printEigenMathematica( _laToWrVsWam.transpose()
		, cout, "laToWrVsWam");	
	printEigenMathematica( _elToWrVsWam.transpose()
		, cout, "elToWrVsWam");	
	printEigenMathematica( _shoToLaVsWam.transpose()
		, cout, "shoToLaVsWam");	
	printEigenMathematica( _shoToWrVsWam.transpose()
		, cout, "shoToWrVsWam");	
	printEigenMathematica( _shoToThVsWam.transpose()
		, cout, "shoToThVsWam");	
	printEigenMathematica( _shoToPiVsWam.transpose()
		, cout, "shoToPiVsWam");	
	printEigenMathematica(_shoToEeVsWam.transpose()
		, cout, "shoToEeVsWam");	
	printEigenMathematica( _elOffsetVsWam.transpose()
		, cout, "elOffsetVsWam");	
	printEigenMathematica( _wOffsetVsWam.transpose()
		, cout, "wOffsetVsWam");	
	//end-effector orientation 
	printEigenMathematica( _EeOrientation
		, cout, "EeOrientation");	


	// (b.1) inPts ROBOT
	printEigenMathematica( _inPtsRPiInChest.transpose(), cout
		, "inPtsRPiInChest");	
	printEigenMathematica( _inPtsRThInChest.transpose(), cout
		, "inPtsRThInChest");	
	printEigenMathematica( _inPtsRWrInChest.transpose(), cout
		, "inPtsRWrInChest");	
	printEigenMathematica( _inPtsRLaInChest.transpose(), cout
		, "inPtsRLaInChest");	
	printEigenMathematica( _inPtsRElInChest.transpose(), cout
		, "inPtsRElInChest");	
	printEigenMathematica( _inPtsRUaInChest.transpose(), cout
		, "inPtsRUaInChest");	
	printEigenMathematica( _inPtsRShInChest.transpose(), cout
		, "inPtsRShInChest");	

	//InVECTORS Robot anchored to chest plane
	printEigenMathematica( _shoToUaVsWamInChest.transpose()
		, cout, "shoToUaVsWamInChest");	
	printEigenMathematica( _shoToElVsWamInChest.transpose()
		, cout, "shoToElVsWamInChest");	
	printEigenMathematica( _laToWrVsWamInChest.transpose()
		, cout, "laToWrVsWamInChest");		
	printEigenMathematica( _elToWrVsWamInChest.transpose()
		, cout, "elToWrVsWamInChest");	
	printEigenMathematica( _shoToLaVsWamInChest.transpose()
		, cout, "shoToLaVsWamInChest");	
	printEigenMathematica( _shoToWrVsWamInChest.transpose()
		, cout, "shoToWrVsWamInChest");
	printEigenMathematica( _shoToThVsWamInChest.transpose()
		, cout, "shoToThVsWamInChest");
	printEigenMathematica( _shoToPiVsWamInChest.transpose()
		, cout, "shoToPiVsWamInChest");
	printEigenMathematica( _shoToEeVsWamInChest.transpose()
		, cout, "shoToEeVsWamInChest");
	//(b.2) inVectors OFFSETS WAM in chest plane 	
	printEigenMathematica( _elOffsetVsWamInChest.transpose()
		, cout, "elOffsetVsWamInChest");	
	printEigenMathematica( _wOffsetVsWamInChest.transpose()
		, cout, "wOffsetVsWamInChest");	
	//end-effector orientation 
	printEigenMathematica( _EeOrientationInChest
		, cout, "EeOrientationInChest");	


}

~WamIkProblem()	{
	delete []_startConfigJoints;
	delete []_startConfigCamera;
}
/*
	Except for the shoulder θs;  org_Ji = F(θ_1,...,θ_i-1)
	(InPtSHOULDER) == org_J3 == gives base position
	(InPtELBOW) == org_J4 == F(θ_1, θ_2, θ_3)
	(InPtWRIST) == org_J7 == F(θ_1,..., θ_6)
	(InPtTHUMB & InPtPINKIE) == P_EE == F(θ_1,..., θ_7)
*/
void mapInPtsHtoVecsR ()	{

		//length of WAM links from DH
		T l1wam_shoToUa = T(std::abs(_fk.getDH(2).getD()));
		T l2wam_elOffset = T(std::abs(_fk.getDH(2).getA()));
		T l3wam_wrOffset = T(std::abs(_fk.getDH(3).getA()));
		T l4wam_elToWr = T(std::abs(_fk.getDH(4).getD()));
		T l5wam_wrToEe = T(std::abs(_fk.getDH(
			_fk.getNJoints()-1).getD()));

		//input vectors to CERES solver
		for (size_t i = 0; i < _nPts; i++) {
			// (b.1) take sholder-to-elbow CROSS elbow-to-wrist
			Eigen::Matrix<T, 3, 1> shoToElVhN = 
				(_inPtsREl.col(i) - _inPtsRSh.col(i)).normalized()
				, elToWrVhN = 
					(_inPtsRW.col(i) - _inPtsREl.col(i)).normalized()
				, wrToThVhN = 
					(_inPtsRTh.col(i) - _inPtsRW.col(i)).normalized()
				, wrToPiVhN = 
					(_inPtsRPi.col(i) - _inPtsRW.col(i)).normalized();

			// (b.2) Elbow offset Vector x3 = -Z2.cross(Z3)
			// Z3 = UA x LA
			Eigen::Matrix<T, 3, 1> shoToElCrossEltoWrVsH = 
				shoToElVhN.cross(elToWrVhN).normalized();
			_elOffsetVsWam.col(i) = l2wam_elOffset *
				-shoToElVhN.cross(shoToElCrossEltoWrVsH).normalized();

			_shoToUaVsWam.col(i)= shoToElVhN * l1wam_shoToUa;
			_shoToElVsWam.col(i) = _shoToUaVsWam.col(i)
				+ _elOffsetVsWam.col(i);

			// (b.3) WRIST offset Vector
			// -x4 = -Y4{Z3}.cross(Z4)
			_wOffsetVsWam.col(i) = l3wam_wrOffset *
				-shoToElCrossEltoWrVsH.cross(elToWrVhN).normalized();
			_laToWrVsWam.col(i) = elToWrVhN * l4wam_elToWr;
			_elToWrVsWam.col(i) = _wOffsetVsWam.col(i) 
				+ _laToWrVsWam.col(i);
			_shoToLaVsWam.col(i) = _shoToElVsWam.col(i) 
				+ _wOffsetVsWam.col(i); 
			_shoToWrVsWam.col(i) = _shoToElVsWam.col(i) 
				+ _elToWrVsWam.col(i); 
			_wrToThVsWam.col(i) = wrToThVhN * l5wam_wrToEe;
			_wrToPiVsWam.col(i) = wrToPiVhN * l5wam_wrToEe;

			_shoToThVsWam.col(i) = _shoToWrVsWam.col(i) 
				+ _wrToThVsWam.col(i);
			_shoToPiVsWam.col(i) = _shoToWrVsWam.col(i) 
				+ _wrToPiVsWam.col(i);
			int iTimesThree = i*3;
			//in vecs along row
			_EeOrientation.block(iTimesThree, 0, 3, 3)
				= buildWrRotReferential(i).transpose();
			//Pee = Pwr + d6*Ree.[0 0 1]T
			_wrToEeVsWam.col(i) = l5wam_wrToEe * 
				buildWrRotReferential(i).col(2);
			_shoToEeVsWam.col(i) = _shoToWrVsWam.col(i) 
				+ _wrToEeVsWam.col(i);

		}
}

void mapInPtsHtoInPtsR ()	{
		//inPts WAM to CERES solver
		for (size_t i = 0; i < _nPts; i++) {
			_pShoWam.col(i)= _inPtsRSh.col(i);
			_pUaWam.col(i)= _pShoWam.col(i) + _shoToUaVsWam.col(i);
			_pElWam.col(i)= _pShoWam.col(i) + _shoToElVsWam.col(i);
			_pLaWam.col(i)= _pShoWam.col(i) + _shoToLaVsWam.col(i);
			_pWrWam.col(i)= _pShoWam.col(i) + _shoToWrVsWam.col(i);
			_pThWam.col(i)= _pShoWam.col(i) + _shoToThVsWam.col(i);
			_pPiWam.col(i)= _pShoWam.col(i) + _shoToPiVsWam.col(i);

		}
}

void rotInPtsToChestPlane() {
	for(size_t i = 0; i < _nPts; i++) {
		//ROTATION from camera to CHEST frame
		Eigen::Matrix<T, 3, 1> RCh = _inPtsRCh.col(i)
			, LCh = _inPtsLCh.col(i), MCh = _inPtsMCh.col(i);
		Eigen::Matrix<T, 3, 3> cameraToChestR =
		 buildRefFramefrom3Pts(RCh, MCh, LCh);
		//InPts InChest		
		_inPtsRShInChest.col(i) = cameraToChestR 
			* _inPtsRSh.col(i);
		_inPtsRUaInChest.col(i) = cameraToChestR 
			* _inPtsRUA.col(i);
		_inPtsRElInChest.col(i) = cameraToChestR 
			* _inPtsREl.col(i);
		_inPtsRLaInChest.col(i) = cameraToChestR 
			* _inPtsRLA.col(i);
		_inPtsRWrInChest.col(i) = cameraToChestR 
			* _inPtsRW.col(i);
		_inPtsRPiInChest.col(i) = cameraToChestR 
			* _inPtsRPi.col(i);
		_inPtsRThInChest.col(i) = cameraToChestR 
			* _inPtsRTh.col(i);

		//InVectors InChest
		_shoToUaVsWamInChest.col(i) = cameraToChestR 
			* _shoToUaVsWam.col(i);
	
		_shoToElVsWamInChest.col(i) = cameraToChestR 
			* _shoToElVsWam.col(i);
		_laToWrVsWamInChest.col(i) = cameraToChestR 
			* _laToWrVsWam.col(i);
		_elToWrVsWamInChest.col(i) = cameraToChestR 
			* _elToWrVsWam.col(i);
		
		_shoToLaVsWamInChest.col(i) = cameraToChestR 
			* _shoToLaVsWam.col(i);
		_shoToWrVsWamInChest.col(i) = cameraToChestR 
			* _shoToWrVsWam.col(i);

		_wrToThVsWamInChest.col(i) = cameraToChestR 
			* _wrToThVsWam.col(i);
		_wrToPiVsWamInChest.col(i) = cameraToChestR 
			* _wrToPiVsWam.col(i);
 		_wrToEeVsWamInChest.col(i) = cameraToChestR 
			* _wrToEeVsWam.col(i);
		_shoToThVsWamInChest.col(i) = cameraToChestR 
			* _shoToThVsWam.col(i);
		_shoToPiVsWamInChest.col(i) = cameraToChestR 
			* _shoToPiVsWam.col(i);
		_shoToEeVsWamInChest.col(i) = cameraToChestR 
			* _shoToEeVsWam.col(i);

		//CHECK
		_elOffsetVsWamInChest.col(i) = cameraToChestR 
			*_elOffsetVsWam.col(i);
		_wOffsetVsWamInChest.col(i) = cameraToChestR 
			* _wOffsetVsWam.col(i);

		int iTimesThree = i*3;
		Eigen::Matrix<T, 3, 3> EeOrientationInChestI =
			cameraToChestR * buildWrRotReferential(i);
		_EeOrientationInChest.block(iTimesThree, 0, 3, 3)
				= EeOrientationInChestI.transpose();

	}
}


const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getInPtsRCh() const {
	return _inPtsRCh;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getInPtsRSh() const {
	return _inPtsRSh;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getInPtsRShInChest() const {
	return _inPtsRShInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getInPtsMCh() const {
	return _inPtsMCh;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getInPtsLCh() const {
	return _inPtsLCh;
}

const Eigen::Matrix<T, 3, 1>
	getInWrToThVhN(int index) const {
	Eigen::Matrix<T, 3, 1> wrToThVhN = 
		(_inPtsRTh.col(index) - _inPtsRW.col(index)).normalized();
	return wrToThVhN;
}


const T* getShoToUaVi (int index)	const {
	T* pt = new T[3];
	pt[0] = _shoToUaVsWam(0,index);
	pt[1] = _shoToUaVsWam(1,index);
	pt[2] = _shoToUaVsWam(2,index);
	return pt;
}

const T* getShoToElVi (int index)	const {
	T* pt = new T[3];
	pt[0] = _shoToElVsWam(0,index);
	pt[1] = _shoToElVsWam(1,index);
	pt[2] = _shoToElVsWam(2,index);
	return pt;
}

const T* getElToWrVi (int index) const {
	T* pt = new T[3];
	pt[0] = _elToWrVsWam(0,index);
	pt[1] = _elToWrVsWam(1,index);
	pt[2] = _elToWrVsWam(2,index);
	return pt;
}

const T* getShoToLaVi (int index)	const {
	T* pt = new T[3];
	pt[0] = _shoToLaVsWam(0,index);
	pt[1] = _shoToLaVsWam(1,index);
	pt[2] = _shoToLaVsWam(2,index);
	return pt;
}

const T* getShoToWrVi (int index)	const {
	T* pt = new T[3];
	pt[0] = _shoToWrVsWam(0,index);
	pt[1] = _shoToWrVsWam(1,index);
	pt[2] = _shoToWrVsWam(2,index);
	return pt;
}
const T* getElOffsetVi (int index)	const {
	T* pt = new T[3];
	pt[0] = _elOffsetVsWam(0,index);
	pt[1] = _elOffsetVsWam(1,index);
	pt[2] = _elOffsetVsWam(2,index);
	return pt;
}

const T* getWrOffsetVi (int index)	const {
	T* pt = new T[3];
	pt[0] = _wOffsetVsWam(0,index);
	pt[1] = _wOffsetVsWam(1,index);
	pt[2] = _wOffsetVsWam(2,index);
	return pt;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getShoToUaVsWam() const {
	return _shoToUaVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToElVsWam() const {
	return _shoToElVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getShoToLaVsWam() const {
	return _shoToLaVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getLaToWrVsWam() const {
	return _laToWrVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getElToWrVsWam() const {
	return _elToWrVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToWrVsWam() const {
	return _shoToWrVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getWrToThVsWam() const {
	return _wrToThVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToThVsWam() const {
	return _shoToThVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getWrToPiVsWam() const {
	return _wrToPiVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToPiVsWam() const {
	return _shoToPiVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getWrToEeVsWam() const {
	return _wrToEeVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToEeVsWam() const {
	return _shoToEeVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getElOffsetVsWam() const {
	return _elOffsetVsWam;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getWrOffsetVsWam() const {
	return _wOffsetVsWam;
}


const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getShoToUaVsWamInChest() const {
	return _shoToUaVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToElVsWamInChest() const {
	return _shoToElVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getShoToLaVsWamInChest() const {
	return _shoToLaVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getElToWrVsWamInChest() const {
	return _elToWrVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getLaToWrVsWamInChest() const {
	return _laToWrVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToWrVsWamInChest() const {
	return _shoToWrVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getWrToThVsWamInChest() const {
	return _wrToThVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToThVsWamInChest() const {
	return _shoToThVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getWrToPiVsWamInChest() const {
	return _wrToPiVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToPiVsWamInChest() const {
	return _shoToPiVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getWrToEeVsWamInChest() const {
	return _wrToEeVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
	getShoToEeVsWamInChest() const {
	return _shoToEeVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getElOffsetVsWamInChest() const {
	return _elOffsetVsWamInChest;
}

const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	getWrOffsetVsWamInChest() const {
	return _wOffsetVsWamInChest;
}

const Eigen::Matrix<T, 3, 1> getLaPlaneNormalVhN(int index) 	
	const {
	// (b.1) take la-to-wr CROSS la-to-el
	Eigen::Matrix<T, 3, 1> laToWrVhN = 
		(_inPtsRLA.col(index) - _inPtsRW.col(index)).normalized()
		, laToElVhN = (_inPtsRLA.col(index) 
			- _inPtsREl.col(index)).normalized();
	const Eigen::Matrix<T, 3, 1> retV 
		= laToWrVhN.cross(laToElVhN).normalized();
	return retV;
}

//take sholder-to-elbow CROSS elbow-to-wrist
const Eigen::Matrix<T, 3, 1> getUaCrossLaVhN(int index) 	
	const {
	Eigen::Matrix<T, 3, 1> shoToElVhN = (_inPtsREl.col(index) 
		- _inPtsRSh.col(index)).normalized()
		, elToWrVhN = (_inPtsRW.col(index) 
			- _inPtsREl.col(index)).normalized();
	const Eigen::Matrix<T, 3, 1> retV = 
		shoToElVhN.cross(elToWrVhN).normalized();
	return retV;
}

const T getNpts()	const {
	return _inPtsAlongRows.rows();
}

const T getPtsBlockSize()	const {
	//Last three joints only contribute to Ee orientation
	const int  nJoints = _fk.getNJoints()-3;
	return  nJoints * (_polynomialOrder + 1);
}

T* getStartConfigJoints()	{
	return _startConfigJoints;
}

const T getCameraBlockSize()	{
	//3 for rotation, 3 for trans
	const int cameraBlockSize = 3; 
	return cameraBlockSize;
}

T* getStartConfigCamera()	{
	return _startConfigCamera;
}

ForwardKin<T>& getWam() {
	return _fk;
}

std::vector<T> getJointMinAngles() {
	return _jointMinAngles;
}

std::vector<T> getJointMaxAngles() {
	return _jointMaxAngles;
}

Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& getJsWam() {
	return _JsWam;
}

void setJointLimits(Problem* problem, double* startConfigJ, int toDof)	{

		//set the joint limits into ceres for Theta_1, Theta_2
		int ctr = 0;

		size_t  nJoints = _fk.getNJoints();
		for (int joint = 0; joint < nJoints && joint < toDof; joint++) {
			problem->SetParameterLowerBound(startConfigJ, ctr
				, _jointMinAngles[joint]);
			problem->SetParameterUpperBound(startConfigJ, ctr
				, _jointMaxAngles[joint]);
			ctr++;
			for(size_t term = 0; term < _polynomialOrder; term++) {
				problem->SetParameterLowerBound(startConfigJ, ctr
					, _jointMinAngles[joint]);
				problem->SetParameterUpperBound(startConfigJ, ctr
					, _jointMaxAngles[joint]);
				ctr++;
			}
		}

}

const T computeJ4FromHumanElbAngle(int index) const {
	/*From Kim Rosen-2015 Eqn. (2)*/
	T lUa = (_inPtsREl.col(index) -
		_inPtsRSh.col(index)).norm();
	T lLa = (_inPtsRW.col(index) -
		_inPtsREl.col(index)).norm();
	T lWrToShV = 	(_inPtsRW.col(index) 
		- _inPtsRSh.col(index)).norm();
	
	T top = (lLa * lLa) + (lUa * lUa) 
		- (lWrToShV * lWrToShV);
	T bot = T(2.0) * lLa * lUa; 

	cout << "lUa " << lUa << endl;
	cout << "lLa " << lLa << endl;
	cout << "lWrToShV " << lWrToShV << endl;
	cout << "top " << top << endl;
	cout << "bot " << bot << endl;

	return T(M_PI) - T(top / bot);
/*
	Eigen::Matrix<T, 3, 1> shoToElVhN = 
		(_inPtsREl.col(index) - _inPtsRSh.col(index)).normalized()
		, elToWrVhN = 
			(_inPtsRW.col(index) - _inPtsREl.col(index)).normalized();

	// (b.1.2) find Human elbow angle, i.e. <)l0.l1
	T elThetaH = T(ceres::acos(shoToElVhN.dot(elToWrVhN)));

	return T(M_PI) - elThetaH;
*/
}

const Eigen::Matrix<T, 3, 3> buildWrRotReferential(int index) 
	const {
	/* First, a vector is defined from the wr marker to the
		 th, and a second vector from the wr marker to the 
			pinkie.*/
	Eigen::Matrix<T, 3, 1> wrToThVhN = 
		(_inPtsRTh.col(index) - _inPtsRW.col(index)).normalized()
		, wrToPiVhN = 
			(_inPtsRPi.col(index) - _inPtsRW.col(index)).normalized();
	Eigen::Matrix<T, 3, 3> retT;
	/* Z	The sum of these two vectors lies on the hand plane
		 and is parallel to the z-axis (i.e. axis of rotation)
		 of the rotated referential.*/
	retT.col(2) = (wrToPiVhN + wrToThVhN).normalized();
	/* The cross-product of these vectors gives the 
			normal to the hand plane (pointing down) and 
			is parallel to the x-axis  */
	retT.col(0) = (wrToPiVhN.cross(wrToThVhN)).normalized();
 	/* The remaining Y axis can easily be computed 
		 as the cross-product of Z cross X axis */
	retT.col(1) = (retT.col(2).cross(retT.col(0))).normalized();
	return retT;
}

const Eigen::Matrix<T, 4, 4> getHwrInCamera(int index) 
	const {

	Eigen::Matrix<T, 4, 4> retT 
		= Eigen::Matrix<T, 4, 4>::Identity(4,4);
	retT.block(0,0,3,3) = buildWrRotReferential(index);
	retT.block(0,3,3,1) =  _inPtsRW.col(index);
	return retT;
}

protected:
	const int  _polynomialOrder = 3;
	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &_JsWam;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		&_inPtsAlongRows;
	int _nPts;
	Eigen::Matrix<T, 3, Eigen::Dynamic> _inPtsRW,
			_inPtsRTh, _inPtsRPi, _inPtsRLA, _inPtsREl, _inPtsRUA
			, _inPtsRSh, _inPtsRCh, _inPtsMCh, _inPtsLCh, _inPtsLSh;

	std::vector<T> _jointMinAngles, _jointMaxAngles;

	T* _startConfigJoints;
	T* _startConfigCamera;

	Eigen::Matrix<T,3, Eigen::Dynamic> 
		_shoToUaVsWam, _shoToElVsWam
		, _laToWrVsWam, _elToWrVsWam
		, _shoToLaVsWam, _shoToWrVsWam
		,	_wrToThVsWam, _wrToPiVsWam
		, _shoToThVsWam, _shoToPiVsWam
		, _wrToEeVsWam,   _shoToEeVsWam;

	Eigen::Matrix<T, 3, Eigen::Dynamic> 
		_elOffsetVsWam, _wOffsetVsWam;

	Eigen::Matrix<T, 3, Eigen::Dynamic> 
		_pShoWam, _pUaWam, _pElWam, _pLaWam
		, _pWrWam, _pThWam, _pPiWam;

	Eigen::Matrix<T, 3, Eigen::Dynamic> 
		_inPtsRWrInChest, _inPtsRThInChest, _inPtsRPiInChest
		, _inPtsRLaInChest, _inPtsRElInChest, _inPtsRUaInChest
		, _inPtsRShInChest;
	
	Eigen::Matrix<T, 3, Eigen::Dynamic> 
		_shoToUaVsWamInChest, _shoToElVsWamInChest
		, _laToWrVsWamInChest, _elToWrVsWamInChest
		, _shoToLaVsWamInChest, _shoToWrVsWamInChest
		,	_wrToThVsWamInChest, _wrToPiVsWamInChest
		, _shoToThVsWamInChest, _shoToPiVsWamInChest
		, _wrToEeVsWamInChest, _shoToEeVsWamInChest;

	Eigen::Matrix<T, 3, Eigen::Dynamic> 
		_elOffsetVsWamInChest, _wOffsetVsWamInChest;

	Eigen::Matrix<T, Eigen::Dynamic,3> 
		_EeOrientation, _EeOrientationInChest;

};

#endif
