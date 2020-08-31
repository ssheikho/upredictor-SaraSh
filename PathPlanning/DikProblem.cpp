#include "DikProblem.h"
#include "BasicFormulas.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Core"
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "glog/logging.h"

using namespace ceres;

DikProblem::DikProblem (ForwardKin<double> &fk
	, std::vector<double> jointMinAngles
	, std::vector<double> jointMaxAngles
	, bool use_quaternions
	, Eigen::MatrixXd inPts	
	) : _fk(fk)
	, _jointMinAngles(jointMinAngles)
	, _jointMaxAngles(jointMaxAngles)	
	, _use_quaternions(use_quaternions)
	, _inPts(inPts)
	, _nPts(inPts.cols())
	, _inPRPi(_inPts.block(Markers::RPi,0,3,_nPts))
	, _inPRTh(_inPts.block(Markers::RTh,0,3,_nPts))
	, _inPRWr(_inPts.block(Markers::RWr,0,3,_nPts))
	, _inPRLA(_inPts.block(Markers::RLA,0,3,_nPts))
	, _inPREl(_inPts.block(Markers::REl,0,3,_nPts))
	, _inPRUA(_inPts.block(Markers::RUA,0,3,_nPts))
	, _inPRSh(_inPts.block(Markers::RSh,0,3,_nPts))
	,	_inPRCh(_inPts.block(Markers::RCh,0,3,_nPts))
	,	_inPMCh(_inPts.block(Markers::MCh,0,3,_nPts))
	,	_inPLCh(_inPts.block(Markers::LCh,0,3,_nPts))
	,	_inPLSh(_inPts.block(Markers::LSh,0,3,_nPts))

	, _shoToUaVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _shoToElVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _laToWrVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _elToWrVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _shoToLaVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _shoToWrVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _shoToThVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _shoToPiVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _wrToThVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _wrToPiVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _wrToEeVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _shoToEeVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _elOffsetVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _wOffsetVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	//Orientation as mathematica matrices:
	// Chest
	, _quaternionsCh(Eigen::MatrixXd::Zero(4,_nPts))
	, _rMatsCh(Eigen::MatrixXd::Zero(3,3*_nPts))
	, _xyzEulersCh(Eigen::MatrixXd::Zero(3,_nPts))
	// end-effectors
	, _quaternionsEe(Eigen::MatrixXd::Zero(4,_nPts))
	, _rMatsEe(Eigen::MatrixXd::Zero(3,3*_nPts))
	, _xyzEulersEe(Eigen::MatrixXd::Zero(3,_nPts))

	,	_inPtsWam(Eigen::MatrixXd::Zero(inPts.rows(),_nPts))
	, _pEeWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _pPiWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _pThWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _pWrWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _pLaWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _pElWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _pUaWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _pShWam(Eigen::MatrixXd::Zero(3,_nPts))

	,	_inLvsWam(Eigen::MatrixXd::Zero(inPts.rows(),_nPts))
	,	_lvelShWam(Eigen::MatrixXd::Zero(3,_nPts))
	,	_lvelUaWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _lvelElWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _lvelLaWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _lvelWrWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _lvelThWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _lvelPiWam(Eigen::MatrixXd::Zero(3,_nPts))
	,	_AngularVWrWam(Eigen::MatrixXd::Zero(3,_nPts))
	//OUT-scaled VECTORS on the robot
	,	_outShoToUaVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outShoToElVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outElToWrVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outLaToWrVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outShoToLaVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outShoToWrVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	,	_outWrToThVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outWrToPiVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outShoToThVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outShoToPiVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outWrToEeVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	,	_outShoToEeVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outElOffsetVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outWOffsetVsWam(Eigen::MatrixXd::Zero(3,_nPts))
	//OUT-scaled positions on the robot
	, _outPtsWam(Eigen::MatrixXd::Zero(inPts.rows(),_nPts))
	, _outPShWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outPUaWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outPElWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outPLaWam(Eigen::MatrixXd::Zero(3,_nPts))	
	, _outPWrWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outPThWam(Eigen::MatrixXd::Zero(3,_nPts))
	, _outPPiWam(Eigen::MatrixXd::Zero(3,_nPts)) {

  _theta_index = new int[_nPts];
  _qBase_index = new int[_nPts];	
	// the angle-axis rotations.
	_nParameters = 3 * _nPts + 7 * _nPts;
  _parameters = new double[_nParameters];
	// initialize qBase to identity
  for (int i = 0; i < _nParameters; ++i) 
		*(_parameters + i) = 0.0;
  
	// Switch the angle-axis rotations to quaternions.
	if(_use_quaternions) {
		_nParameters = 4 * _nPts + 7 * _nPts;
    double* quaternion_parameters = 
			new double[_nParameters];
		double* original_cursor = _parameters;
    double* quaternion_cursor = quaternion_parameters;
    for (int i = 0; i < _nPts; ++i) {
      AngleAxisToQuaternion(
				original_cursor, quaternion_cursor);
      quaternion_cursor += 4;
      original_cursor += 3;
    }
    // Copy the rest of the parameters, i.e. thetas.
    for (int i = 0; i < 7 * _nPts; ++i) {
      *quaternion_cursor++ = *original_cursor++;
    }
    // Swap in the quaternion parameters.
    delete []_parameters;
    _parameters = quaternion_parameters;
  }

	printInPtsH();
	mapInPtsHtoVecsR();
	mapInPtsHtoInPtsR();
	printPtsR();
	printVecsR();
}

DikProblem::~DikProblem () {
	delete []_theta_index;
  delete []_qBase_index;
  delete []_parameters;
}


void DikProblem::mapInPtsHtoVecsR ()	{
	//length of WAM links from DH
	double l1wam_shoToUa = std::abs(_fk.getDH(2).getD())
		, l2wam_elOffset = std::abs(_fk.getDH(2).getA())
		, l3wam_wrOffset = std::abs(_fk.getDH(3).getA())
		, l4wam_elToWr = std::abs(_fk.getDH(4).getD())
		, l5wam_wrToEe = std::abs(_fk.getDH(
			_fk.getNJoints()-1).getD());

	//input vectors to CERES solver
	for (int i = 0; i < _nPts; i++) {
		// (b.1) take sholder-to-elbow CROSS elbow-to-wrist
		Eigen::Matrix<double, 3, 1> shoToElVhN = 
			(_inPREl.col(i) - _inPRSh.col(i)).normalized()
			, elToWrVhN = 
				(_inPRWr.col(i) - _inPREl.col(i)).normalized()
			, wrToThVhN = 
				(_inPRTh.col(i) - _inPRWr.col(i)).normalized()
			, wrToPiVhN = 
				(_inPRPi.col(i) - _inPRWr.col(i)).normalized();

		// (b.2) Elbow offset Vector x3 = -Z2.cross(Z3)
		// Z3 = UA x LA
		Eigen::Matrix<double, 3, 1> shoToElCrossEltoWrVsH = 
			shoToElVhN.cross(elToWrVhN).normalized();
		_elOffsetVsWam.col(i) = 
			l2wam_elOffset * (-shoToElVhN.cross(
				shoToElCrossEltoWrVsH).normalized());
		_shoToUaVsWam.col(i)= shoToElVhN * l1wam_shoToUa;
		_shoToElVsWam.col(i) = _shoToUaVsWam.col(i)
			+ _elOffsetVsWam.col(i);

		// (b.3) WRIST offset Vector
		// -x4 = -Y4{Z3}.cross(Z4)
		_wOffsetVsWam.col(i) = 
			l3wam_wrOffset *	(-shoToElCrossEltoWrVsH.cross(
				elToWrVhN).normalized());
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

		//Orientation as mathematica matrices:
		int iTimesThree = i*3;
		//chest Rot matrix
		Mat3 rMatsCh_i = buildRefFramefrom3Pts(
			RChP(i), MChP(i), LChP(i));
		_rMatsCh.block(0,iTimesThree, 3, 3) =rMatsCh_i;
		//chest XYZ-euler angles
		_xyzEulersCh.col(i) = rMatsCh_i.eulerAngles(0, 1, 2);
		// end-effectors
		Mat3 rMatsEe_i = buildWrRotReferential(i);
	 _rMatsEe.block(0,iTimesThree, 3, 3)
			= buildWrRotReferential(i);
		_xyzEulersEe.col(i) = 
			rMatsEe_i.eulerAngles(0, 1, 2);
	
		//Pee = Pwr + d6*Ree.[0 0 1]T
		_wrToEeVsWam.col(i) = l5wam_wrToEe * 
			buildWrRotReferential(i).col(2);
		_shoToEeVsWam.col(i) = _shoToWrVsWam.col(i) 
			+ _wrToEeVsWam.col(i);

	}
}

void DikProblem::mapInPtsHtoInPtsR ()	{
	//inPts WAM to CERES solver
	for (size_t i = 0; i < _nPts; i++) {
		_inPtsWam.block(Markers::RSh,i,3,1) =
			_inPRSh.col(i);
		_inPtsWam.block(Markers::RUA,i,3,1) =
			_inPRSh.col(i) + _shoToUaVsWam.col(i);
		_inPtsWam.block(Markers::REl,i,3,1) =
			_inPRSh.col(i) + _shoToElVsWam.col(i);
		_inPtsWam.block(Markers::RLA,i,3,1) =
			_inPRSh.col(i) + _shoToLaVsWam.col(i);
		_inPtsWam.block(Markers::RWr,i,3,1) =
			_inPRSh.col(i) + _shoToWrVsWam.col(i);
		_inPtsWam.block(Markers::RTh,i,3,1) =
			_inPRSh.col(i) + _shoToThVsWam.col(i);
		_inPtsWam.block(Markers::RPi,i,3,1) =
			_inPRSh.col(i) + _shoToPiVsWam.col(i);
	}
	_pEeWam = _inPRSh + _shoToEeVsWam;
	_pPiWam = _inPtsWam.block(Markers::RPi,0,3,_nPts);
	_pThWam = _inPtsWam.block(Markers::RTh,0,3,_nPts);
  _pWrWam = _inPtsWam.block(Markers::RWr,0,3,_nPts);
	_pLaWam = _inPtsWam.block(Markers::RLA,0,3,_nPts);
	_pElWam = _inPtsWam.block(Markers::REl,0,3,_nPts);
	_pUaWam = _inPtsWam.block(Markers::RUA,0,3,_nPts);
	_pShWam = _inPtsWam.block(Markers::RSh,0,3,_nPts);

}


Eigen::Matrix<double, 3, 3> DikProblem::buildWrRotReferential
	(int index) {
	/* First, a vector is defined from the wr marker to the
		 th, and a second vector from the wr marker to the 
			pinkie.*/
	Eigen::Matrix<double, 3, 1> 
		wrToThVhN = (_inPRTh.col(index) - 
			_inPRWr.col(index)).normalized()
		, wrToPiVhN = (_inPRPi.col(index) - 
				_inPRWr.col(index)).normalized();
	Eigen::Matrix<double, 3, 3> retT;
	/* Z	The sum of these two vectors lies on the hand 
			plane and is parallel to the z-axis (i.e. axis 
			of rotation) of the rotated referential.*/
	retT.col(2) = (wrToPiVhN + wrToThVhN).normalized();
	/* The cross-product of these vectors gives the 
			normal to the hand plane (pointing down) and 
			is parallel to the x-axis  */
	retT.col(0) = (wrToPiVhN.cross(wrToThVhN)).normalized();
 	/* The remaining Y axis can easily be computed 
		 as the cross-product of Z cross X axis */
	retT.col(1) = (retT.col(2).cross(
		retT.col(0))).normalized();
	return retT;
}

void DikProblem::setOutVecsR () {

}

void DikProblem::setOutPtsR () {

}

// For print in Mathematica 
//(a) for human
void DikProblem::printInPtsH()	{
	printEigenMathematica	(_inPRPi.transpose(),
		cout, "inPtsRPi");	
	printEigenMathematica	(_inPRTh.transpose(),
		cout, "inPtsRTh");	
	printEigenMathematica	(_inPRWr.transpose(),
		cout, "inPtsRW");	
	printEigenMathematica	(_inPRLA.transpose()
		, cout, "inPtsRLA");	
	printEigenMathematica	(_inPREl.transpose()
		, cout, "inPtsREl");	
	printEigenMathematica	(_inPRUA.transpose()
		, cout, "inPtsRUA");	
	printEigenMathematica	(_inPRSh.transpose()
		, cout, "inPtsRSh");	
	printEigenMathematica	(_inPRCh.transpose()
		, cout, "inPtsRCh");	
	printEigenMathematica	(_inPMCh.transpose()
		, cout, "inPtsMCh");	
	printEigenMathematica	(_inPLCh.transpose()
		, cout, "inPtsLCh");	
	printEigenMathematica	(_inPLSh.transpose()
		, cout, "inPtsLSh");	
}

//(b) inPts wam, i.e. scaled to wam
void DikProblem::printPtsR()	{
	printEigenMathematica( _pEeWam.transpose()
		, cout, "inPtsEeWam");	
	printEigenMathematica( _pPiWam.transpose()
		, cout, "inPtsRPiWam");	
	printEigenMathematica( _pThWam.transpose()
		, cout, "inPtsRThWam");	
	printEigenMathematica( _pWrWam.transpose()
		, cout, "inPtsRWWam");	
	printEigenMathematica( _pLaWam.transpose()
		, cout, "inPtsRLAWam");	
	printEigenMathematica( _pElWam.transpose()
		, cout, "inPtsRElWam");	
	printEigenMathematica( _pUaWam.transpose()
		, cout, "inPtsRUAWam");	
	printEigenMathematica( _pShWam.transpose()
		, cout, "inPtsRShWam");	
}

//(c). inVecs Robot
void DikProblem::printVecsR()	{
	printEigenMathematica( _shoToUaVsWam.transpose()
		, cout, "shoToUaVsWam");	
	printEigenMathematica( _elOffsetVsWam.transpose()
		, cout, "elOffsetVsWam");	
	printEigenMathematica( _wOffsetVsWam.transpose()
		, cout, "wOffsetVsWam");	
	printEigenMathematica( _laToWrVsWam.transpose()
		, cout, "laToWrVsWam");	
	printEigenMathematica( _elToWrVsWam.transpose()
		, cout, "elToWrVsWam");	
	printEigenMathematica( _wrToEeVsWam.transpose()
		, cout, "wrToEeVsWam");	

	printEigenMathematica( _shoToElVsWam.transpose()
		, cout, "shoToElVsWam");	
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
	//Orientation as mathematica matrices:
	// Chest
	printEigenMathematica(_rMatsCh
		, cout, "rMatsCh");	
	printEigenMathematica(_xyzEulersCh.transpose()
		, cout, "xyzEulersCh");	
	// end-effectors
	printEigenMathematica(_rMatsEe
		, cout, "rMatsEe");	
	printEigenMathematica(_xyzEulersEe.transpose()
		, cout, "xyzEulersEe");
}

