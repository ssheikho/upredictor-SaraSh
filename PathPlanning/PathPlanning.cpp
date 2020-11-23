#include <iostream>

#include <Plane.h>
#include <UBCUtil.h>
#include <Plane.h>
#include <LinearAlgebra.h>
#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "ForwardKin.h"
#include "DikProblem.h"
#include "DikSolver.h"
#include "MayaAnimation.h"


#include "Ellipse3D.h"
#include "ParseMathematica.h"
#include "ParseCSV.h"
#include <Eigen/SVD>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>       /* pow */
#include "FitFunctions.h"
#include "LinearAlgebra.cpp"
#include "RigidTrans2D.h"

#include "BaysFitFunctions.h"
#include "EllipseConicConstraints.h"
using namespace ceres;
using namespace Eigen;
using namespace std;
/*
//ensure theta is within [-pi,pi], 
//	&& theta(t) is a monotonic function 
MatrixXd fixThetas(MatrixXd a) {
	MatrixXd retVal(a.rows(), a.cols());

	double curTheta = a(0,0);
	while(curTheta >= -M_PI) curTheta -= M_PI;
	while(curTheta <= M_PI) curTheta += M_PI;
	retVal(0,0) = curTheta;

	for(size_t i = 1; i < a.cols(); i++) {
		double curTheta = a(0,i);
		while(curTheta >= -M_PI) curTheta -= 2.0 * M_PI;
		while(curTheta <= M_PI) curTheta += 2.0 * M_PI;

		double testA = curTheta + 2.0 * M_PI;
		double testB = curTheta - 2.0 * M_PI;
		
		double distNull = curTheta - retVal(0, i - 1);
		double distA = testA - retVal(0, i - 1);
		double distB = testB - retVal(0, i - 1);

		retVal(0, i) = distNull < distA ? curTheta : testA;
		retVal(0, i) = distNull < distB ? curTheta : testB;
	}

	return retVal;
}
*/

int main(int argc, char **argv) {
	
	string fileName = argv[1];
	cout << "(* " << fileName << "*)" << endl;

	/************************************************/
	/*** Copy input marker data into Eigen Matrix ***/
	std::vector<std::string> wristMarkers =
 		{"x", "y", "z"};
	std::vector<std::string> wSMarkers =	
	{	
   "RWristX",	"RWristY",	"RWristZ" 
	,	"RThumbX",	"RThumbY",	"RThumbZ" 
	, "RLArmX",	"RLArmY",	"RLArmZ"		
	, "RElbX",	"RElbY",	"RElbZ"				
	, "RUArmX",	"RUArmY",	"RUArmZ"		
	, "RShoX",	"RShoY",	"RShoZ"
	, "RChestX",	"RChestY",	"RChestZ"	
	, "MChestX",	"MChestY",	"MChestZ"
	, "LChestX", "LChestY", "LChestZ"
	, "LShoX", "LShoY", "LShoZ" 
	, "RPinkieX",	"RPinkieY",	"RPinkieZ" };
	int nMarkers = wSMarkers.size();

	// nx3
	std::string metaFileName;
	std::string reachCtr = getPartID(fileName) + "P"
		+fileName.at ((unsigned)(fileName.length()-5));

	int phase = 1;
	//whichPhase(fileName);		
	// Calculating the length of string 
	metaFileName = 
		fileName.substr(0,fileName.find("_clean"))+".csv";	
	string headerLine = 
		getHeaderRowCSV(metaFileName);

	int nPts = countRowsCSV(fileName, wristMarkers[0]);
	MatrixXd inPtsAlongRows = MatrixXd::Zero(nPts,nMarkers);
	int nPtsTotal = countRowsCSV(
		metaFileName, wSMarkers[0]);
	MatrixXd inPtsAlongRowsTotal = 
		MatrixXd::Zero(nPtsTotal,nMarkers);

	int ctr = 0;
	while ((ctr < nMarkers)) {
		int startMarkerDelim = 
			headerLine.find(wSMarkers[ctr]);
		if (startMarkerDelim!=std::string::npos) { 
			fillMatCSV(ctr, metaFileName
				, wSMarkers[ctr], wSMarkers[ctr+1]
				, wSMarkers[ctr+2], inPtsAlongRowsTotal);
		}
		ctr +=3;
	}
/*
	cout << "Phase " << phase << endl;
	if (phase ==1)	
		inPtsAlongRows = inPtsAlongRowsTotal.topRows(nPts);
	
	else	
		inPtsAlongRows = inPtsAlongRowsTotal.bottomRows(nPts);
*/


	/*** END - Reading input CSV Motion Files ***/
	/************************************************/

	/************************************************/
	/*** Wam IK SOLVER 
		(1) CERES solve 
	  given marker inPts, solve for IK of WAM & q_base **/
	// joint limits
	std::vector<double> jointMinAngles, jointMaxAngles;
	jointMinAngles.push_back(-2.6);
	jointMinAngles.push_back(-2.0);
	jointMinAngles.push_back(-2.8);
	jointMinAngles.push_back(0.0/*-0.9*/);
	jointMinAngles.push_back(-4.76);
	jointMinAngles.push_back(-1.6);
	jointMinAngles.push_back(-3.0);

	jointMaxAngles.push_back(2.6);
	jointMaxAngles.push_back(2.0);
	jointMaxAngles.push_back(2.8);
	jointMaxAngles.push_back(3.1);
	jointMaxAngles.push_back(1.24);
	jointMaxAngles.push_back(1.6);
	jointMaxAngles.push_back(3.0);

	//generates the DH T matrix for all the joints	     
	ForwardKin<double> wam =
		ForwardKin<double>::generateBarrettWAM();
	// number of optimization parameters 
	size_t toDof = 7;
	size_t nJoints = 7;
	// positions in meter
	MatrixXd inPts = 0.01 * inPtsAlongRows.transpose();
	MatrixXd inPtsTot = 0.01 * inPtsAlongRowsTotal.transpose();


	/*** (1) Human motion from vicon data**/
	/*********************************************************/
	MatrixXd JsWam = MatrixXd::Zero(nJoints,  nPts); 
	// (a). inverse position
	bool use_quaternions = true;

/*
	DikProblem *dIkProbTot =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			inPtsTot);

	DikProblem *dIkProb =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			inPts);

	Eigen::MatrixXd xyzEulersCh = 
		MayaAnimation::getChEulerAngles (dIkProb);
	printEigenMathematica(xyzEulersCh.transpose()
		, cout, "xyzEulersCh");
	printEigenMathematica(xyzEulersCh.transpose()*(180.0/M_PI)
		, cout, "xyzEulersChInDeg");
	Eigen::MatrixXd xyzEulersEl = 
		MayaAnimation::getElEulerAngles (dIkProb);
	printEigenMathematica(xyzEulersEl.transpose()
		, cout, "xyzEulersEl");
	printEigenMathematica(xyzEulersEl.transpose()*(180.0/M_PI)
		, cout, "xyzEulersElInDeg");

	Eigen::MatrixXd xyzEulersWr = MayaAnimation::
		getWrEulerAngles (dIkProb);
	printEigenMathematica(xyzEulersWr.transpose()
		, cout, "xyzEulersWr");
	printEigenMathematica(xyzEulersWr.transpose()*(180.0/M_PI)
		, cout, "xyzEulersWrInDeg");

	// (a.2). inverse position in Chest frame
	MatrixXd inPtsInBase = MatrixXd::Zero(nMarkers,nPts),
		inPtsInChest = MatrixXd::Zero(nMarkers,nPts),
		inPRCh = inPts.block(Markers::RCh,0,3,nPts),
		inPMCh = inPts.block(Markers::MCh,0,3,nPts),
		inPLCh = inPts.block(Markers::LCh,0,3,nPts);

	Eigen::MatrixXd xyzEulersChF3Pts  = 
		MatrixXd::Zero(3,nPts);
	for (int i = 0; i < nPts; i++) {
		Mat3 rMatsBase_i = MayaAnimation::
			buildBaseRotReferentialAt (dIkProb, i);
		//chest Rot matrix 
		Eigen::Matrix<double, 3, 1> 
			RChP = inPRCh.col(i), 
			MChP = inPMCh.col(i),
			LChP = inPLCh.col(i);
		
		Mat3 rMatsCh_i = buildRefFramefrom3Pts(
			RChP, MChP, LChP);
		xyzEulersChF3Pts.col(i) = rMatsCh_i.eulerAngles(0, 1, 2); 
		//inPts in Chest frame
		for (size_t j = 0; j < nMarkers; ){	
			inPtsInChest.block(j,i,3,1) = 
				rMatsCh_i.transpose() * inPts.block(j,i,3,1);
			j = j+3;
		}
	}

	printEigenMathematica(xyzEulersChF3Pts.transpose()
		, cout, "xyzEulersChF3Pts");

	DikProblem *dIkProbInChest =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			inPtsInChest);

	for (int i = 0; i < nPts; i++)
		ceres::SolveProblemAt(dIkProbInChest, i);

	ceres::SolveProblem(dIkProbInChest);

	Eigen::MatrixXd xyzEulersBase = 
		MayaAnimation::getBaseEulerAngles (dIkProbInChest);
	printEigenMathematica(xyzEulersBase.transpose()
		, cout, "xyzEulersBase");
	printEigenMathematica(xyzEulersBase.transpose()*(180.0/M_PI)
		, cout, "xyzEulersBaseInDeg");

	// (b). inverse orientation
//	for (int i = 0; i < nPts; i++)
//		ceres::SolveOrientationProblemAt(dIkProbInChest, i);

//	ceres::SolveOrientationProblem(dIkProbInChest);
*/


	/*** (2) ELLIPSE fitted model of Human motion **/
	/*********************************************************/

	//Read ellipse fitted data JUST FOR THE WRIST
/*
	string wholeFile = getWholeFileMathematica("outEf.txt");
	MatrixXd efPtsHand = inPts;

	MatrixXd efPtsAlongRows = MatrixXd::Zero(nPts,nMarkers);
	efPtsAlongRows = parseMathematica(wholeFile);
	MatrixXd efPts = efPtsAlongRows.transpose();

	efPtsHand.block(Markers::RPi,0,3,nPts) 
		= efPts.block(Markers::RPi,0,3,nPts);
	efPtsHand.block(Markers::RTh,0,3,nPts) 
		= efPts.block(Markers::RTh,0,3,nPts);
	efPtsHand.block(Markers::RWr,0,3,nPts) 
		= efPts.block(Markers::RWr,0,3,nPts);

	DikProblem *dIkProb =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			efPtsHand);
*/

	//Read ellipse fitted data
	string wholeFile = getWholeFileMathematica("Ef10-2.txt");
	MatrixXd efPtsAlongRows = MatrixXd::Zero(nPts,nMarkers);
	efPtsAlongRows = parseMathematica(wholeFile);


	//get inPts in MotionBuilder frame
	MatrixXd efPtsAlongRowsMB = MotionBuilderAnimation
		::getInPtsInMbFrame(efPtsAlongRows);
	MatrixXd efPtsMB =  100.0*efPtsAlongRowsMB.transpose();
	

	DikProblem *dIkProbEfMB =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			efPtsMB);

	int startFrame = 593;//FIX
	MotionBuilderAnimation::keyAdd(efPtsMB, startFrame);
	
	cout <<" END inPts in MotionBuilder frame"<< endl;

	MatrixXd efPts = 100.0*efPtsAlongRows.transpose();

	DikProblem *dIkProbEf =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			efPts);

	Eigen::MatrixXd xyzEulersChEf = 
		MayaAnimation::getChEulerAngles (dIkProbEf);
	printEigenMathematica(xyzEulersChEf.transpose()
		, cout, "xyzEulersChEf");
	printEigenMathematica(xyzEulersChEf.transpose()*(180.0/M_PI)
		, cout, "xyzEulersChEfInDeg");

	Eigen::MatrixXd xyzEulersElEf = 
		MayaAnimation::getElEulerAngles (dIkProbEf);
	printEigenMathematica(xyzEulersElEf.transpose()
		, cout, "xyzEulersElEf");
	printEigenMathematica(xyzEulersElEf.transpose()*(180.0/M_PI)
		, cout, "xyzEulersElEfInDeg");

	Eigen::MatrixXd xyzEulersWrEf = MayaAnimation::
		getWrEulerAngles (dIkProbEf);
	printEigenMathematica(xyzEulersWrEf.transpose()
		, cout, "xyzEulersWrEf");
	printEigenMathematica(xyzEulersWrEf.transpose()*(180.0/M_PI)
		, cout, "xyzEulersWrEfInDeg");



	// efPts in Chest frame
	MatrixXd 
		efPtsInChest = MatrixXd::Zero(nMarkers,nPts),
		inPRCh = inPts.block(Markers::RCh,0,3,nPts),
		inPMCh = inPts.block(Markers::MCh,0,3,nPts),
		inPLCh = inPts.block(Markers::LCh,0,3,nPts);
	Eigen::MatrixXd xyzEulersChF3Pts  = 
		MatrixXd::Zero(3,nPts);
	for (int i = 0; i < nPts; i++) {
		int iTimesThree = i*3;
		//chest Rot matrix 
		Eigen::Matrix<double, 3, 1> 
			RChP = inPRCh.col(i), 
			MChP = inPMCh.col(i),
			LChP = inPLCh.col(i);
		
		Mat3 rMatsCh_i = buildRefFramefrom3Pts(
			RChP, MChP, LChP);
		xyzEulersChF3Pts.col(i) = rMatsCh_i.eulerAngles(0, 1, 2); 
		for (size_t j = 0; j < nMarkers; ){
			efPtsInChest.block(j,i,3,1) = 
				rMatsCh_i.transpose() * efPts.block(j,i,3,1);

			j = j+3;
		}
	}

	// (a). inverse position of efPts in Chest frame
	DikProblem *dIkProbEfInChest =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			efPtsInChest);


	for (int i = 0; i < nPts; i++)
		ceres::SolveProblemAt(dIkProbEfInChest, i);

	ceres::SolveProblem(dIkProbEfInChest);

	Eigen::MatrixXd xyzEulersBaseEf = 
		MayaAnimation::getBaseEulerAngles (dIkProbEfInChest);
	printEigenMathematica(xyzEulersBaseEf.transpose()
		, cout, "xyzEulersBaseEf");
	printEigenMathematica(xyzEulersBaseEf.transpose()*(180.0/M_PI)
		, cout, "xyzEulersBaseEfInDeg");


/*
MatrixXd 
		inPtsRPi = efPts.block(Markers::RPi,0,3,nPts),
		inPtsRTh = efPts.block(Markers::RTh,0,3,nPts),
		inPtsRWr = efPts.block(Markers::RWr,0,3,nPts),
		inPtsRLA = efPts.block(Markers::RLA,0,3,nPts),
		inPtsREl = efPts.block(Markers::REl,0,3,nPts),
		inPtsRUA = efPts.block(Markers::RUA,0,3,nPts),
		inPtsRSh = efPts.block(Markers::RSh,0,3,nPts);

	//(a) ShoToUa
	Ellipse3D *e3dUA(new Ellipse3D(cartToHom(inPtsRUA)));

	MatrixXd thetasUA = e3dUA.findThetas();
	MatrixXd efPtsRUA = homToCart(
		e3dUA.getPointAtThetasH(fixThetas(thetasUA)));

	//(b) shoToEl
	Ellipse3D e3dEl(cartToHom(inPtsREl));
	MatrixXd thetasEl = e3dEl.findThetas();
	MatrixXd efPtsREl = homToCart(
		e3dEl.getPointAtThetasH(fixThetas(thetasEl)));
	//(c) shoToLa
	Ellipse3D e3dLA(cartToHom(inPtsRLA));
	MatrixXd thetasLA = e3dLA.findThetas();
	MatrixXd efPtsRLA = homToCart(
		e3dLA.getPointAtThetasH(fixThetas(thetasLA)));
	//(d) shoToWr
	Ellipse3D e3dWr(cartToHom(inPtsRWr));
	MatrixXd thetasWr = e3dWr.findThetas();
	MatrixXd efPtsRWr = homToCart(
		e3dWr.getPointAtThetasH(fixThetas(thetasWr)));
	//(e) shoToTh
	Ellipse3D e3dTh(cartToHom(inPtsRTh));
	MatrixXd thetasTh = e3dTh.findThetas();
	MatrixXd efPtsRTh = homToCart(
		e3dTh.getPointAtThetasH(fixThetas(thetasTh)));
	//(e) shoToPi
	Ellipse3D e3dPi(cartToHom(inPtsRPi));
	MatrixXd thetasPi = e3dPi.findThetas();
	MatrixXd efPtsRPi = homToCart(
		e3dPi.getPointAtThetasH(fixThetas(thetasPi)));

*/
	/*** END - WAM Kinematics solver***/
	/************************************************/


	return 0;
}
