#include <iostream>

#include <UBCUtil.h>
#include <LinearAlgebra.h>
#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "ForwardKin.h"
#include "DikProblem.h"
#include "DikSolver.h"
#include "MayaAnimation.h"

#include "ParseMathematica.h"
#include "ParseCSV.h"
#include <Eigen/SVD>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>       /* pow */
#include "FitFunctions.h"
#include "RigidTrans2D.h"

#include "BaysFitFunctions.h"
using namespace ceres;
using namespace Eigen;
using namespace std; 
class Ellipse3D *e3d;
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
//	cout << "(* " << fileName << "*)" << endl;

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

	bool use_quaternions = true;

	// input marker data to copy into Eigen Matrix 
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
	
	/*** (H). Human motion from vicon data**/
	/********************************************/

	

	// motion phase
	int phase = 
		whichPhase(fileName);		
		
	string phaseS = to_string(phase);
	// start and end frame of the motion
	int startF, endF;
	getStartFrame(fileName, startF, endF);
	
	// (H1). Read & copy input CSV data into Eigen
	// nx3
	std::string metaFileName;
	std::string reachCtr = getPartID(fileName) + "P"
		+fileName.at ((unsigned)(fileName.length()-5));

	// Calculating the length of string 
	metaFileName = 
		fileName.substr(0,fileName.find("_clean"))+".csv";	
	string headerLine = 
		getHeaderRowCSV(metaFileName);

	int nPts = countRowsCSV(fileName, wristMarkers[0]);
	MatrixXd inPtsAlongRows = 
		MatrixXd::Zero(nPts,nMarkers);
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


	if (phase == 1)	
		inPtsAlongRows = 
			inPtsAlongRowsTotal.topRows(nPts);
	
	else if (phase == 2)	{
		inPtsAlongRows = 
			inPtsAlongRowsTotal.bottomRows(nPts);
		startF = endF - nPts +1;
	}
	else {
		cout << 
			"Error! can't determine the motion Phase" 
			<< endl;
		}
	// positions in meter
	MatrixXd inPts = 0.01 
		*inPtsAlongRows.transpose();


	// (H1).END- Read & copy input CSV data into Eigen
	/*** (2) ELLIPSE fitted model of Human motion **/
	/***************************************************


	string fileNameH = 
		fileName.substr(0,(unsigned)(
			fileName.length()-7))+".csv";	

	// motion phase
	int phase = 
		whichPhase(fileNameH);		

	string phaseS = to_string(phase);


	// start and end frame of the motion
	int startF, endF;
	getStartFrame(fileNameH, startF, endF);


	string wholeFile = 
		getWholeFileMathematica(fileName);


	MatrixXd inPtsAlongRows = 100.0
		* parseMathematica(wholeFile);
	int nPts = inPtsAlongRows.rows();

	// eliminate first pt - 
	//   due to incorrect ef solution
	if (phase == 1)	
		startF = startF + 1;
 	else if (phase == 2)	
		startF = endF - nPts +1;

	// add EF to the namespace 
	phaseS = "EF"+phaseS;

	// positions in meter
	MatrixXd inPts = 0.01 
		*inPtsAlongRows.transpose();
*/
/***************************************************/
/**
"RWristX",	"RWristY",	"RWristZ" 			0
	,	"RThumbX",	"RThumbY",	"RThumbZ" 	3
	, "RLArmX",	"RLArmY",	"RLArmZ"				6
	, "RElbX",	"RElbY",	"RElbZ"				  9	
	. "RUArmX",	"RUArmY",	"RUArmZ"		   12
	, "RShoX",	"RShoY",	"RShoZ"				 15
	, "RChestX",	"RChestY",	"RChestZ"	 18
	, "MChestX",	"MChestY",	"MChestZ"  21
	, "LChestX", "LChestY", "LChestZ"    24
	, "LShoX", "LShoY", "LShoZ" 				 27
	, "RPinkieX",	"RPinkieY",	"RPinkieZ" 30
	};
**/

	// (H2).  Animation
	// 			- key vicon marker Positions in MB

	// inpts in MotionBuilder frame
	MatrixXd inPtsAlongRowsMB = 
		MotionBuilderAnimation
			::getInPtsInMbFrame(inPtsAlongRows);
	MatrixXd inPtsMB = inPtsAlongRowsMB.transpose();
	// add key to inPts marker positions
//	cout <<" Start inPts in MotionBuilder frame"<< endl;
/*		*/
	int mkr = 30;
	std::string marker_i = wSMarkers[mkr].substr(
		0, ((unsigned)(wSMarkers[mkr].length())));
	//EDITED to eliminate first pt
	MatrixXd inPtsMkr_i = 
		inPtsMB.block(mkr,1,3,nPts-1);
	
	if (phase == 1)	
		cout << "# Frames " << startF-1 << "-" << endF << endl;
		
//	MotionBuilderAnimation::keyAddSingleMarker(
//		inPtsMkr_i, startF, marker_i);
//		cout << "\n" << endl;
	
	MotionBuilderAnimation::keyAdd(
		inPtsMB, startF, wSMarkers);

//	cout <<" END inPts in MotionBuilder frame"<< endl;
	
/*

	// (H3). Print Human motion from vicon (inPts)
	//				& Scaled human motion to WAM (inPtsWam)
	DikProblem *dIkProb =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			inPts, phaseS);

 	// (H3.1). Get human arm joint ORIENTATIONS 
 	//					from vicon inPts in VICON Frame - 
 	//							Rj_i
 	

	// chest orientation
	// in Vicon
	Eigen::MatrixXd rMatsCh = MayaAnimation::
		buildChRotReferentials(dIkProb);
	Eigen::MatrixXd xyzEulersCh = 
		MayaAnimation::getChEulerAngles (dIkProb);

	printEigenMathematica(rMatsCh.transpose()
		, cout, "rMatsCh" + phaseS);
	printEigenMathematica(xyzEulersCh.transpose()
		, cout, "xyzEulersCh" + phaseS);

	// in MAYA ChestFrame_0
	Eigen::MatrixXd rMatsChInChF = MayaAnimation::
		buildChRotReferentialsInChF(dIkProb);
	Eigen::MatrixXd xyzEulersChInChF = 
		MayaAnimation::getChEulerAnglesInChF (dIkProb);

	printEigenMathematica(rMatsChInChF.transpose()
		, cout, "rMatsChInChF" + phaseS);
	printEigenMathematica(xyzEulersChInChF.transpose()
		, cout, "xyzEulersChInChF" + phaseS);
*/


/*
		// elbow
		Eigen::MatrixXd xyzEulersEl = 
			MayaAnimation::getElEulerAngles (dIkProb);
		printEigenMathematica(xyzEulersEl.transpose()
			, cout, "xyzEulersEl");
		printEigenMathematica(xyzEulersEl.transpose()*(180.0/M_PI)
			, cout, "xyzEulersElInDeg");

		//wrist
		Eigen::MatrixXd xyzEulersWr = MayaAnimation::
			getWrEulerAngles (dIkProb);
		printEigenMathematica(xyzEulersWr.transpose()
			, cout, "xyzEulersWr");
		printEigenMathematica(xyzEulersWr.transpose()*(180.0/M_PI)
			, cout, "xyzEulersWrInDeg");
	*/

/*
	// (H4). Transform (inPts) &  (inPtsWam) to CHEST Frame
	//				& Print output		
	MatrixXd
		inPtsInChest = MatrixXd::Zero(nMarkers,nPts);

	for (int i = 0; i < nPts; i++) {	
		//Rchest_in_Vicon
		Mat3 rMatsCh_i = rMatsCh.block(0,i*3,3,3);
		//inPts in Chest frame
		for (size_t j = 0; j < nMarkers; ){	
			inPtsInChest.block(j,i,3,1) = 
				rMatsCh_i.transpose() * inPts.block(j,i,3,1);

			j = j+3;
		}
	}

	DikProblem *dIkProbInChest =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			inPtsInChest, phaseS+"inCh");

*/

	/*** (IK). WAM InverseKinematics from Human motion **/
	/******************************************************

	// (IK1). inverse position in chest frame


	// A. In Human Chest Frame
	cout << "(*IK solution in Human Chest Frame*)" << endl;
	
	cout << "(*(1) SolveProblemAt*)" << endl;
	for (int i = 1; i < nPts; i++)
		ceres::SolveProblemAt(dIkProbInChest, i);
	computeFitErrPtsR(dIkProbInChest);
	Eigen::MatrixXd 
		outThetasWamInChest = 
			ceres::getOutThetasWam(dIkProbInChest),
		outQuatsBaseInChest = 
			ceres::getOutQuatsBase(dIkProbInChest),
		xyzEulersBaseInChest = MayaAnimation::
			getBaseEulerAngles (dIkProbInChest);
	printEigenMathematica(outThetasWamInChest, cout
		, "outThetasWam"+ phaseS + "inCh" + "SolveProbAt");
	printEigenMathematica(outQuatsBaseInChest, cout
		, "outQuatsBase"+ phaseS + "inCh" + "SolveProbAt");
	printEigenMathematica(xyzEulersBaseInChest
		.transpose(), cout, "xyzEulersBase"+ phaseS
		 + "inCh" + "SolveProbAt");
		
	cout << "(*(2) SolveProblem*)" << endl;
	ceres::SolveProblem(dIkProbInChest);
	//	computeFitErrPtsR(dIkProbInChest);
	outThetasWamInChest = 
		ceres::getOutThetasWam(dIkProbInChest);
	outQuatsBaseInChest = 
		ceres::getOutQuatsBase(dIkProbInChest);
	xyzEulersBaseInChest = MayaAnimation::
		getBaseEulerAngles (dIkProbInChest);

	printEigenMathematica(outThetasWamInChest, cout
		, "outThetasWam"+ phaseS + "inCh" + "SolveProb");
	printEigenMathematica(outQuatsBaseInChest, cout
		, "outQuatsBase"+ phaseS + "inCh" + "SolveProb");
	printEigenMathematica(xyzEulersBaseInChest.transpose()
		, cout, "xyzEulersBase"+ phaseS+"inCh"+"SolveProb");


	cout << "(*(3) SolveProblemFPrevPt*)" << endl;
	ceres::SolveProblemFPrevPt(dIkProbInChest);
	computeFitErrPtsR(dIkProbInChest);
//	ceres::SolveProblemFNextPt(dIkProbInChest);
	outThetasWamInChest = 
		ceres::getOutThetasWam(dIkProbInChest);
	outQuatsBaseInChest = 
		ceres::getOutQuatsBase(dIkProbInChest);
	xyzEulersBaseInChest = MayaAnimation::
		getBaseEulerAngles (dIkProbInChest);
	printEigenMathematica(outThetasWamInChest
		, cout, "outThetasWam"+ phaseS+"inCh" 
		+	"SolveProbFPrevPt");
	printEigenMathematica(outQuatsBaseInChest
		, cout, "outQuatsBase"+ phaseS+"inCh"
		+	"SolveProbFPrevPt");
	printEigenMathematica(xyzEulersBaseInChest.transpose()
		, cout, "xyzEulersBaseInChest"+ phaseS+"inCh"
		+	"SolveProbFPrevPt");
		

	cout << "(*(4) SolveProblem2*)" << endl;
	ceres::SolveProblem(dIkProbInChest);
	//	computeFitErrPtsR(dIkProbInChest);
	outThetasWamInChest = 
		ceres::getOutThetasWam(dIkProbInChest);
	outQuatsBaseInChest = 
		ceres::getOutQuatsBase(dIkProbInChest);
	printEigenMathematica(outThetasWamInChest, cout
		, "outThetasWam"+ phaseS + "inCh" + "SolveProb2");
	printEigenMathematica(outQuatsBaseInChest, cout
		, "outQuatsBase"+ phaseS + "inCh" + "SolveProb2");
		

	Eigen::MatrixXd  rMatsBase = MayaAnimation::
		buildBaseRotReferentials(dIkProbInChest);
	printEigenMathematica(rMatsBase.transpose()
		, cout, "rMatsBase"+ phaseS+"inCh");			
	xyzEulersBaseInChest = MayaAnimation:: 
		getBaseEulerAngles (dIkProbInChest);
	printEigenMathematica(xyzEulersBaseInChest.transpose()
		, cout, "xyzEulersBase"+ phaseS+"inCh");
	printEigenMathematica(xyzEulersBaseInChest.transpose()
		*(180.0/M_PI), cout, "xyzEulersBaseInDeg"
			+ phaseS+"inCh");

	ceres::printOutPtsWam(dIkProbInChest,rMatsCh);

**/
/*
	// B. In Vicon Frame
	cout << "solution in Vicon Frame" << endl;
	for (int i = 1; i < nPts; i++)
		ceres::SolveProblemAt(dIkProb, i);

	ceres::SolveProblem(dIkProb);
	//gradually Error-> 0
	ceres::SolveProblemFPrevPt(dIkProb);


	cout  << "END SolveProblemFPrevPt" << endl;

	computeFitErrPtsR(dIkProb);


	//Error starts at 0, but is increasing
	ceres::SolveProblemFNextPt(dIkProb);

	computeFitErrPtsR(dIkProb);

	Eigen::MatrixXd xyzEulersBase = 
		MayaAnimation::getBaseEulerAngles (dIkProb);
	printEigenMathematica(xyzEulersBase.transpose()
		, cout, "xyzEulersBase"+ phaseS);
	printEigenMathematica(xyzEulersBase.transpose()*(180.0/M_PI)
		, cout, "xyzEulersBaseInDeg"+ phaseS);

	Eigen::MatrixXd outThetasWam = 
		ceres::getOutThetasWam(dIkProb);
	printEigenMathematica(outThetasWam
		, cout, "outThetasWam"+ phaseS);
*/
	// (IK2). inverse orientation
//	for (int i = 0; i < nPts; i++)
//		ceres::SolveOrientationProblemAt(dIkProbInChest, i);

//	ceres::SolveOrientationProblem(dIkProbInChest);



/*
	// (IK3).Maya- key WAM joint angles
	
	// spine/chest rotation & translation
	Eigen::MatrixXd inPtsChMean = MatrixXd::Zero(3,nPts);
	for (int i = 0; i < nPts; i++)
		inPtsChMean.col(i) = (dIkProbInChest->RChP(i) 
			+ dIkProbInChest->MChP(i) + dIkProbInChest->LChP(i))/3.0;

 	MayaAnimation::TranskeyAddHjoint(inPtsChMean
		, "Spine", startF);
 	MayaAnimation::RotkeyAddHjoint(xyzEulersChInChF
		, "Spine", startF);

	// A. In Human Chest Frame

//	cout << "Animation of IK solution in Human Chest Frame:" 
//			 << "\n\n\n" << endl;
	MayaAnimation::RotkeyAddHjoint(xyzEulersBaseInChest
		, "Base", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWamInChest.row(0)
			, "Shoulder_Yaw_J1", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWamInChest.row(1)
			, "Shoulder_Pitch_J2", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWamInChest.row(2)
			, "Shoulder_UpperArm_J3", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWamInChest.row(3)
			, "Elbow_ForeArm_J4", startF);
*/
			
			
/*

	cout << "Animation of IK solution in Vicon Frame:" 
			<< "\n\n\n" << endl;

	// B. In Vicon Frame
	MayaAnimation::RotkeyAddHjoint(xyzEulersBase
		, "Base", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWam.row(0)
			, "Shoulder_Yaw_J1", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWam.row(1)
			, "Shoulder_Pitch_J2", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWam.row(2)
			, "Shoulder_UpperArm_J3", startF);

	MayaAnimation::RotkeyAddRjoint(outThetasWam.row(3)
			, "Elbow_ForeArm_J4", startF);
*/

/*
	MayaAnimation::RotkeyAddRjoint(outThetasWam.row(4)
			, "Wrist_Yaw_J5", startF);
	MayaAnimation::RotkeyAddRjoint(outThetasWam.row(5)
			, "Wrist_Pitch_J6", startF);
	MayaAnimation::RotkeyAddRjoint(outThetasWam.row(6)
			, "Wrist_Palm_J7", startF);

*/
		
		





	/*** END - WAM Kinematics solver***/
	/************************************************/


	return 0;
}
