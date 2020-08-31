#include <iostream>
#include <UBCUtil.h>
#include <Plane.h>
#include <LinearAlgebra.h>
#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "ForwardKin.h"
#include "DikProblem.h"
#include "DikSolver.h"


#include "ParseMathematica.h"
#include "ParseCSV.h"
#include <Eigen/SVD>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>       /* pow */
#include "FitFunctions.h"
#include "Ellipse3D.h"
#include "Plane.h"
#include "LinearAlgebra.cpp"
#include "RigidTrans2D.h"

#include "BaysFitFunctions.h"
#include "EllipseConicConstraints.h"
using namespace ceres;
using namespace Eigen;
using namespace std;

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
	int nPts = countRowsCSV(fileName, wristMarkers[0]);
	MatrixXd inPtsAlongRows = MatrixXd::Zero(nPts,nMarkers);
	std::string metaFileName;
	std::string reachCtr = getPartID(fileName) + "P"
		+fileName.at ((unsigned)(fileName.length()-5));
	// Calculating the length of string 
	metaFileName = 
		fileName.substr(0,fileName.find("_clean"))+".csv";	
	string headerLine = 
		getHeaderRowCSV(metaFileName);

	int ctr = 0;
	while ((ctr < nMarkers)) {
		int startMarkerDelim = 
			headerLine.find(wSMarkers[ctr]);
		if (startMarkerDelim!=std::string::npos) { 
			fillMatCSV(ctr, metaFileName
				, wSMarkers[ctr], wSMarkers[ctr+1]
				, wSMarkers[ctr+2], inPtsAlongRows);
		}
		ctr +=3;
	}

	/** TO IMPORT TO MOTION BUILDER 
			X' = -X; Y' = Z; Z'=Y;
	**/
	MatrixXd inPtsAlongRowsMB = MatrixXd::Zero(nPts,nMarkers);
	for (size_t j = 0; j < nMarkers/3; j++){
		size_t jTimesThree = j * 3;	
		inPtsAlongRowsMB.col(jTimesThree) 
			=	-inPtsAlongRows.col(jTimesThree);
		inPtsAlongRowsMB.col(jTimesThree+1) 
			= inPtsAlongRows.col(jTimesThree+2);
		inPtsAlongRowsMB.col(jTimesThree+2) 
			= inPtsAlongRows.col(jTimesThree+1);
	}
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

	/*** (1) Human motion from vicon data**/
	/*********************************************************/
	MatrixXd JsWam = MatrixXd::Zero(nJoints,  nPts); 
	// (a). inverse position
	bool use_quaternions = true;
	DikProblem *dIkProb =
		new DikProblem	(wam, jointMinAngles,
			jointMaxAngles, use_quaternions, 
			inPts);


	for (int i = 0; i < nPts; i++)
		ceres::SolveProblemAt(dIkProb, i);

	ceres::SolveProblem(dIkProb);

	// (b). inverse orientation

	/*** (2) ELLIPSE fitted model of Human motion **/
	/*********************************************************/
	MatrixXd 
		inPtsRPi = inPts.block(Markers::RPi,0,3,nPts),
		inPtsRTh = inPts.block(Markers::RTh,0,3,nPts),
		inPtsRWr = inPts.block(Markers::RWr,0,3,nPts),
		inPtsRLA = inPts.block(Markers::RLA,0,3,nPts),
		inPtsREl = inPts.block(Markers::REl,0,3,nPts),
		inPtsRUA = inPts.block(Markers::RUA,0,3,nPts),
		inPtsRSh = inPts.block(Markers::RSh,0,3,nPts);
/*
	//(a) ShoToUa
	Ellipse3D e3dUA(cartToHom(inPtsRUA));

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
