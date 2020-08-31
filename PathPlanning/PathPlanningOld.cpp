#include <iostream>
#include <UBCUtil.h>
#include <Plane.h>
#include <LinearAlgebra.h>
#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "ForwardKin.h"
#include "BasePoseSolver.h"
#include "WamIkProblem.h"
#include "WamIkPoseSolver.h"
#include "IKfitFromLaToW.h"
#include "IKperInputPt.h"
#include "NaiveIK.h"
#include "ParseMathematica.h"
#include "ParseCSV.h"
#include <Eigen/SVD>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>       /* pow */

using namespace ceres;
using namespace Eigen;
using namespace std;

int main(int argc, char **argv) {
	
	string fileName = argv[1];
	cout << "(* " << fileName << endl;

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
	int nRows = countRowsCSV(fileName, wristMarkers[0]);
	MatrixXd inPtsAlongRows = MatrixXd::Zero(nRows,nMarkers);
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
	MatrixXd inPtsAlongRowsMB = MatrixXd::Zero(nRows,nMarkers);
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
	/*** WAM Kinematics solver ***/
	// joint limits
	std::vector<double> jointMinAngles, jointMaxAngles;
	jointMinAngles.push_back(-2.6);
	jointMinAngles.push_back(-2.0);
	jointMinAngles.push_back(-2.8);
	jointMinAngles.push_back(-0.9);
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

	/*** END - WAM Kinematics solver***/
	/************************************************/


/*** Wam IK SOLVER ---From 3d Position to joint space valuse  
(1) CERE solve from inPts to startConfigs:
			startConfig_j = [a0, a1,a2,a3]_j j = {0,1,...,6}
(2) Calculate thetas from startConfigs and t:  
			t_i = frame rate * #frames
			thetas_i = a0 + a1*t + a2*t^2 + a3*t^3...
***/
	// number of optimization parameters 
	size_t toDof = 7;
	size_t nJoints = 7;
	size_t polynomialOrder = 3;
	size_t nParams = nJoints * (polynomialOrder + 1);

	/*** Wam IK SOLVER 
		(1) CERES solve from inPts to startConfigs:
			startConfig_j = [a0, a1,a2,a3]_j j = {0,1,...,6}
	  given marker inPts, solve for IK of WAM  **/
	/*********************************************************/
	// positions in meter 
	inPtsAlongRows = 0.01 * inPtsAlongRows; 
	MatrixXd JsWam = MatrixXd::Zero(nRows,nJoints); 

	WamIkProblem<double> *wamIkProb =
		new WamIkProblem<double>(wam, JsWam, inPtsAlongRows		
			, jointMinAngles, jointMaxAngles);
	ceres::SolveProblem(wamIkProb);

/*
cout << "θ_4 = " << endl;
	for (int i = 0; i < nRows; i++)
		cout << wamIkProb->computeJ4FromHumanElbAngle(i) << endl;

	ceres::solveJ4FromHumanElbAngle(wamIkProb);
*/

	//** For print in Mathematica *//
	printEigenMathematica( JsWam, cout
		, "JsWam"+ reachCtr);	
/*

	/// inPts ///
	//(a) for human
	printEigenMathematica( inPtsRW.transpose(), cout
		, "inPtsRW"+ reachCtr);	
	printEigenMathematica( inPtsRTh.transpose(), cout
		, "inPtsRTh"+ reachCtr);	
	printEigenMathematica( inPtsRPk.transpose(), cout
		, "inPtsRPk"+ reachCtr);	
	printEigenMathematica( inPtsRLA.transpose(), cout
		, "inPtsRLA"+ reachCtr);	
	printEigenMathematica( inPtsREl.transpose(), cout
		, "inPtsREl"+ reachCtr);	
	printEigenMathematica( inPtsRUA.transpose(), cout
		, "inPtsRUA"+ reachCtr);	
	printEigenMathematica( inPtsRSh.transpose(), cout
		, "inPtsRSh"+ reachCtr);	
	printEigenMathematica( inPtsRCh.transpose(), cout
		, "inPtsRCh"+ reachCtr);	
	printEigenMathematica( inPtsMCh.transpose(), cout
		, "inPtsMCh"+ reachCtr);	
	printEigenMathematica( inPtsLCh.transpose(), cout
		, "inPtsLCh"+ reachCtr);	
	/// (b) for wam, i.e. scaled to wam
		//(b.1) inPts
	printEigenMathematica( inPtsRWrWam.transpose(), cout
		, "inPtsRWrWam"+ reachCtr);	
	printEigenMathematica( inPtsRLaWam.transpose(), cout
		, "inPtsRLaWam"+ reachCtr);	
	printEigenMathematica( inPtsRElWam.transpose(), cout
		, "inPtsRElWam"+ reachCtr);	
	printEigenMathematica( inPtsRUaWam.transpose(), cout
		, "inPtsRUaWam"+ reachCtr);	
	printEigenMathematica( inPtsRShWam.transpose(), cout
		, "inPtsRShWam"+ reachCtr);	
		//(b.2) inVectors
	printEigenMathematica( shoToElVsH.transpose() * shoToElLwam
		, cout, "inUaVsWam"+ reachCtr);	
	printEigenMathematica( elOffsetVsWam.transpose(), cout
		, "inElOffsetVsWam"+ reachCtr);	
	printEigenMathematica( shoToElVsWam.transpose(), cout
		, "inShoToElVsWam"+ reachCtr);	
	printEigenMathematica( wOffsetVsWam.transpose(), cout
		, "inWrOffsetVsWam"+ reachCtr);	
	printEigenMathematica( elToWrVsH.transpose() * elToWrLwam
		, cout, "inLaVsWam"+ reachCtr);	
	printEigenMathematica( elToWrVsWam.transpose(), cout
		, "inElToWrVsWam"+ reachCtr);	

	/// (1.1) Find T_fromCamera^ToBase base reference sys 
	///		via computing the bfp of the chest ///
	//(a) for human
	printEigenMathematica( inPtsRTHinBaseH.transpose(), cout
		, "inPtsRTHinBaseH"+ reachCtr);	
	printEigenMathematica( inPtsRWinBaseH.transpose(), cout
		, "inPtsRWinBaseH"+ reachCtr);	
	printEigenMathematica( inPtsRTHinBaseH.transpose(), cout
		, "inPtsRTHinBaseH"+ reachCtr);	
	printEigenMathematica( inPtsRLAinBaseH.transpose(), cout
		, "inPtsRLAinBaseH"+ reachCtr);	
	printEigenMathematica( inPtsRELinBaseH.transpose(), cout
		, "inPtsRELinBaseH"+ reachCtr);	
	printEigenMathematica( inPtsRUAinBaseH.transpose(), cout
		, "inPtsRUAinBaseH"+ reachCtr);	
	printEigenMathematica( inPtsRSHinBaseH.transpose(), cout
		, "inPtsRSHinBaseH"+ reachCtr);	
	//(b) for robot
	//(b.1) inPts(i.e. SAME AS VECTORS)
	printEigenMathematica( inPtsRTHinBaseWam.transpose(), cout
		, "inPtsRTHinBaseWam"+ reachCtr);	
	printEigenMathematica( inPtsRWinBaseWam.transpose(), cout
		, "inPtsRWinBaseWam"+ reachCtr);	
	printEigenMathematica( inPtsRLAinBaseWam.transpose(), cout
		, "inPtsRLAinBaseWam"+ reachCtr);	
	printEigenMathematica( inPtsRELinBaseWam.transpose(), cout
		, "inPtsRELinBaseWam"+ reachCtr);	
	printEigenMathematica( inPtsRUAinBaseWam.transpose(), cout
		, "inPtsRUAinBaseWam"+ reachCtr);	
	//(b.2) inVectors OFFSETS WAM in base 
	printEigenMathematica( elOffsetVsWaminBaseWam.transpose()
		, cout, "inElOffsetVsWamInBaseWam"+ reachCtr);	
	printEigenMathematica( wOffsetVsWaminBaseWam.transpose()
		, cout, "inWrOffsetVsWamInBaseWam"+ reachCtr);	

	/// outPts ///
	// in Base //
		// (1) output Pts in BASE (i.e. SAME AS VECTORS)
	printEigenMathematica( outPtsRUAinBaseFj2.transpose()
		, cout, "outPtsRUAinBaseFj2" + reachCtr);	
	printEigenMathematica( outPtsRELinBaseFj3.transpose()
		, cout, "outPtsRELinBaseFj3" + reachCtr);	
	printEigenMathematica( outPtsRLAinBaseFj4.transpose()
		, cout, "outPtsRLAinBaseFj4" + reachCtr);	
	printEigenMathematica( outPtsRWinBaseFj5.transpose()
		, cout, "outPtsRWinBaseFj5" + reachCtr);	

		// (2) outVectors OFFSETS WAM in base 
	printEigenMathematica(outElOffsetVsWamInBase.transpose()
		, cout, "outElOffsetVsWamInBaseWam" + reachCtr);	
	printEigenMathematica( outWrOffsetVsWamInBase.transpose()
		, cout, "outWrOffsetVsWamInBase" + reachCtr);	

	// in CAMERA //

		// (1) output Pts in CAMERA
	printEigenMathematica( outPtsRElWam.transpose(), cout
		, "outPtsRElWam"+ reachCtr);					
	printEigenMathematica( outPtsRLaWam.transpose(), cout
		, "outPtsRLaWam"+ reachCtr);		
	printEigenMathematica( outPtsRWrWam.transpose(), cout
		, "outPtsRWrWam"+ reachCtr);	
		// (2) output Vecotors in CAMERA
	printEigenMathematica( outUaVsWam.transpose()
		, cout, "outUaVsWam"+ reachCtr);	
	printEigenMathematica( outElOffsetVsWam.transpose()
		, cout, "outElOffsetVsWam"+ reachCtr);	
	printEigenMathematica( outShoToElVsWam.transpose()
		, cout, "outShoToElVsWam"+ reachCtr);	
	printEigenMathematica( outWrOffsetVsWam.transpose()
		, cout, "outWrOffsetVsWam"+ reachCtr);	
	printEigenMathematica( outLaVsWam.transpose() * elToWrLwam
		, cout, "outLaVsWam"+ reachCtr);	
	printEigenMathematica( outElToWrVsWam.transpose()
		, cout, "outElToWrVsWam"+ reachCtr);	

	printEigenMathematica( JsWam, cout
		, "JsWam"+ reachCtr);	


	//fit error 
	printEigenMathematica( fitErrUAinBaseFj2.transpose()
		, cout, "fitErrUAinBaseFj2"+ reachCtr);
	printEigenMathematica( fitErrElinBaseFj3.transpose()
		, cout, "fitErrElinBaseFj3"+ reachCtr);
	printEigenMathematica( fitErrLAinBaseFj4.transpose()
		, cout, "fitErrLAinBaseFj4"+ reachCtr);
	printEigenMathematica( fitErrWinBaseFj5.transpose()
		, cout, "fitErrWinBaseFj5"+ reachCtr);
*/
	return 0;
}
	/**  TEST (1.2.A) - transform back SCALED pts in base
				to CAMERA frame    
		i.e. T * inPtsInBaseScaled = inPtsinBase 
	Eigen::MatrixXd inPtsRWScaled = MatrixXd::Zero(3,nRows)
		, inPtsRLaScaled = MatrixXd::Zero(3,nRows)
		, inPtsRElScaled = MatrixXd::Zero(3,nRows)
		, inPtsRUaScaled = MatrixXd::Zero(3,nRows)
		, inPtsRShScaled = MatrixXd::Zero(3,nRows);


		for(size_t i = 0; i < nRows; i++) {
			//Transform BACK from base to CAMERA frame
			Vector3d Rsho = inPtsRSh.col(i), RCh = inPtsRCh.col(i)
				, MCh = inPtsMCh.col(i), LCh = inPtsLCh.col(i);
			MatrixXd TfromBaseToCamera = MatrixXd::Identity(4,4);
			TfromBaseToCamera.block(0,0,3,3) = 
				buildRefFramefrom3Pts(RCh, MCh, LCh);
			//sho marker as the origin of the WAM
			TfromBaseToCamera.block(0,3,3,1) = Rsho;
			
			inPtsRWScaled.col(i) = homToCart( 
				TfromBaseToCamera * cartToHom(inPtsRWinBaseScaled.col(i)));
			inPtsRLaScaled.col(i) = homToCart( 
				TfromBaseToCamera * cartToHom(inPtsRLAinBaseScaled.col(i)));
			inPtsRElScaled.col(i) = homToCart( 
				TfromBaseToCamera * cartToHom(inPtsRELinBaseScaled.col(i)));
			inPtsRUaScaled.col(i) = homToCart( 
				TfromBaseToCamera * cartToHom(inPtsRUAinBaseScaled.col(i)));
			inPtsRShScaled.col(i) = homToCart( 
				TfromBaseToCamera * cartToHom(inPtsRSHinBaseScaled.col(i)));
		}
	**/

