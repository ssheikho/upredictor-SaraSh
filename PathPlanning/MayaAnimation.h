#ifndef MAYA_ANIMATION_H
#define MAYA_ANIMATION_H

#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "DifferentialIKErrorTerm.h"
#include "UBCUtil.h"
#include "DikProblem.h"
#include "DikSolver.h"

#include "ceres/ceres.h"
using std::vector;

using namespace Eigen;
using namespace std;
typedef Eigen::Matrix<double, 3, 3> Mat3;
typedef Eigen::Matrix<double, 3, 1> Vec3;
typedef Eigen::Vector4d Vec4;


namespace MotionBuilderAnimation {

	/** Import Vicon Markers to MOTION BUILDER 
			X' = -X; Y' = Z; Z'=Y;
	**/
	MatrixXd getInPtsInMbFrame(MatrixXd inPtsAlongRows)	{
		double nPts = inPtsAlongRows.rows();
		double nMarkers = inPtsAlongRows.cols();

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
	return inPtsAlongRowsMB;
	}
	
	void keyAdd(MatrixXd inPtsMB, int startFrame
		, std::vector<std::string> markers)	{
		int nMarkers = inPtsMB.rows();
		int nPts = inPtsMB.cols();
		int mkr = 0;
		int mkrLabel = 0;
		while ((mkr < nMarkers)) {
		cout << "if marker.Name == marker_labels[" 
			<< mkrLabel <<"]:" << endl;
			cout << "	Translation = marker.PropertyList.Find("
			 	<< "'Lcl Translation' )" << endl;
			cout << "	xTrans, yTrans, zTrans"
				<<" = marker.Translation" << endl;
		 	cout << "	marker.Translation.SetAnimated(True)" << endl;
		 	
			MatrixXd inPtsMkr_i = 
				inPtsMB.block(mkr,0,3,nPts);
			cout << "		#" << markers[mkr].substr(0, ((
				unsigned)(markers[mkr].length()-1))) << endl;
			for (size_t i = 0; i < nPts; i++){
				cout << "	marker.Translation.GetAnimationNode()"	
					<< ".KeyAdd(FBTime(0, 0, 0,"<< startFrame + i 
					<< ",0), [" << inPtsMkr_i(0,i) << "," 
					<<  inPtsMkr_i(1,i)  << "," << inPtsMkr_i(2,i)
					<< "])" << endl;
			}
			mkrLabel++;
			mkr += 3;
		}
	}

	void keyAddSingleMarker(MatrixXd inPtsMkr_i, int startFrame
		, std::string marker_i)	{
		int nPts = inPtsMkr_i.cols();
	
		cout << "#" << marker_i.substr(0, ((
			unsigned)(marker_i.length()-1))) << endl;
		for (size_t i = 0; i < nPts; i++){
			cout << "marker.Translation.GetAnimationNode()"	
				<< ".KeyAdd(FBTime(0, 0, 0,"<< startFrame + i 
				<< ",0), [" << inPtsMkr_i(0,i) << "," 
				<<  inPtsMkr_i(1,i)  << "," << inPtsMkr_i(2,i)
				<< "])" << endl;
		}
	}




}

namespace MayaAnimation {


	
	void RotkeyAddHjoint(Eigen::MatrixXd xyzEul, std::string joint
		, int startFrame)	{

		int nPts = xyzEul.cols();
		//rotAxis 0->x, 1->y, 2->z

		cout << "$label = (" << "\"" << joint 
				 << "\"" << ");" << endl;

		for (int j = 0; j < nPts; j++)
			cout << "setKeyframe -attribute \"rotateZ\" -t "
					 << (j+startFrame)/100.0 << "sec -value " 
					 << xyzEul(2,j) << " $label;" << endl;
				
		for (int j = 0; j < nPts; j++)
				cout << "setKeyframe -attribute \"rotateY\" -t "
						 << (j+startFrame)/100.0 << "sec -value " 
						 << xyzEul(1,j) << " $label;" << endl;


		for (int j = 0; j < nPts; j++)
				cout << "setKeyframe -attribute \"rotateX\" -t "
						 << (j+startFrame)/100.0 << "sec -value " 
						 << xyzEul(0,j) << " $label;" << endl;
		}




	
	void TranskeyAddHjoint(MatrixXd inPts, 
		std::string joint, int startFrame)	{

		int nPts = inPts.cols();
		//rotAxis 0->x, 1->y, 2->z

		cout << "$label = (" << "\"" << joint 
				 << "\"" << ");" << endl;

		for (int j = 0; j < nPts; j++)
			cout << "setKeyframe -attribute \"translateZ\" -t "
					 << (j+startFrame)/100.0 << "sec -value " 
					 << inPts(2,j) << " $label;" << endl;
		

		for (int j = 0; j < nPts; j++)
				cout << "setKeyframe -attribute \"translateX\" -t "
						 << (j+startFrame)/100.0 << "sec -value " 
						 << inPts(0,j) << " $label;" << endl;
			
		for (int j = 0; j < nPts; j++)
				cout << "setKeyframe -attribute \"translateY\" -t "
						 << (j+startFrame)/100.0 << "sec -value " 
						 << inPts(1,j) << " $label;" << endl;
	}



	
	void RotkeyAddRjoint(Eigen::MatrixXd outTheta_i
		, std::string joint, int startFrame)	{

		int nPts = outTheta_i.cols();
		//rotAxis 0->x, 1->y, 2->z

		cout << "$label = (" << "\"" << joint 
				 << "\"" << ");" << endl;

		for (int j = 0; j < nPts; j++)
			cout << "setKeyframe -attribute \"rotateZ\" -t "
					 << (j+startFrame)/100.0 << "sec -value " 
					 << outTheta_i(0,j) << " $label;" << endl;
	}

	void RotkeyAddRjoints(Eigen::MatrixXd outThetasWam
		, int startF)	{
		RotkeyAddRjoint(outThetasWam.row(0)
			, "Shoulder_Yaw_J1", startF);

		RotkeyAddRjoint(outThetasWam.row(1)
			, "Shoulder_Pitch_J2", startF);

		RotkeyAddRjoint(outThetasWam.row(2)
			, "Shoulder_UpperArm_J3", startF);

		RotkeyAddRjoint(outThetasWam.row(3)
			, "Elbow_ForeArm_J4", startF);
/*
		RotkeyAddRjoint(outThetasWam.row(4)
			, "Wrist_Yaw_J5", startF);
		RotkeyAddRjoint(outThetasWam.row(5)
			, "Wrist_Pitch_J6", startF);
		RotkeyAddRjoint(outThetasWam.row(6)
			, "Wrist_Palm_J7", startF);
*/
	}


	Eigen::MatrixXd MapEllipticalMrksToSkeleton (
		DikProblem* dikProb, MatrixXd efPts) {
		
		const int nObservations = dikProb->nObservations();
		int nCols = efPts.cols();
		MatrixXd retSkeletonPts = 
			MatrixXd::Zero(nObservations,nCols);

	
	}

// Checks if a matrix is a valid rotation matrix.
bool isRotationMatrix(Mat3 &R)
{
    Mat3 Rt = R.transpose();
    Mat3 shouldBeIdentity = Rt * R;
    Mat3 I = Mat3::Identity(3,3);
    
    return  (I-shouldBeIdentity).norm() < 1e-6;
    
}

// Calculates rotation matrix to euler angles
// The result is the same as MATLAB except the order
// of the euler angles ( x and z are swapped ).

/*
Note that the inverse sine and cosine functions yield two possible values for the argument. In this geometrical description, only one of the solutions is valid. When Euler angles are defined as a sequence of rotations, all the solutions can be valid, but there will be only one inside the angle ranges. This is because the sequence of rotations to reach the target frame is not unique if the ranges are not previously defined.[2]

For computational purposes, it may be useful to represent the angles using atan2(y, x). For example, in the case of proper Euler angles:

    α = atan2 ⁡ ( Z 1 , − Z 2 ) , {\displaystyle \alpha =\operatorname {atan2} (Z_{1},-Z_{2}),} \alpha = \operatorname{atan2}(Z_1 , -Z_2),
    γ = atan2 ⁡ ( X 3 , Y 3 ) . {\displaystyle \gamma =\operatorname {atan2} (X_{3},Y_{3}).} \gamma =\operatorname {atan2} (X_{3},Y_{3}).
*/
Vec3 rotationMatrixToEulerAnglesC(Mat3 &R) {

    assert(isRotationMatrix(R));
    
    float sy = sqrt(R(0,0) * R(0,0) 
			+  R(1,0) * R(1,0) );

    bool singular = sy < 1e-6; // If

    float x, y, z;
    if (!singular)
    {
        x = atan2(R(2,1) , R(2,2));
        y = atan2(-R(2,0), sy);
        z = atan2(R(1,0), R(0,0));
    }
    else
    {
        x = atan2(-R(1,2), R(1,1));
        y = atan2(-R(2,0), sy);
        z = 0;
    }
    return Vec3(x, y, z);
    
    
    
}

Vec3 rotationMatrixToEulerAngles(Mat3 &R, bool whichSol)
{
    assert(isRotationMatrix(R)); 
		double x1 = 0.0, y1 = 0.0, z1 = 0.0
			, x2 = 0.0, y2 = 0.0, z2 = 0.0;

 if (abs(1.00 - R(2,0)) < 0.001) {
		y1 = -M_PI/2.0
				, x1 = -z1 + atan2(-R(0,1), -R(0,2));
	} else if (abs(R(2,0) + 1.00) < 0,001) {
		y1 = M_PI/2.0
				, x1 = z1 + atan2(R(0,1), R(0,2));		
	} else	{
			 y1 = -asin(R(2,0))
        , y2 = M_PI -y1;

				x1 = atan2(R(2,1)/cos(y1), R(2,2)/cos(y1)),
				x2 = atan2(R(2,1)/cos(y2), R(2,2)/cos(y2));

				z1 = atan2(R(1,0)/cos(y1), R(0,0)/cos(y1)),
				z2 = atan2(R(1,0)/cos(y2), R(0,0)/cos(y2));
		}
/*getDeltaQ*/

/*
double thetaX, thetaY, thetaZ;
if (R(0,2) < +1.0) {
	if (R(0,2) > -1.0) {
		thetaY = asin(R(0,2));
		thetaX = atan2(-R(1,2) , R(2,2)) ;
		thetaZ = atan2(-R(0,1) ,R(0,0)) ;
	}
	else	{
		thetaY = -M_PI/2.0;
		thetaX = -atan2(R(1,0) , R(1,1)) ;
		thetaZ = 0.0;
	}
}
else	{
	
	thetaY = M_PI/2.0;
	thetaX = atan2(R(1,0) , R(1,1)) ;
	thetaZ = 0.0;

}
*/	if (whichSol)
	    return Vec3(x1, y1, z1);
		else 
	    return Vec3(x2, y2, z2);
}



Eigen::MatrixXd fixEulThetas(Eigen::MatrixXd inEulThetas) {

	Eigen::MatrixXd retVal(
		inEulThetas.rows(), inEulThetas.cols());

	for (int row = 0; row < inEulThetas.rows(); row++) {
		/*
		//at t_0
		double curTheta = inEulThetas(row,0);

		while(curTheta >= -M_PI) curTheta -= 2.0 * M_PI;
		while(curTheta <= M_PI) curTheta += 2.0 * M_PI;
		retVal(row,0) = curTheta;

			while(curTheta <= -2.0 * M_PI) curTheta += 2.0 * M_PI;
			while(curTheta >= 2.0 * M_PI) curTheta -= 2.0 * M_PI;
*/		
//		while(curTheta >= -M_PI) curTheta -= 2.0 * M_PI;
//		while(curTheta <= M_PI) curTheta += 2.0 * M_PI;

		for(size_t i = 2; i < inEulThetas.cols(); i++) {
			double curTheta = inEulThetas(row,i);

			double testA = curTheta + 2.0 * M_PI;
			double testB = curTheta - 2.0 * M_PI;

			double distNull = abs( curTheta - retVal(row, i - 1));
			double distA = abs( testA - retVal(row, i - 1));
			double distB = abs( testB - retVal(row, i - 1));
		
			double testC = distA < distB ? testA : testB;
			//double distC = distA < distB ? distA : distB;
		  double distC = testC - retVal(row, i - 1);
/*
cout << "curTheta = " << curTheta << endl;
cout << "testA = " << testA << endl;
cout << "testB = " << testB << endl;

cout << "distNull = " << distNull << endl;
cout << "distA = " << distA << endl;
cout << "distB = " << distB << endl;
*/
		  retVal(row, i) = distNull < distC ? curTheta : testC;
		}
	}
		return retVal;

}

/*
// The returned angles are such that we have the following equality: 
Vector3f ea_i = 
	rMatsChInChF.block(0,i*3,3,3).eulerAngles(0,1 ,2); 

// AngleAxis (const Scalar &angle
// 		, const MatrixBase< Derived > &axis)
 
mat == AngleAxisf(ea_i[0], Vector3f::UnitZ())
     * AngleAxisf(ea_i[1], Vector3f::UnitX())
     * AngleAxisf(ea_i[2], Vector3f::UnitZ()); 

*/

/** Build chest local frame **/
	Mat3 buildChRotReferentialAt	(
		DikProblem* dikProb, int index) {
		
		// (1). chest Rot matrix in Maya 
		//        local chest frame
		Mat3 rMatChInVicon_i = Mat3::Identity(3,3);

		Vec3
			MchToRchVhN = (dikProb->RChP(index) - 
				dikProb->MChP(index)).normalized(),
			MchToLchVhN = (dikProb->LChP(index) - 
				dikProb->MChP(index)).normalized();

		//Xch - Sum of these two vectors, pointing up
		rMatChInVicon_i.col(0) = (MchToRchVhN 
			+ MchToLchVhN).normalized();
		//Zch - Cross-product of these two vectors
		//	 pointing to back
		rMatChInVicon_i.col(2) = ((MchToRchVhN).cross(
			MchToLchVhN)).normalized();
		// in case midCh marker is placed 
		// 	below Lch & rCh markers
		Vec3 rMatChInVicon_0 = Vec3::Zero(3,1);
		rMatChInVicon_0(2,0) = 1.0;
		double dotProduct = 
			rMatChInVicon_i.col(0).dot(rMatChInVicon_0);

		if (dotProduct < 0.0)	{
			rMatChInVicon_i.col(0) = 
				-rMatChInVicon_i.col(0);
			rMatChInVicon_i.col(2) = 
				-rMatChInVicon_i.col(2);
		}
			
		//Ych - Zch X Xch
		rMatChInVicon_i.col(1) = 
			((rMatChInVicon_i.col(2)).cross(
				rMatChInVicon_i.col(0))).normalized();
		
		return rMatChInVicon_i;
	}

	Eigen::MatrixXd buildChRotReferentials	(
		DikProblem* dikProb) {

		const int nObservations = 
			dikProb->nObservations();
		MatrixXd rMatsCh = 
			MatrixXd::Zero(3,3*(nObservations+1));
		for (int i = 0; i < nObservations; i++)	{
			Mat3 rmCh_i = 
				buildChRotReferentialAt (dikProb, i);
			rMatsCh.block(0,i*3,3,3) = rmCh_i;

		}
		return rMatsCh;
	}

	/** Build chest local frame **/
	Eigen::Matrix<double, 3, 3>
	 buildChRotReferentialInChFAt	(
		DikProblem* dikProb, int index) {
		Mat3 rMatChInVicon_i = 
			buildChRotReferentialAt	(dikProb, index);
		// (2). offset - vicon frame in Maya local  
		//				chest frame which is const. over time
		Mat3 rMatViconWrtChInMaya = 
			Eigen::MatrixXd::Zero(3,3);
		rMatViconWrtChInMaya(0,2) = 1.0;
		rMatViconWrtChInMaya(1,1) = -1.0;
		rMatViconWrtChInMaya(2,0) = 1.0;
		
		// (3). chest Rot matrix in Maya 
		//				local chest frame
		Mat3 rMatsChInMaya_i = 
			rMatViconWrtChInMaya * rMatChInVicon_i;
		return rMatsChInMaya_i;
	}

	Eigen::MatrixXd buildChRotReferentialsInChF	(
		DikProblem* dikProb) {

		const int nObservations = 
			dikProb->nObservations();
		MatrixXd rMatsChInChF = 
			MatrixXd::Zero(3,3*nObservations);
		for (int i = 0; i < nObservations; i++)	{
			Mat3 rmChInChF_i = 
				buildChRotReferentialInChFAt (dikProb, i);
			rMatsChInChF.block(0,i*3,3,3) = rmChInChF_i;
		}
		return rMatsChInChF;
	}

	// chest orientation as Eurler angles
	// 		The returned angles are in the  
	//		ranges [0:pi]x[-pi:pi]x[-pi:pi].

	//in Vicon Frame
	Eigen::MatrixXd getChEulerAngles 
		(DikProblem* dikProb) {

		const int nObservations = 
			dikProb->nObservations();
		MatrixXd xyzEulersCh = 
			MatrixXd::Zero(3,nObservations);
		MatrixXd rMatsCh = 
			buildChRotReferentials(dikProb);

		for (int i = 0; i < nObservations; i++)	{
			Mat3 rmCh_i = rMatsCh.block(0,i*3,3,3);
			Vec3 eaCh_i = 
//				rotationMatrixToEulerAngles(rmCh_i, true);
				rmCh_i.eulerAngles(0, 1, 2); 

			xyzEulersCh.col(i) = eaCh_i;

		}

		//ensure theta(t) is a monotonic function 
//		MatrixXd xyzEulersChFixed = 
//			fixEulThetas(xyzEulersCh);

		return xyzEulersCh;
	}

	// In Maya local spine frame at 0
	Eigen::MatrixXd getChEulerAnglesInChF (
		DikProblem* dikProb) {
		const int nObservations = 
			dikProb->nObservations();
		MatrixXd xyzEulersChInChF = 
			MatrixXd::Zero(3,nObservations);
		MatrixXd rMatsChInChF = 
			buildChRotReferentialsInChF(dikProb);


		for (int i = 0; i < nObservations; i++)	{
			Mat3 rmChInChF_i = 
				rMatsChInChF.block(0,i*3,3,3);

			Vec3 eaChInChF_i = 
					rmChInChF_i.eulerAngles(0, 1, 2); 

//				rotationMatrixToEulerAngles
//					(rmChInChF_i,true);
			xyzEulersChInChF.col(i) = eaChInChF_i;		
		}

		//ensure theta(t) is a monotonic function 
//		MatrixXd xyzEulersChFixedInChF =
//			fixEulThetas(xyzEulersChInChF);

		return xyzEulersChInChF;
	}

	
	/** get Base ORIENTATION in local frame **/
	Eigen::MatrixXd buildBaseRotReferentials
		( DikProblem* dikProb) {
		string phase = dikProb->phase();
		int ikFrameBase = phase.find("inCh");
		const int nObservations = 
			dikProb->nObservations();
		MatrixXd rMatBaseInBaseF =
		 MatrixXd::Zero(3,3*nObservations);
		MatrixXd xyzEulersBase =
		 MatrixXd::Zero(3,nObservations);
		const int qBase_block_size =
				dikProb->qBase_block_size();
		double* qsBase = dikProb->mutable_qsBase();

		for (int i = 0; i < nObservations; i++)	{
   		double* qbase_i = qsBase
			+ qBase_block_size * i;
			Eigen::Map< Eigen::Quaternion<double> > 
				q_b_quat(qbase_i);
			//as Rotation Matrix		
			Mat3 rmBase_i = 
				q_b_quat.normalized().toRotationMatrix();


			// offset joint orientations 
			//			 chest and base offset at zero
			Mat3  
				rMatViconInCh0 = Mat3::Zero(3,3)
				, rMatChInBase0 = Mat3::Zero(3,3);

			// (1). chest and base offset at zero
			rMatChInBase0(0,0) = -1.0;
			rMatChInBase0(1,2) = -1.0;
			rMatChInBase0(2,1) = -1.0;

			// (2). Vicon and chest offset at zero
	/*
			//(I) constant
			rMatViconInCh0(0,2) = 1.0;
			rMatViconInCh0(1,1) = -1.0;
			rMatViconInCh0(2,0) = 1.0;
	*/
			//(II) changing wrt time
			Mat3 rMatViconInCh_i = buildChRotReferentialAt(
				dikProb, i).transpose();
			Mat3 rMatBaseInBaseF_i = Mat3::Zero(3,3);
			// A. solver in chest frame
			if (ikFrameBase!=std::string::npos) {
					rMatBaseInBaseF_i = 
				rMatChInBase0 * rmBase_i;
				}
			//B. Vicon frame
			else {
				Mat3 rMatBaseInCh_i = 
					rMatViconInCh_i * rmBase_i;
				rMatBaseInBaseF_i = rMatChInBase0 
					* rMatBaseInCh_i;
					//.transpose();
			}
			rMatBaseInBaseF.block(0,i*3,3,3) 
				= rMatBaseInBaseF_i;
		}
		return rMatBaseInBaseF;
	}

	Eigen::MatrixXd getBaseEulerAngles (
		DikProblem* dikProb) {
		const int nObservations 
			= dikProb->nObservations();
			
		MatrixXd rMatBaseInBaseF = 
			buildBaseRotReferentials(dikProb);
		Mat3 rMatBaseInBaseF_i = Mat3::Zero(3,3); 
		MatrixXd xyzEulersBase = 
		MatrixXd::Zero(3,nObservations),	
			xyzEulersBase2 = 
				MatrixXd::Zero(3,nObservations);
			
		
		for (int i = 0; i < nObservations; i++)	{
		
			rMatBaseInBaseF_i = 
				rMatBaseInBaseF.block(0,i*3,3,3);

			Vec3 eaBase_i = 
//				rotationMatrixToEulerAngles
//					(rMatBaseInBaseF_i,true);
				rMatBaseInBaseF_i.eulerAngles(0, 1, 2); 
			xyzEulersBase.col(i) = eaBase_i;
			
			Vec3 eaBase_i2 = rotationMatrixToEulerAngles
				(rMatBaseInBaseF_i,false);							
			xyzEulersBase2.col(i) = eaBase_i2;
	
		}
//		printEigenMathematica(xyzEulersBase2.transpose()
//			, cout, "xyzEulersBase2");

		return xyzEulersBase;
	}

	Eigen::Matrix<double, 3, 3> buildElRotReferential	
		(DikProblem* dikProb, int index) {
	/** Build elbow local frame **/
	Eigen::Vector3d 
		ElToLaVhN = (dikProb->LaP(index) 
			- dikProb->ElP(index)).normalized(),
		ElToWrVhN = (dikProb->WrP(index) 
			- dikProb->ElP(index)).normalized();
	Mat3 Rel_i;
	//Xel - Sum of these two vectors, pointing up
	Rel_i.col(0) = (ElToLaVhN 
		+ ElToWrVhN).normalized();
	//Yel - Cross-product of these two vectors 
	Rel_i.col(1) = ((ElToLaVhN).cross(
		ElToWrVhN)).normalized();
	//Zel - Xel X Yel
	Rel_i.col(2) = ((Rel_i.col(0)).cross(
		Rel_i.col(1))).normalized();


	return Rel_i;
	}

	Eigen::Matrix<double, 3, 3> buildUaRotReferential	
		(DikProblem* dikProb, int index) {
	/** Build UA local frame **/
	Eigen::Vector3d 
		ShToUaVhN = (dikProb->UaP(index) 
			- dikProb->shoP(index)).normalized(),
		ShToElVhN = (dikProb->ElP(index) 
			- dikProb->shoP(index)).normalized();
	Mat3 RUa_i;
	//Xel - Sum of these two vectors, pointing up
	RUa_i.col(0) = (ShToUaVhN 
		+ ShToElVhN).normalized();
	//Yel - Cross-product of these two vectors 
	RUa_i.col(1) = ((ShToUaVhN).cross(
		ShToElVhN)).normalized();
	//Zel - Xel X Yel
	RUa_i.col(2) = ((RUa_i.col(0)).cross(
		RUa_i.col(1))).normalized();
	return RUa_i;
	}

	Eigen::Matrix<double, 3, 3> buildElRotReferentialInElFAt	
		(DikProblem* dikProb, int index) {
	
	Mat3 Rel_i = buildElRotReferential(dikProb, index);

	Mat3 RUa_i = buildUaRotReferential(dikProb, index);

	/** orientation of UA frame at home position	
		 wrt El frame - CONSTANT over time **/	
	Mat3 RUaInEl0 = Mat3::Identity(3,3);

	/** orientation of El frame at time_i w.r.t 
				Elbows's local frame **/

	Mat3 RElInEl0 = RUaInEl0 
		* RUa_i.transpose() * Rel_i;

	return RElInEl0;
	}

	Eigen::MatrixXd getElEulerAngles (
		DikProblem* dikProb) {

		const int nObservations = 
			dikProb->nObservations();
		MatrixXd xyzEulersEl = 
			MatrixXd::Zero(3,nObservations);

		for (int i = 0; i < nObservations; i++)	{
		Mat3 rMatEl_i = 
			buildElRotReferentialInElFAt(dikProb, i);

			Vec3 eaEl_i = rMatEl_i.eulerAngles(0, 1, 2); 
			xyzEulersEl.col(i) = eaEl_i;
		}
		return xyzEulersEl;
	}

	/** Build Wrist local frame wrt chest frame**/
	Eigen::Matrix<double, 3, 3>
		buildWrRotReferential	(
			DikProblem* dikProb, int index)	{	
		Eigen::Vector3d 
			wrToThVhN  = (dikProb->ThP(index) 
				- dikProb->WrP(index)).normalized(),
			wrToPiVhN  = (dikProb->PiP(index) 
				- dikProb->WrP(index)).normalized();
	Mat3 RWr_i;
	//Xwr - Sum of these two vectors, along the hand
	RWr_i.col(0) = wrToThVhN  + wrToPiVhN ;
	//Ywr - Cross-product of these two vectors pointing up
	RWr_i.col(1) = ((wrToPiVhN).cross(
		wrToThVhN)).normalized();
	//Zwr - Xwr X Ywr
	RWr_i.col(2) = ((RWr_i.col(0)).cross(
		RWr_i.col(1))).normalized();

	return RWr_i;
	}
	/** orientation of wrist frame at time_i w.r.t 
		wrist's local frame **/
	Eigen::Matrix<double, 3, 3> buildWrRotReferentialInWrFAt	(
		DikProblem* dikProb, int index) {
		Mat3 R_ch_i = 
			buildChRotReferentialAt(dikProb, index);

		Mat3 rMatEl_i = 
			buildElRotReferential(dikProb, index);
		Mat3 rMatElInCh_i = R_ch_i.transpose() * rMatEl_i;

		Mat3 rMatWr_i = 
			buildWrRotReferential(dikProb, index);

		Mat3 rMatWrInCh_i = R_ch_i.transpose() * rMatWr_i;

		/** orientation of wrist frame at home position	
			 wrt elbow frame - CONSTANT over time **/	
		Mat3 RelInWr0 = Mat3::Identity(3,3);

		/** orientation of wrist frame at time_i w.r.t 
					wrist's local frame **/

		Mat3 RWrInWr0 = RelInWr0 
			* rMatEl_i.transpose() * rMatWr_i;

//		Mat3 RWrInWr0 = RelInWr0 
//			* rMatElInCh_i.transpose() * rMatWrInCh_i;
	
		return rMatWr_i;
	}

	Eigen::MatrixXd getWrEulerAngles (
		DikProblem* dikProb) {

		const int nObservations = 
			dikProb->nObservations();
		MatrixXd xyzEulersWr = 
			MatrixXd::Zero(3,nObservations);
		MatrixXd xzxEulersWr = 
			MatrixXd::Zero(3,nObservations);

		for (int i = 0; i < nObservations; i++)	{
			Mat3 rmWr_i = 
				buildWrRotReferentialInWrFAt (dikProb, i);

			Vec3 eaWrXZX_i = rmWr_i.eulerAngles(0, 2, 0); 
			xzxEulersWr.col(i) = eaWrXZX_i;
			Vec3 eaWr_i = rmWr_i.eulerAngles(0, 1, 2); 
			xyzEulersWr.col(i) = eaWr_i;
		}

	printEigenMathematica(xzxEulersWr.transpose()
		, cout, "xzxEulersWr");
		return xyzEulersWr;
	}
}
#endif

