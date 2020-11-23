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
	
	void keyAdd(MatrixXd inPtsMB, int startFrame)	{
		int nMarkers = inPtsMB.rows();
		int nPts = inPtsMB.cols();
		int mkr = 0;
		while ((mkr < nMarkers)) {
			MatrixXd inPtsMkr_i = 
				inPtsMB.block(mkr,0,3,nPts);
			for (size_t i = 0; i < nPts; i++){
				cout << "marker.Translation.GetAnimationNode()"	
					<< ".KeyAdd(FBTime(0, 0, 0,"<< startFrame + i 
					<< ",0), [" << inPtsMkr_i(0,i) << "," 
					<<  inPtsMkr_i(1,i)  << "," << inPtsMkr_i(2,i)
					<< "])" << endl;
			}
			mkr += 3;
		}
	}



}

namespace MayaAnimation {

	/** Build chest local frame **/
	Eigen::Matrix<double, 3, 3> buildChRotReferentialAt	(
		DikProblem* dikProb, int index) {
		
		// (1). chest Rot matrix in Maya local chest frame
		Mat3 rMatChInVicon_i = Mat3::Identity(3,3);

		Vec3
			MchToRchVhN = (dikProb->RChP(index) - 
				dikProb->MChP(index)).normalized(),
			MchToLchVhN = (dikProb->LChP(index) - 
				dikProb->MChP(index)).normalized();

		//Xch - Sum of these two vectors, pointing up
		rMatChInVicon_i.col(0) = (MchToRchVhN 
			+ MchToLchVhN).normalized();
		//Zch - Cross-product of these two vectors pointing to back
		rMatChInVicon_i.col(2) = ((MchToRchVhN).cross(
			MchToLchVhN)).normalized();
		//Ych - Zch X Xch
		rMatChInVicon_i.col(1) = ((rMatChInVicon_i.col(2)).cross(
			rMatChInVicon_i.col(0))).normalized();
		
		return rMatChInVicon_i;
	}

	/** Build chest local frame **/
	Eigen::Matrix<double, 3, 3> buildChRotReferentialInChFAt	(
		DikProblem* dikProb, int index) {
		Mat3 rMatChInVicon_i = 
			buildChRotReferentialAt	(dikProb, index);
		// (2). offset - vicon frame in Maya local chest frame
		//				which is const. over time
		Mat3 rMatViconWrtChInMaya = 
			Eigen::MatrixXd::Zero(3,3);
		rMatViconWrtChInMaya(0,2) = 1.0;
		rMatViconWrtChInMaya(1,1) = -1.0;
		rMatViconWrtChInMaya(2,0) = 1.0;
		
		// (3). chest Rot matrix in Maya local chest frame
		Mat3 rMatsChInMaya_i = 
			rMatViconWrtChInMaya * rMatChInVicon_i;
		return rMatsChInMaya_i;
	}

	Eigen::MatrixXd getChEulerAngles (DikProblem* dikProb) {

		const int nObservations = dikProb->nObservations();
		MatrixXd xyzEulersCh = MatrixXd::Zero(3,nObservations);

		for (int i = 0; i < nObservations; i++)	{
			Mat3 rmCh_i = 
				buildChRotReferentialInChFAt (dikProb, i);
			Vec3 eaCh_i = rmCh_i.eulerAngles(0, 1, 2); 
			xyzEulersCh.col(i) = eaCh_i;
		}
		return xyzEulersCh;
	}

	
	/** Build Base local frame **/
	Eigen::Matrix<double, 3, 3> buildBaseRotReferentialAt	(
		DikProblem* dikProb, int index) {
		
		// (1). chest Rot matrix in Maya local chest frame
		Mat3 rMatBaseInVicon_i = Mat3::Identity(3,3);
		Vec3
			MchToRchVhN = (dikProb->RChP(index) - 
				dikProb->MChP(index)).normalized(),
			MchToLchVhN = (dikProb->LChP(index) - 
				dikProb->MChP(index)).normalized();

		//Xch - Sum of these two vectors, pointing up
		rMatBaseInVicon_i.col(0) = (MchToRchVhN 
			+ MchToLchVhN).normalized();
		//Ych - Cross-product of these two vectors pointing to back
		rMatBaseInVicon_i.col(1) = ((MchToRchVhN).cross(
			MchToLchVhN)).normalized();
		//Ych - Zch X Xch
		rMatBaseInVicon_i.col(2) = ((rMatBaseInVicon_i.col(0))
			.cross(rMatBaseInVicon_i.col(1))).normalized();
		
		return rMatBaseInVicon_i;
	}

	Eigen::MatrixXd getBaseEulerAngles (
		DikProblem* dikProb) {

		const int nObservations = dikProb->nObservations();
		MatrixXd xyzEulersBase = MatrixXd::Zero(3,nObservations);

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
			Vec3 eaBase_i = rmBase_i.eulerAngles(0, 1, 2); 
			xyzEulersBase.col(i) = eaBase_i;
		}
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

	Eigen::MatrixXd getElEulerAngles (DikProblem* dikProb) {

		const int nObservations = dikProb->nObservations();
		MatrixXd xyzEulersEl = MatrixXd::Zero(3,nObservations);

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
		buildWrRotReferential	(DikProblem* dikProb, int index)	{	
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
		Mat3 R_ch_i = buildChRotReferentialAt(dikProb, index);

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

	Eigen::MatrixXd getWrEulerAngles (DikProblem* dikProb) {

		const int nObservations = dikProb->nObservations();
		MatrixXd xyzEulersWr = MatrixXd::Zero(3,nObservations);
		MatrixXd xzxEulersWr = MatrixXd::Zero(3,nObservations);

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

