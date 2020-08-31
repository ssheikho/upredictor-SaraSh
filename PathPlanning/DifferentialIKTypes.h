#ifndef DIFFERENTIAL_IK_TYPES_H
#define DIFFERENTIAL_IK_TYPES_H

#include <istream>
#include <map>
#include <string>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "BasicFormulas.h"

typedef Eigen::Matrix<double, 3, 3> Mat3;
typedef Eigen::Matrix<double, 7, 1> Vec7;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector4d Vec4;
using std::vector;

//namespace ceres {

	// all_markers is a vector of all tracked markers existing 
	// in the a reach sample. all_markers to be bundled as one
	// parameter.
	enum Markers {
	 			RWr = 0,
				RTh = 3,
				RPi = 30,
				RLA = 6,
				REl = 9, 
				RUA = 12, 
				RSh = 15, 
				RCh = 18, 
				MCh = 21,
				LCh = 24,
				LSh = 27,
	};

	// A EuclideanCamera is the location and rotation of the 
	// camera. track identifies which time-frame this camera 
	// represents.
	// R is a 3x1 matrix of Euler parameters representing the 
	// rotation of the base of wam w.r.t camera.
	// t is a translation vector representing its positions.
	struct EuclideanCamera {
		EuclideanCamera() : track(-1) {}
		EuclideanCamera(const EuclideanCamera &c) : track(c.track)
		, R(c.R), t(c.t) {}
		int track;
		Vec3 R;
		Vec3 t;
	};

	// A Point is the 3D location of a marker at time track.
	// X represents the 3D position of the track.
	struct EuclideanPoint {
		EuclideanPoint() : track(-1) {}
		EuclideanPoint(const EuclideanPoint &p) 
		: track(p.track)
		, X(p.X) {}
		int track;
		Vec3 X;
	};


	typedef std::vector<EuclideanPoint, 
		Eigen::aligned_allocator<EuclideanPoint> >
	    VectorOfmarkerPts;

	// A Theta is the 7-vec of WAM joint angles
	// at time-frame 'track' 
	struct jointSpaceTheta {
		jointSpaceTheta() : track(-1) {}
		jointSpaceTheta(const jointSpaceTheta &theta) : 
			track(theta.track), thetasI(theta.thetasI) {}
		int track;
		Vec7 thetasI;
	};

	struct CartesianVelocity {
		CartesianVelocity() : track(-1) {}
		CartesianVelocity(const CartesianVelocity &cv) 
		: track(cv.track)
		, linearVelocity(cv.linearVelocity) 
		, anglularVelocity(cv.anglularVelocity)	{}
		int track;
		Vec3 linearVelocity;
		Vec3 anglularVelocity;
	};

	struct Pose3d { 
/*
		Mat3 getRmatfrom3Pts(Vec3 pA, Vec3 pB)	{
			return	buildRefFramefrom3Pts(pA, p, pB);
		}
		
		Eigen::Vector3d getEulAnglesFromRmat(Mat3 Rmat)	{
			return getEulAnglesFromRmat(Rmat);
		}
*/
		Eigen::Vector3d p;
	  Eigen::Quaterniond q;
 		Eigen::Vector3d eulPars;

		/*	Eigen provided macro - If you define a structure 
				having members of fixed-size vectorizable Eigen types, 
				you must overload its "operator new" so that it 
				generates 16-bytes-aligned pointers.*/
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	typedef std::map<int, Pose3d, std::less<int>,
		Eigen::aligned_allocator<std::pair<const int, Pose3d> > >
    	MapOfPoses;

		//what goes into IK error term at each time-frame
	struct Constraint3d {
 		int id_begin;
	  int id_end;
  	// The transformation that represents the pose of the 
		//	wrist frame w.r.t. the camera
  	Pose3d t_cWr;
  	// The inverse of the covariance matrix of inPts
		//The order of the entries are x, y, z, delta orientation.
	  Eigen::Matrix<double, 6, 6> information;

	  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	typedef std::vector<Constraint3d, 
		Eigen::aligned_allocator<Constraint3d> >
	    VectorOfConstraints;
//}  // namespace ceres

#endif  
