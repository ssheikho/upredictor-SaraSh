#ifndef PLANE_H
#define PLANE_H

#include <Eigen/Core>

using namespace Eigen;

class Plane {
public:
	//Points assumed to be Cartesian, 1 point per row
	Plane(MatrixXd inPts);
	~Plane();


	Vector4d getPlane();
	Vector3d getPlaneNormalVect();
	double getNorm();
	Vector4d getPlaneN();
	Vector3d getPlaneNormalVectN();
	Vector3d getPlaneV1N();
	Vector3d getPlaneV2N();
	Vector3d getPlaneOrigin();
	double getPlaneZoff();

	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> 
		getPlaneEigenSys();


	double getPlaneNtoZ();
	double getPlaneYawZ();
	double getPlaneDipX();
	double getPlaneTiltY();


protected:
	Vector4d _plane;
	Vector3d _planeNormalVect;
	double _norm;
	Vector4d _planeN;
	Vector3d _planeNormalVectN, _planeOrigin;
	Vector3d _planeV1N, _planeV2N;
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd>  _eigenSys;
	double _planeZoff; 
	double _planeNtoZ, _planeYawZ, 
		_planeDipX, _planeTiltY;

};

#endif
