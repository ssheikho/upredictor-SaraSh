#ifndef BALL_CENTROIDS_ON_RACK_H
#define BALL_CENTROIDS_ON_RACK_H

#include <Eigen/Core>
#include <string>

#include "ParseCSV.h"
#include <UBCUtil.h>
using namespace Eigen;
using namespace std;

class BallCentroidsOnRack {
public:
	BallCentroidsOnRack(string inLine);
	~BallCentroidsOnRack();

	MatrixXd getMarkerLocs();
	MatrixXd getBallCentroids();

	void setMarkerLocs();
	void setBallCentroids();
	
	template<typename T, size_t N>
	T * end(T (&ra)[N]) {
	    return ra + N;
	}

	static std::vector<std::string> _rackMarkers;	

protected:
	
	string _inLine;
	int _nCols;	
	MatrixXd _rackMarkerLoc;
	//MatrixXd _Offsets;
	MatrixXd _ballCentroids;
};

#endif
