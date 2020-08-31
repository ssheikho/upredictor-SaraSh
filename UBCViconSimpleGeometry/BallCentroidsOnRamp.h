#ifndef BALL_CENTROIDS_ON_RAMP_H
#define BALL_CENTROIDS_ON_RAMP_H

#include <Eigen/Core>
#include <string>

using namespace Eigen;
using namespace std;

class BallCentroidsOnRamp {
public:
	BallCentroidsOnRamp(string inLine);
	~BallCentroidsOnRamp();

	MatrixXd getMarkerLocs();
	MatrixXd getBallCentroids();

	void setMarkerLocs();
	void setBallCentroids();
	
	template<typename T, size_t N>
	T * end(T (&ra)[N]) {
	    return ra + N;
	}

	static std::vector<std::string> _rampMarkers;	

protected:
	
	string _inLine;
	int _nCols;	
	MatrixXd _rampMarkerLoc;
	//MatrixXd _Offsets;
	MatrixXd _rampBallCentroids;
};

#endif
