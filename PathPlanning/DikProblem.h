#ifndef DIK_PROBLEM_H
#define DIK_PROBLEM_H

#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "DifferentialIKErrorTerm.h"
#include "UBCUtil.h"

#include "ceres/ceres.h"

using namespace Eigen;
using namespace std;
typedef Eigen::Matrix<double, 3, 3> Mat3;
typedef Eigen::Matrix<double, 7, 1> Vec7;
typedef Eigen::Matrix<double, 3, 1> Vec3;
typedef Eigen::Vector4d Vec4;
using std::vector;

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::CauchyLoss;

typedef std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >
    VectorOfQuaternionds;

/* Denote the (Scaled/maped) positions of input markers
	 	w.r.t. the inertial frame (camera) by a three-vector 
		O_j (which gives the coordinates of the origin of the 
		joint j frame with respect to the camera frame) */

class DikProblem {

public:
	DikProblem (ForwardKin<double> &fk
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles
		, bool use_quaternions
		, Eigen::MatrixXd inPts	
		);

	~DikProblem ();

	std::vector<double> jointMinAngles()	{
		return _jointMinAngles;
	}

	std::vector<double> jointMaxAngles()	{
		return _jointMaxAngles;
	}

	int qBase_block_size() const {
		return _use_quaternions ? 4 : 3; 
	}

  int theta_block_size()	const {
		return 7;                         
	}

	int nObservations () {
		return _nPts;	
	}

	int nMarkers () {
		return _inPts.rows() / 3;	
	}

  int nParameters()	const { 
		return _nParameters;           
	}

  const int* theta_index()	const {
		return _theta_index;              
	}

  const int* qBase_index()	const { 
		return _qBase_index;             
	}


  const double* parameters()   const { 
		return _parameters;               
	}

  const double* qsBase()	const { 
		return _parameters;               
	}
  double* mutable_qsBase() { 
		return _parameters;               
	}
  double* mutable_thetas() {
    return _parameters  + qBase_block_size() * nObservations();
  }

	Eigen::Matrix<double, Dynamic, 1> markerPtsAt	(
			int index)	{
		return _inPts.col(index).transpose();
	}
	
	Eigen::Matrix<double, Dynamic, 1> markerPtsWamAt	(
			int index)	{
		return _inPtsWam.col(index);
	}

	ForwardKin<double>& getWam()	{
		return _fk;	
	}

	Eigen::Matrix<double, 3, 1> LChP (int i)	{
		return _inPLCh.col(i);
	}

	Eigen::Matrix<double, 3, 1> RChP (int i)	{
		return _inPRCh.col(i);
	}

	Eigen::Matrix<double, 3, 1> MChP (int i)	{
		return _inPMCh.col(i);
	}

	Eigen::Matrix<double, 3, 1> shoP (int i)	{
		return _inPtsWam.block(Markers::RSh,i,3,1);
	}

	Eigen::Matrix<double, 3, 1> UaP (int i)	{
		return _inPtsWam.block(Markers::RUA,i,3,1);
	}

	Eigen::Matrix<double, 3, 1> ElP (int i)	{
		return _inPtsWam.block(Markers::REl,i,3,1);
	}

	Eigen::Matrix<double, 3, 1> LaP (int i)	{
		return _inPtsWam.block(Markers::RLA,i,3,1);
	}

	Eigen::Matrix<double, 3, 1> WrP (int i)	{
		return _inPtsWam.block(Markers::RWr,i,3,1);
	}

	Eigen::Matrix<double, 3, 1> PiP (int i)	{
		return _inPtsWam.block(Markers::RPi,i,3,1);
	}

	Eigen::Matrix<double, 3, 1> ThP (int i)	{
		return _inPtsWam.block(Markers::RTh,i,3,1);
	}
	
	void mapInPtsHtoVecsR ();
	void mapInPtsHtoInPtsR ();
	Eigen::Matrix<double, 3, 3> 
		buildWrRotReferential	(int index);
	Mat3 buildJointsRotReferential(int index);
	//ordered as (x_i y_i z_i)T
	Eigen::MatrixXd pEeWam () {
		return _pEeWam;
	}

	Eigen::MatrixXd pPiWam () {
		return _pPiWam;
	}

	Eigen::MatrixXd pThWam () {
		return _pThWam;
	}

	Eigen::MatrixXd pWrWam () {
		return _pWrWam;
	}

	Eigen::MatrixXd pLaWam (){
		return _pLaWam;
	}

	Eigen::MatrixXd pElWam (){
		return _pElWam;
	}

	Eigen::MatrixXd pUaWam(){
		return _pUaWam;
	}

	Eigen::MatrixXd pShWam(){
		return _pShWam;
	}

	Eigen::MatrixXd shoToUaVsWam(){
		return _shoToUaVsWam;
	}

	Eigen::MatrixXd shoToElVsWam(){
		return _shoToElVsWam;
	}

	Eigen::MatrixXd shoToLaVsWam(){
		return _shoToLaVsWam;
	}

	Eigen::MatrixXd elOffsetVsWam() {
		return _elOffsetVsWam;
	}

	Eigen::MatrixXd wrOffsetVsWam(){
		return _wOffsetVsWam;
	}

	Eigen::MatrixXd laToWrVsWam (){
		return _laToWrVsWam;
	}

	Eigen::MatrixXd shoToWrVsWam (){
		return _shoToWrVsWam;
	}

	Eigen::MatrixXd wrToEeVsWam (){
		return _wrToEeVsWam;
	}

	Eigen::MatrixXd shoToEeVsWam (){
		return _shoToEeVsWam;
	}

	//Orientation as mathematica matrices:
	// Chest
	Eigen::MatrixXd quaternionsCh() {
		return _quaternionsCh; }

	Eigen::MatrixXd rMatsCh() {
		return _rMatsCh; }

	Eigen::MatrixXd xyzEulersCh() {
		return _xyzEulersCh;	}
/*	
	Eigen::MatrixXd rMatsChInMaya() {
		return _rMatsChInMaya; }

	Eigen::MatrixXd xyzEulersChInMaya() {
		return _xyzEulersChInMaya;	}

	// end-effectors
	Eigen::MatrixXd quaternionsEe() {
		return _quaternionsEe;	}

	Eigen::MatrixXd rMatsEe() {
		return _rMatsEe;	}

	Eigen::MatrixXd xyzEulersEe() {
		return _xyzEulersEe;	}
*/
	Eigen::MatrixXd outShoToUaVsWam(){
		return _outShoToUaVsWam;
	}

	Eigen::MatrixXd outShoToElVsWam(){
		return _outShoToElVsWam;
	}

	Eigen::MatrixXd outElOffsetVsWam() {
		return _outElOffsetVsWam;
	}

	Eigen::MatrixXd outWrOffsetVsWam(){
		return _outWOffsetVsWam;
	}

	Eigen::MatrixXd outLaToWrVsWam (){
		return _outLaToWrVsWam;
	}

	Eigen::MatrixXd outShoToWrVsWam (){
		return _outShoToWrVsWam;
	}


	void setOutVecsR ();
	void setOutPtsR ();

	// For print in Mathematica 
	//(a) for human
	void printInPtsH();
	//(b) inPts wam, i.e. scaled to wam
	void printPtsR();
	//(c). inVecs Robot
	void printVecsR();
	
protected:

	//the robot
	ForwardKin<double> &_fk;
	//	Eigen::Matrix <double, 7, Eigen::Dynamic>& _JsWam;
	std::vector<double> _jointMinAngles, _jointMaxAngles;
	//parametrization of qBase 
  bool _use_quaternions;

	//in marker points 
	Eigen::MatrixXd  &_inPts;
	int _nPts;
	Eigen::MatrixXd _inPRWr, _inPRTh
		, _inPRPi , _inPRLA, _inPREl, _inPRUA, _inPRSh
		, _inPRCh, _inPMCh, _inPLCh, _inPLSh;

	//scaled VECTORS on the robot
		Eigen::MatrixXd
			_shoToUaVsWam, _shoToElVsWam
		, _elToWrVsWam, _laToWrVsWam
		, _shoToLaVsWam, _shoToWrVsWam
		,	_wrToThVsWam, _wrToPiVsWam
		, _shoToThVsWam, _shoToPiVsWam
		, _wrToEeVsWam,   _shoToEeVsWam
		, _elOffsetVsWam, _wOffsetVsWam;

	//scaled positions on the robot
	Eigen::MatrixXd _inPtsWam;
	Eigen::MatrixXd _pShWam, _pUaWam, _pElWam, _pLaWam
		, _pWrWam, _pThWam, _pPiWam, _pEeWam;

	//Orientation as mathematica matrices:
	// Chest
	Eigen::MatrixXd 
		_quaternionsCh, _rMatsCh, _xyzEulersCh;
	// chest in Maya
	Eigen::MatrixXd 
		_rMatsChInMaya, _xyzEulersChInMaya;
/*
	// elbow
	Eigen::MatrixXd
		_rMatsEl, _xyzEulersEl;

	// end-effectors
	Eigen::MatrixXd
		_quaternionsEe, _rMatsEe, _xyzEulersEe;
	Eigen::MatrixXd
		_rMatsEeLocal, _xyzEulersEeLocal
			,_rMatsEeLocalFvicon, _xyzEulersEeLocalFvicon;
*/

	// linear velocity of scaled markers on the robot
	Eigen::MatrixXd _inLvsWam;
	Eigen::MatrixXd
	_lvelShWam, _lvelUaWam, _lvelElWam, _lvelLaWam
		, _lvelWrWam, _lvelThWam, _lvelPiWam;

	// angular velocity of robot end-effector
	Eigen::Matrix<double, 3, Eigen::Dynamic> 
		_AngularVWrWam;
	
  int _nParameters;
  int* _theta_index;
  int* _qBase_index;
  // The parameter vector is laid out as follows
  // [camera_1, ..., camera_n, point_1, ..., point_m]
  double* _parameters;

	//*** OUTpoints***/	
	//OUT-scaled VECTORS on the robot
		Eigen::MatrixXd
			_outShoToUaVsWam, _outShoToElVsWam
		, _outElToWrVsWam, _outLaToWrVsWam
		, _outShoToLaVsWam, _outShoToWrVsWam
		,	_outWrToThVsWam, _outWrToPiVsWam
		, _outShoToThVsWam, _outShoToPiVsWam
		, _outWrToEeVsWam,   _outShoToEeVsWam
		, _outElOffsetVsWam, _outWOffsetVsWam;
	//OUT-scaled positions on the robot
	Eigen::MatrixXd _outPtsWam;
	Eigen::MatrixXd _outPShWam, _outPUaWam, _outPElWam,
		_outPLaWam	, _outPWrWam, _outPThWam, _outPPiWam;

};

#endif

