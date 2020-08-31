
	void initializeqsBaseToZero()	{    	
		for (int i = 0; i < _nPts; i++) {
			double* qiBase = new double[4];//quaternion
			_qsBase.push_back(qiBase);
		}
	}


	void initializeThetasWamToZero()	{ 

		for (int i = 0; i < _nPts; i++) {

  	 	double* thetas_i = new double[_fk.getNJoints()];
			for (int i = 0; i < _fk.getNJoints(); i++) 
				thetas_i[i] = 0.0;

			_thetasWamV.push_back(thetas_i);
		}
	}


	// returns a pointer to the memory location of the
	// first entry of the matrix.
	double* qBaseAt	(int i) {
		return _qsBase[i];
	}

	double* thetasWamAt	(int i)	{    	
		return _thetasWamV[i];
	}

	Eigen::Matrix<double, Eigen::Dynamic, 1> 
		diffLvelmarkersWamAt (int index)	{
		//inPts WAM to CERES solver
		if (index > 0)
			return _inPtsWam.col(index) - _inPtsWam.col(index-1);
		else 
			return _inPtsWam.col(index) - _inPtsWam.col(index);
	}

	void setDiffLvelmarkersWam  ()	{
		for (int i = 1; i < _inPtsWam.cols(); i++)
			_inLvsWam.col(i) = diffLvelmarkersWamAt(i);

		_lvelPiWam = _inLvsWam.block(Markers::RPi,0,3,_nPts);
		_lvelThWam = _inPtsWam.block(Markers::RTh,0,3,_nPts);
	  _lvelWrWam = _inPtsWam.block(Markers::RWr,0,3,_nPts);
		_lvelLaWam = _inPtsWam.block(Markers::RLA,0,3,_nPts);
		_lvelElWam = _inPtsWam.block(Markers::REl,0,3,_nPts);
		_lvelUaWam = _inPtsWam.block(Markers::RUA,0,3,_nPts);
		_lvelShWam = _inPtsWam.block(Markers::RSh,0,3,_nPts);
	}


	Eigen::Matrix<double, 3, 3> buildWrRotReferential
		(int index) {
		/* First, a vector is defined from the wr marker to the
			 th, and a second vector from the wr marker to the 
				pinkie.*/
		Eigen::Matrix<double, 3, 1> 
			wrToThVhN = (_inPRTh.col(index) - 
				_inPRWr.col(index)).normalized()
			, wrToPiVhN = (_inPRPi.col(index) - 
					_inPRWr.col(index)).normalized();
		Eigen::Matrix<double, 3, 3> retT;
		/* Z	The sum of these two vectors lies on the hand 
				plane and is parallel to the z-axis (i.e. axis 
				of rotation) of the rotated referential.*/
		retT.col(2) = (wrToPiVhN + wrToThVhN).normalized();
		/* The cross-product of these vectors gives the 
				normal to the hand plane (pointing down) and 
				is parallel to the x-axis  */
		retT.col(0) = (wrToPiVhN.cross(wrToThVhN)).normalized();
	 	/* The remaining Y axis can easily be computed 
			 as the cross-product of Z cross X axis */
		retT.col(1) = (retT.col(2).cross(
			retT.col(0))).normalized();
		return retT;
	}

	Eigen::Matrix<double, 3, 1> getEulAnglesFromRmat	(
		Eigen::Matrix<double, 3, 3> Rmat) {
		Eigen::Matrix<double, 3, 1> retV;
		retV(0,0) = ceres::atan2(Rmat(1,2), Rmat(0,2));
		retV(1,0) = ceres::atan2(ceres::sqrt(Rmat(0,2) *
			Rmat(0,2) + Rmat(1,2) * Rmat(1,2)), Rmat(2,2));
		retV(2,0) = ceres::atan2(Rmat(2,1), -Rmat(2,0));
		return retV;
	}

	//	get RotMat from Euler Parameters 
	Eigen::Matrix<double, 3, 3> getRmatFromEulAngles	
		(double *EulAngles)	{
		MatrixXd Rmat = MatrixXd::Zero(3,3);
		//1st col
		Rmat(0,0) = cos(EulAngles[1]) *
			cos(EulAngles[0]) * cos(EulAngles[2]) - 
				sin(EulAngles[2]) * sin(EulAngles[2]);
		Rmat(1,0) = cos(EulAngles[0]) *
			sin(EulAngles[2]) + cos(EulAngles[1]) *
				cos(EulAngles[2]) * sin(EulAngles[0]);
		Rmat(2,0) = -cos(EulAngles[2]) *
			sin(EulAngles[1]);
		//2nd col
		Rmat(0,1) = -cos(EulAngles[2]) *
			sin(EulAngles[0]) - cos(EulAngles[1]) *
				ceres::cos(EulAngles[0]) * sin(EulAngles[2]);
		Rmat(1,1) = cos(EulAngles[0]) *
			cos(EulAngles[2]) - cos(EulAngles[1]) *
				sin(EulAngles[0]) * sin(EulAngles[2]);
		Rmat(2,1) = sin(EulAngles[1]) *
			sin(EulAngles[2]);
		//3rd col
		Rmat(0,2) = cos(EulAngles[0]) *
			sin(EulAngles[1]);
		Rmat(1,2) = sin(EulAngles[0]) *
			sin(EulAngles[1]);
		Rmat(2,2) = cos(EulAngles[1]);
		return Rmat;
}

	Eigen::Matrix<double, 3, 1> EulAnglesWrAt	(int index)	{
		Eigen::Matrix<double, 3, 1> EulAnglesWr =
		Eigen::Matrix<double, 3, 1> ::Zero(3,1);
		Eigen::Matrix<double, 3, 3> Rmat_i =
			buildWrRotReferential(index);
			EulAnglesWr = getEulAnglesFromRmat(Rmat_i);
			return EulAnglesWr;
	}

	//Spong eqn. (4.103)
	Eigen::Matrix<double, 3, 3> mapFncEulAnglesToAngularVel
		(Eigen::Matrix<double, 3, 1> inEulAngles)	{


	// B(α) - map EulAngles to anglular velocity;	
	//		 if α = [φ, θ, ψ]^T => ω = B(α)α̇
  Eigen::Matrix<double, 3, 3> retMappingMat = 
		Eigen::MatrixXd::Identity(3,3);
		//	1st col
		retMappingMat(0,0) = ceres::cos(inEulAngles(2,0))  
			* ceres::sin(inEulAngles(1,0));
		retMappingMat(1,0) = ceres::sin(inEulAngles(2,0))
			* ceres::sin(inEulAngles(1,0));
		retMappingMat(2,0) = ceres::cos(inEulAngles(1,0));
		//			2nd col
		retMappingMat(0,1) = -ceres::sin(inEulAngles(2,0));
		retMappingMat(1,1) = ceres::cos(inEulAngles(2,0));
		return retMappingMat;
	}

	Eigen::Matrix<double, 3, 1> 
		diffAngularVelWrAt (int index)	{
		int indexMinusOne = 0;
		if (index > 0)
			indexMinusOne = index - 1;
 
		//Euil angles α_i = [φ, θ, ψ]^T 
		Eigen::Matrix<double, 3, 1> eulAngles_i = 
			EulAnglesWrAt(index);
		// B(α) - map EulAngles to anglular velocity;	
		Eigen::Matrix<double, 3, 3> mapFnc = 
			mapFncEulAnglesToAngularVel(eulAngles_i);
		//α̇
		Eigen::Matrix<double, 3, 1> eulAngles_iMinusOne = 
			getEulAnglesFromRmat (buildWrRotReferential
				(indexMinusOne));	
		//ω = B(α)α̇
		Eigen::Matrix<double, 3, 1> retAngVel = 
			mapFnc * (eulAngles_i - eulAngles_iMinusOne);
		return retAngVel;
	}

	Eigen::Matrix<double, 6, Eigen::Dynamic> 
		posesMatWamforMarker (int markerId)	{
		Eigen::Matrix<double, 6, Eigen::Dynamic> posesMat
			=	Eigen::Matrix<double, 6, Eigen::Dynamic>::
				Zero(6, _nPts);
		for (int index = 0; index < _nPts; index++) {
			posesMat.block(0,index,3,1) = 
		  	_inPtsWam.block(markerId,0,3,_nPts);
			posesMat.block(3,index,3,1) = 
		  	EulAnglesWrAt(index); /*ONLY WIRSt right now*/

		}
		return posesMat;
	}

	//set the joint limits into ceres 
	void setJointLimits(Problem* problem
		, double* thetas_i, int toDof)	{
		size_t  nJoints = _fk.getNJoints();
		for (int joint = 0; joint < nJoints 
			&& joint < toDof; joint++) {
			problem->SetParameterLowerBound(thetas_i, joint
				, _jointMinAngles[joint]);
			problem->SetParameterUpperBound(thetas_i, joint
				, _jointMaxAngles[joint]);
		}
	}




	
