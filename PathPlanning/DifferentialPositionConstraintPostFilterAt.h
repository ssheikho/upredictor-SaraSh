struct DifferentialPositionConstraintPostFilterAt	{
	DifferentialPositionConstraint (
		ForwardKin<double> &fk
		// observed scaled Vectors (Cartesian Positions)
		// on the wam at t_i
		//wrist center owr = o7 − d7.R.[0 0 1]^T			
		, const Eigen::Matrix<double,3,1> shToUaVi
		, const Eigen::Matrix<double,3,1> elOffsetVi 
		, const Eigen::Matrix<double,3,1> wrOffsetVi
		, const Eigen::Matrix<double,3,1> laToWrVi
		, const Eigen::Matrix<double,3,1> shoToWrVi
		//inverse of covariance matrix (x y z)
		, const Eigen::Matrix<double, 3, 3>& iMatShToUaWam
		, const Eigen::Matrix<double, 3, 3>& iMatElOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatWrOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatLaToWrWam
		, const Eigen::Matrix<double, 3, 3>& iMatShToWrWam
		, Eigen::Quaternion<double> deltaQbase_i
		) : _fk(fk)
		,_a(_fk.getAvect(7))
		, _alpha(_fk.getAlphaVect(7))
		, _d(_fk.getDvect(7)) 
		,	_shToUaVi(shToUaVi)
		,	_elOffsetVi(elOffsetVi)
		,	_wrOffsetVi(wrOffsetVi)
		,	_laToWrVi(laToWrVi)
		,	_shoToWrVi(shoToWrVi)
		, _iMatShToUaWam(iMatShToUaWam)
		, _iMatElOffsetWam(iMatElOffsetWam)
		, _iMatWrOffsetWam(iMatWrOffsetWam)
		, _iMatLaToWrWam(iMatLaToWrWam)
		, _iMatShToWrWam(iMatShToWrWam)
	{}


		//oUA = o3 - a3.R3.[1 0 0]^T
	template <typename U>
	bool operator()(const U* const qBase_i
		, const U* const thetas_i
		, U* residuals) const {

		// (1). Form the homogeneous transformation matrices
		// 			A_i by substituting the parameters for each 
		//			link into DH-matrix of eqn. (3.10).
		//first 5 joints determine position constraints
		int toDof = 5;   
		//initialize the matrices to zero;
		//avoid calling external eigen::identity function
		Eigen::Matrix<U, 4, 4> tUa, tEl, tLa, tWr;
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				tUa(i,j) = U(0.0);	tEl(i,j) = U(0.0);
				tWr(i,j) = U(0.0); tLa(i,j) = U(0.0);
			}
			tUa(i,i) = U(1.0); tEl(i,i) = U(1.0);
				tWr(i,i) = U(1.0);	tLa(i,i) = U(1.0);
		}		

		Eigen::Matrix<U, 4, 4> A_j;
		for(int i = 0; i < toDof; i++) {
			//1st row 
			A_j(0,0) = 
				U(ceres::cos(thetas_i[i]));
			A_j(0,1) = U(-ceres::cos(_alpha(i,0)) 
				* ceres::sin(thetas_i[i]));
			A_j(0,2) = U(ceres::sin(_alpha(i,0)) 
				* ceres::sin(thetas_i[i]));
			A_j(0,3) = 
				U(_a(i,0) * ceres::cos(thetas_i[i]));
			//2nd row 
			A_j(1,0) = 
				U(ceres::sin(thetas_i[i]));
			A_j(1,1) = U(ceres::cos(_alpha(i,0)) 
				* ceres:: cos(thetas_i[i]));
			A_j(1,2) = U(-ceres::cos(thetas_i[i]) 
				* ceres::sin(_alpha(i,0)));
			A_j(1,3) = 
				U(_a(i,0) * ceres::sin(thetas_i[i]));
			//3rd row 
			A_j(2,0) = U(0.0);
			A_j(2,1) = U(ceres::sin(_alpha(i,0)));
			A_j(2,2) = U(ceres::cos(_alpha(i,0)));
			A_j(2,3) = U(_d(i,0));
			//4th row 
			A_j(3,0) = U(0.0);
			A_j(3,1) = U(0.0);
			A_j(3,2) = U(0.0);
			A_j(3,3) = U(1.0);

			
			// (2). Form T^n_0 = A 1 · · · A n . This then gives 
			//			the position (and orientation for wrist later)
			//			of the tool frame expressed in base .
			//			coordinates

			//	"A1 * A2 * vUA": -> Elbow center, pUa(i)				
			if (i < 2)
				tUa = tUa * A_j;
			//	"A1 * A2 * A3 * Origin": -> pEl(i)
			//	elbowOffset(i) = pUa(i)	- pEl(i)
			if (i < 3)
				tEl = tEl * A_j;
			if (i < 4)
				tLa = tLa * A_j;
			//	"A1 * A2 * ...* A5 * Origin": -> Wrist center
			//	pWr (center) = o5 = o6 OR (ALTERNATIVE)
			//							 = o4 + d5.R4.[0 0 1]^T			
			if (i < 5)
				tWr = tWr * A_j;
		}


		//(3).	estimate marker poitions w.r.t. moving base
		//	"θ1, θ2": -> shoToUa(i)				
		Eigen::Matrix<U, 4, 1> shToUaTrans;
		shToUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);
		Eigen::Matrix	<U, 3, 1> pUaEstimatedInBase = 
			(tUa * shToUaTrans
				).template block<3, 1>(0, 0);//.hnormalized();

		//	"θ1, θ2, θ3": -> elbowOffset(i)				
		Eigen::Matrix	<U, 3, 1> pElEstimatedInBase =
			tEl.template block<3, 1>(0, 3);
		//	"θ1, θ2, θ3, θ4": -> WrOffset(i)				
		Eigen::Matrix	<U, 3, 1> pLaEstimatedInBase =
			tLa.template block<3, 1>(0, 3);
		//	"θ1, θ2, θ3 θ4": -> wrist CENTER pWr(i)
		Eigen::Matrix	<U, 3, 1> pWrEstimatedInBase =
			tWr.template block<3, 1>(0, 3);

		//(4).	represent estimated marker poitions in
		//			the camera frame
		//	void QuaternionRotatePoint(
		//		const T q[4], const T pt[3], T result[3])
    // We use QuaternionRotatePoint as it does not assume 	
		// that the quaternion is normalized, since one of the 
    // ways to run the bundle adjuster is to let Ceres 
    // optimize all 4 quaternion parameters without a
    //  local parameterization.

		//	"shToUa(i)"
		Eigen::Matrix	<U, 3, 1> shToUaViEstimated;
    QuaternionRotatePoint(qBase_i,
			pUaEstimatedInBase.data(), shToUaViEstimated.data());
		//	"shToEl(i)"
		Eigen::Matrix	<U, 3, 1> shToElViEstimated;
    QuaternionRotatePoint(qBase_i,
			pElEstimatedInBase.data(), shToElViEstimated.data());
		//	"shToLa(i)"
		Eigen::Matrix	<U, 3, 1> shToLaViEstimated;
    QuaternionRotatePoint(qBase_i,
			pLaEstimatedInBase.data(), shToLaViEstimated.data());
		//	"shToWr(i)"
		Eigen::Matrix	<U, 3, 1> shToWrViEstimated;
    QuaternionRotatePoint(qBase_i,
			pWrEstimatedInBase.data(), shToWrViEstimated.data());

    //(5).	Compute the residuals.
    // [ position         ]   [ delta_p        ]
    Eigen::Map<Eigen::Matrix<U, 12, 1> >
			residualsMat(residuals);

		//	"shToUa(i)"
		residualsMat.template block<3, 1>(0, 0) =
			shToUaViEstimated - _shToUaVi.template cast<U>();
		//elbowOffset(i)		
		residualsMat.template block<3, 1>(3, 0) =
			shToElViEstimated - shToUaViEstimated -
					_elOffsetVi.template cast<U>();

		//WristOffset(i)		
		residualsMat.template block<3, 1>(6, 0) =
			shToLaViEstimated - shToElViEstimated -
					_wrOffsetVi.template cast<U>();

		//	"LaToWr(i)"
		residualsMat.template block<3, 1>(9, 0) = 
			shToWrViEstimated - shToLaViEstimated - 
				_laToWrVi.template cast<U>();

		//(5). 	Scale the residuals by the
		//			measurement uncertainty.

		//	"iMatShToUa(i)"
    residualsMat.template block<3, 1>(0, 0).applyOnTheLeft(
			_iMatShToUaWam.template cast<U>());

		//	"iMatElOffset(i)"
    residualsMat.template block<3, 1>(3, 0).applyOnTheLeft(
			_iMatElOffsetWam.template cast<U>());

		//	"iMatWrOffset(i)"
    residualsMat.template block<3, 1>(6, 0).applyOnTheLeft(
			_iMatWrOffsetWam.template cast<U>());

		//	"iMatLaToWr(i)"
		residualsMat.template block<3, 1>(9, 0).applyOnTheLeft(
			_iMatLaToWrWam.template cast<U>()); 

    return true;
	}

	static ceres::CostFunction* Create(ForwardKin<double> &fk
		, const Eigen::Matrix<double,3,1> shToUaVi
		, const Eigen::Matrix<double,3,1> elOffsetVi 
		, const Eigen::Matrix<double,3,1> wrOffsetVi
		, const Eigen::Matrix<double,3,1> laToWrVi
		, const Eigen::Matrix<double,3,1> shoToWrVi
		//inverse of covariance matrix (x y z)
		, const Eigen::Matrix<double, 3, 3>& iMatShToUaWam
		, const Eigen::Matrix<double, 3, 3>& iMatElOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatWrOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatLaToWrWam
		, const Eigen::Matrix<double, 3, 3>& iMatShToWrWam)	{

		return (new ceres::AutoDiffCostFunction
			<DifferentialPositionConstraint, 12, 4, 7> (
				new DifferentialPositionConstraint (fk, shToUaVi,
					elOffsetVi, wrOffsetVi,	laToWrVi, shoToWrVi,
					iMatShToUaWam, iMatElOffsetWam, iMatWrOffsetWam, 						iMatLaToWrWam, iMatShToWrWam)));
	}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	ForwardKin<double> &_fk;
	Eigen::Matrix<double, Eigen::Dynamic, 1> 
	 	_a, _alpha, _d;
	//measured/scaled Cartesian vectores on the wam at t_i
	const Eigen::Matrix<double, 3, 1> _shToUaVi, _elOffsetVi,
		_wrOffsetVi,	_laToWrVi, _shoToWrVi;
	//measurements covariance (x y z)
	const Eigen::Matrix<double, 3, 3> _iMatShToUaWam,
		_iMatElOffsetWam, _iMatWrOffsetWam,
	  _iMatLaToWrWam, _iMatShToWrWam;
};


