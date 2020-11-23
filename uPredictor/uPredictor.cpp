

#include <Ellipse3D.h>
#include <Plane.h>
#include <UBCUtil.h>
#include <Plane.h>
#include <LinearAlgebra.h>
#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "ForwardKin.h"
#include "DikProblem.h"
#include "DikSolver.h"


#include "ParseMathematica.h"
#include "ParseCSV.h"
#include <Eigen/SVD>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>       /* pow */
#include "FitFunctions.h"
#include "LinearAlgebra.cpp"
#include "RigidTrans2D.h"

#include "BaysFitFunctions.h"
#include "EllipseConicConstraints.h"
using namespace ceres;
using namespace Eigen;
using namespace std;

#include <iostream>

using namespace std;

class IntFactory : public TypeFactory<int> {
public:
	int *getInstance() {
		return new int();
	}
};

int main()
{
	IntFactory intF;
	SimplePool<int> sp(intF, 10);
	{
		SimplePool<int>::PooledPtr aP = sp.getInstance();
		cout << sp.poolSize() << endl;
		vector<SimplePool<int>::PooledPtr> v;
		v.push_back(aP);
		cout << sp.poolSize() << endl;
	}
	cout << sp.poolSize() << endl;

	getchar();

    return 0;
}

