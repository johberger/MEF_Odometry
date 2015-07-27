#include "src/cpp/common.hpp"
#include "src/cpp/MinEnFilter.hpp"
#include <mex.h>

using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
// 	/* Check for proper number of arguments */
// 	if (nrhs != 2) {
// 		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
// 			"MEF_matlab_interface requires two input arguments.");
// 	} else if (nlhs >= 1) {
// 		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
// 			"MEF_matlab_interface requires no output argument.");
// 	}

	// read command
	char option[128];
	mxGetString(prhs[0],option,128);
	
	static MinEnFilter *F;
	
	if (strcmp(option,"init")==0) {
		// initialize filter
		cout << "Initialization" << endl;
		int arg;
		Matrix3d Camera;
		double q;
		double s1,s2;
		double alpha;
		int numSteps, order;
		
		arg = 1;
			Camera = Eigen::Map<MatrixXd>((double *) mxGetData(prhs[arg]),(int) mxGetM(prhs[arg]),(int) mxGetN(prhs[arg]));
		arg = 2;
			q =  mxGetScalar(prhs[arg]);
		arg = 3;
			s1 = mxGetScalar(prhs[arg]);
		arg = 4;
			s2 = mxGetScalar(prhs[arg]);
		arg = 5;	
			alpha = mxGetScalar(prhs[arg]);
		arg = 6;	
			numSteps = mxGetScalar(prhs[arg]);
		arg = 7;
			order = mxGetScalar(prhs[arg]);
		
		// initialize filter with constructor		
		F = new MinEnFilter(Camera, q, s1, s2, alpha, numSteps, order);
		

	} else if (strcmp(option,"iterate")==0) {
		// iterate filter
		int arg;
		double x,y, *result;
		
		MatrixXd XY;
		VectorXd Depth;
		MatrixXd Flow;
		
		
		arg = 1;
			XY = Eigen::Map<MatrixXd>((double *) mxGetData(prhs[arg]),(int) mxGetM(prhs[arg]),(int) mxGetN(prhs[arg]));
		arg = 2;
			Depth = Eigen::Map<MatrixXd>((double *) mxGetData(prhs[arg]),(int) mxGetM(prhs[arg]),(int) mxGetN(prhs[arg]));
		arg = 3;
			Flow = Eigen::Map<MatrixXd>((double *) mxGetData(prhs[arg]),(int) mxGetM(prhs[arg]),(int) mxGetN(prhs[arg]));
			
		Matrix4d E = F->iterate(XY, Depth, Flow);
  
		plhs[0] = mxCreateDoubleMatrix(4,4, mxREAL);
		result = (double*) mxGetData(plhs[0]);
		memcpy(result, E.data(), 4*4*sizeof(double));
		
			
	} else {
		cout << "Option \"" << option << "\" not available." << endl;
	}
	

  
//   mexcpp(*vin1, *vin2);
  return;
}