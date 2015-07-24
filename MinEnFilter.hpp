#ifndef MINENFILTER
#define MINENFILTER

#include "common.hpp"
#include "LieGroup.hpp"
#include "LieAlgebra.hpp"

class MinEnFilter {
  
private:
	static const size_t NUM_FIXEDPOINT_ITERATIONR = 4;
	int order; // order of kinematic model: order=0: constant velocity, order=1: constant acceleration, order=2: constant jerk, etc.
	int numSteps;
	double delta; // 
	double alpha, s1, s2, q1, q2;
	int nPoints;
	Matrix3d camera;
	LieGroup currG; // current state of Lie Group
	MatrixXd currP; // current second order information matirx
	MatrixXd invW;
  
	
	MatrixXd inhomCoordinates(const MatrixXd& homCoordinates2D, bool isFlow);
	LieAlgebra computeGradient(const LieGroup& E, const MatrixXd& Gk, const MatrixXd& Disps_inhom);
	MatrixXd DAk (const MatrixXd& Gk,const MatrixXd& Disps_inhom, int direction);
	MatrixXd computeHessian(const MatrixXd& Gk, const MatrixXd& Disps_inhom);
	MatrixXd integratePImplicit(const LieAlgebra & grad, const MatrixXd& Hessian);
	MatrixXd integratePExplicit(const LieAlgebra & grad, const MatrixXd& Hessian);
	LieGroup integrateGImplicit(const LieAlgebra& grad, const MatrixXd& Hessian, const MatrixXd& Gk, const MatrixXd& Disps_inhom);
	LieGroup integrateGExplicit(const LieAlgebra& grad, const MatrixXd& Hessian, const MatrixXd& Gk, const MatrixXd& Disps_inhom);
	LieAlgebra dynamicsE(const LieAlgebra & grad, const LieGroup & G);
	
	MatrixXd tildeGammaAst(const Matrix4d& M);
	MatrixXd tildeGamma(const Matrix4d& M);
	LieAlgebra MotionPrior(const LieGroup & G);
	
	/** care - c++ clone of the matlab function care. Solves the continuous time algebraic Riccati equation
	 * 
	 * 
	 */
	MatrixXd care(const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q);
  
public:
	MinEnFilter();
	~MinEnFilter();
	
	/** MinEnFilter Constructor. There are the following paramters
	 * @param Camera: 3x3 matrix with internal camera parameters
	 * @param q: weighting parameter (>0) for the data term
	 * @param s1: weighting parameter for latent variable: rotational component
	 * @param s2: weighting parameter for latent variable: translational component
	 * @param alpha: decay rate (>=0)
	 * @param numSteps: number of integration numSteps
	 * @param order: desired order of filter, i.e. 1: constant velocity, 2: constant acceleration, 3: constant jerk, etc.
	 */
	MinEnFilter(Matrix3d &Camera, double q, double s1, double s2, double alpha, double numSteps, int order);
	
	/** iterates filter from time step t_{k} to time step t_{k+1}
	* @param XY: matrix with n rows and 2 columns with XY positions
	* @param Depth: column vector with depth values
	* @param Flow: matrix with n rows and 2 columns for flow in X-Y-direction
	*/
	Matrix4d iterate(MatrixXd& XY, VectorXd& Depth, MatrixXd& Flow);
	
};



#endif