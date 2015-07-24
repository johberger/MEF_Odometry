#include "MinEnFilter.hpp"

#include <eigen3/Eigen/LU>
// #include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Cholesky>
#include <limits>


MinEnFilter::MinEnFilter() {

}

MinEnFilter::~MinEnFilter() {

}



MinEnFilter::MinEnFilter( Matrix3d& Camera, double q, double s1, double s2, double alpha, double numSteps, int order) {
	this->camera = Camera;
	this->q1 = q / numSteps;
	this->q2 = q / numSteps;
	this->s1 = s1;
	this->s2 = s2;
	this->alpha = alpha;
	this->numSteps  = numSteps;
	this->nPoints = 0;
	delta = 1/numSteps;
	this->order = order;
		
	VectorXd w(6);
	w << s1,s1,s1,s2,s2,s2;
	MatrixXd invW_temp = (w.cwiseInverse()).asDiagonal();
	
	cout << "\n\n delta = " << delta << "\n\n";
	cout << "\n\n q1 = " << q1 << "\n\n";
	cout << "\n\n q2 = " << q2 << "\n\n";
	cout << "\n\n s1 = " << s1 << "\n\n";
	cout << "\n\n s2 = " << s2 << "\n\n";
	cout << "\n\n numSteps = " << numSteps << "\n\n";
	cout << "\n\n alpha = " << alpha << "\n\n";	
	
	if (order == 1) {
		currG = LieGroup(MatrixXd::Identity(4,4));
		currP = MatrixXd::Identity(6,6);
		invW = invW_temp;
	} else if (order >= 2) {
		currG = LieGroup(MatrixXd::Identity(4,4),VectorXd::Zero((order-1)*6), order);
		currP = MatrixXd::Identity(6*order, 6*order);
		invW = MatrixXd::Zero(6*order, 6*order);
		for (size_t i=0; i<order; i++) {
			invW.block(i*6,i*6,6,6) = invW_temp;
		}
	} else {
		cerr << "Parameter order must be positive." << endl;
		exit(0);
	}
	
	
}


Matrix4d MinEnFilter::iterate(MatrixXd& XY, VectorXd& Depth, MatrixXd& Flow) {
	cout << "Iterate:" << endl;
	
	LieAlgebra grad;
	MatrixXd Hessian;
	
	// check dimensions
	nPoints = XY.rows();
	
	if (nPoints!= Depth.rows()) {
		cerr << "Dimensions of XY do not coincide with Dimensions of Depth." << endl;
	}
	if (nPoints!= Flow.rows()) {
		cerr << "Dimensions of XY do not coincide with Dimensions of Flow." << endl;
	}
	
	MatrixXd XYZ_inhom  = inhomCoordinates(XY, false);		// inhomogenous Coordinatestes
	MatrixXd Flow_inhom = inhomCoordinates(Flow, true);		// inhomogenous (relative) flow
	MatrixXd Disps_inhom = XYZ_inhom + Flow_inhom;			// inhomogenous Disparities
	MatrixXd temp(nPoints, 4);
	
	temp << XYZ_inhom, Depth.cwiseInverse();
	MatrixXd Gk = (Depth*MatrixXd::Ones(1,4)).cwiseProduct(temp);

	grad = computeGradient(currG , Gk, Disps_inhom);
// 	cout << "\n\n grad.E: \n\n" << grad.E << endl;
// 	cout << "\n\n grad.V: \n\n" << grad.V << endl;
	
	Hessian = computeHessian(Gk, Disps_inhom);
// 	cout << "\n\n Hessian: \n\n" << Hessian << endl;
	

	
	currG = integrateGImplicit(grad, Hessian, Gk, Disps_inhom);
	cout << "\n\n Result of integrateGImplicit: \n\n" << currG.E << endl;
	
	currP = integratePExplicit(grad,Hessian);
	cout << "\n\n Result of integratePExplicit: \n\n" << currP << endl;
	
	return currG.E;
	
}



MatrixXd MinEnFilter::inhomCoordinates(const MatrixXd& homCoordinates2D, bool isFlow) {
	MatrixXd homCoordinates3d;
	if (isFlow) {
		homCoordinates3d = MatrixXd::Zero(nPoints,3);
	} else {
		homCoordinates3d = MatrixXd::Ones(nPoints,3);
	}
	
	homCoordinates3d.leftCols(2) = homCoordinates2D;
	MatrixXd inhomCoordiantes;
	inhomCoordiantes = camera.lu().solve(homCoordinates3d.transpose());
	
	return inhomCoordiantes.transpose();
}

LieAlgebra MinEnFilter::computeGradient(const LieGroup& S, const MatrixXd& Gk, const MatrixXd& Disps_inhom) {
	/* DO NOT CHANGE THIS FUNCTION */
	LieAlgebra L;
	MatrixXd iEg = ((S.E).lu().solve(Gk.transpose())).transpose();
	VectorXd kappa = iEg.col(2);
	MatrixXd h = iEg.leftCols(2).cwiseQuotient(kappa*MatrixXd::Ones(1,2));
	MatrixXd z = Disps_inhom.leftCols(2) - h;
	MatrixXd Qvec(1,2);
	Qvec << q1,q2;
	MatrixXd qvec = MatrixXd::Ones(nPoints,1) * Qvec;
	MatrixXd b12 = (kappa.cwiseInverse()*MatrixXd::Ones(1,2)).cwiseProduct(qvec.cwiseProduct(z));
	MatrixXd kappa2 = kappa.cwiseProduct(kappa);
	MatrixXd aux1 = (iEg.leftCols(2)).cwiseProduct(z.cwiseProduct(qvec));
	MatrixXd b3 = -((kappa2.cwiseInverse()*MatrixXd::Ones(1,2)).cwiseProduct(aux1)).rowwise().sum();
	MatrixXd B(nPoints,4);
	B << b12, b3, MatrixXd::Zero(nPoints,1);
	MatrixXd Ones = MatrixXd::Ones(1,4);
	MatrixXd A1 = iEg.col(0) * Ones;
	MatrixXd A2 = iEg.col(1) * Ones;
	MatrixXd A3 = iEg.col(2) * Ones;
	MatrixXd A4 = iEg.col(3) * Ones;
	
	MatrixXd G1 = (B.cwiseProduct(A1)).colwise().sum();
	MatrixXd G2 = (B.cwiseProduct(A2)).colwise().sum();
	MatrixXd G3 = (B.cwiseProduct(A3)).colwise().sum();
	MatrixXd G4 = (B.cwiseProduct(A4)).colwise().sum();

	Matrix4d G;
	G << 	G1, 
			G2,
			G3,
			G4;

	G = MEFcommon::prSE(G.transpose());
	if (order == 1) {
		L = LieAlgebra(G);
	} else if (order >=2) {
		L = LieAlgebra(G, VectorXd::Zero((order-1)*6),order);
	} else {
		cerr << "Parameter order must be greater than zero." << endl;
		exit(1);
	}
	
	return L;
}

MatrixXd MinEnFilter::DAk(const MatrixXd& Gk, const MatrixXd& Disps_inhom, int direction) {
	MatrixXd Ones = MatrixXd::Ones(1,4);
	MatrixXd Eaux = MatrixXd::Identity(6,6);
	Matrix4d eta = MEFcommon::matSE(Eaux.col(direction-1));
	MatrixXd iEg = (currG.E.lu().solve(Gk.transpose())).transpose();
	VectorXd kappa = iEg.col(2);
	MatrixXd h = iEg.leftCols(2).cwiseQuotient(kappa*MatrixXd::Ones(1,2));
	MatrixXd z = Disps_inhom.leftCols(2) - h;
	VectorXd kappaM1 = kappa.cwiseInverse();
	VectorXd kappaM2 = (kappa.cwiseProduct(kappa)).cwiseInverse();
	VectorXd kappaM3 = (kappa.cwiseProduct(kappa.cwiseProduct(kappa))).cwiseInverse();	
	MatrixXd lambda = ((currG.E.lu().solve(eta))*(iEg.transpose())).transpose();
	VectorXd xi = (lambda.col(2)).cwiseProduct(kappaM2);
	
	// first part of Hessian
	VectorXd zeta1 = -2* kappaM3.cwiseProduct((lambda.col(2)).cwiseProduct(iEg.col(0)));
	VectorXd zeta2 = -2 * kappaM3.cwiseProduct((lambda.col(2)).cwiseProduct(iEg.col(1)));
	VectorXd eta1 = kappaM2.cwiseProduct(lambda.col(0));
	VectorXd eta2 = kappaM2.cwiseProduct(lambda.col(1));
	VectorXd a31 = zeta1 + eta1;
	VectorXd a32 = zeta2 + eta2;
	MatrixXd b1row = q1 * ((z.col(0))*Ones).cwiseProduct(iEg);
	MatrixXd b2row = q2 * ((z.col(1))*Ones).cwiseProduct(iEg);	
	MatrixXd c1row = ((xi * Ones).cwiseProduct(b1row)).colwise().sum();
	MatrixXd c2row = ((xi * Ones).cwiseProduct(b2row)).colwise().sum();
	MatrixXd aux1(nPoints,4);
	aux1 << (a31.cwiseProduct(b1row.col(0)) + a32.cwiseProduct(b2row.col(0))) , (a31.cwiseProduct(b1row.col(1)) + a32.cwiseProduct(b2row.col(1))), (a31.cwiseProduct(b1row.col(2)) + a32.cwiseProduct(b2row.col(2))), (a31.cwiseProduct(b1row.col(3)) + a32.cwiseProduct(b2row.col(3)));
	MatrixXd c3row = aux1.colwise().sum();
	MatrixXd C(4,4);
	C << 	c1row,
			c2row,
			c3row,
			MatrixXd::Zero(1,4);
	
	// % second part of Hessian
	VectorXd rho1 = -q1 * kappaM2.cwiseProduct(iEg.col(0));
	VectorXd rho2 = -q2 * kappaM2.cwiseProduct(iEg.col(1));
	MatrixXd Frow1 = kappaM1.cwiseProduct(lambda.col(0)) - kappaM2.cwiseProduct((lambda.col(2)).cwiseProduct(iEg.col(0)));
	MatrixXd Frow2 = kappaM1.cwiseProduct(lambda.col(1)) - kappaM2.cwiseProduct((lambda.col(2)).cwiseProduct(iEg.col(1)));
	MatrixXd G1row1 = (Frow1*Ones).cwiseProduct(iEg);
	MatrixXd G1row2 = (Frow2*Ones).cwiseProduct(iEg);
	MatrixXd G2row1 = ((z.col(0))*Ones).cwiseProduct(lambda);
	MatrixXd G2row2 = ((z.col(1))*Ones).cwiseProduct(lambda);
	MatrixXd Grow1 = G1row1 - G2row1;
	MatrixXd Grow2 = G1row2 - G2row2;
	MatrixXd h1row = (q1*(kappaM1*Ones).cwiseProduct(Grow1)).colwise().sum();
	MatrixXd h2row = (q2 * (kappaM1*Ones).cwiseProduct(Grow2)).colwise().sum();
	MatrixXd aux2(nPoints,4);
	aux2 << rho1.cwiseProduct(Grow1.col(0)) + rho2.cwiseProduct(Grow2.col(0)), rho1.cwiseProduct(Grow1.col(1)) + rho2.cwiseProduct(Grow2.col(1)), rho1.cwiseProduct(Grow1.col(2)) + rho2.cwiseProduct(Grow2.col(2)), rho1.cwiseProduct(Grow1.col(3)) + rho2.cwiseProduct(Grow2.col(3));
	MatrixXd h3row = aux2.colwise().sum();
	MatrixXd H(4,4);
	H << 	h1row,
			h2row,
			h3row,
			MatrixXd::Zero(1,4);
			
	return C+H;
}


MatrixXd MinEnFilter::computeHessian(const MatrixXd& Gk, const MatrixXd& Disps_inhom) {
	VectorXd h1 = MEFcommon::vecSE(MEFcommon::prSE(DAk(Gk, Disps_inhom, 1)));
	VectorXd h2 = MEFcommon::vecSE(MEFcommon::prSE(DAk(Gk, Disps_inhom, 2)));
	VectorXd h3 = MEFcommon::vecSE(MEFcommon::prSE(DAk(Gk, Disps_inhom, 3)));
	VectorXd h4 = MEFcommon::vecSE(MEFcommon::prSE(DAk(Gk, Disps_inhom, 4)));
	VectorXd h5 = MEFcommon::vecSE(MEFcommon::prSE(DAk(Gk, Disps_inhom, 5)));
	VectorXd h6 = MEFcommon::vecSE(MEFcommon::prSE(DAk(Gk, Disps_inhom, 6)));

	MatrixXd H(6,6);
	H << h1, h2, h3, h4, h5, h6;
	return H;	
}





MatrixXd MinEnFilter::integratePImplicit(const LieAlgebra & grad, const MatrixXd & Hessian) {
	
	const double eps = 1e-8;
	int dim = order * 6;
	
	MatrixXd A(dim,dim);
	MatrixXd D = MatrixXd::Zero(dim,dim);
	MatrixXd Gamma = MatrixXd::Zero(dim,dim);
	MatrixXd Id = MatrixXd::Identity(dim,dim);
	LieAlgebra FG;
	MatrixXd B(6,6);
	MatrixXd C = invW;

	// define function that applies a function to each row of a matrix
	if (order == 1) {
		A = -tildeGammaAst(grad.E);
// 		cout << "\n\n Matrix A: \n\n" << A << endl;
// 		cout << "\n\n tildeGamma(grad.E): \n\n" << tildeGamma(grad.E) << endl;
		D = (tildeGamma(grad.E) + Hessian);
// 		cout << "\n\n Matrix D: \n\n" << D << endl;
	} else if ( order>=2 ) {
		Gamma.topLeftCorner(6,6) = tildeGammaAst(grad.E);
		FG = MotionPrior(currG);
		B = MEFcommon::adjoint(FG.E); // kronSE(FG.e.mat, eye(4)) - kronSE(eye(4),FG.e.mat);
		A.topLeftCorner(6,6) = -B;
		A.block(0, 5, dim-6 ,dim-6) = MatrixXd::Identity(dim-6,dim-6);
		A = A-Gamma;
		D.topLeftCorner(6,6) = tildeGamma(grad.E) + Hessian;
	} else {
		cerr << "Parameter order must be positive" << endl;
		exit(1);
	}

	// ensure that D is symmetric
	D = 0.5 * (D + D.transpose());
// 	// ensure that D is positive definite
// 	[V,eigs]=eig(D);
// 	d=diag(eigs);
// 	if ( min(eigs)<eps) {
// 		d(d<=eps)=eps;
// 		D = V*diag(d)*V';
// 	}  
	
	MatrixXd tildeQ = currP + delta * C;
// 	cout << "\n\n tildeQ \n\n " << tildeQ << endl;
	MatrixXd tildeA = (delta * A - (1+delta*alpha)/2 * Id).transpose();
// 	cout << "\n\n tildeA  \n\n " << tildeA << endl;	
	
	Eigen::LLT<MatrixXd> lltOfD(D);
	MatrixXd cholD = lltOfD.matrixL();
	MatrixXd tildeB = sqrt(delta)*cholD;
// 	cout << "\n\n tildeB  \n\n " << tildeB << endl;	

   // solve algebraic Riccati equation to find solution of implicit Euler scheme, i.e.
   return care(tildeA, tildeB, tildeQ); 

}

MatrixXd MinEnFilter::integratePExplicit(const LieAlgebra& grad, const MatrixXd& Hessian) {
	const double eps = 1e-8;
	int dim = order * 6;
	
	MatrixXd A = MatrixXd::Zero(dim,dim);
	MatrixXd D = MatrixXd::Zero(dim,dim);
	MatrixXd Gamma = MatrixXd::Zero(dim,dim);
	MatrixXd Id = MatrixXd::Identity(dim,dim);
	LieAlgebra FG;
	MatrixXd B(6,6);
	MatrixXd C = invW;
	MatrixXd R;

	// define function that applies a function to each row of a matrix
	if (order == 1) {
		A = -tildeGammaAst(grad.E);
// 		cout << "\n\n Matrix A: \n\n" << A << endl;
// 		cout << "\n\n tildeGamma(grad.E): \n\n" << tildeGamma(grad.E) << endl;
		D = (tildeGamma(grad.E) + Hessian);
// 		cout << "\n\n Matrix D: \n\n" << D << endl;
	} else if ( order>=2 ) {
		Gamma.topLeftCorner(6,6) = tildeGammaAst(grad.E);
// 		cout << "\n\n Gamma= \n\n" << Gamma << endl;
		FG = MotionPrior(currG);
// 		cout << "\n\n FG.E= \n\n" << FG.E << endl;
// 		cout << "\n\n FG.V= \n\n" << FG.V << endl;
		B = MEFcommon::adjoint(FG.E); // kronSE(FG.e.mat, eye(4)) - kronSE(eye(4),FG.e.mat);
// 		cout << "\n\n B= \n\n" << B << endl;
		A.topLeftCorner(6,6) = -B;
		A.block(0, 6, dim-6 ,dim-6) = MatrixXd::Identity(dim-6,dim-6);
		A = A-Gamma;
// 		cout << "\n\n A= \n\n" << A << endl;
		D.topLeftCorner(6,6) = tildeGamma(grad.E) + Hessian;
// 		cout << "\n\n D= \n\n" << D << endl;
	} else {
		cerr << "Parameter order must be positive" << endl;
		exit(1);
	}
	
	R = currP + delta * ( -alpha * currP + C + A * currP + currP * A.transpose() - currP * D * currP);
	
	return R;
	
}

LieGroup MinEnFilter::integrateGImplicit(const LieAlgebra& grad, const MatrixXd& Hessian, const MatrixXd& Gk, const MatrixXd& Disps_inhom ) {
	
	LieAlgebra K = LieAlgebra(order);
	LieGroup R;
	LieAlgebra S;
	LieGroup U;
	LieGroup M1;
	
	for (size_t i=0; i < NUM_FIXEDPOINT_ITERATIONR; i++) {
		U = (0.5 * K).exp_g();
		S = computeGradient(currG * U, Gk, Disps_inhom);
		K = delta * dynamicsE(S, currG);
	}
	return currG * K.exp_g();
}


LieAlgebra MinEnFilter::dynamicsE(const LieAlgebra& grad, const LieGroup& G) {

	VectorXd v = grad.vec_g();
	VectorXd v1 = currP * v;
	LieAlgebra A = MotionPrior(G);
	LieAlgebra B = LieAlgebra(v1, order);
	return  A-B;
}


LieAlgebra MinEnFilter::MotionPrior(const LieGroup & G) {
	if (G.order > 1) {
		VectorXd v = VectorXd::Zero(G.order*6);
		v.segment(0,6) = G.V.segment(0,6);
		return LieAlgebra(v, G.order);
	} else {
		return LieAlgebra(1);
	}
}


MatrixXd MinEnFilter::tildeGamma(const Matrix4d& M) {
	VectorXd eta = MEFcommon::vecSE(M);
	VectorXd rows(12);
	VectorXd cols(12);
	VectorXd vals(12);
	VectorXd r(36);
	MatrixXd G = MatrixXd::Zero(36,6);
	MatrixXd R(6,6);
	
	rows << 3, 14, 7, 9, 2, 13, 6, 10, 17, 12, 16, 5; 	// Important: Matlab indices
	cols << 2, 1, 3, 1, 3, 2, 5, 6, 4, 4, 5, 6;			// Important: Matlab indices
	vals << 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, 1, 1, 1, -1, -1, -1;
	
	for (size_t i = 0; i<12; i++) {
		G(rows(i)-1, cols(i)-1) = vals(i);
	}
	R << G*eta;
	return R;
}


MatrixXd MinEnFilter::tildeGammaAst(const Matrix4d& M) {
	VectorXd eta = MEFcommon::vecSE(M);
	VectorXd rows(12);
	VectorXd cols(12);
	VectorXd vals(12);
// 	VectorXd r(36);
	MatrixXd G = MatrixXd::Zero(36,6);
	MatrixXd R(6,6);
	
	rows << 9, 2, 13, 3, 14, 7, 30, 34, 23, 24, 28, 35;	// Important: Matlab indices
	cols << 1, 3, 2, 2, 1, 3, 1, 2, 3, 2, 3, 1;			// Important: Matlab indices
	vals << 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, 1., 1., 1., -1., -1., -1.;
	
	for (size_t i = 0; i<12; i++) {
		G(rows(i)-1, cols(i)-1) = vals(i);
	}
	R << G*eta; 
	return R;
}

MatrixXd MinEnFilter::care(const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q) {
	
	MatrixXd X = MatrixXd::Zero(6*order, 6*order);

	return X;
}




