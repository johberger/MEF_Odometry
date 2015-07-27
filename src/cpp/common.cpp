#include "common.hpp"

Vector6d MEFcommon::vecSE(const Matrix4d& m) {
  double s2 = std::sqrt(2);
  Vector6d v(6);
  
  v << s2*m(2,1), s2*m(0,2), s2*m(1,0), m(0,3), m(1,3), m(2,3);
  
  return v;
}


Matrix4d MEFcommon::matSE(const Vector6d& v) {
	double s2 = std::sqrt(2);
	Matrix4d m;
  
	m << 	0, -v(2)/s2, v(1)/s2, v(3),
			v(2)/s2, 0, -v(0)/s2, v(4),
			-v(1)/s2, v(0)/s2, 0, v(5),
			0, 0, 0, 0;
  return m;

}

Vector3d MEFcommon::vecSO(const Matrix3d& m) {
	return std::sqrt(2) * vex(m);
}


Matrix3d MEFcommon::matSO(const Vector3d& v) {
	return mex(v) / sqrt(2);
}

Vector3d MEFcommon::vex(const Matrix3d& m) {
	Vector3d v;
	v <<  m(2,1),m(0,2), m(1,0); 
	return v;
}

Matrix3d MEFcommon::mex(const Vector3d& v) {
	Matrix3d m;
	m <<	0, -v(2), v(1),
			v(2), 0, -v(0),
			-v(1), v(0), 0;
	return m;
}


Matrix4d MEFcommon::expSE(const Matrix4d& m) {
	Vector3d v = m.block(0,3,3,1);
	Matrix3d R = m.block(0,0,3,3);
  
	if (vex(R).norm() <= 1e-20) {
		Matrix4d M = Matrix4d::Identity();
		M.block(0,3,3,1) = v;
		return M;
	}
	Vector3d omega = vex(R);
	double norm_o = omega.norm();
	
	Matrix3d aux = Matrix3d::Identity(3,3) + (1-std::cos(norm_o))/(std::pow(norm_o,2)) * R + (norm_o - std::sin(norm_o))/(std::pow(norm_o,3)) * R * R;
	Matrix4d result = Matrix4d::Identity(4,4);
	
	result.block(0,3,3,1) = aux*v;
	result.block(0,0,3,3) = expSO(R);
	
	return result;
}

Matrix3d MEFcommon::expSO(const Matrix3d& m) {
	// Exponential map on SO(3) with Rodrigues's formula
	Vector3d omega = vex(m);
	double norm_o = omega.norm();
	Matrix3d result =  Matrix3d::Identity(3,3) + (std::sin(norm_o)/norm_o)*m + ((1-std::cos(norm_o))/std::pow(norm_o,2)) * m * m;
	return result;
}

Matrix3d MEFcommon::logSO(const Matrix3d& R) {
	double theta = std::acos(0.5*(R.trace()-1));
	cout << "R*R' \n" << R*R.transpose() << endl;
	Matrix3d result = (1/(2 * std::sin(theta))) * (R - R.transpose());
	return theta*result;
}

Matrix4d MEFcommon::logSE(const Matrix4d& m) {
	Vector3d v = m.block(0,3,3,1);
	Matrix3d R = m.block(0,0,3,3);
	Matrix3d Q = logSO(R);
	Vector3d omega = vex(Q);
	double norm_o = omega.norm();
	
	if (norm_o < 1e-20) {
		return Matrix4d::Identity(4,4);
	}
	double aux1 = (2 * std::sin(norm_o) - norm_o * (1 + std::cos(norm_o))) / (1 * norm_o * norm_o * std::sin(norm_o));
	Matrix3d aux2 = Matrix3d::Identity() - 0.5 * Q + aux1 * Q * Q;
	Matrix4d result = Matrix4d::Zero();
	result.block(0,3,3,1) = aux2*v;
	result.block(0,0,3,3) = Q;
	
	return result;
}

Matrix4d MEFcommon::prSE(const Matrix4d& m) {

	Matrix3d R = m.block(0,0,3,3);	
	Matrix3d Q = 0.5 * (R - R.transpose());
	Matrix4d result = Matrix4d::Zero();
	result.block(0,0,3,3) = Q;
	result.block(0,3,3,1) = m.block(0,3,3,1);

	return result;	
}

Matrix< double, 6, 6 > MEFcommon::adjoint(const Matrix4d& m) {
	VectorXd v = vecSE(m);
	MatrixXd R = MatrixXd::Zero(6,6);
	R.block(0,0,3,3) = matSO(v.segment(0,3));
	R.block(3,0,3,3) = matSO(v.segment(3,3));
	R.block(3,3,3,3) = matSO(v.segment(0,3));
	
	return R;
}




