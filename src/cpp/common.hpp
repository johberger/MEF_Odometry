#ifndef COMMON
#define COMMON

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, 6,1> Vector6d;
typedef Eigen::Matrix<double, 4,1> Vector4d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;
using Eigen:: VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix;
using namespace std;


// Function Definitions

namespace MEFcommon{
	Vector6d vecSE(const Matrix4d & m);
	Matrix4d matSE(const Vector6d & v);
	Matrix4d expSE(const Matrix4d & m);
	Matrix4d logSE(const Matrix4d & m);
	Matrix3d expSO(const Matrix3d & m);
	Matrix3d logSO(const Matrix3d & m);
	Vector3d vecSO(const Matrix3d & m);
	Matrix3d matSO(const Vector3d & v);
	Vector3d vex(const Matrix3d & m);
	Matrix3d mex(const Vector3d & v);
	Matrix4d prSE(const Matrix4d & m);
	Matrix<double,6,6> adjoint(const Matrix4d & m);
}


#endif