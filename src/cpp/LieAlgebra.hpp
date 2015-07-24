#ifndef LIEALGEBRA
#define LIEALGEBRA

#include "common.hpp"

struct LieAlgebra;
#include "LieGroup.hpp"


struct LieAlgebra {
// class LieAlgebra {
  
public:
	LieAlgebra(const Matrix4d& e, const VectorXd& v, int order);
	LieAlgebra(const VectorXd& v, int order);
	LieAlgebra(const Matrix4d& e, int order=1);
	LieAlgebra(int order);
	LieAlgebra();
	~LieAlgebra();

	LieAlgebra operator+(const LieAlgebra& rhs); 
	LieAlgebra operator-(const LieAlgebra& rhs);
	
	friend LieAlgebra operator *(double scalar, const LieAlgebra& rhs);

	LieGroup exp_g() const;
	VectorXd vec_g() const;
	Matrix4d E; 
	VectorXd V;
	int order;
  	
};

#endif // LieAlgebra