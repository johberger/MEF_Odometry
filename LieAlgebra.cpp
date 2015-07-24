#include "LieAlgebra.hpp"


LieAlgebra::LieAlgebra() {

}


LieAlgebra::LieAlgebra(const Matrix4d& e, int order) {
	this->E = e;
	this->order = 1;
}


LieAlgebra::LieAlgebra(int order) {
	
	if (order < 1 ) {
		cerr << "Parameter order must be positive." << endl;
		order = 2;
	}
	this->order = order;
	this->E = Matrix4d::Zero();
	
	if (order >= 2) {
		this->V = VectorXd::Zero((order-1)*6);
	}
	
}


LieAlgebra::LieAlgebra(const Matrix4d& e, const VectorXd& v, int order) {
	if (v.rows()!= (order-1)*6) {
			cerr << "LieAlgebra(const Matrix4d& e, const VectorXd& v, int order): Dimensions do not coincide." << endl;
	}
	this->E = e;
	this->V = v;
	this->order = order;
}

LieAlgebra::LieAlgebra(const VectorXd& v, int order) {
	this->order = order;
	int r = v.rows();
	if (r%6 != 0) {
		cerr << "Wrong size of vector!" << endl;
		exit(1);
	} else if (r/6 != order) {
		cerr << "Order does not fit to size of vector!" << endl;
		exit(1);
	}
	E = MEFcommon::matSE(v.segment(0,6));
	V = v.tail(r-6);
}


LieAlgebra LieAlgebra::operator+(const LieAlgebra& rhs) {
	LieAlgebra temp(*this);
	temp.E += rhs.E;
	if (temp.V.rows()!=rhs.V.rows()){
		cerr << "operator+: dimenions do not coincide";
		exit(1);
	}
	temp.V += rhs.V;
	return temp;
}

LieAlgebra LieAlgebra::operator-(const LieAlgebra& rhs) {
	LieAlgebra temp(*this);
	temp.E -= rhs.E;
	if (temp.V.rows()!=rhs.V.rows()){
		cerr << "operator-: dimenions do not coincide";
		exit(1);
	}
	temp.V -= rhs.V;
	return temp;
}


LieAlgebra operator*(double scalar, const LieAlgebra& rhs){
	return LieAlgebra(scalar * rhs.E, scalar * rhs.V, rhs.order);
}


LieAlgebra::~LieAlgebra() {

}

LieGroup LieAlgebra::exp_g() const {
	 return LieGroup(MEFcommon::expSE(E),V,order);
}


VectorXd LieAlgebra::vec_g() const {;
	VectorXd r = VectorXd::Zero(6 * order);
	r.segment(0,6) = MEFcommon::vecSE(E);
		
	if (order > 1) {
		if ( V.rows() != (order-1)*6) {
			cerr << "vec_g: Dimensions do not coincide" << endl;
			return r;
		}
		r.segment(6, V.rows()) = V;
	} 
	return r;
}
