#include "LieGroup.hpp"

LieGroup::LieGroup() {
  
}

LieGroup::LieGroup(const Matrix4d & E) {
	this->E = E;
	this->order = 1;
}


LieGroup::LieGroup(const Matrix4d & E, const VectorXd &  V, int order){
	this->E = E;
	this->V = V;
	this->order = order;
}


LieGroup::~LieGroup() {

}

LieGroup LieGroup::operator*(const LieGroup & rhs) {
	LieGroup temp(*this);
	
	if (temp.V.rows() != rhs.V.rows()) {
		cerr << "operator*: dimensions do not coincide" << endl;
		return temp;
	}
	temp.E *= rhs.E;
	temp.V += rhs.V;	
	return temp;
}

LieGroup LieGroup::operator*(const LieAlgebra & rhs) {
	LieGroup temp(*this);
	if (temp.V.rows() != rhs.V.rows()) {
		cerr << "operator*: dimensions do not coincide" << endl;
		return temp;
	}
	temp.E *= rhs.E;
	temp.V = rhs.V;
	return temp;
}

LieAlgebra LieGroup::log_g() {
	return LieAlgebra(MEFcommon::logSE(E),V,order);
}





