#include <iostream>
#include "common.hpp"
#include "LieAlgebra.hpp"
#include "LieGroup.hpp"
#include "MinEnFilter.hpp"
// #include "LieAlgebraSE.hpp"
using namespace MEFcommon;

using std::cout;
using std::endl;

int main(int argc, char **argv) {
	
	// Variables
	Matrix4d m_test, m_test1, m_test2;
	Vector6d v_test, v_test1;
	Matrix3d n_test, n_test1, n_test2;
	Vector4d w_test, w_test1;
	
    std::cout << "UNIT TESTS \n ================================================================================" << std::endl;

    
	std::cout << "UNIT TEST OF MATSE / VECSE OPERATRIONS" << std::endl;
	// Create random matrix
	m_test = Matrix4d::Random();   
    
	// Apply vecse operator
	v_test = vecSE(m_test);
	// Apply inverse operatrion
	m_test1 = matSE(v_test);
	// Apply again inverse operation
	v_test1 = vecSE(m_test1);	
	if (v_test == v_test1) {
		std::cout << "Test of vecse / matse successful. \n " << std::endl;
	}
	
	std::cout << "UNIT TEST OF EXPONENTIAL/ LOGARITHMIC MAP ON S0(3)" << std::endl;
	n_test << 0, -3, 2,  3, 0, -1, -2, 1, 0;
	n_test *= 0.1;
    n_test1 = expSO(n_test);
	n_test2 = logSO(n_test1);
	cout << "Exponential of \n" << n_test << "\n is \n" << n_test1 << endl;
	cout << "Logarithm of \n" << n_test1 << "\n is \n" << n_test2 << endl;
	
    
  	std::cout << "UNIT TEST OF EXPONENTIAL/ LOGARITHMIC MAP ON SE(3)" << std::endl;
	m_test << 0, -3, 2, 4, 3, 0, -1, 5, -2, 1, 0, 6, 0, 0, 0, 0;
	m_test *= 0.1;
    m_test1 = expSE(m_test);
	m_test2 = logSE(m_test1);
	cout << "Exponential of \n" << m_test << "\n is \n" << m_test1 << endl;
	cout << "Logarithm of \n" << m_test1 << "\n is \n" << m_test2 << endl;  
	
	std::cout << "UNIT TEST OF PROJECTION ONTO SE(3)" << std::endl;
	m_test = Matrix4d::Random();
	cout << "Projection of m_test defined by \n" << m_test << "\n onto Lie algebra SE(3) is \n" << prSE(m_test) << endl;
	
	cout << "UNIT TEST OF ADJOINT ON SE(3)" << endl;
	
	VectorXd test(6);
	test <<  23,-0.5,15,0.54,34.253,-1;
	cout << "\n\n Test vector: \n\n" << test << endl;
	cout << "\n\n Adjoint vector: \n\n" << adjoint(matSE(test)) << endl;
	
	
	
	
   
	cout << "UNIT TEST OF LIE ALGEBRA / LIE GROUP" << endl;
	m_test << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16;
	v_test << 1,3,4,5,6,2;
	m_test1 << 1,3,2,4,7,5,9,3,7,4,2,5,7,8,9,4;
	v_test1 << 0, 3, 7, 1, 6, 3;
	
	
	LieGroup L1(m_test,v_test,2);
	LieGroup L2(m_test1, v_test1,2);
	
	
	cout << "\n\n L1: \n" << L1.E<< "\n\n" << L1.V << endl;
	cout << "\n\n L2: \n" << L2.E<< "\n\n" << L2.V << endl;
	

	
// 	F = MEF(n_test, 1,2,3,4,5);
	
	
	
// 	LieGroup L = L1*L2;
// 	cout << "Product from L1 and L2:  \n\n" << L.getE() << "\n\n" << L.getV() << endl;

	
	return 0;
}
