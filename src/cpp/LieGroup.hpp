#ifndef LIEGROUP
#define LIEGROUP

#include "common.hpp"

struct LieGroup;
#include "LieAlgebra.hpp"

struct LieGroup {
  
public:
  LieGroup();
  LieGroup(const Matrix4d & E);
  LieGroup(const Matrix4d & E, const VectorXd & V, int order);
  ~LieGroup();  
  
  /** This matix returns the log map of the current Lie group element
   * @return: Element on Lie algebra
   */
  LieAlgebra log_g();
  
  
  /** Lie group group as (matrix multiplication)
   * @param rhs: Lie Group element that is multiplied form the right hand side to 
   * current element
   */
  LieGroup operator*(const LieGroup& rhs);
  
  /** This method is a shorthand for the tangent map of the left translation evaluated at identity
   * @param rhs: Lie Algebra at RHS that shoulb transported to T_{lhs}G, where G denote the Lie Group
   * @return: tangent vector on T_{lhs}G
   */
  LieGroup operator*(const LieAlgebra& rhs);
  
  Matrix4d E;
  VectorXd V;
  int order;

};


#endif