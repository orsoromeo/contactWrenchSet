//   g++ -I /usr/local/include/eigen3/ test_glpk_wrapper.cpp -L /usr/local/lib -l glpk -o test_glpk_wrapper
//  ./test_glpk_wrapper
//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <Eigen/Dense>
//#include <glpk.h>

#include <contact_wrench_set/LP/glpk_wrapper.h>

int main(void)
{

GLPK_wrapper glpk_wrapper;
//glpk_wrapper.call_glpk();

Eigen::MatrixXd A;
Eigen::VectorXd b(2);
Eigen::VectorXd obj(2);
A.resize(2,2);
A(0,0) = -1.0;
A(0,1) = 1.0;
A(1,0) = 1.0;
A(1,1) = 1.0;
b(0) = 1.0;
b(1) = 0.0;
obj.setZero();
Eigen::VectorXd solution;
std::cout<<A<<std::endl;
glpk_wrapper.solve_lp(obj, b, A, solution);
std::cout<<"The solution is: "<<solution.transpose()<<std::endl;

}
