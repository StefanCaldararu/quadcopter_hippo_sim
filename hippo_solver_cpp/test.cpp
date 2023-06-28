#include "solver.h"
#include <iostream>
#include <Eigen/Dense>

int main(int argc, char** argv){
    solver s;

    /**
    NU:  [[ 5.25789609e-01]
[ 0.00000000e+00]
 [ 0.00000000e+00]
 [ 6.07577885e+00]
 [-2.92309708e-03]
 [-6.11626196e-06]]
u-rdot:  -0.025789609296966276   -0.00565409718462398   0.119875532846412
esc vals  [[1752.31740169]
 [1488.24588776]
 [1752.64641186]
 [1743.58768337]]

    **/
    s.computeThrusts(-0.025789609296966276, -0.00565409718462398,0.119875532846412, 0.525789609, -0.00292309708, -0.00000611626196 );
    Eigen::Vector4d result = s.getESC();
    std::cout << result[0] << std::endl;
    std::cout << result[1] << std::endl;
    std::cout << result[2] << std::endl;
    std::cout << result[3] << std::endl;
    return 0;
}