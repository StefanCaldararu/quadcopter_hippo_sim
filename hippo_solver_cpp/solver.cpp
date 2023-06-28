//The input for our solver is our desired thrusts, t_0 (linear x), t_4 (ang. pitch?), and t_5 (ang. yaw)
#include "solver.h"
#include <Eigen/Dense>
//TODO: move this into the constructor.

double solver::f2esc(double f){
    double lower_bound = 1500;
    double upper_bound = 2000;
    double e = 1e-6;
    double x = 1600; //init_guess
    while (std::abs(function(x)-f)>e){
        x = x-((function(x)-f)/deriv(x));
        if(x<lower_bound){
            x = lower_bound;
            break;
        }
        if(x>upper_bound){
            x = upper_bound;
            break;
        }
    }
    return x;
}
double solver::function(double x){
    return x*x*x*p3+x*x*p2+x*p1+p0;
}
double solver::deriv(double x){
    return 3*x*x*p3+2*x*p2+p1;
}
solver::solver(){
    A <<1, 1, 1,
        0.04,0,-0.04,
        0,0.04,0;

    m = 1.43;
    p1 = 1478e-4;
    p2 = -1130e-7;
    p3 = 2838e-11;
    p0 = -6330e-2;
}

Eigen::Vector4d solver::getESC(){
    Eigen::Vector3d b; 
    b << thrust0, thrust4, thrust5;
    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
    Eigen::Vector4d ESCs;
    ESCs << 0,0,0,0;
    for(int i = 0;i<3;i++){
        ESCs[i] = f2esc(x[i]);
    }
    return ESCs;
}
void solver::computeThrusts(double udot, double qdot, double rdot, double nu0, double nu4, double nu5){
    //From M_T^inv, and D_A
    thrust0 = udot*0.39370079+nu0*5.39;
    thrust4 = qdot*37.01373209+0.007*nu4;
    thrust5 = rdot*37.01373209+0.007*nu5;
}
