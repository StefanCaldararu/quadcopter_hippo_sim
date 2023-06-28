//The input for our solver is our desired thrusts, t_0 (linear x), t_4 (ang. pitch?), and t_5 (ang. yaw)
#ifndef solver
#define solver

#include <Eigen/Dense>
//TODO: move this into the constructor.

class solver{
    public:
        solver();
        Eigen::Vector4d getESC();
        void computeThrusts(double udot, double qdot, double rdot, double nu0, double nu4, double nu5);
        
    private:
        double m;
        double p1;
        double p2;
        double p3;
        double p0;
        //our thrusts;
        double thrust0;
        double thrust4;
        double thrust5;
        //assumption that one motor has failed! This will be motor number 4, so f4 = 0.
        Eigen::Matrix3d A;
        
        double f2esc(double f);
        double function(double x);
        double deriv(double x);
};
#endif  