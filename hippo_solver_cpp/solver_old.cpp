//The input for our solver is our desired thrusts, t_0 (linear x), t_4 (ang. pitch?), and t_5 (ang. yaw)
#include<vector>
#include <Eigen/Dense>
//TODO: move this into the constructor.
int MODE = 0;

class solver{
    private:
        double m = 1.43;
        Eigen::Vector3d m1p <<  -0.02,
                                0.0,
                                0.04;
        Eigen::Vector3d m2p <<  -0.02,
                                0.04,
                                0.0;
        Eigen::Vector3d m3p <<  -0.02,
                                0.0,
                                -0.04;
        Eigen::Vector3d m4p <<  -0.02,
                                -0.04,
                                0.0;
        double p1 = 1478e-4;
        double p2 = -1130e-7;
        double p3 = 2838e-11;
        double p0 = -6330e-2;
        //our thrusts;
        double thrust0 = 0;
        double thrust4 = 0;
        double thrust5 = 0;
        //this vector is <t1, t2, t3, t4> the thrusts of the four motors.
        //TODO: there should be a reduction here to solving a system of linear equations, because we are only evaluating the cubic function on a subdomain where there is an inverse, and the rest of the function is linear. implement this...
        std::vector<double,4> evaluate_functions(std::vector<double, 4> ESCs){
            Eigen::Vector3d f1 = eval_force(ESCs[0]);
            Eigen::Vector3d f2 = eval_force(ESCs[1]);
            Eigen::Vector3d f3 = eval_force(ESCs[2]);
            Eigen::Vector3d f4 = eval_force(ESCs[3]);

            std::vector<double,3> f1 = eval_force(ESCs[0]);
            std::vector<double,3> f2 = eval_force(ESCs[1]);
            std::vector<double,3> f3 = eval_force(ESCs[2]);
            std::vector<double,3> f4 = eval_force(ESCs[3]);

            //check if one of the motors is dead!
            if(MODE == 1){
                f4.clear();
                for(int i = 0;i<3;i++)
                    f4.push_back(0);
            }
            //TODO: pretty sure these are linear...
            //FIXME: do we need to incorporate the roll? only thing the moments affect are the roll...
            //Evaluate the equations
            double eq1 = f1[0]+f2[0]+f3[0]+f4[0]-thrust0;
            double eq2 = 0;
            double eq3 = f1[0]*m1p[2]+f3[0]*m3p[2]-thrust4;
            double eq4 = f2[0]*m2p[1]+f4[0]*m2p[1]-thrust5;
            std::vector<double, 4> results{eq1, eq2, eq3, eq4};
            return results;
        }

        //calculate inverse of Jacobian...


        Eigen::Vector3d eval_force(double t){
            Eigen::Vector3d ret <<  p3*(t*t*t)+p2*(t*t)+p1*t+p0,
                                    0,
                                    0;
            return ret;
        }

        void set_thrusts(double t0, double t4, double t5){
            thrust0 = t0;
            thrust4 = t4;
            thrust5 = t5;
        }
}