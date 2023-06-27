//The input for our solver is our desired thrusts, t_0 (linear x), t_4 (ang. pitch?), and t_5 (ang. yaw)
#include<vector>
#include <Eigen/Dense>
//TODO: move this into the constructor.
int MODE = 0;

class solver{
    private:
        double m = 1.43;
        double p1 = 1478e-4;
        double p2 = -1130e-7;
        double p3 = 2838e-11;
        double p0 = -6330e-2;
        //our thrusts;
        double thrust0 = 0;
        double thrust4 = 0;
        double thrust5 = 0;
        //assumption that one motor has failed! This will be motor number 4, so f4 = 0.
        Eigen::Matrix3d A <<    1, 1, 1,
                                0.04,0,-0.04,
                                0,0.04,0;
        Eigen::Vector4d getESC(){
            Eigen::Vector3d b <<t0, t4, t5;
            Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
            Eigen::Vector4d ESCs << 0,0,0,0;
            for(int i = 0;i<3;i++){
                ESCs[i] = f2esc(x[i]);
            }
            return ESCs;
        }
        
        double f2esc(double f){
            double lower_bound = 1500;
            double upper_bound = 2000;
            double e = 1e-6;
            double x = 1600; //init_guess
            while (std::abs(function(x)-y)>e){
                x = x-((function(x)-y)/deriv(x));
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

        double function(double x){
            return x*x*x*p3+x*x*p2+x*p1+p0;
        }
        double deriv(double x){
            return 3*x*x*p3+2*x*p2+p1;
        }

        Eigen::Vector3d eval_force(double t){
            Eigen::Vector3d ret <<  p3*(t*t*t)+p2*(t*t)+p1*t+p0,
                                    0,
                                    0;
            return ret;
        }
        //TODO: constants in here will likely need to be tuned?
        void computeThrusts(double udot, double qdot, double rdot, double nu0, double nu4, double nu5){
            //From M_T^inv, and D_A
            thrust0 = udot*0.39370079+nu0*5.39;
            thrust4 = qdot*37.01373209+0.007*nu4;
            thrust5 = rdot*37.01373209+0.007*nu5;
        }
}