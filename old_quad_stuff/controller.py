#author: Stefan Caldararu
#this solves the system of equations required when we see some motor failure, or just in general to generate control inputs that give us our desired trajectory...
#reference for equations, as well as some example possible solutions: https://ieeexplore.ieee.org/document/6906588

from scipy.optimize import fsolve
import numpy as np

class controller(object):
    def __init__(self):
        self.initial_guess = [0,0.3,0.958,0,0,0, 0.1,437.00751841778504,562.2395063000737,437, 0]
        
    def equations(vars): 
    #constants: 
        k_f = 6.41*0.000001#coefficient for angular velocity of rotors to force
        k_tau = 1.69*0.01
        m = 0.5#KG #mass of vehicle
        g = np.array([
            [0],
            [0],
            [-9.81]
        ])
        gamma = 2.75*0.001#drag coefficient
        l = 0.17#m#distance from center of vehicle to rotors
        I_B_xx = 3.2*0.001#vehicle body inertia... I_B_yy is the same 
        I_B_zz = 1#vehicle body inertia...
        #I_P_xx = 0#vehicle propeller inertia
        I_P_zz = 1.5*0.00001#vehicle propeller inertia
        #total inertias, sum of propeller and body
        I_T_xx = 3.2*0.001#KGm**2
        I_T_zz = 5.5*0.001#KGm**2
        I_B_zz = I_T_zz-4*I_P_zz
        #variables: 
        n_x, n_y, n_z, p, q, r, epsilon, o1, o2, o3, o4  = vars[0], vars[1], vars[2], vars[3], vars[4], vars[5], vars[6], vars[7], vars[8], vars[9], vars[10]

        eq1 = k_f*(o2**2-o4**2)*l-(I_T_zz-I_T_xx)*q*r-I_P_zz*q*(o1+o2+o3+o4) #(12)
        eq2 = k_f*(o3**2-o1**2)*l+(I_T_zz-I_T_xx)*p*r+I_P_zz*p*(o1+o2+o3+o4)#(13)
        eq3 = -gamma*r+k_tau*k_f*(o1**2-o2**2+o3**2-o4**2)#(14)


        eq4 = n_x-epsilon*p#(15)
        eq5 = n_y-epsilon*q#(15)
        eq6 = n_z-epsilon*r#(15)

        eq7 = np.sqrt(n_x**2+n_y**2+n_z**2)-1#(17)
        eq8 = k_f*n_z*(o1**2+o2**2+o3**2+o4**2)-m*np.sqrt(g[0,0]**2+g[1,0]**2+g[2,0]**2)#(18)
        eq9 = o4
        #eq10 = o1-o3
        eq10 = o2
        eq11 = n_x -0.1
        #eq11 = o2/o3
        return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11]


    def solve(self,thetadot, ):
        #when one motor fails, we do have to spin. But it is around the origin, so n = [0,0,1]. now, R = 0, so need to calculate T, what our period is

        self.thetadot = thetadot
        solution = fsolve(self.equations, self.initial_guess)
        self.initial_guess = solution

        # print("Solution:")
        # print("n =", solution[0])
        # print("n =", solution[1])
        # print("n =", solution[2])
        # print("1 =", solution[7])
        # print("2 =", solution[8])
        # print("3 =", solution[9])
        # print("4 =", solution[10])
        return solution#solution[7], solution[8], solution[9], solution[10]
