#author: Stefan Caldararu
#this solves the system of equations required when we see some motor failure, or just in general to generate control inputs that give us our desired trajectory...
#reference for equations, as well as some example possible solutions: https://ieeexplore.ieee.org/document/6906588

from scipy.optimize import fsolve
import numpy as np

class controller(object):
    def __init__(self):
        #constants
        self.k_f = 6.41*0.000001#coefficient for angular velocity of rotors to force
        self.k_tau = 1.69*0.01
        self.m = 0.5#KG #mass of vehicle
        self.g = np.array([
            [0],
            [0],
            [-9.81]
        ])
        self.gamma = 2.75*0.001#drag coefficient
        self.l = 0.17#m#distance from center of vehicle to rotors
        self.I_B_xx = 3.2*0.001#vehicle body inertia... I_B_yy is the same 
        self.I_B_zz = 1#vehicle body inertia...
        #I_P_xx = 0#vehicle propeller inertia
        self.I_P_zz = 1.5*0.00001#vehicle propeller inertia
        #total inertias, sum of propeller and body
        self.I_T_xx = 3.2*0.001#KGm**2
        self.I_T_zz = 5.5*0.001#KGm**2
        self.I_B_zz = self.I_T_zz-4*self.I_P_zz#vehicle body inertia...
        self.initial_guess = [0,0.289,0.958,0,0,0, 0.1,437.00751841778504,562.2395063000737,437, 0]
    def equations(self, vars): 
        n_x, n_y, n_z, pdot, qdot, rdot, epsilon, omega1, omega2, omega3, omega4 = vars[0], vars[1], vars[2], vars[3], vars[4], vars[5], vars[6], vars[7], vars[8], vars[9], vars[10]
        p = self.thetadot[0,0]
        q = self.thetadot[1,0]
        r = self.thetadot[2,0]
        #equations
        eq1 = self.k_f*(omega2**2-omega4**2)*self.l - (self.I_T_zz-self.I_T_xx)*q*r-self.I_P_zz*q*(omega1+omega2+omega3+omega4)-self.I_B_xx*pdot
        eq2 = self.k_f*(omega3**2-omega2**2)*self.l+(self.I_T_zz-self.I_T_xx)*p*r+self.I_P_zz*p*(omega1+omega2+omega3+omega4)-self.I_B_xx*qdot
        eq3 = -self.gamma*r+self.k_tau*self.k_f*(omega1**2-omega2**2+omega3**2-omega4**2)-self.I_B_zz*rdot
        eq4 = epsilon*p-n_x
        eq5 = epsilon*q-n_x
        eq6 = epsilon*r-n_x
        eq7 = np.sqrt(n_x**2+n_y**2+n_z**2)-1
        eq8 = self.k_f*(omega1**2+omega2**2+omega3**2+omega4**2)-self.m*np.sqrt(self.g[0,0]**2+self.g[1,0]**2+self.g[2,0]**2)
        eq9 = omega4
        eq10 = omega1-omega3
        eq11 = omega2/omega3-0.5
        #for now, this is assuming we want linear velocity to be 0? YES, because otherwise eq1-eq3 would have a pdot value!
        return[eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11]


    def solve(self,thetadot):
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
        return solution[7], solution[8], solution[9], solution[10]
