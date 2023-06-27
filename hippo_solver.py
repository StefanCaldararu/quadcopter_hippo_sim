#author: Stefan Caldararu
#a solver to get ESC values for the desired accelerations.

import numpy as np
from scipy.optimize import fsolve
class solver(object):
    def __init__(self):
        self.thrust = np.zeros((6,1), dtype = float)
        self.MODE = 1 #1 if 1 engine dies

    def equations(self, esc):
        e1, e2, e3, e4 = esc[0], esc[1], esc[2], esc[3]
        m1p = np.array([-0.02,0,0.04])
        m3p = np.array([-0.02, 0, -0.04])
        m2p = np.array([-0.02, 0.04, 0])
        m4p = np.array([-0.02, -0.04, 0])
        p1 = 1478e-4
        p2 = -1130e-7
        p3 = 2838e-11
        p0 = -6330e-2
        f1 = np.array([p3*(e1**3)+p2*(e1**2)+p1*(e1)+p0,0,0])
        f2 = np.array([p3*(e2**3)+p2*(e2**2)+p1*(e2)+p0,0,0])
        f3 = np.array([p3*(e3**3)+p2*(e3**2)+p1*(e3)+p0,0,0])
        f4 = np.array([p3*(e4**3)+p2*(e4**2)+p1*(e4)+p0,0,0])
        if(self.MODE == 1):
            f4 = np.array([0,0,0])
        if(self.MODE == 2):
            f4 = np.array([0,0,0])
            f2 = np.array([0,0,0])
        eq1 = f1[0]+f2[0]+f3[0]+f4[0]-self.thrust[0,0]
        #eq1 = 0
        #Dont care about roll... so set this equation equal to 0. moments don't matter.
        eq2 = 0#np.cross(m1p,f1)[0]+np.cross(m2p,f2)[0]+np.cross(m3p,f3)[0]+np.cross(m4p,f4)[0]
        eq3 = np.cross(m1p,f1)[1]+np.cross(m2p,f2)[1]+np.cross(m3p,f3)[1]+np.cross(m4p,f4)[1]-self.thrust[4,0]
        eq4 = np.cross(m1p,f1)[2]+np.cross(m2p,f2)[2]+np.cross(m3p,f3)[2]+np.cross(m4p,f4)[2]-self.thrust[5,0]
        # print(eq4)
        # print(eq3)

        return [eq1, eq2, eq3, eq4]

    def eq2(self, thrust):
        t0,t4,t5 = thrust[0],thrust[1],thrust[2]
        m = 1.43
        M_A = np.diag([-1.11, -2.80, -2.80, -0.00451, -0.0163, -0.0163])
        I_g = np.diag([0.002408, 0.010717, 0.010717]) 
        D_A = np.diag([-5.39, -17.36, -17.36, -0.00114, -0.007, -0.007])

        M_T = np.array([
            [m-M_A[0,0], 0, 0],
            [0, m-M_A[2,2], 0],
            [0, 0, I_g[1,1]-M_A[4,4]]#TODO: make sure I_g correct here...
        ])
        D = np.diag([-D_A[0,0], -D_A[2,2], -D_A[4,4]])

        C = np.array([
            [0,0,0],
            [0,0,-(m-M_A[0,0])*self.nu[0,0]],
            [0, (M_A[2,2]-M_A[0,0])*self.nu[0,0], 0] #beacuse x_g = 0, 3,3 = 0
        ])

        thrust1 = np.array([[t0], [0], [t4]])
        s1 = np.array([[self.nu[0,0]], [self.nu[2,0]], [self.nu[4,0]]])
        RHS1 = thrust1-D@s1-C@s1
        accel1 = np.linalg.inv(M_T)@RHS1
        #want to use accel1[0,0], accel1[2,0]
        M_T = np.array([
            [m-M_A[1,1], 0, 0],
            [0, I_g[0,0]-M_A[3,3], 0],
            [0, 0, I_g[2,2]-M_A[5,5]]#TODO: make sure I_g correct here...
        ])
        D = np.diag([-D_A[1,1], -D_A[3,3], -D_A[5,5]])
        C = np.array([
            [0,0,(m-M_A[0,0])*self.nu[0,0]],
            [0,0,0],
            [(M_A[0,0]-M_A[1,1])*self.nu[0,0], 0, 0]#because x_g = 0, 3,3 = 0
        ])
        thrust2 = np.array([[0], [0], [t5]])
        #state for this system: 
        s2 = np.array([[self.nu[1,0]], [self.nu[3,0]], [self.nu[5,0]]])
        RHS2 = thrust2-D@s2-C@s2
        accel2 = np.linalg.inv(M_T)@RHS2
        #want accel2[2,0]

        return [accel1[0,0]-self.udot, accel1[2,0]-self.qdot, accel2[2,0]-self.rdot]
        
    def solve(self, udot, qdot, rdot, nu):
        self.nu = nu
        m = 1.43
        M_A = np.diag([-1.11, -2.80, -2.80, -0.00451, -0.0163, -0.0163])
        D_A = np.diag([-5.39, -17.36, -17.36, -0.00114, -0.007, -0.007])
        I_g = np.diag([0.002408, 0.010717, 0.010717])
        self.qdot = qdot
        self.rdot = rdot
        self.udot = udot

        thrust_init_guess = np.array([1, 0.1, 0.1])
        thrust = fsolve(self.eq2, thrust_init_guess)
        self.thrust[0,0] = thrust[0]
        self.thrust[4,0] = thrust[1]
        self.thrust[5,0] = thrust[2]



        # self.thrust[0,0] = (m-M_A[0,0])*udot - D_A[0,0]*nu[0,0]
        # self.thrust[4,0] = (I_g[1,0]-M_A[4,4])*qdot-D_A[4,4]*nu[4,0]-(m-M_A[0,0])*nu[0,0]*nu[4,0]
        # self.thrust[5,0] = (I_g[2,0]-M_A[5,5])*rdot - M_A[5,5]*nu[5,0]-(m-M_A[0,0])*nu[0,0]*nu[5,0]
        # print("thrust 0: ",self.thrust[0,0])
        # print("thrust 4: ",self.thrust[4,0])
        # print("thrust 5: ",self.thrust[5,0])

        initial_guess = [1900, 1900, 1650, 1900]
        solution = fsolve(self.equations, initial_guess)
        # print("SOLUTION: ", solution)


        
        
        escVals = np.array([
            [solution[0]],
            [solution[1]],
            [solution[2]],
            [solution[3]]
        ])

        return escVals
