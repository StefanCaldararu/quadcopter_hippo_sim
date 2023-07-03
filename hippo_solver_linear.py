#author: Stefan Caldararu
#a solver to get ESC values for the desired accelerations.

import numpy as np
from scipy.optimize import fsolve
class solver(object):
    def __init__(self):
        self.thrust0 = 0
        self.thurst4 = 0
        self.thrust5 = 0
        self.A = np.array([
            [1,1,1],
            [25, 0, -25],
            [0,25, 0]
        ])

    def compute_thrusts(self, udot, qdot, rdot, nu0, nu4, nu5):
        self.thrust0 = (1.43+1.11)*udot+nu0*5.39
        self.thrust4 = (0.010717+0.0163)*qdot+nu4*(-0.007)
        self.thrust5 = (0.010717+0.0163)*rdot+nu5*(-0.007)
        #print("THRUST BAD: ", self.thrust0, " ",self.thrust4, " ",self.thrust5 )
    def getESC(self):
        print(self.thrust0)
        f2 = 25*self.thrust5
        f1 = (self.thrust0+25*self.thrust4-25*self.thrust5)/2
        f3 = (self.thrust0-25*self.thrust4-25*self.thrust5)/2
        # if(f2<0):
        #     input()
        x = np.array([f1, f2, f3])
        # b = np.array([self.thrust0, self.thrust4, self.thrust5])
        # x = np.linalg.solve(self.A, b)
        #print("X BAD: ", x)
        ESCs = np.zeros(3)
        for i in range(0,3):
            ESCs[i] = self.f2esc(x[i])
        ESCs = np.vstack(np.append(ESCs,0))
        #print("BAD ESC: ", ESCs)
        return ESCs

    def f2esc(self, f):
        neg = 1
        if(f<0):
            f = abs(f)
            neg = -1
        lower_bound = 1500
        upper_bound = 2000
        e = 1e-6
        x = 1600 #init guess
        while(abs(self.function(x)-f)>e):
            x=  x-((self.function(x)-f)/self.deriv(x))
            if(x<lower_bound):
                x = lower_bound
                break
            if(x>upper_bound):
                x = upper_bound
                break
        if(neg == -1):
            return 3000-x
        else:
            return x


    def function(self,x):
        p1 = 1478e-4
        p2 = -1130e-7
        p3 = 2838e-11
        p0 = -6330e-2
        return x*x*x*p3+x*x*p2+x*p1+p0
    def deriv(self,x):
        p1 = 1478e-4
        p2 = -1130e-7
        p3 = 2838e-11
        return 3*x*x*p3+2*x*p2+p1
