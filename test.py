from scipy.optimize import fsolve
import numpy as np
from hippo_solver import solver
from PD_hippo import PD
import hippo_basic_sim
import matplotlib.pyplot as plt



# M_A = np.diag([-1.11, -2.80, -2.80, -0.00451, -0.0163, -0.0163])
# m = 1.43
# I_g = np.diag([0.002408, 0.010717, 0.010717]) 
# M_T = np.array([
#             [m-M_A[0,0], 0, 0],
#             [0, m-M_A[2,2], 0],
#             [0, 0, I_g[1,1]-M_A[4,4]]#TODO: make sure I_g correct here...
#         ])
# print(np.linalg.inv(M_T))
