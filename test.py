from scipy.optimize import fsolve
from scipy.spatial.transform import Rotation
import numpy as np
from hippo_solver import solver
from PD_hippo import PD
import hippo_basic_sim
import matplotlib.pyplot as plt

ang = np.array([3.1415,0,0])
R = Rotation.from_rotvec(ang)
print(R.as_matrix())