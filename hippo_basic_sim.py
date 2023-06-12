#author: Stefan Caldararu
#dynamics from "Marine Craft Hydrodynamics and Motion Control" by Thor I. Fossen
#this uses the general model, and doesn't work well... problems with the C matrices
#sim the dynamics of the hippocampus, given the old BSC066 parameters.
#Ignoring Coriolis, and assuming values from
import numpy as np 
from quad_vis import vis
from scipy.spatial.transform import Rotation
from PD_hippo import PD

M_A = np.diag([-1.11, -2.80, -2.80, -0.00451, -0.0163, -0.0163])
I_g = np.diag([0.002408, 0.010717, 0.010717])
m = 1.43
#FIXME: make sure we need the lower half of this matrix
M_RB = np.diag([m, m, m, 0.002408, 0.010717, 0.010717])

#for now, assume Damping is diagonal constant...
D_A = np.diag([-5.39, -17.36, -17.36, -0.00114, -0.007, -0.007])
# def calcD_A(nu):
#     D_A = np.zeros((6,6), dtype = float)
#     return D_A




def main():
    dt = 0.05
    eta = np.array([[0],[0],[3],[0.0],[0],[0]])
    nu = np.zeros((6,1), dtype = float)
    nu[4,0] = 0.1

    esc = 0
    ESC = np.array([[esc],[esc],[esc],[esc]])
    start = 0
    end = 20
    time = np.arange(start, end, dt)
    plt = vis(0.5)#FIXME: make this so that it can vis horizontal
    
    for t in time:
        # if(t>1):
        #     esc = 0
        #     ESC = np.array([[esc],[esc],[esc],[esc]])
        dest = np.array([[3],[1],[4]])
        #get errors 
        e_d, e_theta, e_psi = PD.getErrors(eta, dest)
        #print(e_theta)
        udot, qdot, rdot = PD.getAccels(e_d, e_theta, e_psi, nu)
        # print(e_theta)
        # print(qdot)
        thrust = computeThrust(ESC)
        #print(thrust)
        eta, nu = motion_model(eta, nu, thrust, dt)
        plt.plot(eta[:3], eta[-3:])
        print(nu)



#from the ESC values of the 4 motors, return the thrust.
def computeThrust(ESC):
    #the force for each motor...
    f = np.array([np.zeros((3,1), dtype = float), np.zeros((3,1), dtype = float), np.zeros((3,1), dtype = float), np.zeros((3,1), dtype = float)])
    #calibrated based off of Timo paper, 8v eq 5.1
    p1 = 1478e-4
    p2 = -1130e-7
    p3 = 2838e-11
    p0 = -6330e-2
    #for now, assuming thrusters can only go forward. TODO: fix this...
    for i in range(0,4):
        escval = ESC[i,0]
        f[i,0,0] = p3*escval**3+p2*escval**2+p1*escval+p0
        if(f[i,0,0]<0):
            f[i,0,0] = 0
    #motor positions, around the body reference frame.
    #assuming each motor is positioned 2cm behind com, and 4 cm outside the body.
    m1p = np.array([
        [-0.02],
        [0],
        [0.04]
    ])
    m2p = np.array([
        [-0.02],
        [0],
        [-0.04]
    ])
    m3p = np.array([
        [-0.02],
        [0.04],
        [0]
    ])
    m4p = np.array([
        [-0.02],
        [-0.04],
        [0]
    ])
    locs = np.array([m1p, m2p, m3p, m4p])
    #now using lennarts model, assuming 0 moment produced by each motor.
    thrust = np.zeros((6,1), dtype = float)
    for i in range(0,4):
        fi = f[i]
        thrust[0,0] = thrust[0,0]+fi[0,0]
        moment = np.cross(locs[i].flatten(), fi.flatten())
        moment = np.reshape(moment,(3,1))
        thrust[3,0] = thrust[3,0]+moment[0,0]
        thrust[4,0] = thrust[4,0]+moment[1,0]
        thrust[5,0] = thrust[5,0]+moment[2,0]
        #TODO: here assuming moment about motors is 0...
        
    
    return thrust
def motion_model(eta, nu, thrust, dt):


    #Split int two secitons, longitudinal dynamics (using u,w,q)

    M_T = np.array([
        [m-M_A[0,0], 0, 0],
        [0, m-M_A[2,2], 0],
        [0, 0, I_g[1,1]-M_A[4,4]]#TODO: make sure I_g correct here...
    ])
    D = np.diag([-D_A[0,0], -D_A[2,2], -D_A[4,4]])

    C = np.array([
        [0,0,0],
        [0,0,-(m-M_A[0,0])*nu[0,0]],
        [0, (M_A[2,2]-M_A[0,0])*nu[0,0], 0] #beacuse x_g = 0, 3,3 = 0
    ])

    #assuming hydrostatics is 0.
    #thrust for this system: 
    t1 = np.array([[thrust[0,0]], [thrust[2,0]], [thrust[4,0]]])
    #state for this system: 
    s1 = np.array([[nu[0,0]], [nu[2,0]], [nu[4,0]]])

    RHS1 = t1-D@s1-C@s1
    accel1 = np.linalg.inv(M_T)@RHS1
    nu[0,0] = nu[0,0]+accel1[0,0]*dt
    nu[2,0] = nu[2,0]+accel1[1,0]*dt
    nu[4,0] = nu[4,0]+accel1[2,0]*dt
    
    #Now we are doing lateral dynamics (v,p,r). For us, the linear dynamics of  v should be 0, but this is because the thrust in this dimension will be 0.

    M_T = np.array([
        [m-M_A[1,1], 0, 0],
        [0, I_g[0,0]-M_A[3,3], 0],
        [0, 0, I_g[2,2]-M_A[5,5]]#TODO: make sure I_g correct here...
    ])
    D = np.diag([-D_A[1,1], -D_A[3,3], -D_A[5,5]])
    C = np.array([
        [0,0,(m*M_A[0,0])*nu[0,0]],
        [0,0,0],
        [(M_A[0,0]-M_A[1,1])*nu[0,0], 0, 0]#because x_g = 0, 3,3 = 0
    ])
    t2 = np.array([[thrust[1,0]], [thrust[3,0]], [thrust[5,0]]])
    #state for this system: 
    s2 = np.array([[nu[1,0]], [nu[3,0]], [nu[5,0]]])
    RHS2 = t2-D@s2-C@s2
    accel2 = np.linalg.inv(M_T)@RHS2
    nu[1,0] = nu[1,0]+accel2[0,0]*dt
    nu[3,0] = nu[3,0]+accel2[1,0]*dt
    nu[5,0] = nu[5,0]+accel2[2,0]*dt

    #nu is now propogated. 
    body2world = np.zeros((6,6), dtype=float)
    theta = np.array([[eta[3,0]],[eta[4,0]],[eta[5,0]]])

    body2world[:3, :3] = body2worldRot(theta)
    body2world[3:,3:] = body2worldTrans(theta)
    globalChange = body2world@nu
    eta = eta+dt*globalChange#FIXME: do we need to include nudot here too?
    return eta, nu
#Rotation and transformation matrices...
def body2worldTrans(Theta):
    phi = Theta[0,0]
    theta = Theta[1,0]
    T = np.array([
        [1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
        [0, np.cos(phi), -np.sin(phi)],
        [0, 0,1]#FIXME: for some reason this doesn't work, div by 0... np.sin(phi)/np.tan(theta), np.cos(phi)/np.cos(theta)]
    ])
    return T
def world2bodyTrans(theta):
    T = body2worldTrans(theta)
    return np.linalg.inv(T)
def body2worldRot(Theta):
    phi = Theta[0,0]
    theta = Theta[1,0]
    psi = Theta[2,0]
    R = np.array([  
        [np.cos(psi)*np.cos(theta), np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi), np.sin(psi)*np.sin(phi)+ np.cos(psi)*np.cos(phi)*np.sin(theta)],

        [np.sin(psi)*np.cos(theta), np.cos(psi)*np.cos(phi)+np.sin(phi)*np.sin(theta)+np.sin(psi), np.sin(theta)*np.sin(psi)*np.cos(phi)-np.cos(psi)*np.sin(phi)],

        [-np.sin(theta), np.cos(theta)*np.sin(phi), np.cos(theta)*np.cos(phi)]
    ])
    return R
def world2bodyRot(theta):
    R = body2worldRot(theta)
    return np.linalg.inv(R)
if __name__ == '__main__':
    main()