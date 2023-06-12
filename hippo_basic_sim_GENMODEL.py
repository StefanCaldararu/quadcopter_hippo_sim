#author: Stefan Caldararu
#dynamics from "Marine Craft Hydrodynamics and Motion Control" by Thor I. Fossen
#this uses the general model, and doesn't work well... problems with the C matrices
#sim the dynamics of the hippocampus, given the old BSC066 parameters.
#Ignoring Coriolis, and assuming values from
import numpy as np 
from quad_vis import vis

M_A = np.diag([-1.11, -2.80, -2.80, -0.00451, -0.0163, -0.0163])
I_g = np.diag([0.002408, 0.010717, 0.010717])
m = 1.43
#FIXME: make sure we need the lower half of this matrix
M_RB = np.diag([m, m, m, 0.002408, 0.010717, 0.010717])

#total mass load
M_T = M_RB+M_A
M_Tinv = np.linalg.inv(M_T)

def calcS(lamb):
    l1 = lamb[0,0]
    l2 = lamb[1,0]
    l3 = lamb[2,0]
    S = np.array([
        [0,-l3, l2],
        [l3, 0, -l1],
        [-l2, l1, 0]
    ])
    return S

def calcC_RB(nu):
    o_b = np.array([[nu[3,0]],[nu[4,0]],[nu[5,0]]])
    #S(o_b)
    S1 = calcS(o_b)
    #S(I_g@o_b)
    S2 = calcS(I_g@o_b)
    C_RB = np.zeros((6, 6), dtype=float)
    C_RB[:3, :3] = m*S1
    C_RB[3:, 3:] = S2
    # C_RB = np.array([
    #     [m*S1[0,0], m*S1[0,1], m*S1[0,2], 0,0,0],
    #     [m*S1[1,0], m*S1[1,1], m*S1[1,2], 0,0,0],
    #     [m*S1[2,0], m*S1[2,1], m*S1[2,2], 0,0,0],
    #     [0,0,0,S2[2,0], S2[2,1], S2[2,2]],
    #     [0,0,0,S2[2,0], S2[2,1], S2[2,2]],
    #     [0,0,0,S2[2,0], S2[2,1],S2[2,2]]
    # ])
    return C_RB

def calcC_A(nu):
    v1 = np.array([[nu[0,0]], [nu[1,0]], [nu[2,0]]])
    v2 = np.array([[nu[3,0]], [nu[4,0]], [nu[5,0]]])
    A11 = M_A[:3, :3]
    A12 = M_A[:3, -3:]
    A21 = M_A[-3:, :3]
    A22 = M_A[-3:, -3:]
    #S(A11@v1+A12@v2)
    S1 = calcS(A11@v1+A12@v2)
    #S(A21@v1+A22@v2)
    S2 = calcS(A21@v1+A22@v2)
    C_A = np.zeros((6,6), dtype=float)
    C_A[:3,3:] = -1*S1
    C_A[3:,:3] = -1*S1
    C_A[3:,3:] = -1*S2
    return C_A

#for now, assume Damping is diagonal constant...
D_A = np.diag([5.39, 17.36, 17.36, -0.00114, -0.007, -0.007])
# def calcD_A(nu):
#     D_A = np.zeros((6,6), dtype = float)
#     return D_A




def main():
    dt = 0.05
    eta = np.array([[0],[0],[3],[0.0],[0],[0]])
    nu = np.zeros((6,1), dtype = float)
    #nu[3,0] = 0.1

    esc = 1650
    ESC = np.array([[esc],[esc],[esc],[0]])
    start = 0
    end = 20
    time = np.arange(start, end, dt)
    plt = vis(0.5)#FIXME: make this so that it can vis horizontal
    
    for t in time:
        # if(t>0.1):
        #     esc = 1650
        #     ESC = np.array([[esc],[esc],[esc],[esc]])

        thrust = computeThrust(ESC)
        #print(thrust)
        eta, nu = motion_model(eta, nu, thrust, dt)
        plt.plot(eta[:3], eta[-3:])
        # print(eta)



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
    
    C_RB = calcC_RB(nu)
    C_A = calcC_A(nu)
    
    #Right hand side of our equation...
    #TODO: why do we need to add damping? this is not the way it should work...
    RHS = thrust-D_A@nu-C_RB@nu-C_A@nu
    nudot = M_Tinv@RHS
    print(nudot)
    nu = nu+dt*nudot
    #nu was in the body reference frame. eta is with respect to the global reference frame. So, need to transfer both...
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