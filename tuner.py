#author: Stefan Caldararu
#same as hippo_basic_sim, just set up to tune the high level Stanley controller
import numpy as np 

from scipy.spatial.transform import Rotation
from PD_hippo import PD
from hippo_solver import solver
from skopt import gp_minimize

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




def main(inp):
    a1 = inp[0]
    b1 = inp[1]
    c1 = inp[2]
    d1 = inp[3]
    e1 = inp[4]
    a2 = inp[5]
    b2 = inp[6]
    c2 = inp[7]
    d2 = inp[8]
    e2 = inp[9]

    dt = 0.01
    eta = np.array([[0],[0],[0],[1],[0.0],[0.0],[0.0]])
    nu = np.zeros((6,1), dtype = float)
    #nu[0,0] = 1.0
    #nu[4,0] = 1
    #nu[5,0] = 1

    esc = 1600
    ESC = np.array([[esc],[esc],[esc],[esc]])
    start = 0
    end = 5
    time = np.arange(start, end, dt)
    #plot = vis(0.5)#FIXME: make this so that it can vis horizontal
    sol = solver()
    count = 0
    udot, qdot, rdot = 0,0,0
    dhist = []

    for t in time:
        # if(t<1):
        #     udot, qdot, rdot = 0.3, 0, 0.1
        if(count%10 == 0):
            udot = -(nu[0,0]-0.5)

            orientation = np.array([eta[3,0], eta[4,0], eta[5,0], eta[6,0]])
            #orientation = q2e(orientation)
            qdot = a1*eta[2,0]+b1*orientation[1]+c1*orientation[2]+d1*orientation[3]+e1*nu[4,0]#1000*compAngles(0,eta[4,0])-1*eta[3,0]
            rdot = a2*eta[1,0]+b2*orientation[1]+c2*orientation[2]+d2*orientation[3]+e2*nu[5,0]#1000*compAngles(0,eta[5,0])-1*eta[2,0]
            ESC = np.clip(sol.solve(udot, qdot, rdot, nu), 1500, 2000)

        thrust = computeThrust(ESC)
        # print(ESC)
        eta, nu = motion_model(eta, nu, thrust, dt)

        eta[3,0] = fix_angle(eta[3,0])
        eta[4,0] = fix_angle(eta[4,0])
        eta[5,0] = fix_angle(eta[5,0])
        count = count+1
        dhist.append(np.sqrt(eta[1,0]**2+eta[2,0]**2))
    return np.max(np.abs(dhist))

def q2e(quaternion):
    if(quaternion[0] == 0):
        return [0,0,0]
    # Normalize the quaternion to ensure unit magnitude
    quaternion /= np.linalg.norm(quaternion)

    # Set the roll angle to 0
    euler_angles = [0, 0, 0]

    # Compute the yaw and pitch angles
    quaternion_rotation = Rotation.from_quat(quaternion)
    yaw_pitch_rotation = quaternion_rotation.inv() * Rotation.from_euler('xyz', euler_angles, degrees=True)
    yaw_pitch_angles = yaw_pitch_rotation.as_euler('xyz', degrees=True)

    # Update the yaw and pitch angles in the Euler angles
    euler_angles[1] = yaw_pitch_angles[1]  # pitch
    euler_angles[2] = yaw_pitch_angles[0]  # yaw

    return euler_angles

def fix_angle(theta):
    if(theta>np.pi):
        theta = theta-2*np.pi
    if(theta<-np.pi):
        theta = theta+2*np.pi
    return theta

def compAngles(ref, act):
    if( (ref>0 and act>0) or (ref<=0 and act <=0)):
        err_theta = ref-act
    elif( ref<=0 and act > 0):
        if(abs(ref-act)<abs(2*np.pi+ref-act)):
            err_theta = -abs(act-ref)
        else:
            err_theta = abs(2*np.pi + ref- act)
    else:
        if(abs(ref-act)<abs(2*np.pi-ref+act)):
            err_theta = abs(act-ref)
        else: 
            err_theta = -abs(2*np.pi-ref+act)
    return err_theta


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
    m3p = np.array([
        [-0.02],
        [0],
        [-0.04]
    ])
    m2p = np.array([
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
    #nu[2,0] = nu[2,0]+accel1[1,0]*dt #no lateral accel...
    nu[4,0] = nu[4,0]+accel1[2,0]*dt
    # print("actual qdot is: ", accel1[2,0])
    # print("actual udot is: ", accel1[0,0])
    
    #Now we are doing lateral dynamics (v,p,r). For us, the linear dynamics of  v should be 0, but this is because the thrust in this dimension will be 0.

    M_T = np.array([
        [m-M_A[1,1], 0, 0],
        [0, I_g[0,0]-M_A[3,3], 0],
        [0, 0, I_g[2,2]-M_A[5,5]]#TODO: make sure I_g correct here...
    ])
    D = np.diag([-D_A[1,1], -D_A[3,3], -D_A[5,5]])
    C = np.array([
        [0,0,(m-M_A[0,0])*nu[0,0]],
        [0,0,0],
        [(M_A[0,0]-M_A[1,1])*nu[0,0], 0, 0]#because x_g = 0, 3,3 = 0
    ])
    t2 = np.array([[thrust[1,0]], [thrust[3,0]], [thrust[5,0]]])
    #state for this system: 
    s2 = np.array([[nu[1,0]], [nu[3,0]], [nu[5,0]]])
    RHS2 = t2-D@s2-C@s2
    accel2 = np.linalg.inv(M_T)@RHS2
    #nu[1,0] = nu[1,0]+accel2[0,0]*dt #no lateral accel...
    nu[3,0] = nu[3,0]+accel2[1,0]*dt
    nu[5,0] = nu[5,0]+accel2[2,0]*dt

    #nu is now propogated. 
    #want to get linear velocities in world frame (using rotaiton), and apply rotation to our quaternions...
    orientation = np.array([eta[3,0],eta[4,0], eta[5,0], eta[6,0]])
    if(orientation[0] == 0):
        orientation = np.array([1,0,0,0])
    R = Rotation.from_quat(orientation)
    R = R.inv()
    linear_vel = np.array([nu[0,0], nu[1,0], nu[2,0]])
    world_linear_vel = R.apply(linear_vel)
    ang_vel = np.array([nu[3,0], nu[4,0], nu[5,0]])
    mag_ang_vel = np.sqrt(ang_vel[0]**2+ang_vel[1]**2+ang_vel[2]**2)
    new_orientation = np.zeros(4)

    if(mag_ang_vel!=0):
        #angular_velocity_global = np.dot(R, ang_vel)
        qdelta = Rotation.from_rotvec(ang_vel * dt)#.as_matrix()
        #qdelta = np.array([dt*mag_ang_vel, ang_vel[0]/mag_ang_vel,ang_vel[1]/mag_ang_vel,ang_vel[2]/mag_ang_vel])
        #rot = Rotation.from_quat(qdelta)
        #print(rot.as_quat())
        new_R = qdelta*R


        new_orientation = new_R.as_quat()
        magnitude = np.linalg.norm(new_orientation)

        # Renormalize the quaternion
        #print(new_orientation)
        new_orientation = new_orientation / magnitude

    eta[0,0] = eta[0,0]+dt*world_linear_vel[0]
    eta[1,0] = eta[1,0]+dt*world_linear_vel[1]
    eta[2,0] = eta[2,0]+dt*world_linear_vel[2]
    eta[3,0] = new_orientation[0]
    eta[4,0] = new_orientation[1]
    eta[5,0] = new_orientation[2]
    eta[6,0] = new_orientation[3]
    # body2world = np.zeros((6,6), dtype=float)
    # theta = np.array([[eta[3,0]],[eta[4,0]],[eta[5,0]]])

    # body2world[:3, :3] = body2worldRot(theta)
    # body2world[-3:,-3:] = body2worldTrans(theta)
    # globalChange = body2world@nu#compGlobalChange(theta, nu)
    # # globalChange = local_to_global_velocity(eta, nu)
    # print("eta",eta)
    # print(nu)
    # print(globalChange)
    # # if(globalChange[4,0]!=globalChange[5,0]):
    #     # input("ERROR!!")
    # eta = eta+dt*globalChange#FIXME: do we need to include nudot here too?
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
    bounds = [(-50,50),(-50,50),(-50,50),(-50,50),(-50,50),(-50,50),(-50,50),(-50,50),(-50,50),(-50,50)]

    # Perform Bayesian optimization
    result = gp_minimize(
        func=main,
        dimensions=bounds,
        acq_func='EI',  # Expected Improvement
    )

    # Retrieve the best parameters and the minimum value
    best_params = result.x
    min_value = result.fun
    print(best_params)
        # main(0.4, 0.9, 0.9)

