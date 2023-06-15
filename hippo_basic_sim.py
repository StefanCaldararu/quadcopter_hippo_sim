#author: Stefan Caldararu
#dynamics from "Marine Craft Hydrodynamics and Motion Control" by Thor I. Fossen
#this uses the general model, and doesn't work well... problems with the C matrices
#sim the dynamics of the hippocampus, given the old BSC066 parameters.
#Ignoring Coriolis, and assuming values from
from wsgiref.simple_server import WSGIRequestHandler
import numpy as np 
from quad_vis import vis
from scipy.spatial.transform import Rotation
from PD_hippo import PD
from hippo_solver import solver
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

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
    dt = 0.001
    #orientation represented as quaternion
    eta = np.array([[0],[0],[0],[1],[0.0],[0.0],[0.0]])
    #angular velocities no quaternions...
    nu = np.zeros((6,1), dtype = float)
    #nu[0,0] = 1.0
    #nu[4,0] = 0.0
    #nu[5,0] = 1

    esc = 1600
    ESC = np.array([[esc],[esc],[esc],[esc]])
    start = 0
    end = 20
    time = np.arange(start, end, dt)
    #plot = vis(0.5)#FIXME: make this so that it can vis horizontal
    sol = solver()
    count = 0
    udot, qdot, rdot = 0,0,0
    xhist = []
    yhist = []
    zhist = []

    for t in time:
        # if(t<1):
        #     udot, qdot, rdot = 0.3, 0, 0.1
        if(count%100 == 0):
            udot = -(nu[0,0]-0.5)
            a1, b1, c1, d1, e1, a2, b2, c2, d2, e2 = 50, 26, 25, 27, -10, -46, -41, 49, 18, -19
            orientation = np.array([eta[3,0], eta[4,0], eta[5,0], eta[6,0]])
            #orientation = q2e(orientation)
            qdot = a1*eta[2,0]+b1*orientation[1]+c1*orientation[2]+d1*orientation[3]+e1*nu[4,0]#1000*compAngles(0,eta[4,0])-1*eta[3,0]
            rdot = a2*eta[1,0]+b2*orientation[1]+c2*orientation[2]+d2*orientation[3]+e2*nu[5,0]#1000*compAngles(0,eta[5,0])-1*eta[2,0]
            ESC = np.clip(sol.solve(udot, qdot, rdot, nu), 1500, 2000)
            #ESC[3,0] = 0
            R = Rotation.from_quat(orientation).as_matrix()
            #plot.plot(eta[:3], R)

        #ESC = np.array([[1850],[1850],[1900],[1900]])
        
        # print("Angle diff: ",compAngles(0,eta[4,0]))
        # print("speed is: ", nu[0,0])
        # print("desired qdot is: ", qdot)
        # print("rdot is: ", rdot)
        
        thrust = computeThrust(ESC)
        # print(ESC)
        eta, nu = motion_model(eta, nu, thrust, dt)
        #print("nu is: ",nu)
        # print("eta is: ", eta)
        # print("nu is: ", nu)
        #print("and r is: ", nu[5,0])
        eta[3,0] = fix_angle(eta[3,0])
        eta[4,0] = fix_angle(eta[4,0])
        eta[5,0] = fix_angle(eta[5,0])
        count = count+1
        xhist.append(eta[0,0])
        yhist.append(eta[1,0])
        zhist.append(eta[2,0])
        
        print(t)
            #input("Press Enter to continue...")
    
        
        #print(eta)
    fig = plt.figure(figsize = (7.5, 7.5))
    plt.plot(xhist, yhist, label = 'XY', color = 'g', linewidth = 0.3)
    plt.plot(xhist, zhist, label = 'XZ', color = 'b', linewidth = 0.3)
    plt.legend()
    plt.xlabel('Position x', fontsize=20)
    plt.ylabel('Position Y-Z', fontsize=20)
    plt.show()

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

#compute body2world
def e2q(Theta):
    roll = Theta[0,0]
    pitch = Theta[1,0]
    yaw = Theta[2,0]

    cy = np.cos(yaw * 0.5)
    sy = np.sin(yaw * 0.5)
    cp = np.cos(pitch * 0.5)
    sp = np.sin(pitch * 0.5)
    cr = np.cos(roll * 0.5)
    sr = np.sin(roll * 0.5)

    qw = cr * cp * cy + sr * sp * sy
    qx = sr * cp * cy - cr * sp * sy
    qy = cr * sp * cy + sr * cp * sy
    qz = cr * cp * sy - sr * sp * cy

    return qw, qx,qy,qz
def q2e(quaternion):
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
# def q2e(e0, e1, e2, e3):
#     R = np.array([
#         [1-2*(e2**2+e3**2), 2*(e1*e2-e3*e0), 2*(e1*e3+e2*e0)],
#         [2*(e1*e2+e3*e0), 1-2*(e1**2+e3**2), 2*(e2*e3-e1*e0)],
#         [2*(e1*e3-e2*e0), 2*(e2*e3+e1*e0), 1-2*(e1**2+e2**2)]
#     ])
#     roll = np.arctan2(R[2, 1], R[2, 2])
#     pitch = np.arctan2(-R[2, 0], np.sqrt(R[2, 1]**2 + R[2, 2]**2))
#     yaw = np.arctan2(R[1, 0], R[0, 0])
#     return roll, pitch, yaw

def closest_euler_angle(reference_euler, quaternion):
    # Convert the reference Euler angles to a rotation matrix
    reference_rotation = Rotation.from_euler('xyz', reference_euler, degrees=True)
    reference_matrix = reference_rotation.as_matrix()

    # Convert the quaternion to a rotation matrix
    quaternion_rotation = Rotation.from_quat(quaternion)
    quaternion_matrix = quaternion_rotation.as_matrix()

    # Find the rotation matrix that minimizes the Frobenius norm difference
    optimal_matrix = np.dot(np.linalg.inv(reference_matrix), quaternion_matrix)

    # Convert the optimal rotation matrix back to Euler angles
    optimal_euler = Rotation.from_matrix(optimal_matrix).as_euler('xyz', degrees=True)

    return optimal_euler

def compGlobalChange(Theta, nu):
    e0, e1, e2, e3 = e2q(Theta)

    R = np.array([
        [1-2*(e2**2+e3**2), 2*(e1*e2-e3*e0), 2*(e1*e3+e2*e0)],
        [2*(e1*e2+e3*e0), 1-2*(e1**2+e3**2), 2*(e2*e3-e1*e0)],
        [2*(e1*e3-e2*e0), 2*(e2*e3+e1*e0), 1-2*(e1**2+e2**2)]
    ])
    T=  0.5*np.array([
        [-e1, -e2, -e3],
        [e0, -e3, e2],
        [e3, e0, -e1],
        [-e2, e1, e0]
    ])

    J = np.zeros([7,6])
    J[:3, :3] = R
    J[3:, 3:] = T
    print(J)
    comp = J@nu

    # phi, theta, psi = q2e(comp[3,0], comp[4,0],comp[5,0],comp[6,0])
    euler = q2e([e0, e1, e2, e3])
    new_nu = np.array([[nu[0,0]],[nu[1,0]], [nu[2,0]], [euler[0]],[euler[1]],[euler[2]]])
    return new_nu



def local_to_global_velocity(eta, nu):
    # Extract location and orientation from eta
    x, y, z, phi, theta, psi = eta.flatten()

    # Calculate the rotation matrix
    rotation_matrix = np.array([
        [np.cos(psi) * np.cos(theta), np.sin(psi) * np.cos(theta), -np.sin(theta)],
        [np.cos(psi) * np.sin(theta) * np.sin(phi) - np.sin(psi) * np.cos(phi), 
         np.sin(psi) * np.sin(theta) * np.sin(phi) + np.cos(psi) * np.cos(phi), 
         np.cos(theta) * np.sin(phi)],
        [np.cos(psi) * np.sin(theta) * np.cos(phi) + np.sin(psi) * np.sin(phi), 
         np.sin(psi) * np.sin(theta) * np.cos(phi) - np.cos(psi) * np.sin(phi), 
         np.cos(theta) * np.cos(phi)]
    ])

    # Extract linear and angular velocities from nu
    u, v, w, p, q, r = nu.flatten()

    # Create the local velocity vector
    local_velocity = np.array([u, v, w])

    # Create the angular velocity vector
    angular_velocity = np.array([p, q, r])

    # Transform velocities to global reference frame
    global_linear_velocity = np.dot(rotation_matrix, local_velocity)
    global_angular_velocity = np.dot(rotation_matrix, angular_velocity)
    global_vel = np.array([[global_linear_velocity[0]],[global_linear_velocity[1]],[global_linear_velocity[2]],[global_angular_velocity[0]],[global_angular_velocity[1]],[global_angular_velocity[2]]])
    return global_vel


def body2worldTrans(Theta):
    phi = Theta[0,0]
    theta = Theta[1,0]

    if(np.tan(theta) == 0 or np.cos(theta) == 0):
        T = np.array([[1,0,-np.sin(theta)],
            [0, np.cos(phi), np.cos(theta)*np.sin(phi)],
            [0, -np.sin(phi),np.cos(theta)*np.cos(phi)]])
    else:

        T = np.array([
        [1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
        [0, np.cos(phi), -np.sin(phi)],
        [0, np.sin(phi)/np.tan(theta), np.cos(phi)/np.cos(theta)]
        ])
    # T = np.array([
    #     [1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
    #     [0, np.cos(phi), -np.sin(phi)],
    #     [0, np.sin(phi)/np.tan(theta), np.cos(phi)/np.cos(theta)]#FIXME: for some reason this doesn't work, div by 0... np.sin(phi)/np.tan(theta), np.cos(phi)/np.cos(theta)]
    # ])
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