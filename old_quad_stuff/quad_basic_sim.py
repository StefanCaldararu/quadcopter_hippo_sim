#author: Stefan Caldararu
#sim of quadcopter dynamics. possible to modify later to reflect hippocampus dynamics. Later sims will revolve around motor failure and maintaining the control of the vehicle. 
#ref: https://ieeexplore.ieee.org/document/6906588

from quad_vis import vis
import numpy as np
from old_quad_stuff.controller import controller


def main():
    dt = 0.05
    x = np.array([[0],[0.008],[3]])
    xdot = np.array([[-0.158],[0],[0]])
    theta = np.array([[1.15],[0],[0]])
    thetadot = np.array([[0],[5.69],[18.89]])
    # x = np.array([[0],[0],[3]])
    # xdot = np.array([[0],[0],[0]])
    # theta = np.array([[0],[0],[0]])
    # thetadot = np.array([[3],[0],[-15]])
    angv = 430.3815
    u = np.array([[angv],[angv],[angv],[angv]])
    start = 0
    end = 20
    time = np.arange(start, end, dt)
    xhist = [x]
    plt = vis(0.5)
    cont = controller()
    for t in time:
        # if(x[2,0]>3.5):
        #     angv = 400
        # elif(x[2,0]<1.0):
        #     angv = 500
        o1, o2, o3, o4 = 556.961, 393.218, 556.961, 0#cont.solve(thetadot)
        #o1, o2, o3, o4 = 0, 712.409, 0, 511.409#618.55
        u = np.array([[o1],[o2],[o3],[o4]])
        #print(angv)
        x, xdot, theta, thetadot = motion_model(x, xdot, theta, thetadot, u, dt)

        plt.plot(x, theta)
        #xhist.append(x)
        #print(x)
        print(t)
    

    
#x is our location (global), x = [x, y, z]
#theta is our orientation (global), theta = [phi (roll), theta(pitch), psi (yaw)]
#thetadot is angular velocity, [p,q,r]
#u is our control input, [omega_1, omega_2, omega_3, omega_4] which for now just relate to the angular velocities of the motors
def motion_model(x, xdot, theta, thetadot, u, dt):
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
    I_B_zz = 1#vehicle body inertia...#TODO
    #I_P_xx = 0#vehicle propeller inertia
    I_P_zz = 1.5*0.00001#vehicle propeller inertia
    #total inertias, sum of propeller and body
    I_T_xx = 3.2*0.001#KGm**2
    #I_T_yy = I_B_xx+I_P_xx*4#FIXME: make sure I_B_xx = I_B_yy
    I_T_zz = 5.5*0.001#KGm**2
    I_B_zz = I_T_zz-4*I_P_zz#vehicle body inertia...#TODO
    #TODO: define constants above appropriately.

    #given u, comopute how much force each of the motors is generating.
    R = calcR(theta)
    f = np.zeros((4,1))
    f_tot = 0
    for i in range(0,4):
        f[i,0] = k_f*u[i,0]**2
        f_tot = f_tot+f[i,0]
    d_ddot = R@np.array([
        [0],
        [0],
        [f_tot]#our change in center of mass...
    ])/m+g #this is our acceleration... our change in center of mass...global reference frame
    omega2_1 = u[0,0]**2
    omega2_2 = u[1,0]**2
    omega2_3 = u[2,0]**2
    omega2_4 = u[3,0]**2
    #change in angular velocity...
    pdot = (k_f*(omega2_2-omega2_4)*l-(I_T_zz-I_T_xx)*thetadot[1,0]*thetadot[2,0]-I_P_zz*thetadot[1,0]*(u[0,0]+u[1,0]+u[2,0]+u[3,0]))/I_B_xx
    qdot = (k_f*(omega2_3-omega2_1)*l+(I_T_zz-I_T_xx)*thetadot[0,0]*thetadot[2,0]+I_P_zz*thetadot[0,0]*(u[0,0]+u[1,0]+u[2,0]+u[3,0]))/I_B_xx
    rdot = (-gamma*thetadot[2,0]+k_tau*k_f*(omega2_1-omega2_2+omega2_3-omega2_3))/I_B_zz

    thetaddot = np.array([
        [pdot],
        [qdot],
        [rdot]
    ])
    print(thetaddot)
    thetadot = thetadot+dt*thetaddot
    theta = theta+thetadot*dt
    xdot = xdot+d_ddot*dt
    x = x+xdot*dt

    return x, xdot, theta, thetadot





#FROM THE PREV ARTICLE....
    # m = 1 #Mass of vehicle??
    # g =  9.81 #Gravity
    # k_D = 
    # #TODO: define constants above...
    # R = calcR(theta)
    # omega = np.array([  
    #     [1,0,-np.sin(theta[1,0])],
    #     [0, np.cos(theta[0,0]), np.cos(theta[1,0])*np.sin(theta[0,0])],
    #     [0, -np.sin(theta[0,0]), np.cos(theta[1,0])*np.sin(theta[0,0])]
    # ])@thetadot
    # #the friction force...
    # F_D = 
    # T_B = 
    # #linear translation...
    # xdotdot = (1/m)*(np.array([
    #     [0],
    #     [0],
    #     [-m*g]
    # ])+R@T_B + F_D)


#function to take vectors from the body reference frame to the inertial/global reference frame.
#theta is our orientation (global), theta = [phi (roll), theta(pitch), psi (yaw)]
def calcR(theta):
    R = np.array([  
        [np.cos(theta[0,0])*np.cos(theta[2,0])-np.cos(theta[1,0])*np.sin(theta[0,0])*np.sin(theta[2,0]), -np.cos(theta[2,0])*np.sin(theta[1,0])-np.cos(theta[0,0])*np.cos(theta[1,0])*np.sin(theta[2,0]), np.sin(theta[1,0])*np.sin(theta[2,0])],
        [np.cos(theta[1,0])*np.cos(theta[2,0])*np.sin(theta[0,0])+np.cos(theta[0,0])*np.sin(theta[2,0]), np.cos(theta[0,0])*np.cos(theta[2,0])*np.cos(theta[1,0])-np.sin(theta[0,0])*np.sin(theta[2,0]), -np.cos(theta[2,0])*np.sin(theta[1,0])],
        [np.sin(theta[0,0])*np.sin(theta[1,0]), np.cos(theta[0,0])*np.sin(theta[1,0]), np.cos(theta[1,0])]
    ])
    return R
if __name__ == '__main__':
    main()