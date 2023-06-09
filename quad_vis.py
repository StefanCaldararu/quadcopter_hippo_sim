#author: Stefan Caldararu
#Visualization of the quadcopter in 3D space...

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np

class vis(object):
    def __init__(self, l):
        self.l = l
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        
    #plot vehicle with locaiton x [x, y, z]
    #and orientation theta = [phi (roll), theta(pitch), psi (yaw)]
    def plot(self,x, R):
        
        #create our points
        p1 = np.array([[self.l], [0], [0]])#location of propeller 1
        p2 = np.array([[0], [self.l], [0]])
        p3 = np.array([[-self.l], [0], [0]])
        p4 = np.array([[0], [-self.l], [0]])
        up = np.array([[0],[0],[self.l]])
        c = np.array([[0],[0],[0]])

        #translate the points by x, the center of the vehicle.
        p1 = p1+x
        p2 = p2+x
        p3 = p3+x
        p4 = p4+x
        up = up+x
        c = c+x
        #rotate the points by theta
        p1 = R@p1
        p2 = R@p2
        p3 = R@p3
        p4 = R@p4
        up = R@up
        c = R@c
        #plot the vehicle location... wait some time
        #create the arms of the drone
        arm1 = [[p1[0,0],p3[0,0]], [p1[1,0],p3[1,0]], [p1[2,0],p3[2,0]]]
        arm2 = [[p2[0,0],p4[0,0]], [p2[1,0],p4[1,0]], [p2[2,0],p4[2,0]]]
        u = [[c[0,0],up[0,0]],[c[1,0],up[1,0]],[c[2,0],up[2,0]]]
        self.ax.cla()
        self.ax.set_xlim(-20, 20)
        self.ax.set_ylim(-20, 20)
        self.ax.set_zlim(-5, 5)
        self.ax.set_xlabel('X-axis')
        self.ax.set_ylabel('Y-axis')
        self.ax.set_zlabel('Z-axis')
        self.ax.set_title('Quadcopter in 3D Space')
        self.ax.grid(True)
        plt.gcf().canvas.mpl_connect('key_release_event', lambda event: [exit(0) if event.key == 'escape' else None])
        self.ax.plot(arm1[0], arm1[1],arm1[2], color = 'red')
        self.ax.plot(arm2[0],arm2[1],arm2[2], color = 'blue')
        self.ax.plot(u[0],u[1],u[2], color = 'green')
        plt.pause(0.1)




# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xlim(-5, 5)
# ax.set_ylim(-5, 5)
# ax.set_zlim(-5, 5)
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis')
# ax.set_zlabel('Z-axis')
# ax.set_title('Line Segment in 3D Space')
# ax.grid(True)

# time = np.arange(0,5,0.05)
# for i in time: 

#     arm1 = [[0, 0],[-0.5,0.5],[i,i]]
#     arm2 = [[-0.5, 0.5],[0,0],[i,i]]
#     ax.cla()
#     ax.set_xlim(-5, 5)
#     ax.set_ylim(-5, 5)
#     ax.set_zlim(-5, 5)
#     ax.set_xlabel('X-axis')
#     ax.set_ylabel('Y-axis')
#     ax.set_zlabel('Z-axis')
#     ax.set_title('Line Segment in 3D Space')
#     ax.grid(True)
#     plt.gcf().canvas.mpl_connect('key_release_event', lambda event: [exit(0) if event.key == 'escape' else None])
#     ax.plot(arm1[0], arm1[1],arm1[2], color = 'red')
#     ax.plot(arm2[0],arm2[1],arm2[2], color = 'blue')
#     plt.pause(0.002)

# plt.show()