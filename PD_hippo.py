import numpy as np

class PD():
    def getErrors(eta, dest):
        e_d = np.sqrt((eta[0,0]-dest[0,0])**2+(eta[1,0]-dest[1,0])**2+(eta[2,0]-dest[2,0])**2)
        #the vector we want to get the angle between with...
        vecComp = np.array([1,0,0])
        #projection of my dest onto x-y
        vec = np.array([dest[0,0]-eta[0,0],dest[1,0]-eta[1,0],0])
        cos = np.dot(vecComp, vec) / (np.linalg.norm(vecComp) * np.linalg.norm(vec))
        
        e_psi = np.arccos(np.clip(cos,-1,1)) - eta[5,0]
        
        #project onto the x-z plane
        vec = np.array([dest[0,0]-eta[0,0], 0, dest[2,0]-eta[2,0]])
        cos = np.dot(vecComp, vec) / (np.linalg.norm(vecComp) * np.linalg.norm(vec))

        e_theta = np.arccos(np.clip(cos,-1,1))-eta[4,0]


        return e_d, e_theta, e_psi


    def getAccels(e_d, e_theta, e_psi, nu):
        c1 = 0.1
        c2 = 0.1
        c3 = 0.03
        if(nu[0,0] == 0):
            udot = 1
        else:
            udot = np.clip(e_d*c3/nu[0,0],-1, 1)
        #if(nu[4,0])
        if(nu[4,0] == 0 and (e_theta > 0)):
            qdot = 1
        elif(nu[4,0] == 0 and (e_theta <0)):
            qdot = -1
        elif(e_theta ==0):
            qdot = 0
        else:
            qdot = e_theta*c1/(nu[4,0]*e_d)
        
        if(nu[5,0] == 0 and (e_psi > 0)):
            rdot = 1
        elif(nu[5,0] == 0 and (e_psi <0)):
            rdot = -1
        elif(e_psi ==0):
            rdot = 0
        else:
             rdot = e_psi*c2/(nu[5,0]*e_d)
        #udot = e_d*c3/nu[0,0]
        return udot, qdot, rdot