import numpy as np
from sympy.matrices import Matrix
import sympy as sp

cos = np.cos
sin = np.sin
pi = np.pi

#Delete -0 values
def del_minz(R):
    for i in range(R.shape[0]):  
        for j in range(R.shape[1]):
            if np.round(R[i,j],3) == -0.000:
                R[i,j] = 0
    return R

#Unitary vector from vector
def unit_v(v):
    return v / np.linalg.norm(v)

#Euler angles from rotation matrix R
def R_eulerangzyz(R,s):
    
    R = del_minz(R)

    if s == '+':
        fi2 = np.arctan2(np.sqrt(R[0,2]**2+R[1,2]**2),R[2,2])
    elif s == '-':
        fi2 = np.arctan2(-np.sqrt(R[0,2]**2+R[1,2]**2),R[2,2])

    if sin(fi2) != 0:
        fi1 = np.arctan2(R[1,2]/sin(fi2),R[0,2]/sin(fi2))
        fi3 = np.arctan2(R[2,1]/sin(fi2),-R[2,0]/sin(fi2))
    else:
        print("Singularity")
        if (fi2/pi)%2 == 0:
            f13 = np.arctan2(R[1,0],R[0,0])
            return fi2, f13, 'fi1+fi3'
        else: 
            f31 = np.arctan2(R[0,1],R[1,1])
            return fi2, f31, 'fi3-fi1'

#Roll, pitch, yaw angles from rotation matrix R
def R_rpy(R,s):
    
    R = del_minz(R)

    if R[2,0] == 0:
        if s == '+':
            fi_p = np.arctan2(R[2,0],np.sqrt(R[2,1]**2+R[2,2]**2))
        elif s == '-':
            fi_p = np.arctan2(R[2,0],-np.sqrt(R[2,1]**2+R[2,2]**2))
    else: 
        if s == '+':
            fi_p = np.arctan2(-R[2,0],np.sqrt(R[2,1]**2+R[2,2]**2))
        elif s == '-':
            fi_p = np.arctan2(-R[2,0],-np.sqrt(R[2,1]**2+R[2,2]**2))

    if cos(fi_p) != 0:
        fi_r = np.arctan2(R[1,0]/cos(fi_p),R[0,0]/cos(fi_p))
        fi_y = np.arctan2(R[2,1]/cos(fi_p),R[2,2]/cos(fi_p))
        return fi_r, fi_p, fi_y
    else:
        print("Singularity")
        if fi_p == pi/2:
            fry = np.arctan2(R[1,2],R[0,2])
            return fi_p, fry, 'fi_r-fi_y'
        elif fi_p == -pi/2: 
            f31 = np.arctan2(-R[0,1],R[1,1])
            return fi_p, fyr, 'fi_r+fi_y'

#Axis/Angle
#Rotation matrix from U, th
def Rf_uth(u,th,au='r'):
    if au == 'd':
        th = np.deg2rad(th)

    u = unit_v(u)

    u_a = np.array([[0, -u[0,2], u[0,1]],
                    [u[0,2], 0, -u[0,0]],
                    [-u[0,1], u[0,0], 0]])

    R = np.eye(3) + u_a*sin(th) + np.dot(u_a,u_a)*(1-cos(th))
    R = del_minz(R)
    return R

#U, th from rotation matrix
def u_th_R(R,s):
    
    R = del_minz(R)

    if s == '+':
        th = np.arctan2(np.sqrt((R[1,0]-R[0,1])**2+(R[2,0]-R[0,2])**2+(R[2,1]-R[1,2])**2)/2,(np.trace(R)-1)/2)
    elif s == '-':
        th = np.arctan2(-np.sqrt((R[1,0]-R[0,1])**2+(R[2,0]-R[0,2])**2+(R[2,1]-R[1,2])**2)/2,(np.trace(R)-1)/2)
    
    if th == 0 or th == -pi or th == pi:
        print("Singularity")
    else:
        u = (1/(2*sin(th)))*np.array([[R[2,1]-R[1,2],R[0,2]-R[2,0],R[1,0]-R[0,1]]])
        return u,th       

#Unit Cuaternions
class  cuaternion:

    def __init__(self, w=0, ex=0, ey=0, ez=0):
        self.w = w
        self.ex = ex
        self.ey = ey
        self.ez = ez
        self.e_v = np.array([[self.ex],[self.ey],[self.ez]])

    def del_minz_c(self):
    
        if np.round(self.w,3) == -0.000:
            self.w = 0
        if np.round(self.ex,3) == -0.000:
            self.ex = 0
        if np.round(self.ey,3) == -0.000:
            self.ey = 0
        if np.round(self.ez,3) == -0.000:
            self.ez = 0
        return cuaternion(self.w, self.ex, self.ey, self.ez)
    
    def __str__(self):
        return del_minz(np.array([[self.w, self.ex, self.ey, self.ez]])).__str__()

    def __mul__(self,other):    
        cmul_w = float(self.w*other.w - np.dot(self.e_v.T,other.e_v))
        print(cmul_w)
        cmul_e = self.w*other.e_v + other.w*self.e_v + np.cross(self.e_v.T,other.e_v.T).T
        print(cmul_e)
        cmul_ex = cmul_e[0,0]
        cmul_ey = cmul_e[1,0]
        cmul_ez = cmul_e[2,0]      
        return cuaternion(cmul_w, cmul_ex, cmul_ey, cmul_ez).del_minz_c()

    def inv(self):
        return cuaternion(self.w,-self.ex,-self.ey,-self.ez).del_minz_c()

    def eye():
        return cuaternion(1,0,0,0)

    def to_nparray(self):
        return np.array([[self.w, self.ex, self.ey, self.ez]])

#Cuaternion from u,th
def cuat_uth(u,th,au='r'):
    
    if au == 'd':
        th = np.deg2rad(th)

    u = unit_v(u)

    w = cos(th/2)
    ex = u[0,0]*sin(th/2)
    ey = u[0,1]*sin(th/2)
    ez = u[0,2]*sin(th/2) 

    return cuaternion(w,ex,ey,ez).del_minz_c()

#Cuaternion from angles of x,y,z
def cuat_rotxyz(th,c,au='r'):

    if au == 'd':
        th = np.deg2rad(th)

    w = cos(pi/2)

    if c == 'x':
       ex = sin(pi/2)
       ey = 0
       ez = 0
    elif c == 'x':
       ex = 0
       ey = sin(pi/2)
       ez = 0
    elif c == 'x':
       ex = 0
       ey = 0
       ez = sin(pi/2)

    return cuaternion(w,ex,ey,ez).del_minz_c()

#Axis/angle from cuaternion
def u_th_C(C):
    
    th = 2*np.arctan2(np.linalg.norm(C.e_v),C.w)
    u = (C.e_v/np.linalg.norm(C.e_v)).T
    
    if (np.round(np.linalg.norm(C.e_v),3) == np.round(sin(th/2),3) and np.round(sin(th/2),3) == 0):
        print("Rotation angle = 0")
    else: 
        return u, th

#Rotation from cuaternion
def Rf_C(C):

    [u,th] = u_th_C(C)

    R = Rf_uth(u,th)
    R = del_minz(R)
    return R

#Cuaternion from Rotation matrix
def C_R(R):
    
    w = 0.5*np.sqrt(1 + R[0,0] + R[1,1] + R[2,2])
    ex = (1/(4*w))*(R[2,1]-R[1,2])
    ey = (1/(4*w))*(R[0,2]-R[2,0])
    ez = (1/(4*w))*(R[1,0]-R[0,1])
    
    C = cuaternion(w,ex,ey,ez).del_minz_c()
    return C

#Vector rotations using Parameterizations
#Vector rotation given u,th
def rv_uth(v,u,th,au='r'):
    
    if au == 'd':
        th = np.deg2rad(th)

    u = unit_v(u)
    
    v_rot = v*cos(th) + np.cross(u,v.T).T*sin(th) + u.T*(np.dot(u,v))*(1-cos(th))

    return v_rot

#Vector rotation given C
def rv_C(v,C):

    w_v = 0
    ex_v = v[0,0]
    ey_v = v[1,0]
    ez_v = v[2,0]
    c_v = cuaternion(w_v,ex_v,ey_v,ez_v)
    
    c_v_rot = (C*c_v)*C.inv()

    v_rot = np.array([[c_v_rot.ex],[c_v_rot.ey],[c_v_rot.ez]])

    return v_rot 

