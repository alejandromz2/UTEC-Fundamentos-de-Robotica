import numpy as np
from sympy.matrices import Matrix
import sympy as sp

#Redefine numpy functions
cos = np.cos
sin = np.sin
pi = np.pi

#Redefine sympy functions
coss = sp.cos
sins = sp.sin

na = np.array

#T given t,R
def T_tR(t,R):
    T = np.vstack((np.hstack((R,t)),np.array([[0, 0, 0, 1]])))
    return T

#R from T
def R_T(T):
    R = T[:3,:3]
    return R

#Pure rotation
def tf_rot(ang,c,ua='r'):

    if ua == 'd':
        ang = np.deg2rad(ang)

    if c == 'x':
        Trotx = np.array([[1, 0, 0, 0],
                         [0, cos(ang), -sin(ang), 0],
                         [0, sin(ang), cos(ang), 0],
                         [0, 0, 0, 1]])
        return Trotx

    elif c == 'y':
        Troty = np.array([[cos(ang), 0, sin(ang), 0],
                         [0, 1 ,0, 0],
                         [-sin(ang), 0, cos(ang), 0],
                         [0, 0, 0, 1]])
        return Troty

    elif c == 'z':
        Trotz = np.array([[cos(ang), -sin(ang), 0, 0],
                         [sin(ang), cos(ang), 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        return Trotz

#Symbolic pure rotation
def stf_rot(ang,c):
   
    if c == 'x':
        Trotx = sp.Matrix([[1, 0, 0, 0],
                         [0, coss(ang), -sins(ang), 0],
                         [0, sins(ang), coss(ang), 0],
                         [0, 0, 0, 1]])
        return Trotx

    elif c == 'y':
        Troty = sp.Matrix([[coss(ang), 0, sins(ang), 0],
                         [0, 1 ,0, 0],
                         [-sins(ang), 0, coss(ang), 0],
                         [0, 0, 0, 1]])
        return Troty

    elif c == 'z':
        Trotz = sp.Matrix([[coss(ang), -sins(ang), 0, 0],
                         [sins(ang), coss(ang), 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        return Trotz

#Pure traslation
def tf_tr(d,c):

    if c == 'x':
        Ttrx = np.array([[1, 0, 0, d],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        return Ttrx

    elif c == 'y':
        Ttry = np.array([[1, 0, 0, 0],
                         [0, 1, 0, d],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        return Ttry

    elif c == 'z':
        Ttrz = np.array([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, d],
                         [0, 0, 0, 1]])
        return Ttrz

#Symbolic Pure traslation
def stf_tr(d,c):

    if c == 'x':
        Ttrx = sp.Matrix([[1, 0, 0, d],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        return Ttrx

    elif c == 'y':
        Ttry = sp.Matrix([[1, 0, 0, 0],
                         [0, 1, 0, d],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        return Ttry

    elif c == 'z':
        Ttrz = sp.Matrix([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, d],
                         [0, 0, 0, 1]])
        return Ttrz

#Complete pure traslation
def tf_trc(dx,dy,dz):
    Ttr = np.array([[1, 0, 0, dx],
                     [0, 1, 0, dy],
                     [0, 0, 1, dz],
                     [0, 0, 0, 1]])
    return Ttr

#Symbolic complete complete pure traslation
def stf_trc(dx,dy,dz):
    Ttr = sp.Matrix([[1, 0, 0, dx],
                     [0, 1, 0, dy],
                     [0, 0, 1, dz],
                     [0, 0, 0, 1]])
    return Ttr

#Define homogeneus coordinate
def def_hc(x,y,z):
    return np.array([[x],
                     [y],
                     [z],
                     [1]])

#New vector from vector and Homogeneus matrix rotation
def np_pT(p1,T):
    p1 = np.array([[p1[0,0]],[p1[1,0]],[p1[2,0]],[1]])
    p2 = np.dot(T,p1)
    p2 = p2[:3,:]
    return p2
