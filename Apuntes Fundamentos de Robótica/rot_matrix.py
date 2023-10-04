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

#Round and select type
def printr(x,t='n',a=3):
    if t != 'n':
        if t == 'd':
            print(np.round(np.rad2deg(x),a))
        elif t == 'r':
            print(np.round(np.deg2rad(x),a))
    else:
         print(np.round(x,a))

#Canonic rotation matrixes
def rot(ang,c,ua='r'):

    if ua == 'd':
        ang = np.deg2rad(ang)
    
    if c == 'x':
        Rx = np.array([[1, 0, 0],
                       [0, cos(ang) ,-sin(ang)],
                       [0, sin(ang), cos(ang)]])
        return Rx

    elif c == 'y':
        Ry = np.array([[cos(ang), 0, sin(ang)],
                       [0, 1 ,0],
                       [-sin(ang), 0, cos(ang)]])
        return Ry

    elif c == 'z':
        Rz = np.array([[cos(ang), -sin(ang), 0],
                       [sin(ang), cos(ang), 0],
                       [0, 0, 1]])
        return Rz

#Symbolic canonic rotation matrixes
def srot(ang,c):

    if c == 'x':
        Rx = sp.Matrix([[1, 0, 0],
                       [0, coss(ang) ,-sins(ang)],
                       [0, sins(ang), coss(ang)]])
        return Rx

    elif c == 'y':
        Ry = sp.Matrix([[coss(ang), 0, sins(ang)],
                        [0, 1, 0],
                        [-sins(ang), 0, coss(ang)]])
        return Ry

    elif c == 'z':
        Rz = sp.Matrix([[coss(ang), -sins(ang), 0],
                       [sins(ang), coss(ang), 0],
                       [0,0,1]])
        return Rz

#Checks if a matrix is a rotation matrix
def check_rotm(R,p=3):
    if ((np.round(np.dot(R,R.T),p) == np.eye(3)).all() and np.round(np.linalg.det(R),p) == 1):
        return True
    else:
        return False

#Define point in the space
def def_point(x,y,z):
    return np.array([[x],[y],[z]])
