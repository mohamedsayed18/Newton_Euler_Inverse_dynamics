#newton euler approach
# RRP robot
import sympy as sp
import numpy as np
import matplotlib.py as plt


#symbols
g, a1, a2, a3 = sp.symbols('g, a1, a2, a3')
k0, k1, k2, k3 = sp.symbols('k0,k1,k2,k3')
kr = [k0, k1, k2, k3]

#joints types
joints = ['R','R','P']

#TODO define q, q_dot, q_ddot symbols
q0,     q1,     q2,     q3 = sp.symbols('q0, q1, q2, q3')
q0_dot, q1_dot, q2_dot, q3_dot = sp.symbols('q0_dot, q1_dot, q2_dot, q3_dot')
q0_ddot, q1_ddot, q2_ddot, q3_ddot = sp.symbols('q0_ddot, q1_ddot, q2_ddot, q3_ddot')

q = [q0,q1,q2,q3]
q_dot = [q0_dot, q1_dot, q2_dot, q3_dot]
q_ddot = [q0_ddot, q1_ddot, q2_ddot, q3_ddot]

Radius_center = []
Radius = []
Forces = []
Moments = []
Acc_center = []
Omega_dot_mass = []
m1 = 2
m2 = 3
m3 = 1
M = [m1, m2, m3]

#distance to centre of mass
lc = [0,0,0]

#initial conditions for the velocities and accelerations
acc = np.array([0, g, 0]).reshape(-1,1)
Acc = [acc]

omega_0 = np.zeros(3,1)
omega_dot_0 = np.zeros(3,1)
Omega=[omega_0]
Omega_dot = [omega_dot_0]

#TODO constant vectors
Radius.append([-lc[1],0,0], [0,0,lc[2]], [0,0,a2+a3+q3])

#Rotations
R_1 = np.array([[sp.cos(q1),-sp.sin(q1), 0],
                [sp.cos(q1),-sp.sin(q1), 0]])

R_2 = np.array([1,0,0],
               [0,0,1],
               [0,-1,0])

R_3 = np.identity(3)
R_4 = np.identity(3)
Rototions = [R_1, R_2, R_3, R_4]

#axis of Rotations
z0 = [0,0,1]
z1 = [0,1,0]
z2 = [0,0,1]
z3 = [0,0,1]
Z = [z0, z1, z2, z3]

#TODO forward recursion
#compute velocities, acceleration, from link 1 to link 3
#eq 7.107 --- 7.111
#velocities
for i in range(1,4):
    if(joints[i-1]=='R')    #Revolute
        Omega.append(R[i].dot(Omega[i-1]+q1_dot.dot(z0))) #Revolute
        Omega_dot.append(R[i].dot(Omega_dot[i-1]+q1_dot.dot(z0))) #Revolute
        Acc.append(R[i].dot(Acc[i-1]) + np.cross(Omega_dot[i],Radius[i-1]) + np.cross(Omega_dot[i],np.cross(Omega_dot[i],Radius[i-1])))
        Acc_center.append(Acc[i] + np.cross(omega_dot[i],Radius_center[i]) + np.cross(Omega[i],np.cross(Omega[i],Radius_center[i])))
        Omega_dot_mass.append(omega_dot[i-1] + kr[i]*q_ddot[i].dot(z[i]))
    else    #prismatic
        Omega.append(R[i].dot(Omega[i-1]))
        Omega_dot.append(R[i].dot(Omega_dot[i-1]))
        Acc.append(R[i].dot(Acc[i-1]+q_ddot[i]*Z[0]) + np.cross(2*q_dot[i]*omega[i],Rotations[i].dot(Z[0])) + np.cross(Omega_dot[i],Radius[i-1]) + np.cross(Omega[i],np.cross(Omega[i],Radius[i-1])))
        Acc_center.append(Acc[i] + np.cross(omega_dot[i],Radius_center[i]) + np.cross(Omega[i],np.cross(Omega[i],Radius_center[i])))
        Omega_dot_mass.append(omega_dot[i-1] + kr[i]*q_ddot[i].dot(z[i]))

# TODO: back recursion
#compute forces, Moments
for i in range(4,1,-1):
    Forces.append(Rotations[i+1].dot(Forces[i+1])+M[i].dot(Acc_center[i]))
    Moments.append(np.cross(-Forces[i],Radius[i+1]+Radius_center[i])+Radius[i+1].dot(Forces[i+1]))
