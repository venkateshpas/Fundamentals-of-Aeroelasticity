import numpy as np
from numpy.linalg import inv
import math
import matplotlib.pyplot as plt
import sympy as sp

alpha_o = sp.symbols('alpha_o', real = True)
# Mass of air foil and flap
m_airfoil = 0.75 * 16
m_flap = 0
m = m_airfoil + m_flap

#Properties of Typical Section
L_wing = 16
b = 0.5 #half chord length
a = 0 #from the center to the elastic axis
c = 0 #from the center to the hinge axis
X_theta = 0 #between elastic axis and the center of gravity of air foil
X_beta = 0 #between hinge axis and the center of gravity of flap
S = 2*b* 1

#Defining the stiffness matrix
L = 0.63 * L_wing
K_h = 2e4 *3/ L**3
K_theta = 1e4/L
Structural_stiffness_matrix = np.array([[K_h,0],[0,K_theta]])

#Defining the aerodynamic stiffness matrix
rho = 1.2
V = 30
q = 0.5 * rho * V**2
CL_alpha  = 2*np.pi
Aerodynamic_stiffness_matrix = np.array([[0, -q * S * CL_alpha],[0, q*S*CL_alpha*(0.5 + a) * b]])

Equivalent_stiffness_matrix = Structural_stiffness_matrix - Aerodynamic_stiffness_matrix

#Defining the aerodynamic Forces matrix
#alpha_o = 10 * np.pi / 180
beta = 0
CL_beta = 0
CM_AC = 0
CM_AC_beta = 0
F_aerodynamic_zeroDOF = np.array([[-q * S * CL_alpha * alpha_o], [q * S * CL_alpha * (0.5 + a) * b * alpha_o]]) + np.array([[-q * S * CL_beta * beta], [q * S * CL_beta * (0.5 + a) * b * beta + q * S * CM_AC_beta * (2*b) * beta]]) + np.array([[0], [q*S*CM_AC*(2*b)]])

d_alpha_F_aero = np.zeros([2,1])
d_alpha_F_aero[0] = -1*sp.diff(-q * S * CL_alpha * alpha_o, alpha_o)
d_alpha_F_aero[1] = -1*sp.diff(q * S * CL_alpha * (0.5 + a) * b * alpha_o, alpha_o)

Equivalent_stiffness = np.vstack([np.hstack([Equivalent_stiffness_matrix,d_alpha_F_aero]),[0,q*S*CL_alpha,q*S*CL_alpha]])
Resultant_forces = np.array([[0],[0],[m_airfoil*9.81]])

X_monolithic = np.dot(inv(Equivalent_stiffness), Resultant_forces) #Monolithic Solution Code
print(X_monolithic)