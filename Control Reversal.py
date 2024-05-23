import numpy as np
import math
from scipy.linalg import eig
import matplotlib.pyplot as plt

# Mass of air foil and flap
m_airfoil = 0.75 * 16
m_flap = 0
m = m_airfoil + m_flap

#Properties of Typical Section
L_wing = 16
b = 0.5 #half chord length
a = -1 #from the center to the elastic axis
c = 0.4 #from the center to the hinge axis
X_theta = 0 #between elastic axis and the center of gravity of air foil
X_beta = 0 #between hinge axis and the center of gravity of flap
S = 2*b* 1

#Defining the stiffness matrix
L = 0.63 * L_wing
K_h = 2e4 *3/ L**3
K_theta = 1e4/L
K_structural_matrix = np.array([[K_h, 0],[0, K_theta]])

#Defining the aerodynamic stiffness matrix
rho = 1.2
V = 100
q = 0.5 * rho * V**2
CL_alpha  = 2*np.pi
CL_beta = CL_alpha/20
CM_AC_beta = -0.1
K_aerodynamic_stiffness = np.array([[0, -q*S*CL_alpha], [0, q*S*CL_alpha*(0.5+a)*b]])

#control_reversal effectiveness
q_divergence = K_theta/(S*CL_alpha*(0.5 + a)*b)
q_reversal = -CL_beta*K_theta/(CL_alpha*S*(2*b)*CM_AC_beta)
V = 0.1
q = 0.5 * rho * V**2
q_plot = []
effectiveness_plot = []
Lift_flexible_plot =  []
Lift_rigid_plot = []

while q<1.5*q_reversal:
    th_be = (q*S*CL_beta*(0.5+a)*b+CM_AC_beta*(2*b)*q*S)/(K_theta-q*S*CL_alpha*(0.5+a)*b)
    Lift_flexible = q*S*(CL_alpha*th_be+CL_beta)
    Lift_rigid = q*S*CL_beta
    effectiveness = 1 + th_be * CL_alpha/CL_beta
    q_plot.append(q)
    effectiveness_plot.append(effectiveness)
    Lift_flexible_plot.append(Lift_flexible)
    Lift_rigid_plot.append(Lift_rigid)
    V = V + 1
    q = 0.5 * rho * V**2

    
#Plot for Effectiveness Curve
plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(q_plot, effectiveness_plot, color='blue',label='effectiveness')
#plt.plot(q_plot, Lift_flexible_plot, color='orange',label='Lift_flexible')
#plt.plot(q_plot, Lift_rigid_plot, color='pink',label='Lift_rigid')
plt.xlabel('divergence pressure')
plt.ylabel('effectiveness')
plt.title('Plot for control reversal')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()
