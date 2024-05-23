import numpy as np
from numpy.linalg import inv
import math
import matplotlib.pyplot as plt

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

#Defining the aerodynamic Forces matrix
alpha_o = 10 * np.pi / 180
beta = 0
CL_beta = 0
CM_AC = 0
CM_AC_beta = 0
F_aerodynamic_zeroDOF = np.array([[-q * S * CL_alpha * alpha_o], [q * S * CL_alpha * (0.5 + a) * b * alpha_o]]) + np.array([[-q * S * CL_beta * beta], [q * S * CL_beta * (0.5 + a) * b * beta + q * S * CM_AC_beta * (2*b) * beta]]) + np.array([[0], [q*S*CM_AC*(2*b)]])


#Monolithic Equation Solver

X_monolithic = np.dot(inv(Structural_stiffness_matrix - Aerodynamic_stiffness_matrix), F_aerodynamic_zeroDOF) #Monolithic Solution Code
F_aerodynamic_final = F_aerodynamic_zeroDOF + np.dot(Aerodynamic_stiffness_matrix, X_monolithic)


#Structural Solver
X = np.zeros((2,1))
error = 5
alpha = [0]
h = [0]
counter = [0]
i = 0
while(error>1e-6):
    F_aerodynamic_structure = np.array([[-q * S * CL_alpha * X[1][0]], [q*S*CL_alpha* (0.5 + a)*b * X[1][0]]]) + np.array([[-q * S * CL_alpha * alpha_o], [q * S * CL_alpha * (0.5 + a) * b * alpha_o]]) + np.array([[-q * S * CL_beta * beta], [q * S * CL_beta * (0.5 + a) * b * beta + q * S * CM_AC_beta * (2*b) * beta]]) + np.array([[0], [q*S*CM_AC*(2*b)]])
    X_new = np.dot(inv(Structural_stiffness_matrix), F_aerodynamic_structure)
    error = math.sqrt((X_new[0][0]-X[0][0])**2+(X_new[1][0]-X[1][0])**2)
    i += 1
    h.append(X[0][0])
    alpha.append(X[1][0])
    counter.append(i)
    X = X_new
    
#Plot for structural solver & the monolithic equation
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(counter,h, label='From structural solver')
ax1.axhline(X_monolithic[0][0], color='r', label='Monolithic Solution')
ax1.set_xlabel('Iteration number')
ax1.set_ylabel('heave displacement')
ax1.set_title('Monolithic vs heave displacement')
ax1.legend()

ax2.plot(counter, alpha, label='Alpha change')
ax2.axhline(X_monolithic[1][0], color='r', label='Monolithic Solution')
ax2.set_xlabel('Iteration number')
ax2.set_ylabel('Angle of attack')
ax2.set_title('Monolithic vs angle of attack')
ax2.legend()

plt.tight_layout()
plt.show()

#Trimming aeroelasticity for a typical section
X = np.zeros((2,1))
alpha = [0]
h = [0]
i = 0
Lift = q * S * CL_alpha * alpha_o
Weight = m_airfoil * 9.81 
while(Lift - Weight > 1):
    counter = 0
    while(counter < 10):
        F_aerodynamic_structure = np.array([[-q * S * CL_alpha * X[1][0]], [q*S*CL_alpha* (0.5 + a)*b * X[1][0]]]) + np.array([[-q * S * CL_alpha * alpha_o], [q * S * CL_alpha * (0.5 + a) * b * alpha_o]]) + np.array([[-q * S * CL_beta * beta], [q * S * CL_beta * (0.5 + a) * b * beta + q * S * CM_AC_beta * (2*b) * beta]]) + np.array([[0], [q*S*CM_AC*(2*b)]])
        X_new = np.dot(inv(Structural_stiffness_matrix), F_aerodynamic_structure)
        X = X_new
        counter += 1
    Lift = q * S * CL_alpha * (alpha_o + X[1][0])
    if(Lift < Weight):
        alpha_o = alpha_o + alpha_o * 0.01
    elif(Lift > Weight):
        alpha_o = alpha_o - alpha_o * 0.01  

print(Lift)
print(Weight)
print(X)