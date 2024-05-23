import numpy as np
import math
from scipy.linalg import eig

L_wing = 16

# Mass of air foil and flap
m_airfoil = 0.75 * 16
m_flap = 0
m = m_airfoil + m_flap

#Typical section parameters
b = 0.5 #half chord length
a = 0 #from the center to the elastic axis
c = 0 #from the center to the hinge axis
X_theta = 0 #between elastic axis and the center of gravity of air foil
X_beta = 0 #between hinge axis and the center of gravity of flap
S = 2*b* 1

#moment of intertias
I_cg_theta = 0.1 * L_wing
I_cg_beta = 0
W = 30 #probably the omega
S_theta = m_airfoil * X_theta * b + m_flap * (c-a+X_beta) * b
S_beta = m_flap * X_beta * b
I_theta = I_cg_theta + I_cg_beta + m_airfoil * (X_theta * b) ** 2 + m_flap * (c-a+X_beta)**2 * b**2 
I_beta = I_cg_beta + m_flap * (X_beta * b) ** 2


mass_matrix = np.array([[m,S_theta,S_beta],[S_theta,I_theta,(c-a)*b*S_beta+I_beta],[S_beta,(c-a)*b*S_beta+I_beta,I_beta]])

#stiffness parameters
K_h = 2e4 * 3/(0.63*L_wing)**3
K_theta = 1e4 /(0.63*L_wing)
K_beta = 0

stiffness_matrix = np.array([[K_h,0,0],[0,K_theta,0],[0,0,K_beta]])

#eigen vectors and eigen values of coupled
eigvals, eigvecs = eig(stiffness_matrix, mass_matrix)

#eigen frequencies
eigen_frequencies_coupled =  np.sqrt(eigvals)
print(eigen_frequencies_coupled)

#uncoupled eigen vectors and eigen values 
mass_matrix_uncoupled = np.array([[m,0,0],[0,I_theta,0],[0,0,I_beta]])
eigvals_un, eigvecs_un = eig(stiffness_matrix, mass_matrix_uncoupled)

#uncoupled eigen frequencies
eigen_frequencies_uncoupled =  np.sqrt(eigvals_un)

#Aerodynamic stiffness parameters
rho = 1.2
V = 30
q = 0.5 * rho * V**2
CLalpha = 2 *  math.pi
mass_aero_matrix = np.array([[m,S_theta],[S_theta,I_theta]])
Stiffness_aero_matrix = np.array([[0, -q * S * CLalpha],[0, q * S * CLalpha * (0.5 + a) * b]])
Stiffness_stru_matrix = np.array([[K_h,0],[0,K_theta]])
eff_stiffness_matrix = Stiffness_stru_matrix - Stiffness_aero_matrix
eigvals_aero, eigvecs_aero = eig(eff_stiffness_matrix, mass_aero_matrix)

eigen_frequencies_aero =  np.sqrt(eigvals_aero)
print(eigen_frequencies_aero) 
