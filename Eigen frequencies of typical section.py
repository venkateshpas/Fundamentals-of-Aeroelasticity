import numpy as np
import math
from scipy.linalg import eig

L_wing = 32

# Mass of air foil and flap
m_airfoil = 20000/L_wing
m_flap = 0
m = m_airfoil + m_flap

#Typical section parameters
b = 4 #half chord length
a = -0.2 #from the center to the elastic axis
c = 1 #from the center to the hinge axis
X_theta = 0.2 #between elastic axis and the center of gravity of air foil
X_beta = 0 #between hinge axis and the center of gravity of flap
S = 2*b* 1

#moment of intertias
I_cg_theta = m_airfoil
I_cg_beta = 0
W = 30 #probably the omega
S_theta = m_airfoil * X_theta * b + m_flap * (c-a+X_beta) * b
S_beta = m_flap * X_beta * b
I_theta = I_cg_theta + I_cg_beta + m_airfoil * (X_theta * b) ** 2 + m_flap * (c-a+X_beta)**2 * b**2 
I_beta = I_cg_beta + m_flap * (X_beta * b) ** 2


mass_matrix = np.array([[m,S_theta,S_beta],[S_theta,I_theta,(c-a)*b*S_beta+I_beta],[S_beta,(c-a)*b*S_beta+I_beta,I_beta]])

#stiffness parameters
K_h = 4 * 10**5
K_theta = 30 * 10**6
K_beta = 0

stiffness_matrix = np.array([[K_h,0,0],[0,K_theta,0],[0,0,K_beta]])

#eigen vectors and eigen values of coupled
eigvals, eigvecs = eig(stiffness_matrix, mass_matrix)

#eigen frequencies
eigen_frequencies_coupled =  np.sqrt(eigvals) * 1/(2*math.pi)
print(eigen_frequencies_coupled)

#uncoupled eigen vectors and eigen values 
mass_matrix_uncoupled = np.array([[m,0,0],[0,I_theta,0],[0,0,I_beta]])
eigvals_un, eigvecs_un = eig(stiffness_matrix, mass_matrix_uncoupled)

#uncoupled eigen frequencies
eigen_frequencies_uncoupled =  np.sqrt(eigvals_un) * 1/(2*math.pi)
print(eigen_frequencies_uncoupled)