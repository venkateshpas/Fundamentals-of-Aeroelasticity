import numpy as np
from scipy.linalg import eig
import matplotlib.pyplot as plt

# Mass of air foil and flap
m_airfoil = 0.75 * 16
m_flap = 0
m = m_airfoil + m_flap

#Properties of Typical Section
L_wing = 16
b = 0.5 #half chord length
a = 0.2 #from the center to the elastic axis
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
K_aerodynamic_stiffness = np.array([[0, -q*S*CL_alpha], [0, q*S*CL_alpha*(0.5+a)*b]])

#Divergence using normal formula: Finding Divergence pressure

q_divergence = K_theta/(S*CL_alpha*(0.5 + a)*b)

#Divergence using eigen value analysis of Kh-Ka matrix
K_effective = K_structural_matrix - K_aerodynamic_stiffness

eigvals, eigvecs = eig(K_effective)

K_gamma = 2e4 *3/ L**3
K_theta = 1e4/L
y_cp = 0.8*L
e_cp = (0.5+a)*b
isoclinic_angle = np.arctan(K_gamma * e_cp / (K_theta * y_cp))#*180/np.pi
print(isoclinic_angle)


#q divergence change in semirigid wing with sweep angle
sweep_angle = isoclinic_angle#*np.pi/180
K_gamma_t = K_gamma/(S*CL_alpha*y_cp*np.cos(sweep_angle))
K_theta_t = K_theta/(S*CL_alpha*y_cp*np.cos(sweep_angle))
e_t = e_cp/y_cp
q_divergence_sweep = K_gamma_t*K_theta_t/(K_gamma_t*e_t-K_theta_t*np.tan(sweep_angle))
print(q_divergence_sweep)

# q_plot = []
# lamb = 0
# q_div = 0
# lamb_plot = []
# while(q_div < 10e3):
#     K_gamma1 = K_gamma/(y_cp* np.cos(lamb)*S*CL_alpha)
#     print(np.cos(lamb))
#     K_theta1 = K_theta/(y_cp* np.cos(lamb)*S*CL_alpha)
#     e1 = e_cp/y_cp
#     q_div = K_gamma1 * K_theta1 / (+K_gamma1 * e1 - K_theta1 * np.tan(lamb))
#     print(q_div)
#     q_plot.append(q_div)
#     lamb_plot.append(lamb)
#     lamb = lamb + np.pi/180 * 5 
# print(q_plot)   
# fig, (ax1, ax2) = plt.subplots(1, 2)
# ax1.plot(lamb_plot,q_plot, label='q_divergence')
# ax1.set_xlabel('Iteration number')
# ax1.set_ylabel('heave displacement')
# ax1.set_title('Monolithic vs heave displacement')
# ax1.legend()
# plt.tight_layout()
# plt.show()