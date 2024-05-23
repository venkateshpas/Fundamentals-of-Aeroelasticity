import numpy as np
import sympy as sp
from TheodoersenEquationsinqAXValues import Theodoersen_matrix
from Theodoersenfunction import Theodoersen_function
import matplotlib.pyplot as plt

rho, a, b, V, theta, h, C, k, t = sp.symbols('rho a b V theta h C k t',real = True)

rho = 1.225 #This
alpha_0 = 5*np.pi/180
CM_AC = 0
CL_alpha = 2*np.pi


# Mass of air foil and flap
m_airfoil = 1.567
m_flap = 0
m = m_airfoil + m_flap

#Geometric parameters
b = 0.127 #This
a = -0.5 #This
c = .5
X_theta = -.5
X_beta = .1
S = 2*b*1

#moment of intertias
I_cg_theta = 1
I_cg_beta = 0.01
W = 30 #probably the omega
S_theta = 0.08587 #m_airfoil * X_theta * b + m_flap * (c-a+X_beta) * b
S_beta = m_flap * X_beta * b
I_theta = 0.01347 #I_cg_theta + I_cg_beta + m_airfoil * (X_theta * b) ** 2 + m_flap * (c-a+X_beta)**2 * b**2 
I_beta = I_cg_beta + m_flap * (X_beta * b) ** 2


mass_matrix = np.array([[m,S_theta],[S_theta,I_theta]])

#Stiffness Matrix
K_h = 2818.8
K_theta = 37.7

Structural_stiffness_matrix = np.array([[K_h,0],[0,K_theta]])

V = np.linspace(10,20,20)

sigma = []

for i in range(len(V)):
    k_1 = 3
    k_2 = 3
    p_1 = 1j *k_1 * V[i]/b
    p_2 = (-0.01 + 1j * k_2)  *V[i]/b
    sigma_diff = np.real(p_2)
    det_matrix_solve3 = 1
    while abs(det_matrix_solve3) > 1e-6:
        C = Theodoersen_function(k_1)
        matrix_solve = p_1**2 * mass_matrix + Structural_stiffness_matrix - 0.5 * rho * V[i]**2 * Theodoersen_matrix(rho,a,b,C,k_1)
        matrix_solve = np.array(matrix_solve,dtype = complex)
        det_matrix_solve1 = np.linalg.det(matrix_solve)
        C = Theodoersen_function(k_2)
        matrix_solve = p_2**2 * mass_matrix + Structural_stiffness_matrix - 0.5 * rho * V[i]**2 * Theodoersen_matrix(rho,a,b,C,k_2)
        matrix_solve = np.array(matrix_solve,dtype = complex)
        det_matrix_solve2 = np.linalg.det(matrix_solve)
        p_3 = (p_2*det_matrix_solve1-p_1*det_matrix_solve2)/(det_matrix_solve1-det_matrix_solve2)
        k_3 = np.imag(p_3)*b/V[i]
        C = Theodoersen_function(k_3)
        matrix_solve = p_3**2 * mass_matrix + Structural_stiffness_matrix - 0.5 * rho * V[i]**2 * Theodoersen_matrix(rho,a,b,C,k_3)
        matrix_solve = np.array(matrix_solve,dtype = complex)
        det_matrix_solve3 = np.linalg.det(matrix_solve)
        sigma_2 = np.real(p_3)
        sigma_1 = np.real(p_2)
        sigma_diff = sigma_2-sigma_1
        p_1 = p_2
        k_1 = np.imag(p_1)*b/V[i]
        # del_1 = np.real(p_1)*b/V[i]
        # p_1 = del_1 + 1j * k_1
        p_2 = p_3
        k_2 = np.imag(p_2)*b/V[i]
        # del_2 = np.real(p_2)*b/V[i]
        # p_2 = del_2 + 1j * k_2

    sigma.append(np.real(p_3))
    if np.real(p_3) >= 0:
        break

print(sigma)

plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(V[:len(sigma)], sigma, '-', color='blue', label = 'Velocity with k = 3')
plt.xlabel('Velocity')
plt.ylabel('Damping')
plt.title('Flutter Plot: PK-method')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()

# k = 2
# sigma = 0

# p_1 = sigma + 1j * k #sp.symbols('p')
# C = Theodoersen_function(k)
# matrix_solve = p_1**2 * mass_matrix + Structural_stiffness_matrix - 0.5 * rho * V[0]**2 * Theodoersen_matrix(rho,a,b,C,k)
# matrix_solve = np.array(matrix_solve,dtype = complex)
# det_matrix_solve1 = np.linalg.det(matrix_solve)
# #print(det_matrix_solve1)
# sigma = -0.01
# p_2 = sigma + 1j * k
# matrix_solve = p_2**2 * mass_matrix + Structural_stiffness_matrix - 0.5 * rho * V[0]**2 * Theodoersen_matrix(rho,a,b,C,k)
# matrix_solve = np.array(matrix_solve,dtype = complex)
# det_matrix_solve2 = np.linalg.det(matrix_solve)
# #print(det_matrix_solve2)

# p_3 = (p_2*det_matrix_solve1-p_1*det_matrix_solve2)/(det_matrix_solve1-det_matrix_solve2)
# print(p_3)

# k = np.imag(p_0)
# print(k)
# C = Theodoersen_function(k)
# matrix_solve = p_0**2 * mass_matrix + Structural_stiffness_matrix - 0.5 * rho * V[0]**2 * Theodoersen_matrix(rho,a,b,C,k)
# matrix_solve = np.array(matrix_solve,dtype = complex)
# det_matrix_solve0 = np.linalg.det(matrix_solve)
# print(det_matrix_solve0)
