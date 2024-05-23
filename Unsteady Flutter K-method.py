#Theodoersen Forces and Moments
import numpy as np
import sympy as sp
from scipy.linalg import eig
from numpy.linalg import inv
from Theodoersenfunction import Theodoersen_function
from TheodoersenEquationsinqAXValues import Theodoersen_matrix
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
 
k = np.linspace(5,0.1,100) #This

omega_1 = []
omega_2 = []
g_1 = []
g_2 = []
V_1 = []
V_2 = []
change_k = 0
before_k = 0
for i in range(len(k)):
    C = Theodoersen_function(k[i])
    A = Theodoersen_matrix(rho, a, b, C, k[i])
# Display the results
    B = np.dot(inv(Structural_stiffness_matrix),(0.5*rho*(b/k[i])**2*A+mass_matrix))
    B = np.array(B, dtype=complex)
    eigenvalues = np.linalg.eigvals(B)
    omega1 = 1/np.sqrt(np.real(eigenvalues)[0])
    omega2 = 1/np.sqrt(np.real(eigenvalues)[1])
    omega_1.append(omega1)
    omega_2.append(omega2)
    g1 = omega1**2 * np.imag(eigenvalues)[0]
    g2 = omega2**2 * np.imag(eigenvalues)[1]
    g_1.append(g1)
    g_2.append(g2)
    V1 = omega1*b/k[i]
    V2 = omega2*b/k[i]
    V_1.append(V1)
    V_2.append(V2)
    if g2 >= 0:
        change_k = k[i]
        before_k = k[i-1]
        k_new = np.linspace(before_k,change_k,1000)
        error = 1
        g2prev = 5
        while error>10e-6:
            C = Theodoersen_function(k_new[i])
            A = Theodoersen_matrix(rho, a, b, C, k_new[i])
        # Display the results
            B = np.dot(inv(Structural_stiffness_matrix),(0.5*rho*(b/k_new[i])**2*A+mass_matrix))
            B = np.array(B, dtype=complex)
            eigenvalues = np.linalg.eigvals(B)
            omega1 = 1/np.sqrt(np.real(eigenvalues)[0])
            omega2 = 1/np.sqrt(np.real(eigenvalues)[1])
            g1 = omega1**2 * np.imag(eigenvalues)[0]
            g2 = omega2**2 * np.imag(eigenvalues)[1]
            V1 = omega1*b/k_new[i]
            V2 = omega2*b/k_new[i]
            error = abs(g2-g2prev)
            g2prev = g2
        print("The velocity at which flutter will occur is "+str(V2))
        break

# omega = sorted(omega)
# g = sorted(g)
# V = sorted(V)

plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(V_1, g_1, '-', color='blue', label = 'Velocity 1')
plt.plot(V_2, g_2, '-', color='orange', label = 'velocity 2')
plt.xlabel('Velocity')
plt.ylabel('Artificial Damping g')
plt.title('Flutter Plot: K-method')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()

plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(V_1, omega_1, '-', color='blue', label = "omega_theta")
plt.plot(V_2, omega_2, '-', color='orange', label = "omega_heave")
plt.xlabel('Velocity')
plt.ylabel('Frequencies')
plt.title('Flutter Plot: K-method')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()