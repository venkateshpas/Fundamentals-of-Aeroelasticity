import numpy as np
import cmath
from scipy.linalg import eig
import matplotlib.pyplot as plt
from sympy import symbols, Matrix, symarray, det, roots, symarray
from sympy.polys.polyfuncs import symmetrize

# Mass of air foil and flap
m_airfoil = 1.567
m_flap = 0
m = m_airfoil + m_flap

#Properties of Typical Section
L_wing = 16
b = 0.127 #half chord length
a = -0.1 #from the center to the elastic axis
c = 0.5 #from the center to the hinge axis
X_theta = -0.5 #between elastic axis and the center of gravity of air foil
X_beta = 0.1 #between hinge axis and the center of gravity of flap
S = 2*b* 1 #change later back to 2b*1

alpha_0 = 5*np.pi/180
CM_AC = 0

#Defining the stiffness matrix
L = 0.63 * L_wing
K_h = 2818.8
K_theta = 37.3
structural_matrix = np.array([[K_h, 0],[0, K_theta]])

#Defining the aerodynamic stiffness matrix
rho = 1.225
V = 23
q = 0.5 * rho * V**2
CL_alpha  = 2*np.pi
CL_beta = CL_alpha/20
CM_AC_beta = -0.1
K_aerodynamic_stiffness = np.array([[0, -q*S*CL_alpha], [0, q*S*CL_alpha*(0.5+a)*b]])
K_aerodynamic_stiffness_tilda = np.array([[0, -S*CL_alpha], [0, S*CL_alpha*(0.5+a)*b]])

#Defining the aerodynamic damping matrix
CL_alphadot = 0.1
CM_alphadot = 0.01
K_aerodynamic_damping = np.array([[-q*S*CL_alpha/V, -q*S*CL_alphadot*(b/V)],[q*S*CL_alphadot*(0.5+a)*b/V, q*S*CL_alphadot*(b/V)*(0.5+a)*b+q*S*CM_alphadot*b*b/V]])
K_aerodynamic_damping_tilda = np.array([[-S*CL_alpha, -S*CL_alphadot*(b)],[S*CL_alpha*(0.5+a)*b, S*CL_alphadot*(b)*(0.5+a)*b+S*CM_alphadot*b*b]])
print(K_aerodynamic_damping)

#moment of intertias
I_cg_theta = 1
I_cg_beta = 0.01
W = 30 #probably the omega
S_theta = 0.08587 #m_airfoil * X_theta * b + m_flap * (c-a+X_beta) * b
S_beta = m_flap * X_beta * b
I_theta = 0.01347 #I_cg_theta + I_cg_beta + m_airfoil * (X_theta * b) ** 2 + m_flap * (c-a+X_beta)**2 * b**2 
I_beta = I_cg_beta + m_flap * (X_beta * b) ** 2
mass_matrix = np.array([[m,S_theta,S_beta],[S_theta,I_theta,(c-a)*b*S_beta+I_beta],[S_beta,(c-a)*b*S_beta+I_beta,I_beta]])

q_divergence = K_theta/(S*CL_alpha*(0.5 + a)*b)

q = np.linspace(0.5,500,200)
V = np.sqrt(2*q/rho)
p1=[]
p2=[]
p3=[]
p4=[]
#Flutter omegas:
ps = symbols('ps')
p_save = []
for i in range(len(q)):
    Kae = structural_matrix[:2, :2] - q[i] * K_aerodynamic_stiffness_tilda 
    A = ps**2 * Matrix(mass_matrix[:2, :2]) + Matrix(Kae) - ps*q[i]/V[i]*Matrix(K_aerodynamic_damping_tilda)
    # Alternatively, if Ca is used in the equation
    # A = ps**2 * Matrix(Ms[:2, :2]) - ps * q[i] / V[i] * Matrix(Ca) + Matrix(Kae)
    DA = det(A)
    DA = symmetrize(DA,ps)
    characteristic_equation = np.poly(DA)
    p = roots(characteristic_equation[1])
    p_save.extend(p.keys())


#print("Roots:", np.real(p_save))
real_parts = []
imaginary_parts = []

for root in p_save:
    real_parts.append(root.as_real_imag()[0])  # Extract real part
    imaginary_parts.append(root.as_real_imag()[1])  # Extract imaginary part

q = [item for item in q for _ in range(4)]
print(q)

plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(q, imaginary_parts, '.', color='blue',label='effectiveness')
plt.xlabel('Dynamic pressure')
plt.ylabel('Frequency')
plt.title('Flutter Plot: Pressure vs Frequency')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()

plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(q, real_parts, '.', color='orange',label='effectiveness')
plt.xlabel('Dynamic pressure')
plt.ylabel('Damping')
plt.title('Flutter Plot: Pressure vs Damping')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()

plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(real_parts, imaginary_parts, '.', color='red',label='effectiveness')
plt.xlabel('Damping')
plt.ylabel('Frequency')
plt.title('Flutter Plot: Damping vs Frequency')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()

# for i in q:
#     a4 = m*I_theta - S_theta**2
#     a2 = K_h*I_theta + m*K_theta - m*i*S*CL_alpha*(0.5+a)*b-S_theta*i*S*CL_alpha
#     a0 = K_h*(K_theta-i*S*CL_alpha*(0.5+a)*b)
#     p_squared = (-a2+cmath.sqrt(a2**2-4*a4*a0))/(2*a4)
#     p_squared_2 = (-a2-cmath.sqrt(a2**2-4*a4*a0))/(2*a4)
#     p1.append(cmath.sqrt(p_squared))
#     p2.append(-cmath.sqrt(p_squared))
#     p3.append(cmath.sqrt(p_squared_2))
#     p4.append(-cmath.sqrt(p_squared_2))


# plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
# plt.plot(q, np.imag(p1), '.', color='blue',label='effectiveness')
# plt.plot(q, np.imag(p3), '.', color='orange',label='effectiveness')
# plt.xlabel('Dynamic pressure')
# plt.ylabel('Frequency')
# plt.title('Plot for control reversal')
# plt.legend()  # Show legend
# plt.grid(True)  # Show grid
# plt.show()

# plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
# plt.plot(q, np.real(p1), color='blue',label='effectiveness')
# plt.plot(q, np.real(p2), color='orange',label='effectiveness')
# plt.xlabel('Dynamic pressure')
# plt.ylabel('Damping')
# plt.title('Plot for control reversal')
# plt.legend()  # Show legend
# plt.grid(True)  # Show grid
# plt.show()

# plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
# plt.plot(np.real(p1), np.imag(p1), '.', color='blue',label='effectiveness')
# plt.plot(np.real(p2), np.imag(p2), '.', color='blue',label='effectiveness')
# plt.plot(np.real(p3), np.imag(p3), '.',color='blue',label='effectiveness')
# plt.plot(np.real(p4), np.imag(p4), '.', color='blue',label='effectiveness')
# plt.xlabel('divergence pressure')
# plt.ylabel('Frequency')
# plt.title('Plot for control reversal')
# plt.legend()  # Show legend
# plt.grid(True)  # Show grid
# plt.show()
