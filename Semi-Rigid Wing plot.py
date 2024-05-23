import numpy as np
import matplotlib.pyplot as plt

L = 16 * 0.63
a = 0.1
b = 0.5
K_gamma = 2818.8
K_theta = 37.3
y_cp = 0.8*L
e_cp = (0.5+a)*b
S = 2*b*1
CL_alpha = 2*np.pi
sweep_angle = np.linspace(0,60,100)
e_t = e_cp/y_cp
K_gamma_t = K_gamma /(S*CL_alpha*np.cos(np.pi*sweep_angle/180)*y_cp)
K_theta_t = K_theta /(S*CL_alpha*np.cos(np.pi*sweep_angle/180)*y_cp)

qdiv = K_gamma_t * K_theta_t/(K_gamma_t*e_t-K_theta_t*np.tan(np.pi*sweep_angle/180))

plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
plt.plot(sweep_angle, qdiv, '.', color='blue',label='divergence pressure vs Sweep angle')
plt.xlabel('Sweep angle in Degrees')
plt.ylabel('divergence pressure')
plt.title('Plot for Divergence Pressure vs Sweep Angle')
plt.legend()  # Show legend
plt.grid(True)  # Show grid
plt.show()

isoclinic_angle = K_gamma*e_t/K_theta
print("The Isoclinic angle is " +str(isoclinic_angle))