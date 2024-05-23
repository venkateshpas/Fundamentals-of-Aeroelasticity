import numpy as np
import sympy as sp
from Theodoersenfunction import Theodoersen_function

rho, a, b, V, theta, h, C, k, t = sp.symbols('rho a b V theta h C k t',real = True)


h_t = h * sp.exp(1j * V/b*k*t)
theta_t = theta * sp.exp(1j * V/b*k*t)
h1 = sp.diff(h_t,t)
h2 = sp.diff(h1,t)
theta1 = sp.diff(theta_t,t)
theta2 = sp.diff(theta1,t)
Lift_h = (-rho*b**2*(np.pi*h2)-2*np.pi*rho*V*b*C*(h1))/(1/2*rho*V**2*h_t)
Lift_h = sp.simplify(Lift_h)
Lift_theta = (-rho*b**2*(V*np.pi*theta1-np.pi*b*a*theta2)-2*np.pi*rho*V*b*C*(V*theta_t+b*(1/2-a)*theta1))/(1/2*rho*V**2*theta_t)
Lift_theta = sp.simplify(Lift_theta)
Moment_h = (-rho*b**2*(-a*np.pi*b*h2)+2*np.pi*rho*V*b**2*(a+1/2)*C*(h1))/(1/2*rho*V**2*h_t)
Moment_h = sp.simplify(Moment_h)
Moment_theta = (-rho*b**2*(np.pi*(1/2-a)*V*b*theta1+np.pi*b**2*(1/8+a**2)*theta2)+2*np.pi*rho*V*b**2*(a+1/2)*C*(V*theta_t+b*(1/2-a)*theta1))/(1/2*rho*V**2*theta_t)
Moment_theta = sp.simplify(Moment_theta)
Theodoersen_matrix = np.array([[Lift_h, Lift_theta],[Moment_h, Moment_theta]])