import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


# k = np.linspace(0, 1000, 1000)
# F = 1 - (0.165 * k**2/(k**2 + 0.0455**2)) - (0.335 * k**2/(k**2 + 0.3**2))
# G = -((0.165*0.0455*k)/(k**2 + 0.0455**2) + (0.335*0.3*k/(k**2 + 0.3**2)))

# fig, ax = plt.subplots(figsize=(8, 6))
# line, = ax.plot([], [], '.', color='orange')

# def init():
#     ax.set_xlim(0.5, 1)
#     ax.set_ylim(-0.2, 0)
#     return line,

# def update(frame):
#     xdata = F[:frame]
#     ydata = G[:frame]
#     line.set_data(xdata, ydata)
#     return line,

# ani = animation.FuncAnimation(fig, update, frames=len(F), init_func=init, blit=True)

# plt.xlabel('F(k)')
# plt.ylabel('G(k)')
# plt.title('Plot of Theodoersen function')
# plt.grid(True)



# plt.show()


def Theodoersen_function(k):
    F = 1 - (0.165 * k**2/(k**2 + 0.0455**2)) - (0.335 * k**2/(k**2 + 0.3**2))
    G = -((0.165*0.0455*k)/(k**2 + 0.0455**2) + (0.335*0.3*k/(k**2 + 0.3**2)))
    C = F + 1j * G
    return C