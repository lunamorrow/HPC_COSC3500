# Visualise optional trajectory dump from Lennard-Jones model
# To access/obtain trajectory, commented-out sections of 'lj_model_optimised.c'
# will need to be uncommented. Additionally, the number of frames and number of
# particles constants set below in lines 11-12 will need to be adjusted.

import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

num_frames = 999
num_particles = 60
num_values = num_frames*num_particles*2*8

# Read binary file
x_values = np.empty(num_frames*num_particles)
y_values = np.empty(num_frames*num_particles)
with open('trajectory.data', 'rb') as file:
    fileContent = file.read()
    index_x = 0
    index_y = 0
    for i in range(0, num_values, 8):
        float_value = struct.unpack('d', fileContent[i:i+8])[0]
        if (int(i/8)%2 == 0):
            x_values[index_x] = float_value
            index_x += 1
        else:
            y_values[index_y] = float_value
            index_y += 1

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = ax.plot([], [], 'ro')

def init():
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    return ln,

def update(frame):
    print(frame)
    lower = int((frame-1)*num_particles)
    upper = int(frame*num_particles)
    xdata = x_values[lower:upper]
    ydata = y_values[lower:upper]
    print(f"x: {xdata}")
    print(f"y: {ydata}")
    print()
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, frames=np.linspace(1, num_frames, num_frames//10),
                    init_func=init, blit=True)
plt.show()