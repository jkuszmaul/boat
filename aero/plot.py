#!/usr/bin/python3
from matplotlib import pyplot as plt
import numpy as np

data = np.genfromtxt("sim.csv")

ts = data[:, 0]
pos = data[:, 1:4]
vel = data[:, 4:7]
quat = data[:, 7:11]
omega = data[:, 11:14]
acc = data[:, 14:17]
trajpos = data[:, 17:20]
trajvel = data[:, 20:23]
Fbz = data[:, 23]
moments = data[:, 24:27]
qds = data[:, 27:31]

plt.figure()
plt.plot(ts, pos)
plt.plot(ts, trajpos)
plt.legend(['X', 'Y', 'Z', 'X_d', 'Y_d', 'Z_d'])
plt.title("Position")

plt.figure()
plt.plot(ts, vel)
plt.plot(ts, trajvel)
plt.legend(['X', 'Y', 'Z', 'X_d', 'Y_d', 'Z_d'])
plt.title("Velocity")

plt.figure()
plt.plot(ts, acc)
plt.legend(['X', 'Y', 'Z'])
plt.title('Acceleration')

(fig, ax) = plt.subplots()
ax.set_prop_cycle('color', ['blue', 'green', 'red', 'cyan'])
plt.plot(ts, quat)
plt.plot(ts, qds, linestyle='dashed')
plt.legend(['X', 'Y', 'Z', 'W', 'X_d', 'Y_d', 'Z_d', 'W_d'])
plt.title('Quaternion')

plt.figure()
plt.subplot(211)
plt.plot(ts, omega)
plt.legend(['X', 'Y', 'Z'])
plt.title('Angular Velocities')
plt.subplot(212)
plt.plot(ts, moments)
plt.legend(['X', 'Y', 'Z'])
plt.title('Moments')

plt.figure()
plt.plot(ts, Fbz)
plt.title('Thrust')

plt.show()
