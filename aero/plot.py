#!/usr/bin/python3
from matplotlib import pyplot as plt
import numpy as np
import sys
prefix = '~/ae-wpi/ae5224/hw3/'
if len(sys.argv) > 1:
  prefix = sys.argv[1]


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
poshat = data[:, 31:34]
velhat = data[:, 34:37]
quathat = data[:, 37:41]
omegahat = data[:, 41:44]

(fig, ax) = plt.subplots()
ax.set_prop_cycle('color', ['blue', 'green', 'red'])
plt.plot(ts, pos)
plt.plot(ts, trajpos, linestyle='dashed')
plt.plot(ts, poshat, linestyle='dotted')
plt.legend(['X', 'Y', 'Z', 'X_d', 'Y_d', 'Z_d', 'X_hat', 'Y_hat', 'Z_hat'])
plt.title("Position")
plt.savefig(prefix + '-pos.eps')

(fig, ax) = plt.subplots()
ax.set_prop_cycle('color', ['blue', 'green', 'red'])
plt.plot(ts, vel)
plt.plot(ts, trajvel, linestyle='dashed')
plt.plot(ts, velhat, linestyle='dotted')
plt.legend(['X', 'Y', 'Z', 'X_d', 'Y_d', 'Z_d', 'X_hat', 'Y_hat', 'Z_hat'])
plt.title("Velocity")
plt.savefig(prefix + '-vel.eps')

plt.figure()
plt.plot(ts, acc)
plt.legend(['X', 'Y', 'Z'])
plt.title('Acceleration')
plt.savefig(prefix + '-acceleration.eps')

(fig, ax) = plt.subplots()
ax.set_prop_cycle('color', ['blue', 'green', 'red', 'cyan'])
plt.plot(ts, quat)
plt.plot(ts, qds, linestyle='dashed')
plt.plot(ts, quathat, linestyle='dotted')
plt.legend(['X', 'Y', 'Z', 'W', 'X_d', 'Y_d', 'Z_d', 'W_d', 'X_hat', 'Y_hat', 'Z_hat', 'W_hat'])
plt.title('Quaternion')
plt.savefig(prefix + '-quaternion.eps')

plt.figure()
ax = plt.subplot(211)
ax.set_prop_cycle('color', ['blue', 'green', 'red'])
plt.plot(ts, omega)
if True:
  plt.plot(ts, omegahat, linestyle='dotted')
  plt.legend(['X', 'Y', 'Z', 'X_hat', 'Y_hat', 'Z_hat'])
else:
  plt.legend(['X', 'Y', 'Z'])#, 'X_hat', 'Y_hat', 'Z_hat'])
plt.title('Angular Velocities')
plt.subplot(212)
plt.plot(ts, moments)
plt.legend(['X', 'Y', 'Z'])
plt.title('Moments')
plt.savefig(prefix + '-moments.eps')

plt.figure()
plt.plot(ts, Fbz)
plt.title('Thrust')
plt.savefig(prefix + '-thrust.eps')

plt.show()
