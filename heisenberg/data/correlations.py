from matplotlib import pyplot as plt
import numpy as np

M = 200
N = 128

data = np.fromfile("cooling" + str(N) + "_1.000000_0.000000_20.bin")
data = np.reshape(data, (M, 3*N + 1))
data = np.hsplit(data, [1, 3*N + 1])
betas = data[0].reshape(M)
data = data[1]
data = np.reshape(data, (3*M*N,))
data = np.reshape(data, (3, M*N), 'F')
data = np.reshape(data, (3, M, N))
correlationX = data[0]
correlationY = data[1]
correlationZ = data[2]

plt.rcParams.update({'font.size': 20})
fig = plt.figure()
ax = fig.subplots(ncols=4, sharey=True, subplot_kw=dict(aspect='auto'))

x = ax[0].pcolormesh(np.arange(N), betas, correlationX, cmap='RdYlGn')
y = ax[1].pcolormesh(np.arange(N), betas, correlationY, cmap='RdYlGn')
z = ax[2].pcolormesh(np.arange(N), betas, correlationZ, cmap='RdYlGn')
s = ax[3].pcolormesh(np.arange(N), betas, correlationX + correlationY + correlationZ, cmap='RdYlGn')

fig.subplots_adjust(right=0.88, left=0.1)
cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.8])
fig.colorbar(s, cax=cbar_ax)

ax[0].set_ylabel('$\\beta$')
ax[0].set_xlabel('$r$')
ax[0].set_title('$\\langle \\sigma_0^x \\sigma_r^x \\rangle$')
ax[1].set_xlabel('$r$')
ax[1].set_title('$\\langle \\sigma_0^y \\sigma_r^y \\rangle$')
ax[2].set_xlabel('$r$')
ax[2].set_title('$\\langle \\sigma_0^z \\sigma_r^z \\rangle$')
ax[3].set_xlabel('$r$')
ax[3].set_title('$\\langle \\vec{\\sigma_0} \\cdot \\vec{\\sigma_r} \\rangle$')

plt.show()
