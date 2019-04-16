from matplotlib import pyplot as plt
import numpy as np

M = 100
N = 1024

data = np.fromfile("heating" + str(N) + "_1.000000_0.000000.bin")
data = np.reshape(data, (M, N + 1))
data = np.hsplit(data, [1, N + 1])
betas = data[0].reshape(M)
print(betas)
correlations = data[1]

fig = plt.figure()
ax = fig.subplots(subplot_kw=dict(aspect='auto'))

print(correlations.shape)
print(np.arange(N).shape)
print((1/betas).shape)

image = ax.pcolormesh(np.arange(N), betas, correlations, cmap='RdYlGn')
bar = fig.colorbar(image)

plt.show()