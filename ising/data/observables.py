from matplotlib import pyplot as plt
import numpy as np

observables = np.fromfile("observables.bin")
observables = np.reshape(observables, (5, int(len(observables) / 5)), 'F')

betas = observables[0]
mags = observables[1]
magsS = observables[2]
chis = betas * (magsS - np.square(mags))
energs = observables[3]
energsS = observables[4]
specs = np.square(betas) * (energsS - np.square(energs))
yes = specs

plt.rcParams.update({'font.size': 20})
fig = plt.figure()
yMax = max(yes)
yMin = min(yes)
xMin = min(betas)
xMax = max(betas)
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()

line1, = ax.plot(betas, yes, '-')
ax.set_ylabel('$M$')
ax.set_xlabel('$\\beta$')

plt.show()
