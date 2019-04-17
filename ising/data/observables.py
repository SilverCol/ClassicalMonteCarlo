from matplotlib import pyplot as plt
import numpy as np
from scipy import optimize as opt

L = 256
observables = np.fromfile('heating' + str(L*L) + '_1.000000_0.000000.bin')
observables = np.reshape(observables, (5, int(len(observables) / 5)), 'F')

betas = observables[0]
mags = observables[1]
magsS = observables[2]
chis = betas * (magsS - np.square(mags))
energs = observables[3]
energsS = observables[4]
specs = np.square(betas) * (energsS - np.square(energs))

# TODO this function isnt finished
def magnet(x, betac):
    y = 1 - np.power(np.sinh(0.8813735870195429 * x / betac), -4)
    return np.power(np.abs(y), .125) * (1 + np.sign(y))/2


popt, pcov = opt.curve_fit(magnet, betas, mags, .4)
betaC = popt[0]

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots(nrows=3, sharex=True)
ax[0].grid()
ax[0].axvline(betaC, color='red')
ax[0].text(.02, .8, '$\\beta_c = %.2f$' % betaC, transform=ax[0].transAxes)
ax[1].grid()
ax[1].axvline(betaC, color='red')
ax[2].grid()
ax[2].axvline(betaC, color='red')

ax[0].plot(betas, mags, '.')
ax[0].plot(betas, magnet(betas, popt[0]), '-')
ax[0].set_ylabel('$M$')
ax[0].set_title('Mreža $%d \\times %d$' % (L, L))
ax[1].plot(betas, chis, '-')
ax[1].set_ylabel('$\\chi$')
ax[2].plot(betas, specs, '-')
ax[2].set_ylabel('$C_V$')
ax[2].set_xlabel('$\\beta$')

plt.show()
