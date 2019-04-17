from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.animation as animation

frames = 1000
L = 127
title = 'Segrevanje Isingove mreže ($%d \\times %d$)' % (L, L)

lattice = np.fromfile("lattice" + str(L*L) + "_-1.000000_0.000000.bin", dtype=np.dtype('b'))
betas = np.fromfile("betas" + str(L*L) + "_-1.000000_0.000000.bin")
temps = 1/betas;
lattice = np.reshape(lattice, (frames, L, L))
lattice -= 48

fig = plt.figure()
ax = fig.subplots(subplot_kw=dict(aspect='equal', autoscale_on=False, xlim=(0, L - 1), ylim=(0, L - 1)))
ax.axis('off')
ax.set_title(title)

image = ax.imshow(lattice[500], mpl.colors.ListedColormap(['#ff7f0e', '#1f77b4']))
text = ax.text(.0, -.05, '', transform=ax.transAxes)


def init():
    image.set_data(lattice[0])
    text.set_text('T %.2f' % temps[0])
    return image, text


def animate(i):
    image.set_array(lattice[4*i])
    text.set_text('T %.2f' % temps[int(4*i/10)])
    return image, text


ani = animation.FuncAnimation(fig, animate, frames=250, interval=40, blit=True, init_func=init)
ani.save('../RENAME.gif', writer='imagemagick', fps=25)
plt.show()
