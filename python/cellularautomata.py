"""

Things I could add:
    vector fields for adhesion (?)
    re-create the anderson05 exp w/o adhesion (only chemotaxis)
    Change chemo-dispersion to? O2?

"""

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
from cancer import *
from multiprocessing import Pool

title = sys.argv[-1]
#title = "twotumors"
images_folder = "../images/"


def plotit(frame=-1):
    c_array = np.zeros((N, N))
    for cell in tumor.cancer_cells():
        c_array[cell.x, cell.y] = 0.5 if cell.dead else 1

    plt.subplot(211)
    plt.imshow(c_array)
    plt.subplot(212)

    phi_rounded = phi.copy()
    phi_range = np.max(phi_rounded)
    for i in range(phi_rounded.shape[0]):
        for j in range(phi_rounded.shape[1]):
            if phi_rounded[i, j] > phi_range*0.75:
                phi_rounded[i, j] = 1.0
            elif phi_rounded[i, j] > phi_range * 0.5:
                phi_rounded[i, j] = 0.75
            elif phi_rounded[i, j] > phi_range*0.25:
                phi_rounded[i, j] = 0.5
            elif phi_rounded[i, j] > phi_range*0.125:
                phi_rounded[i, j] = 0.25

    plt.imshow(phi_rounded)
    plt.colorbar()

    print np.linalg.norm(c_array-phi)

    if frame == -1:
        plt.show()
    else:
        plt.savefig(images_folder + "frame_" + title + "_" + str(t) + ".png")
        plt.clf()


N = 100     # diameter of the grid.

phi = np.zeros((N, N))  # The chemo-attractant field
o2 = np.ones((N, N))

tumor = Tumor()  # Collection of cancer cells.

cancer_cells = [CancerCell(N / 2-25, N / 2-25, tumor), CancerCell(N / 2-25, N / 2 - 1-25, tumor),
                CancerCell(N / 2 - 1-25, N / 2 - 1-25, tumor), CancerCell(N / 2 - 1-25, N / 2-25, tumor)]

cancer_cells += [CancerCell(N / 2+25, N / 2+25, tumor), CancerCell(N / 2+25, N / 2 - 1+25, tumor),
                 CancerCell(N / 2 - 1+25, N / 2 - 1+25, tumor), CancerCell(N / 2 - 1+25, N / 2+25, tumor)]

for c in cancer_cells:
    phi[c.x, c.y] = 1.0


def dispersion(arg):
    print "Starting",arg
    if arg == "ca":
        for _ in range(int(1 / float(dt))):
            phi[1:-1, 1:-1] += Dp * lap(phi) * (dt / dx ** 2)
            for c in tumor.cancer_cells():
                phi[c.x, c.y] = 1.0
        return phi
    elif arg == "o2":
        for _ in range(int(1 /float(dt))):
            o2[1:-1, 1:-1] += Dc * lap(o2) * (dt / dx ** 2)
            for c in tumor.cancer_cells():
                o2[c.x, c.y] = np.max(o2[c.x, c.y] - c.consumption() * dt, 0)
        return o2

p = Pool()

for t in range(50):
    st = time.time()
    print "Iteration:", t

    print np.linalg.norm(o2)
    """ Dispersion of the chemo-attractant. """
    o2, phi = p.map(dispersion, ["o2", "ca"])
    print np.linalg.norm(o2)

    print "Updated fields"

    for c in tumor.cancer_cells():
        c.move(phi)
        phi[c.x, c.y] = 1.0

    for c in tumor.cancer_cells():
        c.life_cycle(o2)

    plotit(t)

    print "Tumor is", tumor.size(), "cells big"

    print "Iteration took", (time.time()-st), "seconds"

plotit()
