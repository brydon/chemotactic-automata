"""

Every outer loop represents one cell cycle
inner loop represents dispersion of the chemo-attractant.

Every cell cycle the cell should choose whether to reproduce or to move. (Right?)

The GWN get's included by allowing the cell to move counter to the chemo-attractant.
It randomly picks a direction to move, chemo-attractant influenes but doesn't determine. (w8d random)


"""

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
from cancer import *

title = sys.argv[-1]
title = "twotumors"


def plotit(frame=-1):
    C = np.zeros((N, N))
    for c in tumor.cancer_cells():
        C[c.x, c.y] = 1

    plt.subplot(211)
    plt.imshow(C)
    plt.subplot(212)

    phi_rounded = phi.copy()
    for i in range(phi_rounded.shape[0]):
        for j in range(phi_rounded.shape[1]):
            if phi_rounded[i,j] > 10**(-3):
                phi_rounded[i, j] = 0.25
            if phi_rounded[i,j] > 10**(-2):
                phi_rounded[i, j] = 0.5
            if phi_rounded[i, j] > 10 ** (-1):
                phi_rounded[i, j] = 0.75
            if phi_rounded[i , j] > 0.5:
                phi_rounded[i, j] = 1.0

    plt.imshow(phi_rounded)

    print np.linalg.norm(C-phi)

    if frame == -1:
        plt.show()
    else:
        plt.savefig("frame_"+title+"_"+str(t)+".png")
        plt.clf()


N = 100     # diameter of the grid.
dx = 1e-2   # dx and dt are for the numerical solutions of the PDE for dispersion of chemoattractant
dt = 1e-5   #

D = 1       # diffusion coeff of the cancer-cells
Dp = 0.01   # diffusion coeff of the chemo-attractant
mu = 1      # uptake parameter of the chemo-attractant by the individual cancer-cells

phi = np.zeros((N, N))  # The chemo-attractant field

tumor = Tumor()  # Collection of cancer cells.

cancer_cells = [CancerCell(N / 2-25, N / 2-25, tumor), CancerCell(N / 2-25, N / 2 - 1-25, tumor),
                CancerCell(N / 2 - 1-25, N / 2 - 1-25, tumor), CancerCell(N / 2 - 1-25, N / 2-25, tumor)]

cancer_cells += [CancerCell(N / 2+25, N / 2+25, tumor), CancerCell(N / 2+25, N / 2 - 1+25, tumor),
                 CancerCell(N / 2 - 1+25, N / 2 - 1+25, tumor), CancerCell(N / 2 - 1+25, N / 2+25, tumor)]

for c in cancer_cells:
    phi[c.x, c.y] = 1.0

for t in range(50):
    st = time.time()
    print "Iteration:", t

    """ Dispersion of the chemo-attractant. """
    for _ in range(int(1/float(dt))):
        phi[1:-1, 1:-1] += Dp * lap(phi) * (dt / dx ** 2)
        for c in tumor.cancer_cells():
            phi[c.x, c.y] = 1.0

    print "Updated phi"

    for c in tumor.cancer_cells():
        c.move(phi)
        phi[c.x, c.y] = 1.0

    for c in tumor.cancer_cells():
        c.life_cycle()

    plotit(t)

    print "Iteration took", (time.time()-st), "seconds"

plotit()