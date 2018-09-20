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
title = "ark_anch"
images_folder = "../images/"


def plotit(frame=-1):
    c_array = np.zeros((N, N))
    for cell in tumor.cancer_cells():
        c_array[cell.x, cell.y] = 0.5 if cell.dead else 1

    plt.subplot(221)
    plt.imshow(c_array)
    plt.subplot(222)
    plt.imshow(o2, cmap="gray")
    plt.colorbar()
    plt.subplot(223)

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


N = 400     # diameter of the grid.

phi = np.zeros((N, N))  # The chemo-attractant field
o2 = np.ones((N, N))

tumor = Tumor()  # Collection of cancer cells.

##############################
#### Starting Connditions ####
##############################

# Single, central solid tumour site
starting_positions = [(N/2, N/2), (N/2 - 1, N/2), (N/2, N/2 - 1), (N/2 - 1, N/2 - 1)]

# Two solid tumour sites
starting_positions = [(N/4, N/4), (N/4 - 1, N/4), (N/4, N/4 - 1), (N/4 - 1, N/4 - 1), (3*N/4, 3*N/4),
                      (3*N/4 - 1, 3*N/4), (3*N/4, 3*N/4 - 1), (3*N/4 - 1, 3*N/4 - 1)]

# Random Tumor Sites, non-solid
starting_positions = [(i, j) for i in range(N) for j in range(N) if np.random.random() > 0.9]

print starting_positions

cancer_cells = [CancerCell(x, y, tumor) for (x, y) in starting_positions]

for c in cancer_cells:
    phi[c.x, c.y] = 1.0


def dispersion(args):
    s = args[0]
    arr = args[1].copy()
    tum = args[2]
    print "\tStarting", s, " -- ", tum.size()
    if s == "ca":
        _const = Dp * (dt / dx ** 2)
        for c in tum.cancer_cells():
            arr[c.x, c.y] = 1.0
        for _ in range(int(1 / float(dt))):
            arr[1:-1, 1:-1] += _const * (arr[2:, 1:-1] + arr[1:-1, 2:] + arr[:-2, 1:-1] + arr[1:-1, :-2] - 4 * arr[1:-1, 1:-1])
        return arr
    elif s == "o2":
        _const = Dc * (dt / dx ** 2)
        for c in tum.cancer_cells():
            arr[c.x, c.y] = np.max(arr[c.x, c.y] - c.consumption(), 0)
        for _ in range(int(1 / float(dt))):
            arr[1:-1, 1:-1] += _const * (arr[2:, 1:-1] + arr[1:-1, 2:] + arr[:-2, 1:-1] + arr[1:-1, :-2] - 4 * arr[1:-1, 1:-1])
        return arr
    print "\t<<Field " + s + " returning>>"


p = Pool()

for t in range(50):
    st = time.time()
    print "Iteration:", t

    print "\t||o2|| =", np.linalg.norm(o2), " ||phi||", np.linalg.norm(phi)
    print "\tHi/Lo =", np.max(o2), np.min(o2)
    """ Dispersion of the chemo-attractant. """
    o2, phi = p.map(dispersion, [("o2", o2, tumor), ("ca", phi, tumor)])
    print "\t||o2|| =", np.linalg.norm(o2), " ||phi|| =", np.linalg.norm(phi)
    print "\tHi/Lo =", np.max(o2), np.min(o2)

    print "\tUpdated fields"

    for c in tumor.cancer_cells():
        c.move(phi)
        phi[c.x, c.y] = 1.0

    print "\tDone moving"

    for c in tumor.cancer_cells():
        c.life_cycle(o2)

    print "\tDone life-cycling"

    plotit(t)

    print "Tumor is", tumor.size(), "cells big"

    print "Iteration took", (time.time()-st), "seconds"

plotit()
