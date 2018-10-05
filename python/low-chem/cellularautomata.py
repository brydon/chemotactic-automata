#!/software/.admin/bins/bin/python2.7

"""

Things I could add:
    vector fields for adhesion (?)
    re-create the anderson05 exp w/o adhesion (only chemotaxis)
    Change chemo-dispersion to? O2?

"""

import numpy as np
import time
#import matplotlib.pyplot as plt
import sys
from cancer import *
from multiprocessing import Pool
import pickle as pik
import os

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

phi = np.zeros((N, N))  # The chemo-attractant field

tumor = Tumor()  # Collection of cancer cells.

##############################
#### Starting Connditions ####
##############################

# Single, central solid tumour site
#starting_positions = [(N/2, N/2), (N/2 - 1, N/2), (N/2, N/2 - 1), (N/2 - 1, N/2 - 1)]

# Two solid tumour sites
starting_positions = [(N/4, N/4), (N/4 - 1, N/4), (N/4, N/4 - 1), (N/4 - 1, N/4 - 1), (3*N/4, 3*N/4),
                      (3*N/4 - 1, 3*N/4), (3*N/4, 3*N/4 - 1), (3*N/4 - 1, 3*N/4 - 1)]

# Random Tumor Sites, non-solid
#starting_positions = [(i, j) for i in range(N) for j in range(N) if np.random.random() > 0.9]

#starting_positions = [(i, j) for i in range(10,40) for j in range(10, 40)]

print starting_positions

cancer_cells = [CancerCell(x, y, tumor) for (x, y) in starting_positions]

for c in cancer_cells:
    phi[c.x, c.y] = 1.0

def scale_down(a, fact=5):
    n = a.shape[0]
    b = np.empty((n/fact, n/fact))
    for i in range(n/fact):
        for j in range(n/fact):
            b[i, j] = np.mean(a[i*fact:(i+1)*fact, j*fact:(j+1)*fact])
    return b

np.set_printoptions(threshold=40*40+1)

num_steps = 100

dat_files = [x for x in os.listdir(os.getcwd()) if x.endswith(".dat")]
dat_files.sort()

cur_done = -1

if len(dat_files) > 0:

    fn = dat_files[-1]
    print "Opening run",fn

    underscore = fn[::-1].find("_")
    dot = fn.find(".")

    cur_done = int(fn[-underscore:dot])

    with open(fn,'rb') as f:
        [tumor, phi] = pik.load(f)

print scale_down(phi, 10)

for t in range(cur_done + 1, num_steps):
    print "Frame",t

    for cell in tumor.cancer_cells():
        phi[cell.x, cell.y] = 1.

    tm = time.time()
    """ Disperse Chemo-attractant """
    Nt = int(1/float(dt))

    _const = Dp * dt / (dx ** 2)

    for x in [Nt/4, Nt/4, Nt/4, Nt-3*(Nt/4)]:
        for _ in range(x):
            phi[1:-1, 1:-1] += _const * (phi[2:, 1:-1] + phi[1:-1, 2:] + phi[:-2, 1:-1] + phi[1:-1, :-2] - 4 * phi[1:-1, 1:-1])
        print "Check-point"

    print "\tdispersed in",time.time() - tm

    tm = time.time()

    """ Cell lifecycle """
    for cell in tumor.cancer_cells():
        cell.move(phi)

    for cell in tumor.cancer_cells():
        cell.life_cycle()

    for cell in tumor.cancer_cells():
        phi[cell.x, cell.y] = 1.

    print "\tlife cycled in", time.time() - tm
    tm = time.time()

    with open('out_dat_frame_' + str(t) + '.dat', 'wb') as f:
        pik.dump([tumor, phi], f)

    print "\tdumped in", time.time() - tm

    break
