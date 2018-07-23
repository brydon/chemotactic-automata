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

title = sys.argv[-1]
title = "twotumors"

class Tumor:
    def __init__(self):
        self.cells = {}

    def add(self, cell):
        self.cells[(cell.x, cell.y)] = cell
        cell.tumor = self

    def cancer_at(self, x, y):
        return (x, y) in self.cells.keys()

    def cancer_cells(self):
        return self.cells.values()

    def move(self, old, new):
        self.cells[new] = self.cells.pop(old)


class CancerCell:
    def __init__(self, x, y):
        self.x, self.y = x, y
        self.tumor = None
        self.stayed = False

    def move(self, phi):
        """
        TBD: Should there be a dx/dt factor in each/most of the Cim terms??
        :param phi:
        :return:
        """
        Cim1 = D - mu * 0.25 * (phi[self.x + 1, self.y] - phi[self.x - 1, self.y])
        Cip1 = D + mu * 0.25 * (phi[self.x + 1, self.y] - phi[self.x - 1, self.y])
        Cjm1 = D - mu * 0.25 * (phi[self.x, self.y + 1] - phi[self.x, self.y - 1])
        Cjp1 = D + mu * 0.25 * (phi[self.x, self.y + 1] - phi[self.x, self.y - 1])
        Cs = dx ** 2 / dt - 4 * D + mu * single_lap(phi, self.x, self.y)

        """
        if self.stayed:
            Cs = 0
        """

        w8s = [Cim1, Cip1, Cjm1, Cjp1, Cs]

        directions = (
        (self.x - 1, self.y), (self.x + 1, self.y), (self.x, self.y - 1), (self.x, self.y + 1), (self.x, self.y))

        for i, d in enumerate(directions[:-1]):
            if self.tumor.cancer_at(d[0], d[1]):
                w8s[i] = 0 # Can't move there if someone exists there

        ox, oy = self.x, self.y

        if np.sum(w8s) == 0:
            return self.x, self.y # only occurs if it's surrounded, so stay still cause ya can't move.

        w8s /= np.sum(w8s)

        self.x, self.y = directions[np.random.choice((0, 1, 2, 3, 4), p=w8s)] # Actually do the move.

        self.tumor.move((ox, oy), (self.x, self.y))

        self.stayed = (self.x == ox and self.y == oy)

        return self.x, self.y


def id(x):
    return x.id


def single_lap(x, i, j, dx=1):
    return (x[i + 1, j] + x[i - 1, j] + x[i, j + 1] + x[i, j - 1] - 4 * x[i, j]) / (dx ** 2)


def lap(x, dx=1):
    return (x[2:, 1:-1] + x[1:-1, 2:] + x[:-2, 1:-1] + x[1:-1, :-2] - 4 * x[1:-1, 1:-1]) / (dx ** 2)


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

phi = np.zeros((N, N)) # The chemo-attractant field

tumor = Tumor() # Collection of cancer cells.

cancer_cells = [CancerCell(N / 2-25, N / 2-25), CancerCell(N / 2-25, N / 2 - 1-25), CancerCell(N / 2 - 1-25, N / 2 - 1-25),
                CancerCell(N / 2 - 1-25, N / 2-25)]

cancer_cells += [CancerCell(N / 2+25, N / 2+25), CancerCell(N / 2+25, N / 2 - 1+25), CancerCell(N / 2 - 1+25, N / 2 - 1+25),
                CancerCell(N / 2 - 1+25, N / 2+25)]

for c in cancer_cells:
    tumor.add(c)
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

    plotit(t)

    print "Iteration took", (time.time()-st), "seconds"

plotit()