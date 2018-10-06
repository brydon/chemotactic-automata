#!/software/.admin/bins/bin/python2.7

"""

Things I could add:
    vector fields for adhesion (?)
    re-create the anderson05 exp w/o adhesion (only chemotaxis)
    Change chemo-dispersion to? O2?

"""

import numpy as np
import time
import sys
import cancer
from multiprocessing import Pool
import pickle as pik
import os

log_str = ""

if os.path.isdir(sys.argv[-2]):
    os.chdir(sys.argv[-2])
else:
    os.mkdir(sys.argv[-2])
    os.chdir(sys.argv[-2])


def initialise_tumor(_tumor, _phi):
    # Single, central solid tumour site
    # starting_positions = [(N/2, N/2), (N/2 - 1, N/2), (N/2, N/2 - 1), (N/2 - 1, N/2 - 1)]

    # Two solid tumour sites
    starting_positions = [(N / 4, N / 4), (N / 4 - 1, N / 4), (N / 4, N / 4 - 1), (N / 4 - 1, N / 4 - 1),
                          (3 * N / 4, 3 * N / 4),
                          (3 * N / 4 - 1, 3 * N / 4), (3 * N / 4, 3 * N / 4 - 1), (3 * N / 4 - 1, 3 * N / 4 - 1)]

    # Random Tumor Sites, non-solid
    # starting_positions = [(i, j) for i in range(N) for j in range(N) if np.random.random() > 0.9]
    # starting_positions = [(i, j) for i in range(10,40) for j in range(10, 40)]

    cancer_cells = [cancer.CancerCell(x, y, _tumor) for (x, y) in starting_positions]

    for c in cancer_cells:
        _phi[c.x, c.y] = 1.0


def scale_down(a, fact=5):
    n = a.shape[0]
    b = np.empty((n/fact, n/fact))
    for i in range(n/fact):
        for j in range(n/fact):
            b[i, j] = np.mean(a[i*fact:(i+1)*fact, j*fact:(j+1)*fact])
    return b


def read_last_run():
    dat_files = [x for x in os.listdir(os.getcwd()) if x.endswith(".dat") and x.startswith(sys.argv[-1])]
    dat_files.sort()

    _cur_done = -1

    if len(dat_files) > 0:
        fn = dat_files[-1]
        print "Opening run", fn

        underscore = fn[::-1].find("_")
        dot = fn.find(".")

        _cur_done = int(fn[-underscore:dot])

        with open(fn, 'rb') as _f:
            _tum, _ph = pik.load(_f)

        return (_cur_done, _tum, _ph)

    return _cur_done


def log_print(*args):
    global log_str
    print args
    log_str += " ".join([str(x) for x in args])
    log_str += "\n"


if __name__ == "__main__":
    N = cancer.N

    chem_lo = mig_lo = 0.1
    chem_hi = mig_hi = 10.
    chem_reg = mig_reg = 1.

    np.set_printoptions(threshold=40*40+1)

    if sys.argv[-1][-3:] == "-nn":
        cancer.gamma = 0.
    else:
        cancer.gamma = 1.

    if sys.argv[-1][:2] == "cl" or sys.argv[-1][:2] == "mh":
        print "Chemotaxis-lo, Migration-hi run"
        cancer.D = mig_hi
        cancer.mu = chem_lo
    elif sys.argv[-1][:2] == "ch" or sys.argv[-1][:2] == "ml":
        print "Chemotaxis-hi, Migration-lo run"
        cancer.D = mig_lo
        cancer.mu = chem_hi
    else:
        print "Chemotaxis-regular, Migration-regular run"
        cancer.D = mig_reg
        cancer.mu = chem_reg

    num_steps = 100

    phi = np.zeros((N, N))  # The chemo-attractant field
    tumor = cancer.Tumor()  # Collection of cancer cells.

    initialise_tumor(tumor, phi)

    cur_done = read_last_run()

    if cur_done != -1:
        cur_done, tumor, phi = cur_done

    print scale_down(phi, 10)

    for t in range(cur_done + 1, num_steps):
        log_print("Frame", t, len(tumor.cancer_cells()))

        if len(tumor.cancer_cells()) == 0:
            break

        for cell in tumor.cancer_cells():
            phi[cell.x, cell.y] = 1.

        tm = time.time()
        """ Disperse Chemo-attractant """
        Nt = int(1/float(cancer.dt))

        _const = cancer.Dp * cancer.dt / (cancer.dx ** 2)

        for _x in [Nt/4, Nt/4, Nt/4, Nt-3*(Nt/4)]:
            for _ in range(_x):
                phi[1:-1, 1:-1] += _const * (phi[2:, 1:-1] + phi[1:-1, 2:] + phi[:-2, 1:-1] + phi[1:-1, :-2]
                                             - 4 * phi[1:-1, 1:-1])
            print "Check-point"

        log_print("\tdispersed in", time.time() - tm)

        tm = time.time()

        """ Cell lifecycle """
        for cell in tumor.cancer_cells():
            cell.move(phi)

        for cell in tumor.cancer_cells():
            cell.life_cycle()

        for cell in tumor.cancer_cells():
            phi[cell.x, cell.y] = 1.

        log_print("\tlife cycled in", time.time() - tm)
        tm = time.time()

        with open('%s_out_dat_frame_%03d.dat' % (sys.argv[-1], t), 'wb') as f:
            pik.dump([tumor, phi], f)

        log_print("\tdumped in", time.time() - tm)

        log_str += str(time.asctime()) + "\n"

        with open('%s_log.txt' % sys.argv[-1], 'a') as f:
            f.write(log_str)

        break
