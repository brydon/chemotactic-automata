import numpy as np
import time

D = 1       # migration strength
mu = 1      # chemo-tactic strengh
gamma = 1   # 1 or 0, turns on/off the GWN term.

N = 400
dx = 1./N
L = N*0.0025  # [#cells]*[cm/cell]

tau_by_l = 16 * 3600 / (L ** 2)
Dp = tau_by_l * (10 ** (-5))
dt = 1. * (dx ** 2 / (4*Dp))


def single_lap(x, i, j, _dx=1):
    return (x[i + 1, j] + x[i - 1, j] + x[i, j + 1] + x[i, j - 1] - 4 * x[i, j]) / (_dx ** 2)


def lap(x, _dx=1):
    return (x[2:, 1:-1] + x[1:-1, 2:] + x[:-2, 1:-1] + x[1:-1, :-2] - 4 * x[1:-1, 1:-1]) / (_dx ** 2)


class Tumor:
    def __init__(self):
        self.cells = {}
        self.next_id = 0
        self.tid = time.time()

    def age(self):
        for cell in self.cells:
            cell.age += 1

    def add(self, cell):
        self.cells[(cell.x, cell.y)] = cell
        cell.tumor = self
        cell.id = self.next_id
        self.next_id += 1

    def cancer_at(self, x, y):
        if (x, y) in self.cells.keys() and not self.cells[(x, y)].dead:
                return True
        return False

    def cancer_cells(self):
        return sorted(self.cells.values(), key=lambda cel: cel.id)

    def move(self, old, new):
        self.cells[new] = self.cells.pop(old)

    def size(self):
        return len(self.cells)


class CancerCell:
    PROLIF_AGE = 1

    def __init__(self, x, y, tum=None):
        self.x, self.y = int(x), int(y)
        self.tumor = tum
        self.stayed = False
        self.age = 0
        self.tumor.add(self)
        self.dead = False
        self.proliferating = False

    @staticmethod
    def consumption():
        return 0.57

    def mitosis(self):
        if self.age >= self.PROLIF_AGE:
            free_spots = []
            potential_spots = ((self.x+1, self.y), (self.x-1, self.y), (self.x, self.y+1), (self.x, self.y-1))

            for spot in potential_spots:
                if not self.tumor.cancer_at(spot[0], spot[1]) and 0 <= spot[0] < N and 0 <= spot[1] < N:
                    free_spots.append(spot)

            if len(free_spots) == 0:
                return False

            new_spot = free_spots[np.random.randint(0, len(free_spots))]

            CancerCell(new_spot[0], new_spot[1], self.tumor)

            return True
        return False

    def life_cycle(self):
        if self.dead:
            return False
        self.age += 1
        self.proliferating = self.mitosis()
        return self.proliferating

    def die(self):
        self.dead = True

    def move(self, phi):
        """

        :param phi:
        :return:
        """

        if self.x == 0 or self.x == phi.shape[0] - 1 or self.y == 0 or self.y == phi.shape[1] - 1:
            return self.x, self.y  # Bdary conditions

        _const = dt/(dx**2)

        _Gamma = np.random.randn()*gamma
        noise_term = _Gamma*dt/(2.*dx*np.sqrt(2*D))

        c_left = (D - mu * 0.25 * (phi[self.x + 1, self.y] - phi[self.x - 1, self.y])) * (dt/dx**2) - noise_term
        c_right = (D + mu * 0.25 * (phi[self.x + 1, self.y] - phi[self.x - 1, self.y])) * (dt/dx**2) + noise_term
        c_down = (D - mu * 0.25 * (phi[self.x, self.y + 1] - phi[self.x, self.y - 1])) * (dt/dx**2) - noise_term
        c_up = (D + mu * 0.25 * (phi[self.x, self.y + 1] - phi[self.x, self.y - 1])) * (dt/dx**2) + noise_term
        c_stay = 1 - dt/(dx ** 2) * (4 * D + mu * single_lap(phi, self.x, self.y))

        w8s = [c_left, c_right, c_down, c_up, c_stay]

        directions = ((self.x - 1, self.y), (self.x + 1, self.y), (self.x, self.y - 1),
                      (self.x, self.y + 1), (self.x, self.y))

        for i, d in enumerate(directions[:-1]):
            if self.tumor.cancer_at(d[0], d[1]):
                w8s[i] = 0  # Can't move there if someone exists there

        ox, oy = self.x, self.y

        if np.sum(w8s) == 0:
            return self.x, self.y  # only occurs if it's surrounded, so stay still cause ya can't move.

        w8s /= np.sum(w8s)

        r = np.random.random()

        tot = 0
        for i in range(0, len(w8s)):
            if tot < r <= w8s[i]:
                self.x, self.y = directions[i]
                break
            tot += w8s[i]

        self.tumor.move((ox, oy), (self.x, self.y))

        self.stayed = (self.x == ox and self.y == oy)

        return self.x, self.y
