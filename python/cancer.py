import numpy as np

D = 1       # diffusion coeff of the cancer-cells
Dp = 0.01   # diffusion coeff of the chemo-attractant
mu = 1      # uptake parameter of the chemo-attractant by the individual cancer-cells
dx = 1e-2   # dx and dt are for the numerical solutions of the PDE for dispersion of chemoattractant
dt = 1e-5   #
Dc = 0.01


def single_lap(x, i, j, dx=1):
    return (x[i + 1, j] + x[i - 1, j] + x[i, j + 1] + x[i, j - 1] - 4 * x[i, j]) / (dx ** 2)


def lap(x, dx=1):
    return (x[2:, 1:-1] + x[1:-1, 2:] + x[:-2, 1:-1] + x[1:-1, :-2] - 4 * x[1:-1, 1:-1]) / (dx ** 2)


class Tumor:
    def __init__(self):
        self.cells = {}
        self.next_id = 0

    def age(self):
        for cell in self.cells:
            cell.age += 1

    def add(self, cell):
        self.cells[(cell.x, cell.y)] = cell
        cell.tumor = self
        cell.id = self.next_id
        self.next_id += 1

    def cancer_at(self, x, y):
        if (x, y) in self.cells.keys() and not cells[(x,y)].dead:
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
        self.x, self.y = x, y
        self.tumor = tum
        self.stayed = False
        self.age = 0
        self.tumor.add(self)
        self.dead = False

    def consumption(self):
        return 0.57

    def mitosis(self):
        if self.age >= self.PROLIF_AGE:
            free_spots = []
            potential_spots = ((self.x+1, self.y), (self.x-1, self.y), (self.x, self.y+1), (self.x, self.y-1))

            for spot in potential_spots:
                if not self.tumor.cancer_at(spot[0], spot[1]):
                    free_spots.append(spot)

            if len(free_spots) == 0:
                return False

            new_spot = free_spots[np.random.randint(0, len(free_spots))]

            CancerCell(new_spot[0], new_spot[1], self.tumor)
            self.tumor.add(self)
            return True

    def life_cycle(self, o2):
        if o2[c.x, c.y] == 0:
            self.die()
        else:
            self.age += 1
            reproduced = self.mitosis()

    def die(self):
        self.dead = True

    def move(self, phi):
        """
        TBD: Should there be a dx/dt factor in each/most of the Cim terms??
        :param phi:
        :return:
        """
        c_left = (D - mu * 0.25 * (phi[self.x + 1, self.y] - phi[self.x - 1, self.y])) * (dt/dx**2)
        c_right = (D + mu * 0.25 * (phi[self.x + 1, self.y] - phi[self.x - 1, self.y])) * (dt/dx**2)
        c_down = (D - mu * 0.25 * (phi[self.x, self.y + 1] - phi[self.x, self.y - 1])) * (dt/dx**2)
        c_up = (D + mu * 0.25 * (phi[self.x, self.y + 1] - phi[self.x, self.y - 1])) * (dt/dx**2)
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

        self.x, self.y = directions[np.random.choice((0, 1, 2, 3, 4), p=w8s)]  # Actually do the move.

        self.tumor.move((ox, oy), (self.x, self.y))

        self.stayed = (self.x == ox and self.y == oy)

        return self.x, self.y
