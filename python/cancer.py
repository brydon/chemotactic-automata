import numpy as np


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
        return (x, y) in self.cells.keys()

    def cancer_cells(self):
        return sorted(self.cells.values(), key=lambda cel: cel.id)

    def move(self, old, new):
        self.cells[new] = self.cells.pop(old)


class CancerCell:
    PROLIF_AGE = 1

    def __init__(self, x, y, tum=None):
        self.x, self.y = x, y
        self.tumor = tum
        self.stayed = False
        self.age = 0
        self.tumor.add(self)

    def mitosis(self):
        if self.age >= self.PROLIF_AGE:
            free_spots = []
            potential_spots = ((self.x+1, self.y), (self.x-1, self.y), (self.x, self.y+1),(self.x, self.y-1))

            for spot in potential_spots:
                if not self.tumor.cancer_at(spot[0], spot[1]):
                    free_spots.append(spot)

            if len(free_spots) == 0:
                return False

            new_spot = np.random.choice(potential_spots)

            CancerCell(new_spot[0], new_spot[1], self.tumor)
            self.tumor.add(self)
            return True

    def life_cycle(self):
        self.age += 1
        reproduced = self.mitosis()

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
