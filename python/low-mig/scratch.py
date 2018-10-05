import numpy as np
import time

# intervals in x-, y- directions, mm
dx = dy = 0.0025
N = 400
L = N*(0.0025) # [#cells]*[cm/cell]

tau_by_l = 16 * (3600) / ((L) ** 2)

dc = tau_by_l * (10 ** (-5))

D = dc

print dc
print L

nx, ny = int(N), int(N)
print nx, ny

dx2, dy2 = dx*dx, dy*dy
dt = dx2 * dy2 / (2 * D * (dx2 + dy2))
print dt

u0 = np.empty((nx, ny))
u = np.empty((nx, ny))
cancer = []


# Initial conditions - ring of inner radius r, width dr centred at (cx,cy) (mm)
cancer

u = u0.copy()

def scale_down(a, fact=5):
    n = a.shape[0]
    b = np.empty((n/fact, n/fact))
    for i in range(n/fact):
        for j in range(n/fact):
            b[i, j] = np.mean(a[i*fact:(i+1)*fact, j*fact:(j+1)*fact])
    return b

np.set_printoptions(threshold=40*40+1)

def do_timestep(u0, u):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )

    u0 = u.copy()
    return u0, u

# Number of timesteps
nsteps = 1

print scale_down(u0, 10)
#fig = plt.figure()
for m in range(nsteps):
    st = time.time()
    u0[10:40, 10:40] = np.max(u0[10:40, 10:40] - 0.57, 0)
    for _ in range(int(1/dt)):
        u0, u = do_timestep(u0, u)
    print time.time()-st
    print m


print scale_down(u, 10)
