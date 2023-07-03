import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import random
from pprint import pprint


# Simulation parameters
L = 10  # Size of the habitat
Nr = 900  # Number of rabbits
Nw = 100  # Number of wolves
sigma = 0.5   # Step size for rabbits and wolves
trd = 100  # Rabbit lifespan
twd = 50  # Wolf lifespan without food
rc = 0.5  # Rabbit detection range
pr = 0.02  # Rabbit reproduction rate
pwe = 0.02  # Wolf success rate in hunting
pwr = 0.02  # Wolf reproduction rate

# Rabbits: [x, y, age]
X = np.random.uniform(low=0.0, high=L, size=(Nr, 1))
Y = np.random.uniform(low=0.0, high=L, size=(Nr, 1))
age = np.random.randint(low=1, high=trd, size=(Nr, 1))

Rabbits = np.concatenate([X, Y, age], axis=1)

# Wolves: [x, y, age]
X = np.random.uniform(low=0.0, high=L, size=(Nw, 1))
Y = np.random.uniform(low=0.0, high=L, size=(Nw, 1))
age = np.ones((Nw, 1))

Wolves = np.concatenate([X, Y, age], axis=1)

num_steps = 3000
list_rab = []
list_wol = []
for t in range(num_steps):
    cell_size = rc
    num_cell = int(L / (2 * cell_size))
    cell = [[] for i in range(num_cell * num_cell)]

    num_rab = Rabbits.shape[0]
    num_wol = Wolves.shape[0]
    # print(num_rab)
################################################################################
    # rabbit behaviours
################################################################################
    # The moves of rabbits
    coor_rab = Rabbits[:, 0: 2]
    dir_x = np.random.randint(low=0, high=2, size=(num_rab, 1)) 
    dir_y = np.ones((num_rab, 1)) - dir_x
    direction = np.concatenate([dir_x, dir_y], axis=1) # move direction
    step = np.random.normal(loc=0, scale=sigma, size=(num_rab, 1))
    move = np.multiply(direction, step)
    new_coor = coor_rab + move
    # The upper bound is connect to the lower bound, 
    #   lower to npper, left to right, right to left. 
    # The rabbit will appear on the other side of the habitat
    new_coor = np.mod(new_coor, L)
    Rabbits[:, 0: 2] = new_coor

    # rabbits die
    age_rab = Rabbits[:, 2] + 1
    Rabbits[:, 2] = age_rab
    mask = age_rab <= trd
    Rabbits = Rabbits[mask, :]
    num_rab = Rabbits.shape[0] # Update the total number of rabbits

    # The replications of rabbits
    mask = np.random.choice([False, True], size=(num_rab), p=[1-pr, pr])
    new_rab = Rabbits[mask, :]
    new_rab[:, 2] = 1
    Rabbits = np.concatenate([Rabbits, new_rab], axis=0)
    num_rab = Rabbits.shape[0]

    # create a cell list for rabbits
    for i in range(num_rab):
        x = Rabbits[i, 0]
        y = Rabbits[i, 1]
        cell_x = int(x / 1)
        cell_y = int(y / 1)
        cell_id = cell_x + cell_y * L
        cell[cell_id] = cell[cell_id] + [i]
    # print(num_rab)
################################################################################
    # wolf behaviours
################################################################################
    # The moves of wolf
    coor_wol = Wolves[:, 0: 2]
    dir_x = np.random.randint(low=0, high=2, size=(num_wol, 1))
    dir_y = np.ones((num_wol, 1)) - dir_x
    direction = np.concatenate([dir_x, dir_y], axis=1) # move direction
    step = np.random.normal(loc=0, scale=sigma, size=(num_wol, 1))
    move = np.multiply(direction, step) # move distance
    new_coor = coor_wol + move # new coordinates
    new_coor = np.mod(new_coor, L)
    Wolves[:, 0: 2] = new_coor

    # # wolf die
    age_wol = Wolves[:, 2] + 1
    Wolves[:, 2] = age_wol
    mask = age_wol <= twd
    Wolves = Wolves[mask, :]
    num_wol = Rabbits.shape[0] # Update the total number of wolves
    
    # # hunting
    eatten_list = np.ones(num_rab)
    for ind_wol, wol in enumerate(Wolves):
        x, y = wol[0], wol[1]
        cw_x = int(x / 1)
        cw_y = int(y / 1)
        num_cell = cw_x + cw_y * L
        # (cw_x - 1) + (cw_y - 1) * 10       cw_x + (cw_y - 1) * 10       (cw_x + 1) + (cw_y - 1) * 10 
        # (cw_x - 1) + cw_y * 10             cw_x + cw_y * 10             (cw_x + 1) + cw_y * 10
        # (cw_x - 1) + (cw_y + 1) * 10       cw_x + (cw_y + 1) * 10       (cw_x + 1) + (cw_y + 1) * 10
        n0 = ((cw_x - 1) % L) + ((cw_y - 1) % L) * L
        n1 = cw_x             + ((cw_y - 1) % L) * L 
        n2 = ((cw_x + 1) % L) + ((cw_y - 1) % L) * L

        n3 = ((cw_x - 1) % L) + cw_y * L 
        n4 = cw_x             + cw_y * L
        n5 = ((cw_x + 1) % L) + cw_y * L

        n6 = ((cw_x - 1) % L) + ((cw_y + 1) % L) * L   
        n7 = cw_x             + ((cw_y + 1) % L) * L    
        n8 = ((cw_x + 1) % L) + ((cw_y + 1) % L) * L

        # check rabbits in the neighbour cell lists
        bound = cell[n0] + cell[n1] + cell[n2] + cell[n3] + cell[n4] + cell[n5] + cell[n6] + cell[n7] + cell[n8]

        for ind_rab in bound:
            if eatten_list[ind_rab] == 1:
                rab = Rabbits[ind_rab, :]
                # distance of the rabbit and the wolf
                d_x = abs(Rabbits[ind_rab, 0] - x)
                d_y = abs(Rabbits[ind_rab, 1] - y)
                if d_x > L/2:
                    d_x = L - d_x
                if d_y > L/2:
                    d_y = L - d_y
                dist = d_x**2 + d_y**2
                if dist < rc**2 and np.random.uniform() < pwe:
                    # wolf eatting
                    Wolves[ind_wol, 2] = 0
                    eatten_list[ind_rab] = 0
                    if np.random.uniform() < pwr:
                        new_wolf = np.array([[x, y, 1]])
                        Wolves = np.concatenate([Wolves, new_wolf], axis=0)


    eatten_list = np.array(eatten_list).astype(bool)   
    Rabbits = Rabbits[eatten_list, :]
    num_rab = Rabbits.shape[0]
    num_wol = Wolves.shape[0]
    list_rab.append(num_rab)
    list_wol.append(num_wol)

    if t % 300 == 0:
        print(t)

plt.plot(list_rab, color='red', label='Rabbits')
plt.plot(list_wol, color='blue', label='Wolves')
plt.legend()
plt.show()