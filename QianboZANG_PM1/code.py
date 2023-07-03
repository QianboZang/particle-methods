import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.animation import FuncAnimation
import random


def get_neighbors(grid, x, y):
    """
    0 1 2
    3 x 4
    5 6 7
    """
    neighbors = [
        grid[(x+100-1)%100, (y+100-1)%100], grid[(x+100-1)%100, (y+100)%100], grid[(x+100-1)%100, (y+100+1)%100],
        grid[(x+100)%100, (y+100-1)%100], grid[(x+100)%100, (y+100+1)%100],
        grid[(x+100+1)%100, (y+100-1)%100], grid[(x+100+1)%100, (y+100)%100], grid[(x+100+1)%100, (y+100+1)%100]
    ]
    return neighbors


def nei_index(grid, x, y):
    index = [
        [(x+100-1)%100, (y+100-1)%100], [(x+100-1)%100, (y+100)%100], [(x+100-1)%100, (y+100+1)%100],
        [(x+100)%100, (y+100-1)%100], [(x+100)%100, (y+100+1)%100],
        [(x+100+1)%100, (y+100-1)%100], [(x+100+1)%100, (y+100)%100], [(x+100+1)%100, (y+100+1)%100]
    ]
    return index


def check_happiness(grid, x, y, H):
    neighbors = get_neighbors(grid, x, y)
    if grid[x, y] == 0:
        return "empty"
    if neighbors.count(grid[x, y]) >= H:
        return "happy"
    return "unhappy"


def get_happy_unhappy_empty(grid, H):
    happy_agents = []
    unhappy_agents = []
    empty_agents = []

    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            if check_happiness(grid, x, y, H) == "happy":
                happy_agents.append([x, y])
            elif check_happiness(grid, x, y, H) == "unhappy":
                unhappy_agents.append([x, y])
            else:
                empty_agents.append([x, y])
    return empty_agents, unhappy_agents


def find_empty(grid, x, y, empty_agents):
    return random.choice(empty_agents)


# 2. find nearest empty cell in 4 horizontal directions
# def find_empty(grid, x, y, empty_agents, unhappy_agents):
#     for coord in unhappy_agents:
#         x = coord[0]
#         y = coord[1]
#         for i in range(1, n+1):
#             if grid[x-i, y] == 0:
#                 return [x-i, y]
#             elif grid[x+i, y] == 0:
#                 return [x+i, y]
#             elif grid[x, y-i] == 0:
#                 return [x, y-i]
#             elif grid[x, y+i] == 0:
#                 return [x, y+i]
# 
#   
# 3. find nearest empty cell using norm
# def find_empty(grid, x, y, empty_agents, unhappy_agents):
#     for coord in unhappy_agents:
#         x = coord[0]
#         y = coord[1]
#         tmp = [[x, y]] * len(empty_agents)
#         dist = np.linalg.norm(np.array(empty_agents) - np.array(tmp), axis=1) 
#     return empty_agents[np.argmin(dist)]


def move(grid, empty_agents, unhappy_agents, H):
    # move unhappy agents to the nearest empty cell
    # 1. iterate all unhappy agents
    for x, y in unhappy_agents:
        # 2. find the nearest empty cell
        nearest_empty = find_empty(grid, x, y, empty_agents)
        grid[nearest_empty[0], nearest_empty[1]] = grid[x, y]
        grid[x, y] = 0
        # 3. delete the unhappy agent [x, y] from the unhappy_agents
        # print(len(unhappy_agents))
        unhappy_agents.remove([x, y])
        # print(len(unhappy_agents))
        # 4. add the [x, y] to the empty_agents
        empty_agents.remove(nearest_empty)
        empty_agents.append([x, y])

        # 5. move the unhappy agent to the nearest empty cell, and update the grid
        if check_happiness(grid, nearest_empty[0], nearest_empty[1], H) == "unhappy":
            unhappy_agents.append(nearest_empty)
        
        # 6. check happy or unhappy for the original neighbour of (x, y)
        index = nei_index(grid, x, y)
        for n in index:
            if check_happiness(grid, n[0], n[1], H) == "unhappy":
                if n not in unhappy_agents:
                    unhappy_agents.append(n)


def main():
    grid_size = 100
    empty_cells = 0.1
    agent_types = 2 # empty: 0, red: 1, and blue: 2
    emp_distribution, red_distribution, blue_distribution = 0.1, 0.45, 0.45
    H = 4 # number of neighbors of the same type for the agent to be "happy"
    cmap = colors.ListedColormap(['white', 'red', 'blue'])
    max_iter = 100
    
    # initail grid
    grid = np.random.choice(
        [0, 1, 2], size=(grid_size, grid_size), 
        p=[emp_distribution, red_distribution, blue_distribution]
    )
    print(grid)
    empty_agents, unhappy_agents = get_happy_unhappy_empty(grid, H)
    imgs = []
    iter = 0
    while(len(unhappy_agents) > 0):
        print(f"iter: {iter}, unhappy_agents: {len(unhappy_agents)}, empty_agents: {len(empty_agents)}")
        # move unhappy agents to the nearest empty cell
        imgs.append(grid.copy())
        move(grid, empty_agents, unhappy_agents, H)
        # update grid
        iter += 1

    fig = plt.figure(figsize=(10, 10))
    # Function to plot each picture in the list
    def animate(i):
        plt.imshow(imgs[i], cmap=cmap)

    # Create the animation
    ani = FuncAnimation(fig, animate, frames=len(imgs), repeat=False)
    plt.show()
    # ani.save('animation.mp4', writer='ffmpeg')


if __name__ == "__main__":
    main()