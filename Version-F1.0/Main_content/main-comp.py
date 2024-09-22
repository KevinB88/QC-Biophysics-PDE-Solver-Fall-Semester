import func_comp as comp
import func_tab
import func_plot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from datetime import datetime


def plt_multiple_csvs(filepaths, x_column, y_column, plot_title, labels=None):
    plt.figure(figsize=(10, 6))
    for i, filepath in enumerate(filepaths):
        data = pd.read_csv(filepath)

        x_data = data[x_column]
        y_data = data[y_column]

        label = labels[i] if labels else f"File {i+1}"

        plt.plot(x_data, y_data, label=label)

        plt.title(plot_title)
        plt.xlabel(x_column)
        plt.ylabel(y_column)

        plt.legend()
        plt.grid(True)
        plt.show()


def quantify(domain_rings, domain_rays, r, d, t, domain_tube_placements):
    d_radius = r / domain_rings
    d_theta = ((2 * math.pi) / domain_rays)
    d_time = (0.1 * min(d_radius * d_radius, d_theta * d_theta * d_radius * d_radius)) / (2 * d)
    time_steps = math.ceil(t / d_time)
    phi_center = 1 / (math.pi * (d_radius * d_radius))

    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f'Current date of operation: {current_time}')
    print(f'Radial Curves: {rings}')
    print(f'Angular Rays: {rays}')
    print(f'Domain Radius: {r}')
    print(f'Diffusion Coefficient: {d}')
    print(f'Run Time: {t}')
    print(f'Delta Radius: {d_radius}')
    print(f'Delta Theta: {d_theta}')
    print(f'Delta Time: {d_time}')
    print(f'Time steps: {time_steps}')
    print(f'Initial Central Density: {phi_center}')
    print(f'Velocity: {v}')
    print(f'Microtubule placements: {domain_tube_placements}')


if __name__ == "__main__":
    rings = 10
    rays = 10
    sim_time = 1
    v = -20
    N = [0, 2, 3]
    N = np.unique(np.sort(N))
    a = 10
    b = 10
    radius = 1
    d_c = 1

    output, td_time, d_rad, diffusive_layer = comp.solve_short(rings, rays, radius, d_c, sim_time, True, 1, a, b, v, N, ['global-max'])


