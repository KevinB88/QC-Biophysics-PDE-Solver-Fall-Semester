import multiprocessing as mp
import func_comp_peak_time as comp
from functools import partial
import time


def solve(N, w):
    rings = 15
    rays = 15
    sim_time = 1
    v = -10
    d_c = 1
    radius = 1

    # This function could also be used to calculate the MFPT more efficiently as well
    peak_time = comp.solve(rings, rays, radius, d_c, sim_time, True, 1, w, w, v, N, True)
    return {f'W={w}', f'Peak-Time={peak_time}'}


if __name__ == "__main__":

    w_list = [1, 10, 100, 1000, 10000]

    N_list = [
        [0, 1, 3],
        [0, 8, 16, 24],
        [0, 4, 8, 12, 16, 20, 24, 28],
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
    ]

    '''
        corresponding-index:            microtubule-configuration

        0                                [0, 8, 16, 24]
        1                                [0, 4, 8, 12, 16, 20, 24, 28]
        2                                [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
        .                                []
        .                                []
        N
    '''
    config_index = 0
    #
    # comp.solve(10, 10, 1, 1, 1, True, 1, 10000, 10000, -10, [0, 1, 3], True)

    with mp.Pool(processes=5) as pool:
        results = pool.map(partial(solve, N_list[config_index]), w_list)

    print(results)
    # # The following prints a dictionary of {w, peak-time} pairings
    #
