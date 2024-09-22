
import func_calc_mfpt as calc
import func_tab as tb
import func_plots as plt

import multiprocessing as mp
import filepaths as fp
from datetime import datetime
from functools import partial
import time


def solve(N_param, rg_param, ry_param, v_param, ext_param, w_param):
    M = rg_param
    N = ry_param
    d_c = 1
    radius = 1
    sim_time = 1
    v = v_param
    ext_fact = ext_param

    mfpt = calc.solve(M, N, radius, d_c, sim_time, True, 1, w_param, w_param, v, N_param, True, ext_fact)
    print(f'Microtubule Configuration: {N_param}')
    print()
    return {f'W: {w_param}', f'MFPT: {mfpt}'}


def parallel_process(rg_param, ry_param, v_param, ext_param, N_list, w_low_bound, w_high_bound):
    w_list = []
    lower_bound = w_low_bound
    upper_bound = w_high_bound
    for x in range(lower_bound, upper_bound+1):
        w_list.append(10 ** x)

    print(w_list)

    if len(w_list) > 6:
        process_count = 6
    else:
        process_count = len(w_list)

    current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    data_filepath = tb.create_directory(fp.mfpt_data_fp, current_time)

    for n in range(len(N_list)):
        with mp.Pool(processes=process_count) as pool:
            mfpt_results = pool.map(partial(solve, N_list[n], rg_param, ry_param, v_param, ext_param), w_list)
        print(mfpt_results)

        tb.data_extraction_pandas(mfpt_results, data_filepath, f'MFPT_Results_N={len(N_list[n])}_v={velocity}_')
        time.sleep(1)

    plt.plot_all_csv_in_directory(data_filepath, N_list, data_filepath, f'MFPT_versus_W_v={velocity}_', True)


if __name__ == "__main__":

    rings = 32
    rays = 32
    velocity = -20
    extension_factor = 50

    microtubule_configs = [
        [0, 8, 16, 24],
        [0, 4, 8, 12, 16, 20, 24, 28],
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
    ]

    w_low = -2
    w_high = 4

    parallel_process(rings, rays, velocity, extension_factor, microtubule_configs, w_low, w_high)

'''
    Common Microtubule configurations

    [0, 8, 16, 24]
    [0, 4, 8, 12, 16, 20, 24, 28]
    [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
'''


