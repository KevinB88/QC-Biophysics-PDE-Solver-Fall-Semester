
import func_calc_mfpt as calc
import func_tab as tb
import func_plots as plt

import multiprocessing as mp
import filepaths as fp
from datetime import datetime
from functools import partial
import time


def set_cores(core_input):
    fp.default_cpu_core_count_fp = 8
    if core_input > fp.default_cpu_core_count_fp:
        core_input = fp.default_cpu_core_count_fp
    return core_input


def solve(N_param, rg_param, ry_param, v_param, ext_param, w_param):
    M = rg_param
    N = ry_param
    d_c = 1
    radius = 1
    sim_time = 1
    v = v_param
    ext_fact = ext_param

    # mfpt = calc.solve_mass_decay(M, N, radius, d_c, sim_time, True, 1, w_param, w_param, v, N_param, True, ext_fact)

    mfpt = calc.solve_mass_decay(M, N, radius, d_c, sim_time, True, 1, w_param, w_param, v, N_param, True, ext_fact)

    print(f'Microtubule Configuration: {N_param}')
    print()
    return {f'W: {w_param}', f'MFPT: {mfpt}'}


def parallel_process(rg_param, ry_param, v_param, ext_param, N_list, w_low_bound, w_high_bound, cores):
    w_list = []
    lower_bound = w_low_bound
    upper_bound = w_high_bound
    for x in range(lower_bound, upper_bound+1):
        w_list.append(10 ** x)

    print(w_list)

    if len(w_list) > 3:
        process_count = 3
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
    velocity = -10
    extension_factor = 20
    # input the number of CPU cores

    microtubule_configs = [
        [0, 4, 8, 12, 16, 20, 24, 28]
    ]

    w_low = -2
    w_high = 4

    path = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/data-package-10-5/NxN_MFPT_N=8_Results'
    data_path = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data'

    first = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/data-package-10-5/NxN_MFPT_N=8_Results/MFPT_Results_N=8_v=-10__2024-10-03_08-53-48_16x16.csv'
    second = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/data-package-10-5/NxN_MFPT_N=8_Results/MFPT_Results_N=8_v=-10__2024-10-04_11-37-23_24x24.csv'
    third = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/data-package-10-5/NxN_MFPT_N=8_Results/MFPT_Results_N=8_v=-10__2024-10-03_04-26-05_32x32.csv'
    fourth = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/data-package-10-5/NxN_MFPT_N=8_Results/MFPT_Results_N=8_v=-10__2024-10-04_07-43-13_40x40.csv'
    fifth = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/data-package-10-5/NxN_MFPT_N=8_Results/MFPT_Results_N=8_v=-10__2024-10-05_22-14-57_48x48.csv'

    plt.plot_all_csv_in_directory_manual([first, second, third, fourth, fifth], ['16x16', '24x24', '32x32', '40x40', '48x48'], data_path, 'Test of convergence, N=8', True, custom_label=True)

