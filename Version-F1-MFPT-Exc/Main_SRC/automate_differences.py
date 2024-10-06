import pandas as pd
from datetime import datetime
import os


def calculate_diff(csv_file_1, csv_file_2, file_name, file_path ):
    df_1 = pd.read_csv(csv_file_1)
    df_2 = pd.read_csv(csv_file_2)

    MFPT_diff = []

    for val_1, val_2 in zip(df_1['MFPT'], df_2['MFPT']):
        MFPT_diff.append(((val_1 - val_2)/val_1) * 100)
    print(MFPT_diff)

    df = pd.DataFrame({
        'W': df_1['W'],
        'MFPT (40x40)': df_1['MFPT'],
        'MFPT (48x48)': df_2['MFPT'],
        '% Diff': MFPT_diff
    })

    current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    full_file_name = f'{file_name}_{current_time}.csv'

    if not os.path.exists(file_path):
        os.makedirs(file_path)

    df.to_csv(os.path.join(file_path, full_file_name), sep=',', index=False)


if __name__ == "__main__":
    MFPT_16 = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/10-4-Saturation-test/MFPT_Results_N=8_v=-10__2024-10-03_08-53-48_16x16.csv'
    MFPT_24 = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/10-4-Saturation-test/MFPT_Results_N=8_v=-10__2024-10-04_11-37-23_24x24.csv'
    MFPT_32 = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/10-4-Saturation-test/MFPT_Results_N=8_v=-10__2024-10-03_04-26-05_32x32.csv'
    MFPT_40 = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/10-4-Saturation-test/MFPT_Results_N=8_v=-10__2024-10-04_07-43-13_40x40.csv'
    MFPT_48 = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/data-package-10-5/NxN_MFPT_N=8_Results/5_MFPT_Results_N=8_v=-10__2024-10-05_22-14-57_48x48.csv'

    destination_path = '/Users/kbedoya88/Desktop/QC24-Fall-Semester/Computer-Science-Research/Biophysics/Kogan-Bedoya-Algorithm/Biophysics_Numerical_PDE_Solver_Fall_Semester/Version-F1-MFPT-Exc/Data/10-4-Saturation-test/Saturation-Report'

    calculate_diff(MFPT_40, MFPT_48, 'diff_test_N_8_', destination_path)
