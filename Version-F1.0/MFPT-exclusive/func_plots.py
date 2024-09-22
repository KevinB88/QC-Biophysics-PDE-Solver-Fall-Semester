import matplotlib.pyplot as plt
import os
import pandas as pd
from datetime import datetime


def plot_all_csv_in_directory(directory_path, N_labels, filepath, title, save_png=False, show_plt=True, transparent=False):
    plt.figure(figsize=(12, 8))

    i = 0
    for filename in os.listdir(directory_path):

        if filename.endswith('.csv'):
            file_path = os.path.join(directory_path, filename)
            df = pd.read_csv(file_path)

            if 'W' in df.columns and 'MFPT' in df.columns:
                plt.scatter(df['W'], df['MFPT'], label=f'{len(N_labels[i])}', linewidth=(10/(i+1)))
                plt.yscale('log')
                plt.xscale('log')
                i += 1

    plt.xlabel('W', fontsize=30, fontname='Times New Roman')
    plt.ylabel('MFPT', fontsize=30, fontname='Times New Roman')
    plt.title(title, fontsize=40, fontname='Times New Roman')

    plt.xticks(fontsize=30, fontname='Times New Roman')
    plt.yticks(fontsize=25, fontname='Times New Roman')

    plt.legend(fontsize=22, frameon=True, edgecolor='black', loc='best')
    plt.ylim(10**-1, 10**1)

    if save_png:
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        if filepath:
            if not os.path.exists(filepath):
                os.makedirs(filepath)
            file = os.path.join(filepath, f'{title}_date{current_time}.png')
            plt.savefig(file, bbox_inches='tight', transparent=transparent)
            print(f'Plot saved to {filepath}')
    plt.show()