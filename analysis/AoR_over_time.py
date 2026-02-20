import numpy as np
import pandas as pd
import os
import re
import sys

# Add parent directory to path to import avalanching_analysis
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from avalanching_analysis.plotter import *
from avalanching_analysis.read_find_functions import *
from avalanching_analysis.analysis_functions import *
    
if __name__ == '__main__':
    # Inputs:

    mixture = "single_material"

    # mixture = "mixture"
    # sim_types = ["low_susc"]
    sim_types = ["high_susc"]

    # sim_types = ["75_high"]#, "50_high", "75_high"]

    fields = ["0000","0010","0100","1000"]
    # fields = ["0000","1000"]

    # root_directory_plots = "/Users/sikka-mac/Obsidian/Second_Brain/01 - Projects/01 - Research/Writing/Paper 2 (Granular Matter)/figures/final_frames/"

    output_csv_path = "/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/values/aor_rms_values_overtime.csv"  # Path to save CSV

    all_results = []  # List to store results from each field

    y_bin_num = 5 # Number of y bins to create for slicing the data

    for sim_type in sim_types:
        print("--------------------------------")
        print(f"\nProcessing sim_type: {sim_type}\n")
        print("--------------------------------")
        for field in fields:
            print("--------------------------------")
            print(f"\nProcessing field: {field}\n")
            print("--------------------------------")
            folder_path = "/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/" + f"{mixture}/{sim_type}/field_{field}uT/"
            frame_i = find_initial_frame(folder_path)
            frame_f = find_end_timestep_from_dumps(folder_path + "results")

            print("--------------------------------")
            print(f"\nProcessing multiple frames for continuos AoR calculation\n")
            print("--------------------------------")

            frames = np.arange(frame_i + 100000, frame_f, step=5000)  # Adjust step as needed

            for frame in frames:
                df_initial, df_final = read_dump_files(folder_path + "results/", frame_i, frame)

                if df_final is None or df_initial is None:
                    print(f"Skipping frame {frame} due to missing data.")
                    continue

                top_par_final = find_top_particles_with_epsilon(df_final)


                angleRepose, SE_AoR, best_r_square, best_upper_limit, par_fit, lower_limit, optimal_model = calculate_AoR(df_initial, df_final)


                # Create a dictionary with the results for the current field
                results = {
                    'mixture' : mixture,
                    'simtype': sim_type,
                    'field': field,
                    'frame': frame,
                    'AoR (deg)': angleRepose,
                    'SE_AoR': SE_AoR,
                    'R2': best_r_square,
                    'lower_limit': lower_limit,

                }
                all_results.append(results)
                # print("\nData in separate lines:")
                # for key, value in results.items():
                #     # Check if the value is a numpy float and convert to regular float
                #     if isinstance(value, np.float64) or isinstance(value, np.float32):
                #         print(f"{key}: {float(value)}")
                #     else:
                #         print(f"{key}: {value}")

    
    # Convert the list of results to a DataFrame
    results_df = pd.DataFrame(all_results)

    # Check if the CSV file exists
    file_exists = os.path.isfile(output_csv_path)

    # Save the results to the CSV file
    if file_exists:
        # If the file exists, append the results
        results_df.to_csv(output_csv_path, mode='a', header=False, index=False)
    else:
        # If the file doesn't exist, create it with headers
        results_df.to_csv(output_csv_path, mode='w', header=True, index=False)

    print(f"\nResults saved to: {output_csv_path}")

