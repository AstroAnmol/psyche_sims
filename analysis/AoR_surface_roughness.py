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
    # sim_types = ["high_susc"]
    sim_types = ["no_susc"]

    # sim_types = ["25_high", "50_high"]#, "75_high"]

    # fields = ["0000","0010","0100","1000"]
    fields = ["0000"]#,"1000"]

    root_directory_plots = "/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/plots/"

    # root_directory_plots = "/Users/sikka-mac/Obsidian/Second_Brain/01 - Projects/01 - Research/Writing/Paper 2 (Granular Matter)/figures/final_frames/"

    output_csv_path = "/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/values/aor_rms_values.csv"  # Path to save CSV

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

            df_initial, df_final = read_dump_files(folder_path + "results/", frame_i, frame_f)

            particle_diameter = df_initial['radius'].iloc[0] * 2  # Assuming all particles have the same radius
            # plot_filename = root_directory_plots + f"{sim_type}_" + f"seed_{seed}_" + f"susc_{susc}_" + f"field_{field}" + ".png"
            
            if df_final is not None and df_initial is not None:

                plot_filename = root_directory_plots + f"{mixture}_{sim_type}_" +  f"field_{field}" + ".png" 
                top_par_final = find_top_particles_with_epsilon(df_final)

                # plot_both_frames(df_final, top_par_final, plot_filename)

                angleRepose, SE_AoR, best_r_square, best_upper_limit, par_fit, lower_limit, optimal_model = calculate_AoR(df_initial, df_final)
                
                plot_both_frames_aor(df_final, par_fit, optimal_model, plot_filename)

                rms_height, se_rms_height, rms_slope, se_rms_slope = calculate_surface_roughness(top_par_final, particle_diameter/2)
                
                # plot_filename = root_directory_plots + f"{sim_type}_" + f"seed_{seed}_" + f"field_{field}" + ".png"
                # Assuming plot_both_frames function is defined elsewhere
                # plot_both_frames_aor(df_final, par_fit, optimal_model, plot_filename)

                # Create a dictionary with the results for the current field
                results = {
                    'mixture' : mixture,
                    'simtype': sim_type,
                    'field': field,
                    'y_bin': 0,
                    'AoR (deg)': angleRepose,
                    'SE_AoR': SE_AoR,
                    'R2': best_r_square,
                    'lower_limit': lower_limit,
                    'RMS_height (mm)': rms_height * 1000,
                    'SE_RMS_height': se_rms_height * 1000,
                    'RMS_slope (deg)': rms_slope,
                    'SE_RMS_slope': se_rms_slope
                }
                all_results.append(results)
                print("\nData in separate lines:")
                for key, value in results.items():
                    # Check if the value is a numpy float and convert to regular float
                    if isinstance(value, np.float64) or isinstance(value, np.float32):
                        print(f"{key}: {float(value)}")
                    else:
                        print(f"{key}: {value}")

                # Print statement to indicate slicing is starting
                print(f"\nSlicing data into {y_bin_num} y bins and plotting each bin separately...")

                df_final = y_bin_slicing(df_final, y_bin_num)
                df_initial = y_bin_slicing(df_initial, y_bin_num)
                particle_diameter = df_initial['radius'].iloc[0] * 2  # Assuming all particles have the same radius
                
                print(df_final['y_bin'].unique())

                for y_bin_label in df_final['y_bin'].unique():
                    print("--------------------------------")
                    print(f"\nProcessing y_bin: {y_bin_label}\n")
                    print("--------------------------------")
                    # Process final DataFrame
                    data_final_y_bin = df_final[df_final['y_bin'] == y_bin_label]
                    top_par_y_bin_final = find_top_particles(data_final_y_bin)

                    data_initial_y_bin = df_initial[df_initial['y_bin'] == y_bin_label]
                    top_par_y_bin_initial = find_top_particles(data_initial_y_bin)


                    angleRepose, SE_AoR, best_r_square, best_upper_limit, par_fit, lower_limit, optimal_model = calculate_AoR(data_initial_y_bin, data_final_y_bin)

                    # rms_height, se_rms_height, rms_slope, se_rms_slope = calculate_surface_roughness(top_par_y_bin_final, particle_diameter/2)

                    plot_filename = root_directory_plots +  f"{mixture}_{sim_type}_"  + f"field_{field}_" + f"y_bin_{y_bin_label}" + ".png"
                    
                    # plot_both_frames(data_final_y_bin, top_par_y_bin_final, plot_filename)
                    plot_both_frames_aor(data_final_y_bin, par_fit, optimal_model, plot_filename)

                    # Create a dictionary with the results for the current field
                    results = {
                        'mixture' : mixture,
                        'simtype': sim_type,
                        'field': field,
                        'y_bin': y_bin_label+1,
                        'AoR (deg)': angleRepose,
                        'SE_AoR': SE_AoR,
                        'R2': best_r_square,
                        'lower_limit': lower_limit,
                        'RMS_height (mm)': rms_height * 1000,
                        'SE_RMS_height': se_rms_height * 1000,
                        'RMS_slope (deg)': rms_slope,
                        'SE_RMS_slope': se_rms_slope
                    }
                    print("\nData in separate lines:")
                    for key, value in results.items():
                        # Check if the value is a numpy float and convert to regular float
                        if isinstance(value, np.float64) or isinstance(value, np.float32):
                            print(f"{key}: {float(value)}")
                        else:
                            print(f"{key}: {value}")
                    all_results.append(results)

            else:
                print(f"Skipping field {field} due to missing final frame data.")


    
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