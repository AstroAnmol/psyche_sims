import matplotlib.pyplot as plt
import pandas as pd

# Read the CSV file containing AoR and RMS values over time
results_df = pd.read_csv("/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/values/aor_rms_values_overtime_mixture.csv")

sim_types = results_df['simtype'].unique()
fields = results_df['field'].unique()
mixture = results_df['mixture'].unique()[0]  # Assuming only one mixture type in the data

# Plot frame vs AoR
for sim_type in sim_types:
    plt.figure(figsize=(10, 6))

    for field in fields:
        data = results_df[(results_df['simtype'] == sim_type) & (results_df['field'] == field)].copy()
        first_frame = data['frame'].iloc[0]
        data['frame_scaled'] = data['frame'] - first_frame
        plt.plot(data['frame_scaled'], -data['AoR (deg)'], label=f'{sim_type} - {field}')

    plt.xlabel('Frame')
    plt.ylabel('Angle of Repose (deg)')
    plt.title(f'Angle of Repose Over Time - {mixture}')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/plots/AoR_over_time_{mixture}_{sim_type}.png")
    plt.show()

# Read the CSV file containing AoR and RMS values over time
results_df = pd.read_csv("/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/values/aor_rms_values_overtime.csv")

sim_types = results_df['simtype'].unique()
fields = results_df['field'].unique()
mixture = results_df['mixture'].unique()[0]  # Assuming only one mixture type in the data

# Plot frame vs AoR
for sim_type in sim_types:
    plt.figure(figsize=(10, 6))

    for field in fields:
        data = results_df[(results_df['simtype'] == sim_type) & (results_df['field'] == field)].copy()
        first_frame = data['frame'].iloc[0]
        data['frame_scaled'] = data['frame'] - first_frame
        plt.plot(data['frame_scaled'], -data['AoR (deg)'], label=f'{sim_type} - {field}')

    plt.xlabel('Frame')
    plt.ylabel('Angle of Repose (deg)')
    plt.title(f'Angle of Repose Over Time - {mixture}')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/plots/AoR_over_time_{mixture}_{sim_type}.png")
    plt.show()