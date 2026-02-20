import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})  # Set global font size for all plots

# Read the CSV file
df = pd.read_csv('/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/values/aor_rms_values.csv')

# Group by mixture, simtype, and field
grouped = df.groupby(['mixture', 'simtype', 'field'])

# Calculate mean and std deviation for AoR, RMS_height, and RMS_slope
stats = grouped[['AoR (deg)', 'RMS_height (mm)', 'RMS_slope (deg)']].agg(['mean', 'std'])
# Calculate weighted mean using SE values as weights (different weights per metric)
weights_AoR = 1 / (df['SE_AoR']**2)
weights_RMS = 1 / (df['SE_RMS_height']**2 + df['SE_RMS_slope']**2)

weighted_stats = grouped.apply(lambda x: pd.Series({
    'AoR (deg)': (x['AoR (deg)'] * weights_AoR[x.index]).sum() / weights_AoR[x.index].sum(),
    'RMS_height (mm)': (x['RMS_height (mm)'] * weights_RMS[x.index]).sum() / weights_RMS[x.index].sum(),
    'RMS_slope (deg)': (x['RMS_slope (deg)'] * weights_RMS[x.index]).sum() / weights_RMS[x.index].sum(),
}))
stats = weighted_stats
print(stats)
# Filter for y_bin = 0
df_y0 = df[df['y_bin'] == 0]

# Group by mixture, simtype, and field
grouped_y0 = df_y0.groupby(['mixture', 'simtype', 'field'])

# Calculate weighted mean for y_bin = 0
weights_AoR_y0 = 1 / (df_y0['SE_AoR']**2)
weights_RMS_y0 = 1 / (df_y0['SE_RMS_height']**2 + df_y0['SE_RMS_slope']**2)

weighted_stats_y0 = grouped_y0.apply(lambda x: pd.Series({
    'AoR (deg)': (x['AoR (deg)'] * weights_AoR_y0[x.index]).sum() / weights_AoR_y0[x.index].sum(),
    'RMS_height (mm)': (x['RMS_height (mm)'] * weights_RMS_y0[x.index]).sum() / weights_RMS_y0[x.index].sum(),
    'RMS_slope (deg)': (x['RMS_slope (deg)'] * weights_RMS_y0[x.index]).sum() / weights_RMS_y0[x.index].sum(),
}))

print("\nValues for y_bin = 0:")
print(weighted_stats_y0)
