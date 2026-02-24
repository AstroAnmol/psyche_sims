import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})  # Set global font size for all plots

# Read the CSV file
df = pd.read_csv('/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/results/values/aor_rms_values_final.csv')

plot_folder = "/Users/sikka-mac/Obsidian/Second_Brain/01 - Projects/01 - Research/Writing/Paper 3/Planetary_Science_Journal/figures/trends/"

# Group by mixture, simtype, and field
grouped = df.groupby(['mixture', 'simtype', 'field'])

# Calculate mean and std deviation for AoR, RMS_height, and RMS_slope
stats = grouped[['AoR (deg)', 'RMS_height (mm)', 'RMS_slope (deg)']].agg(['mean', 'std'])

# Calculate weighted mean using R2 values as weights (different weights per metric)
weights_AoR = (df['R2'])  # Assuming R2 is the measure of confidence for AoR
weights_RMS_height = 1 / (df['SE_RMS_height']**2)
weights_RMS_slope = 1 / (df['SE_RMS_slope']**2)

weighted_stats = grouped.apply(lambda x: pd.Series({
    'AoR (deg)': (x['AoR (deg)'] * weights_AoR[x.index]).sum() / weights_AoR[x.index].sum(),
    'AoR SE': (x['SE_AoR'] * weights_AoR[x.index]).sum() / weights_AoR[x.index].sum(),
    'RMS_height (mm)': (x['RMS_height (mm)'] * weights_RMS_height[x.index]).sum() / weights_RMS_height[x.index].sum(),
    'RMS_height SE': (x['SE_RMS_height'] * weights_RMS_height[x.index]).sum() / weights_RMS_height[x.index].sum(),
    'RMS_slope (deg)': (x['RMS_slope (deg)'] * weights_RMS_slope[x.index]).sum() / weights_RMS_slope[x.index].sum(),
    'RMS_slope SE': (x['SE_RMS_slope'] * weights_RMS_slope[x.index]).sum() / weights_RMS_slope[x.index].sum()
}))
stats = weighted_stats
print(stats.round(2))

# Plot Angle of Repose vs field for each simtype on a single axis
fig, ax = plt.subplots(figsize=(10, 6))
simtypes = stats.index.get_level_values('simtype').unique()

# Define labels and markers for each case
markers = {'25_high': '^', '50_high': 's', '75_high': 'v', 'high_susc': '<', 'low_susc': '>'}
colors = {'basalt': 'C0', 'olivine': 'C1', 'magnetite': 'C2'}
labels = {'25_high': 'KC 1:3', '50_high': 'KC 1:1', '75_high': 'KC 3:1', 'high_susc': 'KC', 'low_susc': 'TR'}

for simtype in simtypes:
    data = stats.xs(simtype, level='simtype')
    
    if simtype in ['no_susc']:
        continue
    
    for mixture in data.index.get_level_values('mixture').unique():
        mixture_data = data.xs(mixture, level='mixture')
        
        # Remove field 1000 for specific mixtures
        if simtype in ['50_high', '75_high', 'high_susc']:
            mixture_data = mixture_data[mixture_data.index != 1000]

        ax.errorbar(
            mixture_data.index,
            -mixture_data['AoR (deg)'],
            yerr=mixture_data['AoR SE'],
            marker=markers.get(simtype, 'o'),
            label=labels.get(simtype, simtype),
            linestyle='None',  # This removes the connecting lines
            capsize=5,
            markersize=12
        )

ax.set_xlabel(r'$B (\mu T)$')
ax.set_ylabel(r'$\bar{\theta}_{B}$ (deg)')
# ax.set_title('AoR vs Field')
ax.set_xscale('symlog')
ax.set_ylim(25, 45)
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(plot_folder + 'aor_vs_field.png', dpi=300)
# plt.show()

# Plot RMS_height vs field for each simtype on a single axis
fig, ax = plt.subplots(figsize=(10, 6))

for simtype in simtypes:
    data = stats.xs(simtype, level='simtype')
    
    if simtype in ['no_susc']:
        continue
    
    for mixture in data.index.get_level_values('mixture').unique():
        mixture_data = data.xs(mixture, level='mixture')
        
        # if simtype in ['50_high', '75_high', 'high_susc']:
        #     mixture_data = mixture_data[mixture_data.index != 1000]

        ax.errorbar(
            mixture_data.index,
            mixture_data['RMS_height (mm)'],
            yerr=mixture_data['RMS_height SE'],
            marker=markers.get(simtype, 'o'),
            label=labels.get(simtype, simtype),
            linestyle='None',
            capsize=5,
            markersize=12
        )

ax.set_xlabel(r'$B (\mu T)$')
ax.set_ylabel(r'$\bar{\xi}$ (mm)')
ax.set_xscale('symlog')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(plot_folder + 'rms_height_vs_field.png', dpi=300)
# plt.show()

# Plot RMS_slope vs field for each simtype on a single axis
fig, ax = plt.subplots(figsize=(10, 6))

for simtype in simtypes:
    data = stats.xs(simtype, level='simtype')
    
    if simtype in ['no_susc']:
        continue
    
    for mixture in data.index.get_level_values('mixture').unique():
        mixture_data = data.xs(mixture, level='mixture')
        
        # if simtype in ['50_high', '75_high', 'high_susc']:
        #     mixture_data = mixture_data[mixture_data.index != 1000]

        ax.errorbar(
            mixture_data.index,
            mixture_data['RMS_slope (deg)'],
            yerr=mixture_data['RMS_slope SE'],
            marker=markers.get(simtype, 'o'),
            label=labels.get(simtype, simtype),
            linestyle='None',
            capsize=5,
            markersize=12
        )

ax.set_xlabel(r'$B (\mu T)$')
ax.set_ylabel(r'$\bar{\theta}_{rms} (deg)$')
ax.set_xscale('symlog')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(plot_folder + 'rms_slope_vs_field.png', dpi=300)
# plt.show()

# Create a figure with 3 subplots in a single row
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Plot 1: AoR vs field
for simtype in simtypes:
    data = stats.xs(simtype, level='simtype')
    if simtype in ['no_susc']:
        continue
    for mixture in data.index.get_level_values('mixture').unique():
        mixture_data = data.xs(mixture, level='mixture')
        if simtype in ['50_high', '75_high', 'high_susc']:
            mixture_data = mixture_data[mixture_data.index != 1000]
        axes[0].errorbar(
            mixture_data.index,
            -mixture_data['AoR (deg)'],
            yerr=mixture_data['AoR SE'],
            marker=markers.get(simtype, 'o'),
            label=labels.get(simtype, simtype),
            linestyle='None',
            capsize=5,
            markersize=12
        )

axes[0].set_xlabel(r'$B (\mu T)$')
axes[0].set_ylabel(r'$\bar{\theta}_{B}$ (deg)')
axes[0].set_xscale('symlog')
axes[0].set_ylim(25, 45)
axes[0].grid(True, alpha=0.3)

# Plot 2: RMS_height vs field
for simtype in simtypes:
    data = stats.xs(simtype, level='simtype')
    if simtype in ['no_susc']:
        continue
    for mixture in data.index.get_level_values('mixture').unique():
        mixture_data = data.xs(mixture, level='mixture')
        axes[1].errorbar(
            mixture_data.index,
            mixture_data['RMS_height (mm)'],
            yerr=mixture_data['RMS_height SE'],
            marker=markers.get(simtype, 'o'),
            label=labels.get(simtype, simtype),
            linestyle='None',
            capsize=5,
            markersize=12
        )

axes[1].set_xlabel(r'$B (\mu T)$')
axes[1].set_ylabel(r'$\bar{\xi}$ (mm)')
axes[1].set_xscale('symlog')
axes[1].grid(True, alpha=0.3)

# Plot 3: RMS_slope vs field
for simtype in simtypes:
    data = stats.xs(simtype, level='simtype')
    if simtype in ['no_susc']:
        continue
    for mixture in data.index.get_level_values('mixture').unique():
        mixture_data = data.xs(mixture, level='mixture')
        axes[2].errorbar(
            mixture_data.index,
            mixture_data['RMS_slope (deg)'],
            yerr=mixture_data['RMS_slope SE'],
            marker=markers.get(simtype, 'o'),
            label=labels.get(simtype, simtype),
            linestyle='None',
            capsize=5,
            markersize=12
        )

axes[2].set_xlabel(r'$B (\mu T)$')
axes[2].set_ylabel(r'$\bar{\theta}_{rms}$ (deg)')
axes[2].set_xscale('symlog')
axes[2].grid(True, alpha=0.3)

# Create a common legend below all subplots
handles, labels_list = axes[0].get_legend_handles_labels()
fig.legend(handles, labels_list, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.05))

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(plot_folder + 'trends_combined.png', dpi=300, bbox_inches='tight')
# plt.show()

# print("\nValues for y_bin = 0:")
# # Filter for y_bin = 0
# df_y0 = df[df['y_bin'] == 0]

# # Group by mixture, simtype, and field
# grouped_y0 = df_y0.groupby(['mixture', 'simtype', 'field'])

# # Calculate weighted mean for y_bin = 0
# weights_AoR_y0 = 1 / (df_y0['SE_AoR']**2)
# weights_RMS_height_y0 = 1 / (df_y0['SE_RMS_height']**2)
# weights_RMS_slope_y0 = 1 / (df_y0['SE_RMS_slope']**2)

# weighted_stats_y0 = grouped_y0.apply(lambda x: pd.Series({
#     'AoR (deg)': (x['AoR (deg)'] * weights_AoR_y0[x.index]).sum() / weights_AoR_y0[x.index].sum(),
#     'AoR SE': x['SE_AoR'].iloc[0],
#     'RMS_height (mm)': (x['RMS_height (mm)'] * weights_RMS_height_y0[x.index]).sum() / weights_RMS_height_y0[x.index].sum(),
#     'RMS_height SE': x['SE_RMS_height'].iloc[0],
#     'RMS_slope (deg)': (x['RMS_slope (deg)'] * weights_RMS_slope_y0[x.index]).sum() / weights_RMS_slope_y0[x.index].sum(),
#     'RMS_slope SE': x['SE_RMS_slope'].iloc[0]
# }))

# print("\nValues for y_bin = 0:")
# print(weighted_stats_y0)
