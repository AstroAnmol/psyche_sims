import networkx as nx
from scipy.spatial.distance import cdist
from scipy.spatial import Voronoi
from collections import Counter

import numpy as np
import pandas as pd
import os
import sys
import tqdm

# Add parent directory to path to import avalanching_analysis
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# from avalanching_analysis.plotter import *
from avalanching_analysis.read_find_functions import *
from avalanching_analysis.analysis_functions import *
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

def find_force_chains(df_step, force_threshold_multiplier=1.5):
    """
    Identify force chains in a given timestep's DataFrame.
    Parameters:
    - df_step: DataFrame containing particle data for a single timestep.
    - force_threshold_multiplier: Multiplier for mean force to define "strong" particles.
    Returns:
    - DataFrame with an additional 'chain_label' column indicating force chain membership.
    """
    # # 1. Filter for high-force particles
    # mean_f = df_step['force_mag'].mean()
    # strong_particles = df_step[df_step['force_mag'] > (mean_f * force_threshold_multiplier)].copy()
    
    # if strong_particles.empty:
    #     return df_step
    strong_particles = df_step.copy()

    # print(f"Identified {len(strong_particles)} strong particles out of {len(df_step)} total particles.")

    # 2. Build a distance matrix
    pos = strong_particles[['x', 'y', 'z']].values
    dist_mat = cdist(pos, pos)
    
    # 3. Create Graph: Link particles within contact distance (e.g., 2 * radius)
    tol = 1e-7  # Small tolerance to account for numerical precision
    adj = dist_mat < 0.01 + tol # Adjust based on your particle diameter
    G = nx.from_numpy_array(adj)
    
    # 4. Map components back to Particle IDs
    chains = list(nx.connected_components(G))
    # Assign labels...
    strong_particles['chain_label'] = -1
    # Only assign a label if the chain has more than 1 particle
    chain_id_counter = 0
    for chain in chains:
        if len(chain) > 1:  # <--- This is the key filter
            strong_particles.iloc[list(chain), strong_particles.columns.get_loc('chain_label')] = chain_id_counter
            chain_id_counter += 1

    return strong_particles

def get_voronoi_bonds(group):
    """Returns a set of frozensets representing particle ID pairs (bonds)."""
    points = group[['x', 'y', 'z']].values
    ids = group['id'].values
    
    # Voronoi requires at least 4 points for 3D
    if len(points) < 4:
        return set()
    
    vor = Voronoi(points)
    bonds = set()
    
    # ridge_points contains indices of the 'points' array
    for p1_idx, p2_idx in vor.ridge_points:
        # Map array indices back to actual Particle IDs
        p1, p2 = ids[p1_idx], ids[p2_idx]
        bonds.add(frozenset([p1, p2]))
    
    return bonds

def get_voronoi_chain_labels(df_step):
    """
    Uses Voronoi ridge points to group particles into topological chains.
    """
    points = df_step[['x', 'y', 'z']].values
    ids = df_step['id'].values
    
    if len(points) < 4:
        df_step['voronoi_chain_id'] = -1
        return df_step

    vor = Voronoi(points)
    
    # 1. Create a Graph where nodes are the actual Particle IDs
    G = nx.Graph()
    G.add_nodes_from(ids)
    
    # 2. Add edges based on Voronoi ridges
    for p1_idx, p2_idx in vor.ridge_points:
        p1, p2 = ids[p1_idx], ids[p2_idx]
        # Only add edge if particles are within contact distance
        particle_distance = np.linalg.norm(points[p1_idx] - points[p2_idx])
        if particle_distance < 0.01 + 1e-07:  # Adjust threshold as needed
            G.add_edge(p1, p2)
    
    # 3. Find connected components (clumps/chains)
    components = list(nx.connected_components(G))
    
    # 4. Map back to the DataFrame
    # Initialize with -1
    df_step['voronoi_chain_id'] = -1
    
    chain_id_counter = 0
    for component in components:
        if len(component) > 1:  # Only label actual groups
            # We use .isin for the particle_id column to map the set to the rows
            df_step.loc[df_step['id'].isin(component), 'voronoi_chain_id'] = chain_id_counter
            chain_id_counter += 1
            
    return df_step

def get_persistent_chain_labels(df_step, current_bond_age_dict, threshold):
    """
    Creates chain IDs based only on Voronoi bonds that have persisted.
    """
    points = df_step[['x', 'y', 'z']].values
    ids = df_step['id'].values
    
    if len(points) < 4:
        df_step['persistent_chain_id'] = -1
        return df_step, {}

    vor = Voronoi(points)
    current_frame_bonds = set()
    
    # Identify all current Voronoi neighbors
    for p1_idx, p2_idx in vor.ridge_points:
        p1, p2 = ids[p1_idx], ids[p2_idx]
        current_frame_bonds.add(frozenset([p1, p2]))

    # Update bond ages
    new_age_dict = {}
    stable_edges = []
    
    for bond in current_frame_bonds:
        # Increment age if it existed before, else start at 1
        age = current_bond_age_dict.get(bond, 0) + 1
        new_age_dict[bond] = age
        
        # Only use this bond for clustering if it's "old enough"
        if age >= threshold:
            stable_edges.append(list(bond))

    # Build Graph from stable edges only
    G = nx.Graph()
    G.add_nodes_from(ids)
    G.add_edges_from(stable_edges)
    
    # Map to DataFrame
    df_step['persistent_chain_id'] = -1
    chain_id_counter = 0
    for component in nx.connected_components(G):
        if len(component) > 1:
            df_step.loc[df_step['id'].isin(component), 'persistent_chain_id'] = chain_id_counter
            chain_id_counter += 1
            
    return df_step, new_age_dict


if __name__ == '__main__':

    # mixture = "single_material"
    # sim_type = "low_susc"
    # field = "1000"
    # folder_path = "/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/" + f"{mixture}/{sim_type}/field_{field}uT/"

    # frame_i = find_initial_frame(folder_path)
    # frame_f = find_end_timestep_from_dumps(folder_path + "results")

    # if frame_i is None or frame_f is None:
    #     print("Could not determine frame range. Exiting.")
    #     sys.exit(1)
    
    # frames_list = list(range(frame_i, frame_f + 1, 1000))
    
    # # --- BOND LIFETIME TRACKER INITIALIZATION ---
    # # We use a Counter to store frozensets of (id1, id2)
    # bond_persistence_counter = Counter()
    # total_analyzed_frames = 0

    # print("Processing frames and tracking bond lifetimes...")
    # for frame in tqdm.tqdm(frames_list):
    #     df_step = read_single_dump_file_csv(folder_path + f"results/dump{frame}.post")
    #     if df_step.empty:
    #         continue
        
    #     # # Calculate force magnitude
    #     # df_step['force_mag'] = np.sqrt(df_step['fx']**2 + df_step['fy']**2 + df_step['fz']**2)
        
    #     # # 1. Standard Force Chain Analysis (Existing)
    #     # df_chains = find_force_chains(df_step)
        
    #     # 2. Voronoi Bond Tracking (New)
    #     current_bonds = get_voronoi_bonds(df_step)
    #     bond_persistence_counter.update(current_bonds)
    #     total_analyzed_frames += 1

    #     # ... [Optional: keep your file saving logic here] ...

    # # --- POST-PROCESSING: CALCULATE LIFETIMES ---
    # print("\nCalculating lifetime statistics...")
    
    # # Convert Counter to a DataFrame
    # lifetime_data = []
    # for bond, count in bond_persistence_counter.items():
    #     p_ids = list(bond)
    #     lifetime_data.append({
    #         'p1': p_ids[0],
    #         'p2': p_ids[1],
    #         'lifetime_count': count,
    #         'persistence_ratio': count / total_analyzed_frames
    #     })

    # df_lifetimes = pd.DataFrame(lifetime_data)
    
    # # Filter for bonds that are "High Persistence" (e.g., present in > 80% of frames)
    # strong_bonds = df_lifetimes[df_lifetimes['persistence_ratio'] > 0.8].sort_values(by='lifetime_count', ascending=False)

    # print(f"Analysis Complete.")
    # print(f"Total Unique Bonds Detected: {len(df_lifetimes)}")
    # print(f"Strong/Persistent Bonds (>80% of time): {len(strong_bonds)}")
    
    # # Save lifetime results
    # df_lifetimes.to_csv(folder_path + "bond_lifetimes.csv", index=False)
    # print(f"Lifetime data saved to {folder_path}bond_lifetimes.csv")

    mixture = "single_material"
    sim_type = "low_susc"
    field = "1000"

    folder_path = "/Users/sikka-mac/Research/Code/SIMS-LIGGGHTS/psyche_sims/" + f"{mixture}/{sim_type}/field_{field}uT/"


    frame_i = find_initial_frame(folder_path)
    frame_f = find_end_timestep_from_dumps(folder_path + "results")


    # --- 1. Global tracker (Put this before your frame loop) ---
    bond_age_dict = {}  # Keys: frozenset(ID1, ID2), Values: number of consecutive frames
    STABILITY_THRESHOLD = 10  # Only "chains" if bond lived for + frames

    if frame_i is None or frame_f is None:
        print("Could not determine frame range. Exiting.")
        sys.exit(1)
    else:
        # # Analyze force chains for a series of frames between frame_i and frame_f
        # for frame in range(frame_i, frame_f + 1, 10):  #
        #     df_step = read_single_dump_file_csv(folder_path + f"results/dump{frame}.post")
        #     if df_step.empty:
        #         continue
        #     df_chains = find_force_chains(df_step)
        #     print(f"Frame {frame}: Found {df_chains['chain_label'].nunique()} force chains.")

        # Create animation of force chains
        # import matplotlib.pyplot as plt
        
        frames_data = []
        frames_list = list(range(frame_i, frame_f + 1, 1000))

        bond_age_tracker = {} # Initialize empty
        
        for frame in tqdm.tqdm(frames_list):
            df_step = read_single_dump_file_csv(folder_path + f"results/dump{frame}.post")
            if df_step.empty:
                continue
            df_step['force_mag'] = np.sqrt(df_step['fx']**2 + df_step['fy']**2 + df_step['fz']**2)
            df_chains = find_force_chains(df_step)
            frames_data.append((frame, df_chains))


            df_chains = get_voronoi_chain_labels(df_chains)
            

            # # Get labels based on persistence
            # df_chains, bond_age_tracker = get_persistent_chain_labels(
            #     df_chains, 
            #     bond_age_tracker, 
            #     threshold=STABILITY_THRESHOLD
            # )
    

            save_path = folder_path + f"results_chains/force_chains_frame_{frame}.post"
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            # df_chains.to_csv(save_path, sep=' ', index=False)
            # Write LAMMPS-style header
            with open(save_path, 'w') as f:
                f.write(f"ITEM: TIMESTEP\n")
                f.write(f"{frame}\n")
                f.write(f"ITEM: NUMBER OF ATOMS\n")
                f.write(f"{len(df_chains)}\n")
                f.write(f"ITEM: BOX BOUNDS ff ff fs\n")
                f.write(f"0 0.3\n0 0.1\n0 0.108043\n")
                f.write(f"ITEM: ATOMS {' '.join(df_chains.columns[df_chains.notna().any()])}\n")
                df_chains.to_csv(f, sep=' ', index=False, header=False)
        
        # fig = plt.figure(figsize=(10, 8))
        # ax = fig.add_subplot(111, projection='3d')
        
        # def update(idx):
        #     ax.clear()
        #     frame, df_chains = frames_data[idx]
            
        #     scatter = ax.scatter(df_chains['x'], df_chains['y'], df_chains['z'], 
        #            c=df_chains['chain_label'], cmap='tab20', s=df_chains['radius']*10000)
        #     ax.set_xlabel('X')
        #     ax.set_ylabel('Y')
        #     ax.set_zlabel('Z')
        #     ax.set_title(f'Force Chains - Frame {frame}')
        #     return scatter,
        
        # anim = FuncAnimation(fig, update, frames=len(frames_data), interval=100, blit=False)
        # frame, df_chains = frames_data[10]
        # scatter = ax.scatter(df_chains['x'], df_chains['y'], df_chains['z'], 
        #            c=df_chains['chain_label'], cmap='tab20', s=df_chains['radius']*10000)
        # ax.set_xlabel('X')
        # ax.set_ylabel('Y')
        # ax.set_zlabel('Z')
        # ax.set_title(f'Force Chains - Frame {frame}')
        # anim.save(folder_path + 'force_chains_animation.mp4')
        # plt.show()
        # print("Animation saved to force_chains_animation.mp4")