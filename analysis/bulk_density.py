import numpy as np
import pandas as pd



def read_single_dump_file_csv(file_path):
    """Reads a single LIGGGHTS dump file into a pandas DataFrame using pd.read_csv."""
    # def check_read_permission(file_path):
    #     """Checks if the given file path has read permission for the current user."""
    #     if os.access(file_path, os.R_OK):
    #         print(f"'{file_path}' has read permission.")
    #         return True
    #     else:
    #         print(f"'{file_path}' does NOT have read permission.")
    #         return False
    # check_read_permission(file_path)
    try:    # id type x y z vx vy vz radius fx fy fz magfx magfy magfz mux muy muz mass 
        column_names = ['id', 'type', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'radius', 'fx', 'fy', 'fz', "mass", 'NaN'] # Add or adjust columns as needed

        # Read the data part of the file
        df = pd.read_csv(file_path, skiprows=9, sep=' ', header=None, names=column_names)
        return df
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None

def total_particles(df):
    """Returns the total number of particles in the DataFrame without floor particles, 
    also remove particles above the fill-height of 100 mm."""
    fill_height = 100e-3  # 100 mm in meters
    # find the lowest particle height
    lowest_height = df['z'].min()
    # define the maximum height for particles to be considered
    max_height = lowest_height + fill_height
    # filter the DataFrame to include only particles below the max height
    df = df[df['z'] <= max_height]
    #find the number of unique types
    unique_types = df['type'].nunique()
    if unique_types > 1:
        floor_type = df['type'].max()  # Assuming the floor particles have the highest type number
        df = df[df['type'] != floor_type]
    return len(df)

def percentage_of_each_type(df):
    """Returns a dictionary with the count of each particle type in the DataFrame 
    without floor particles, also remove particles above the fill-height of 100 mm."""
    fill_height = 100e-3  # 100 mm in meters
    # find the lowest particle height
    lowest_height = df['z'].min()
    # define the maximum height for particles to be considered
    max_height = lowest_height + fill_height
    # filter the DataFrame to include only particles below the max height
    df = df[df['z'] <= max_height]
    unique_types = df['type'].nunique()
    if unique_types > 1:
        floor_type = df['type'].max()  # Assuming the floor particles have the highest type number
        df = df[df['type'] != floor_type]
    total_particles = len(df)
    type_counts = df['type'].value_counts().to_dict()
    for t in type_counts:
        type_counts[t] = (type_counts[t] / total_particles) * 100  # Convert to percentage
    return type_counts

def total_mass(df):
    """Calculates the total mass from the DataFrame removing the floor particles and particles above the fill-height of 100 mm."""
    fill_height = 100e-3  # 100 mm in meters
    # find the lowest particle height
    lowest_height = df['z'].min()
    # define the maximum height for particles to be considered
    max_height = lowest_height + fill_height
    # filter the DataFrame to include only particles below the max height
    df = df[df['z'] <= max_height]
    unique_types = df['type'].nunique()
    # print(f"Number of unique types: {unique_types}")
    if unique_types > 1:
        floor_type = df['type'].max()  # Assuming the floor particles have the highest type number
        df = df[df['type'] != floor_type]
    # Floor particles already removed above
    if df is not None:
        return df['mass'].sum()

def total_par_volume(df):
    """Calculates the total particle volume from the DataFrame removing the floor particles
    and particles above the fill-height of 100 mm."""
    fill_height = 100e-3  # 100 mm in meters
    # find the lowest particle height
    lowest_height = df['z'].min()
    # define the maximum height for particles to be considered
    max_height = lowest_height + fill_height
    # filter the DataFrame to include only particles below the max height
    df = df[df['z'] <= max_height]
    unique_types = df['type'].nunique()
    # print(f"Number of unique types: {unique_types}")
    if unique_types > 1:
        floor_type = df['type'].max()  # Assuming the floor particles have the highest type number
        df = df[df['type'] != floor_type]
    # Floor particles already removed above
    if df is not None:
        radii = df['radius']
        volumes = (4/3) * np.pi * (radii ** 3)
        return volumes.sum()

def avg_par_density(df):
    """Calculates the bulk density from the DataFrame removing the floor particles 
    and particles above the fill-height of 100 mm."""
    fill_height = 100e-3  # 100 mm in meters
    # find the lowest particle height
    lowest_height = df['z'].min()
    # define the maximum height for particles to be considered
    max_height = lowest_height + fill_height
    # filter the DataFrame to include only particles below the max height
    df = df[df['z'] <= max_height]
    particle_mass = total_mass(df)
    particle_volume = total_par_volume(df)
    if particle_volume > 0:
        return particle_mass / particle_volume
    else:
        return 0
    
def bulk_density(df, box_volume):
    """Calculates the bulk density from the DataFrame removing the floor particles 
    and particles above the fill-height of 100 mm."""
    fill_height = 100e-3  # 100 mm in meters
    # find the lowest particle height
    lowest_height = df['z'].min()
    # define the maximum height for particles to be considered
    max_height = lowest_height + fill_height
    # filter the DataFrame to include only particles below the max height
    df = df[df['z'] <= max_height]
    particle_mass = total_mass(df)
    if box_volume > 0:
        return particle_mass / box_volume
    else:
        return 0

def porosity(df, box_volume):
    """Calculates the porosity from the DataFrame removing the floor particles 
    and particles above the fill-height of 100 mm."""
    fill_height = 100e-3  # 100 mm in meters
    # find the lowest particle height
    lowest_height = df['z'].min()
    # define the maximum height for particles to be considered
    max_height = lowest_height + fill_height
    # filter the DataFrame to include only particles below the max height
    df = df[df['z'] <= max_height]
    particle_volume = total_par_volume(df)
    if box_volume > 0:
        return 1 - (particle_volume / box_volume)
    else:
        return 0

if __name__ == "__main__":
    # Example usage
    folder_paths = ["mixture/25_high/", "mixture/50_high/", "mixture/75_high/",
                     "single_material/low_susc/", "single_material/no_susc/"]
    frames = [750000, 5000000]
    for path in folder_paths:
        print(f"Simulation type: {path}")
        for frame in frames:
            print(f"Data at frame: {frame}")
            file_name = f"setup/results/dump{frame}.post"
            file_path = path + file_name
            df = read_single_dump_file_csv(file_path)
            
            box_volume = 100 * 100 * 150 * 1e-9 # in cubic meters'
            print("Total Particles:", total_particles(df))
            print("Percentage of Each Type:", percentage_of_each_type(df))
            print("Average Particle Density (kg/m^3):", avg_par_density(df))
            print("Bulk Density (kg/m^3):", bulk_density(df, box_volume))
            print("Porosity:", porosity(df, box_volume))
            print("-" * 40 )


