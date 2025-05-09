import numpy as np
import matplotlib.pyplot as plt # type: ignore
from mpl_toolkits.mplot3d import Axes3D # type: ignore
from matplotlib.patches import FancyArrowPatch # type: ignore
from mpl_toolkits.mplot3d import proj3d # type: ignore
from scipy.spatial import ConvexHull # type: ignore
import matplotlib.colors as colors # type: ignore

# Constants
G = 9.81  # acceleration due to gravity in m/s²
MACH_3 = 1029/3  # Mach 3 in m/s at sea level
MAX_ACCELERATION = 3 * G  # 3G acceleration
M_TO_MI = 0.000621371  # Conversion from meters to miles
M_TO_KM = 0.001  # Conversion from meters to kilometers

class Arrow3D(FancyArrowPatch):
    """Custom 3D arrow for visualization"""
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)

def calculate_nominal_trajectory(initial_position, initial_velocity, times):
    """Calculate the nominal trajectory without error or maneuvering"""
    trajectories = []
    
    for t in times:
        position = initial_position + initial_velocity * t
        trajectories.append(position)
        
    return np.array(trajectories)

def generate_error_sphere(center, radius, resolution=10000):
    """Generate a sphere mesh for visualization"""
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    
    return x, y, z

def generate_initial_error_samples(initial_position, initial_velocity, error_percentage, num_samples=10000):
    """Generate samples within the initial error bounds for position and velocity"""
    pos_error = initial_position * error_percentage
    vel_error = initial_velocity * error_percentage
    
    # For more uniform distribution of points in 3D space, use spherical sampling
    # Generate random directions for position error
    phi_pos = np.random.uniform(0, 2*np.pi, num_samples)
    theta_pos = np.arccos(2*np.random.uniform(0, 1, num_samples) - 1)  # Uniform on sphere
    r_pos = np.random.uniform(0, 1, num_samples)**(1/3)  # Uniform in volume, cube root for proper density
    
    # Scale by the maximum error and create position samples
    max_pos_error = np.linalg.norm(pos_error)
    x_pos = r_pos * np.sin(theta_pos) * np.cos(phi_pos) * max_pos_error
    y_pos = r_pos * np.sin(theta_pos) * np.sin(phi_pos) * max_pos_error
    z_pos = r_pos * np.cos(theta_pos) * max_pos_error
    
    pos_samples = np.column_stack((
        initial_position[0] + x_pos,
        initial_position[1] + y_pos,
        initial_position[2] + z_pos
    ))
    
    # For each position sample, create several velocity samples
    all_pos_samples = []
    all_vel_samples = []
    
    # Number of velocity samples per position sample
    vel_samples_per_pos = 5
    
    for i in range(num_samples):
        pos = pos_samples[i]
        
        # Generate multiple velocity samples for this position
        for _ in range(vel_samples_per_pos):
            # Generate random direction for velocity error
            phi_vel = np.random.uniform(0, 2*np.pi)
            theta_vel = np.arccos(2*np.random.uniform(0, 1) - 1)
            r_vel = np.random.uniform(0, 1)**(1/3)  # Cube root for uniform density
            
            # Scale by the maximum velocity error
            max_vel_error = np.linalg.norm(vel_error)
            x_vel = r_vel * np.sin(theta_vel) * np.cos(phi_vel) * max_vel_error
            y_vel = r_vel * np.sin(theta_vel) * np.sin(phi_vel) * max_vel_error
            z_vel = r_vel * np.cos(theta_vel) * max_vel_error
            
            vel = np.array([
                initial_velocity[0] + x_vel,
                initial_velocity[1] + y_vel,
                initial_velocity[2] + z_vel
            ])
            
            all_pos_samples.append(pos)
            all_vel_samples.append(vel)
    
    return np.array(all_pos_samples), np.array(all_vel_samples)

def calculate_no_maneuver_error(initial_position, initial_velocity, error_percentage, t_lock, num_samples=1000):
    """Calculate error region without maneuvering"""
    pos_samples, vel_samples = generate_initial_error_samples(
        initial_position, initial_velocity, error_percentage, num_samples
    )
    
    # Propagate each sample to t_lock without maneuvering
    final_positions = []
    for i in range(num_samples):
        final_pos = pos_samples[i] + vel_samples[i] * t_lock
        final_positions.append(final_pos)
    
    return np.array(final_positions)

def generate_random_unit_vectors(num_samples):
    """Generate random unit vectors uniformly distributed on a sphere"""
    # Generate random directions
    vectors = np.random.randn(num_samples, 3)
    
    # Normalize to unit vectors
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    unit_vectors = vectors / norms
    
    return unit_vectors

def calculate_maneuver_error(initial_position, initial_velocity, error_percentage, max_accel_g, t_lock, 
                            num_initial_samples=20000, num_accel_directions=1000):
    """Calculate error region with maximum maneuvering in all directions"""
    # Get initial position and velocity samples within error bounds
    # With the improved sampling, we get more combinations of initial states
    pos_samples, vel_samples = generate_initial_error_samples(
        initial_position, initial_velocity, error_percentage, num_initial_samples
    )
    
    # For each initial state, calculate possible positions with different acceleration directions
    final_positions = []
    final_velocities = []
    
    # Process each initial state (position + velocity combo)
    total_samples = len(pos_samples)
    
    for i in range(total_samples):
        # Generate random acceleration directions for this initial state
        accel_directions = generate_random_unit_vectors(num_accel_directions)
        
        for direction in accel_directions:
            # Apply maximum acceleration in this direction
            accel_vector = direction * max_accel_g * G
            
            # Propagate with constant acceleration
            pos = pos_samples[i].copy()
            vel = vel_samples[i].copy()
            
            # Use smaller time steps for integration
            dt = 0.1  # seconds
            num_steps = int(t_lock / dt)
            
            # Initialize maneuvering path for this sample
            path = [pos.copy()]
            
            for _ in range(num_steps):
                # Update velocity with acceleration
                vel_new = vel + accel_vector * dt
                
                # Check if velocity exceeds Mach 3
                vel_magnitude = np.linalg.norm(vel_new)
                if vel_magnitude > MACH_3:
                    # Scale velocity back to Mach 3
                    vel_new = vel_new * (MACH_3 / vel_magnitude)
                
                # Update position with new velocity (trapezoidal integration)
                pos = pos + 0.5 * (vel + vel_new) * dt
                vel = vel_new
                
                # Append to path
                path.append(pos.copy())
            
            final_positions.append(pos)
            final_velocities.append(vel)
    
    return np.array(final_positions), np.array(final_velocities)

def calculate_error_summary(no_maneuver_positions, maneuver_positions, maneuver_velocities, initial_position):
    """Calculate numerical summary of the error growth"""
    # The first position in no_maneuver_positions should be the nominal trajectory endpoint
    nominal_endpoint = no_maneuver_positions[0]
    
    # Calculate distance from nominal endpoint to each point in the no-maneuver volume
    no_maneuver_distances_from_nominal = np.linalg.norm(no_maneuver_positions - nominal_endpoint, axis=1)
    
    # Calculate distance from nominal endpoint to each point in the maneuver volume
    maneuver_distances_from_nominal = np.linalg.norm(maneuver_positions - nominal_endpoint, axis=1)
    
    # Calculate velocity magnitudes for maneuver case
    maneuver_velocity_magnitudes = np.linalg.norm(maneuver_velocities, axis=1)
    
    # Calculate nominal trajectory error (distance from initial position)
    nominal_trajectory_error = np.linalg.norm(nominal_endpoint - initial_position)
    
    summary = {
        "nominal_trajectory_error": nominal_trajectory_error,
        "max_distance_from_nominal_no_maneuver": np.max(no_maneuver_distances_from_nominal),
        "mean_distance_from_nominal_no_maneuver": np.mean(no_maneuver_distances_from_nominal),
        "max_distance_from_nominal_with_maneuver": np.max(maneuver_distances_from_nominal),
        "mean_distance_from_nominal_with_maneuver": np.mean(maneuver_distances_from_nominal),
        "max_velocity_with_maneuver": np.max(maneuver_velocity_magnitudes),
        "mean_velocity_with_maneuver": np.mean(maneuver_velocity_magnitudes)
    }
    
    return summary

def plot_convex_hull(ax, points, alpha=0.2, color='b', label=None, subsample=True):
    """Plot the convex hull of a set of points"""
    if len(points) < 4:
        # Not enough points for convex hull
        return
    
    try:
        # For very large point sets, subsample to improve performance
        if subsample and len(points) > 100000:
            # Take a random subsample of points for the convex hull calculation
            # This helps with performance while maintaining the approximate shape
            indices = np.random.choice(len(points), min(100000, len(points)), replace=False)
            hull_points = points[indices]
        else:
            hull_points = points
        
        hull = ConvexHull(hull_points)
        
        # Get simplices (triangles) from the convex hull
        simplices = hull.simplices
        
        # Plot each face of the convex hull
        for simplex in simplices:
            x = hull_points[simplex, 0]
            y = hull_points[simplex, 1]
            z = hull_points[simplex, 2]
            ax.plot_trisurf(x, y, z, color=color, alpha=alpha)
        
        # Also plot a scatter of a subset of points inside the hull
        # if len(points) > 500:
        #     # Only plot a subset of points for performance
        #     scatter_indices = np.random.choice(len(points), 500, replace=False)
        #     scatter_points = points[scatter_indices]
        #     ax.scatter(scatter_points[:, 0], scatter_points[:, 1], scatter_points[:, 2], 
        #               color=color, alpha=0.3, s=2)
        
        # Add a proxy artist for the legend
        if label:
            ax.plot([0], [0], [0], color=color, label=label)
    
    except Exception as e:
        print(f"Could not create convex hull: {e}")
        # Fallback: just plot the points as a scatter
        if len(points) > 5000:
            # Subsample for performance
            indices = np.random.choice(len(points), 5000, replace=False)
            scatter_points = points[indices]
        else:
            scatter_points = points
            
        ax.scatter(scatter_points[:, 0], scatter_points[:, 1], scatter_points[:, 2], 
                  color=color, alpha=0.5, s=5, label=label)

def plot_3d_visualization(initial_position, initial_velocity, nominal_trajectory=None, 
                         no_maneuver_positions=None, maneuver_positions=None, maneuver_velocities=None,
                         error_percentage=0.05, t_lock=3.0, existing_ax=None, reuse_summary=None):
    """Create 3D visualization of all components"""
    if existing_ax is None:
        # Create a new figure
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
    else:
        # Use existing axis
        ax = existing_ax
        fig = ax.figure
    
    # If we need to recalculate components
    if nominal_trajectory is None:
        times = np.linspace(0, t_lock, 100)
        nominal_trajectory = calculate_nominal_trajectory(initial_position, initial_velocity, times)
    
    # 1. Plot initial position - convert to miles for display
    ax.scatter(initial_position[0] * M_TO_MI, 
              initial_position[1] * M_TO_MI, 
              initial_position[2] * M_TO_MI, 
              color='green', s=100, label='Initial Position Estimate')
    
    # 2. Plot initial error sphere - convert to miles
    initial_error_radius = np.linalg.norm(initial_position * error_percentage)
    x_sphere, y_sphere, z_sphere = generate_error_sphere(initial_position, initial_error_radius)
    ax.plot_surface(x_sphere * M_TO_MI, 
                   y_sphere * M_TO_MI, 
                   z_sphere * M_TO_MI, 
                   color='green', alpha=0.2, 
                   label='Initial Position Error (±5%)')
    
    # 3. Plot nominal trajectory - convert to miles
    ax.plot(nominal_trajectory[:, 0] * M_TO_MI, 
           nominal_trajectory[:, 1] * M_TO_MI, 
           nominal_trajectory[:, 2] * M_TO_MI, 
           'r--', linewidth=2, label='Nominal Trajectory')
    
    # 4 & 5. Plot error regions if provided - convert to miles
    if no_maneuver_positions is not None:
        # Convert to miles for display
        no_maneuver_miles = no_maneuver_positions * M_TO_MI
        plot_convex_hull(ax, no_maneuver_miles, color='blue', alpha=0.3, 
                        label='Error Region without Maneuvering')
    
    if maneuver_positions is not None and maneuver_velocities is not None:
        # Convert to miles for display
        maneuver_miles = maneuver_positions * M_TO_MI
        plot_convex_hull(ax, maneuver_miles, color='red', alpha=0.2, 
                        label='Error Region with 3G Maneuvering')
    
    # Add velocity arrow for initial velocity - adjust for miles
    vel_magnitude = np.linalg.norm(initial_velocity)
    normalized_vel = initial_velocity / vel_magnitude
    arrow_length = vel_magnitude * 0.8 * M_TO_MI  # Adjust and convert scale factor for miles
    
    arrow = Arrow3D(
        [initial_position[0] * M_TO_MI, 
         initial_position[0] * M_TO_MI + arrow_length * normalized_vel[0]],
        [initial_position[1] * M_TO_MI, 
         initial_position[1] * M_TO_MI + arrow_length * normalized_vel[1]],
        [initial_position[2] * M_TO_MI, 
         initial_position[2] * M_TO_MI + arrow_length * normalized_vel[2]],
        mutation_scale=20, lw=3, arrowstyle='-|>', color='k'
    )
    ax.add_artist(arrow)
    
    # Set equal aspect ratio
    max_range = np.array([
        ax.get_xlim()[1] - ax.get_xlim()[0],
        ax.get_ylim()[1] - ax.get_ylim()[0],
        ax.get_zlim()[1] - ax.get_zlim()[0]
    ]).max() / 2.0
    
    mid_x = (ax.get_xlim()[1] + ax.get_xlim()[0]) / 2
    mid_y = (ax.get_ylim()[1] + ax.get_ylim()[0]) / 2
    mid_z = (ax.get_zlim()[1] + ax.get_zlim()[0]) / 2
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # Labels and title - update with miles units
    ax.set_xlabel('X (mi)')
    ax.set_ylabel('Y (mi)')
    ax.set_zlabel('Z (mi)')
    ax.set_title(f'Target Error Growth (t_lock = {t_lock:.1f} seconds)')
    
    # Create custom legend
    handles, labels = ax.get_legend_handles_labels()
    if existing_ax is None:
        fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.95, 0.95))
    else:
        ax.legend(loc='upper right')
    
    # Add summary text if we have the regions to calculate it
    if reuse_summary is not None:
        summary = reuse_summary
    elif no_maneuver_positions is not None and maneuver_positions is not None and maneuver_velocities is not None:
        summary = calculate_error_summary(no_maneuver_positions, maneuver_positions, maneuver_velocities, initial_position)
    else:
        summary = None
    
    if summary is not None and existing_ax is None:
        # Convert measurements to miles and km for display
        nominal_mi = summary['nominal_trajectory_error'] * M_TO_MI
        nominal_km = summary['nominal_trajectory_error'] * M_TO_KM
        
        max_no_maneuver_mi = summary['max_distance_from_nominal_no_maneuver'] * M_TO_MI
        max_no_maneuver_km = summary['max_distance_from_nominal_no_maneuver'] * M_TO_KM
        
        max_with_maneuver_mi = summary['max_distance_from_nominal_with_maneuver'] * M_TO_MI
        max_with_maneuver_km = summary['max_distance_from_nominal_with_maneuver'] * M_TO_KM
        
        summary_text = (
            f"Range on nominal trajectory: {nominal_mi:.2f} mi ({nominal_km:.2f} km)\n"
            f"Max error without maneuvering: {max_no_maneuver_mi:.2f} mi ({max_no_maneuver_km:.2f} km)\n"
            f"Max error with 3G maneuvering: {max_with_maneuver_mi:.2f} mi ({max_with_maneuver_km:.2f} km)\n"
            f"Max velocity with maneuvering: {summary['max_velocity_with_maneuver']:.2f} m/s"
        )
        
        plt.figtext(0.05, 0.05, summary_text, fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
    
    # Adjust the view angle for better visualization
    ax.view_init(elev=30, azim=45)
    
    if existing_ax is None:
        plt.tight_layout()
        return fig
    else:
        return None

def visualize_error_growth(initial_position, initial_velocity, t_lock, error_percentage=0.05):
    """Main function to visualize error growth"""
    # Convert inputs to numpy arrays
    initial_position = np.array(initial_position)
    initial_velocity = np.array(initial_velocity)
    
    # Calculate and print initial position error radius in miles
    initial_error_radius = np.linalg.norm(initial_position * error_percentage)
    initial_error_radius_mi = initial_error_radius * M_TO_MI
    initial_error_radius_km = initial_error_radius * M_TO_KM
    print(f"Initial position error radius: {initial_error_radius_mi:.2f} mi ({initial_error_radius_km:.2f} km)")
    
    print("Calculating nominal trajectory...")
    # 1. Calculate nominal trajectory
    times = np.linspace(0, t_lock, 100)
    nominal_trajectory = calculate_nominal_trajectory(initial_position, initial_velocity, times)
    
    print("Calculating error region without maneuvering...")
    # 2. Calculate error region without maneuvering
    # Increase number of samples for better coverage
    no_maneuver_positions = calculate_no_maneuver_error(
        initial_position, initial_velocity, error_percentage, t_lock, num_samples=800
    )
    
    print("Calculating error region with maneuvering (this may take a moment)...")
    # 3. Calculate error region with maneuvering
    # This is the most computationally intensive part
    maneuver_positions, maneuver_velocities = calculate_maneuver_error(
        initial_position, initial_velocity, error_percentage, 3.0, t_lock, 
        num_initial_samples=500, num_accel_directions=50
    )
    
    print("Creating visualization...")
    # 4. Visualize
    fig = plot_3d_visualization(
        initial_position, initial_velocity, nominal_trajectory,
        no_maneuver_positions, maneuver_positions, maneuver_velocities,
        error_percentage, t_lock
    )
    
    # 5. Calculate summary statistics
    summary = calculate_error_summary(
        no_maneuver_positions, maneuver_positions, maneuver_velocities, initial_position
    )
    
    return fig, summary

def run_target_simulation():
    """Run the target simulation with user inputs"""
    print("Target Error Visualization Tool")
    print("-------------------------------------------")
    
    # Get user inputs
    px = float(input("Enter initial X position (m): "))
    py = float(input("Enter initial Y position (m): "))
    pz = float(input("Enter initial Z position (m): "))
    
    vx = float(input("Enter initial X velocity (m/s): "))
    vy = float(input("Enter initial Y velocity (m/s): "))
    vz = float(input("Enter initial Z velocity (m/s): "))
    
    t_lock = float(input("Enter lock-on time (seconds): "))
    
    initial_position = np.array([px, py, pz])
    initial_velocity = np.array([vx, vy, vz])
    
    # Run visualization
    fig, summary = visualize_error_growth(initial_position, initial_velocity, t_lock)
    
    # Print summary
    print("\nNumerical Summary:")
    print(f"Position error on nominal trajectory: {summary['nominal_trajectory_error'] * M_TO_MI:.2f} mi ({summary['nominal_trajectory_error'] * M_TO_KM:.2f} km)")
    print(f"Maximum position error without maneuvering: {summary['max_position_error_no_maneuver'] * M_TO_MI:.2f} mi ({summary['max_position_error_no_maneuver'] * M_TO_KM:.2f} km)")
    print(f"Mean position error without maneuvering: {summary['mean_position_error_no_maneuver'] * M_TO_MI:.2f} mi ({summary['mean_position_error_no_maneuver'] * M_TO_KM:.2f} km)")
    print(f"Maximum position error with 3G maneuvering: {summary['max_position_error_with_maneuver'] * M_TO_MI:.2f} mi ({summary['max_position_error_with_maneuver'] * M_TO_KM:.2f} km)")
    print(f"Mean position error with 3G maneuvering: {summary['mean_distance_from_nominal_with_maneuver'] * M_TO_MI:.2f} mi ({summary['mean_distance_from_nominal_with_maneuver'] * M_TO_KM:.2f} km)")
    print(f"Maximum velocity with maneuvering: {summary['max_velocity_with_maneuver']:.2f} m/s")
    print(f"Mean velocity with maneuvering: {summary['mean_velocity_with_maneuver']:.2f} m/s")
    
    plt.show()

# Example usage
if __name__ == "__main__":
    # Example parameters (can be replaced with user inputs)
    initial_position = np.array([0, 0, 10000])  # 10km altitude
    initial_velocity = np.array([500, 0, -100])  # Initial velocity
    t_lock = 5.0  # 5 seconds lock-on time
    
    # For command-line interface, uncomment the following line:
    # run_target_simulation()
    
    # For direct execution with example parameters:
    fig, summary = visualize_error_growth(initial_position, initial_velocity, t_lock)
    plt.show()