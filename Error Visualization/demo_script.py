import numpy as np
import matplotlib.pyplot as plt # type: ignore
from target_vis import visualize_error_growth
import datetime
import os

def demonstrate_target_scenarios():
    """Run multiple target scenarios to demonstrate error growth patterns"""
    startime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    # Define several scenarios
    scenarios = [
        # {
        #     "name": "High-Altitude Intercept",
        #     "position": np.array([5000, 2000, 20000]),  # High altitude (20km)
        #     "velocity": np.array([-300, -250, -100]),   # Incoming, medium speed
        #     "t_lock": 2.0,                             # Short lock-on time
        # },
        {
            "name": "Low-Altitude Cruise Missile",
            "position": np.array([8000, 3000, 800]),   # Low altitude (500m)
            "velocity": np.array([225, -225, 0]),      # Fast, shallow approach
            "t_lock": 2.0,                             # Medium lock-on time
        },
        {
            "name": "Low-Altitude Cruise Missile",
            "position": np.array([8000, 3000, 800]),   # Low altitude (500m)
            "velocity": np.array([225, -225, 0]),      # Fast, shallow approach
            "t_lock": 5.0,                             # Medium lock-on time
        }
        # {
        #     "name": "Low-Altitude Subsonic Cruise Missile",
        #     "position": np.array([8000, 3000, 500]),   # Low altitude (500m)
        #     "velocity": np.array([-150, 150, 0]),      # Subsonic, shallow approach
        #     "t_lock": 5.0,                             # Medium lock-on time
        # }
        # {
        #     "name": "High-Speed Ballistic Target",
        #     "position": np.array([0, 0, 15000]),   # Medium-high altitude
        #     "velocity": np.array([0, 0, -1029]),     # Very fast, steep approach
        #     "t_lock": 1,                             # Very short lock-on time
        # }
        # {
        #     "name": "Maneuvering Target with Long Lock-on",
        #     "position": np.array([2000, 1000, 5000]),  # Medium altitude
        #     "velocity": np.array([-200, -150, -50]),   # Slower speed
        #     "t_lock": 8.0,                             # Long lock-on time - more error growth
        # }
    ]
    
    # Create a figure for the summary plot
    fig_summary, axes = plt.subplots(2, 2, figsize=(18, 14))
    axes = axes.flatten()
    
    # Process each scenario
    for i, scenario in enumerate(scenarios):
        print(f"\nProcessing Scenario {i+1}: {scenario['name']}")
        print("-" * 60)
        print(f"Initial position: {scenario['position']} meters")
        print(f"Initial velocity: {scenario['velocity']} m/s")
        print(f"Initial speed: {np.linalg.norm(scenario['velocity']):.1f} m/s")
        print(f"Lock-on time: {scenario['t_lock']} seconds")
        
        # Run the visualization
        fig, summary = visualize_error_growth(
            scenario['position'], 
            scenario['velocity'], 
            scenario['t_lock']
        )
        
        # Print summary statistics
        # print("\nNumerical Summary:")
        # print(f"Maximum position error without maneuvering: {summary['max_position_error_no_maneuver']:.2f} meters")
        # print(f"Mean position error without maneuvering: {summary['mean_position_error_no_maneuver']:.2f} meters")
        # print(f"Maximum position error with 3G maneuvering: {summary['max_position_error_with_maneuver']:.2f} meters")
        # print(f"Mean position error with 3G maneuvering: {summary['mean_position_error_with_maneuver']:.2f} meters")
        # print(f"Maximum velocity with maneuvering: {summary['max_velocity_with_maneuver']:.2f} m/s")
        # print(f"Mean velocity with maneuvering: {summary['mean_velocity_with_maneuver']:.2f} m/s")
        
        # Save the individual figure
        fpath = f"exports/{startime}"
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        fname = f"{fpath}/scenario_{i+1}_{scenario['name'].replace(' ', '_')}.png"

        fig.savefig(fname, 
                   bbox_inches='tight', dpi=150)
        
        # Add to summary plot
        ax = axes[i]
        
        # Create bar plot of error growth
    #     bar_data = [
    #         # summary['max_position_error_no_maneuver'],
    #         summary['max_position_error_with_maneuver']
    #     ]
        
    #     bars = ax.bar(['No Maneuver', 'With 3G Maneuver'], bar_data, 
    #                  color=['blue', 'red'], alpha=0.7)
        
    #     # Add labels
    #     ax.set_ylabel('Maximum Error (meters)')
    #     ax.set_title(f"Scenario {i+1}: {scenario['name']}\n"
    #                 f"Lock Time: {scenario['t_lock']}s, "
    #                 f"Initial Speed: {np.linalg.norm(scenario['velocity']):.0f} m/s")
        
    #     # Add value labels on bars
    #     for bar in bars:
    #         height = bar.get_height()
    #         ax.text(bar.get_x() + bar.get_width()/2., height + 100,
    #                f'{height:.0f}m',
    #                ha='center', va='bottom', rotation=0)
        
    #     # Add a statement about error growth factor
    #     growth_factor = summary['max_position_error_with_maneuver'] / summary['max_position_error_no_maneuver']
    #     ax.text(0.5, 0.05, 
    #            f"3G maneuvering increases error by {growth_factor:.1f}x",
    #            transform=ax.transAxes, ha='center', 
    #            bbox=dict(facecolor='yellow', alpha=0.5))
    
    # # Add a main title
    # fig_summary.suptitle('Missile Interceptor Error Growth Comparison', 
    #                     fontsize=20, y=0.98)
    
    # # Adjust layout and save
    # plt.tight_layout(rect=[0, 0, 1, 0.96])
    # fig_summary.savefig("scenario_summary_comparison.png", bbox_inches='tight', dpi=150)
    
    print("\nAll scenarios processed. Displaying results...")
    plt.show()

# Run the demonstration
if __name__ == "__main__":
    demonstrate_target_scenarios()