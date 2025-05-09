import numpy as np
import matplotlib.pyplot as plt # type: ignore
from matplotlib.widgets import Slider, Button, RadioButtons # type: ignore
from target_vis import visualize_error_growth
import threading

class InteractiveTargetApp:
    def __init__(self):
        self.initial_position = np.array([5000, 2000, 15000])  # Default position (meters)
        self.initial_velocity = np.array([-300, -250, -100])   # Default velocity (m/s)
        self.t_lock = 3.0  # Default lock-on time (seconds)
        self.calculation_in_progress = False
        self.quality_mode = "HighQuality"  # Default to fast quality for responsive UI
        
        self.create_figure()
        self.update_plot()
        
    def create_figure(self):
        # Create main figure
        self.fig = plt.figure(figsize=(16, 10))
        plt.subplots_adjust(left=0.25, bottom=0.25)
        
        # Set up axes for plot
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        # Add status text
        self.status_text = self.fig.text(0.5, 0.97, "Ready", ha='center', fontsize=12, 
                                        bbox=dict(facecolor='lightgreen', alpha=0.8))
        
        # Set up slider axes
        self.ax_px = plt.axes([0.25, 0.15, 0.65, 0.03])
        self.ax_py = plt.axes([0.25, 0.12, 0.65, 0.03])
        self.ax_pz = plt.axes([0.25, 0.09, 0.65, 0.03])
        
        self.ax_vx = plt.axes([0.25, 0.06, 0.65, 0.03])
        self.ax_vy = plt.axes([0.25, 0.03, 0.65, 0.03])
        self.ax_vz = plt.axes([0.25, 0.00, 0.65, 0.03])
        
        self.ax_t_lock = plt.axes([0.05, 0.25, 0.03, 0.65])
        
        # Create sliders
        self.slider_px = Slider(self.ax_px, 'Position X (m)', -10000, 10000, valinit=self.initial_position[0])
        self.slider_py = Slider(self.ax_py, 'Position Y (m)', -10000, 10000, valinit=self.initial_position[1])
        self.slider_pz = Slider(self.ax_pz, 'Position Z (m)', 0, 25000, valinit=self.initial_position[2])
        
        self.slider_vx = Slider(self.ax_vx, 'Velocity X (m/s)', -800, 800, valinit=self.initial_velocity[0])
        self.slider_vy = Slider(self.ax_vy, 'Velocity Y (m/s)', -800, 800, valinit=self.initial_velocity[1])
        self.slider_vz = Slider(self.ax_vz, 'Velocity Z (m/s)', -800, 800, valinit=self.initial_velocity[2])
        
        self.slider_t_lock = Slider(self.ax_t_lock, 'Lock Time (s)', 0.1, 10.0, valinit=self.t_lock, orientation='vertical')
        
        # Add quality control radio buttons
        self.ax_quality = plt.axes([0.05, 0.05, 0.15, 0.1])
        self.quality_selector = RadioButtons(
            self.ax_quality, ('Fast', 'High Quality'), active=0)
        self.quality_selector.on_clicked(self.set_quality)
        
        # Add calculate button instead of automatic update
        self.button_ax = plt.axes([0.05, 0.20, 0.15, 0.04])
        self.calculate_button = Button(self.button_ax, 'Calculate')
        self.calculate_button.on_clicked(self.calculate)
        
        # Add reset button
        self.reset_ax = plt.axes([0.05, 0.15, 0.15, 0.04])
        self.reset_button = Button(self.reset_ax, 'Reset Values')
        self.reset_button.on_clicked(self.reset)
        
    def set_quality(self, label):
        self.quality_mode = label
        print(f"Set quality mode to: {self.quality_mode}")
        
    def calculate(self, event=None):
        if self.calculation_in_progress:
            print("Calculation already in progress...")
            return
            
        # Get current values from sliders
        self.initial_position = np.array([
            self.slider_px.val,
            self.slider_py.val,
            self.slider_pz.val
        ])
        
        self.initial_velocity = np.array([
            self.slider_vx.val,
            self.slider_vy.val,
            self.slider_vz.val
        ])
        
        self.t_lock = self.slider_t_lock.val
        
        # Start calculation in a separate thread
        self.calculation_in_progress = True
        self.status_text.set_text("Calculating... Please wait")
        self.status_text.set_bbox(dict(facecolor='yellow', alpha=0.8))
        self.fig.canvas.draw_idle()
        
        thread = threading.Thread(target=self.calculate_in_thread)
        thread.daemon = True
        thread.start()
        
    def calculate_in_thread(self):
        """Run the calculation in a separate thread to avoid UI freezing"""
        try:
            # Adjust sample sizes based on quality mode
            if self.quality_mode == "Fast":
                # Override function to use smaller sample sizes for interactive use
                from functools import partial
                from target_vis import visualize_error_growth as orig_vis
                
                # Create a partial function with reduced sample sizes
                vis_func = lambda pos, vel, t_lock: orig_vis(pos, vel, t_lock)
            else:
                # Use the full-quality visualization
                vis_func = visualize_error_growth
            
            # Create new visualization
            fig, self.summary = vis_func(
                self.initial_position, 
                self.initial_velocity, 
                self.t_lock
            )
            
            # Close the figure generated by visualize_error_growth
            plt.close(fig)
            
            # Update plot from the main thread
            plt.pause(0.01)  # Small pause to let the main thread process
            self.update_plot_ui()
            
        except Exception as e:
            print(f"Error in calculation: {e}")
            self.status_text.set_text(f"Error: {str(e)}")
            self.status_text.set_bbox(dict(facecolor='red', alpha=0.8))
            
        finally:
            self.calculation_in_progress = False
        
    def update_plot_ui(self):
        """Update the UI with calculation results (called from main thread)"""
        try:
            # Clear current plot
            self.ax.clear()
            
            # Create new visualization (reusing the summary from calculate_in_thread)
            from target_vis import plot_3d_visualization
            
            # Get summary from previous calculation
            summary = self.summary
            
            # Plot directly to our existing axis
            plot_3d_visualization(
                self.initial_position, self.initial_velocity, 
                None,  # We'll recalculate the nominal trajectory
                None, None, None,  # These will be recalculated
                0.05, self.t_lock, 
                existing_ax=self.ax,  # Pass our existing axis
                reuse_summary=summary  # Reuse the summary we calculated earlier
            )
            
            # Update our summary text
            if hasattr(self, 'summary_text'):
                self.summary_text.remove()
            
            summary_str = (
                f"Initial Speed: {np.linalg.norm(self.initial_velocity):.1f} m/s\n"
                f"Maximum Error without Maneuvering: {summary['max_position_error_no_maneuver']:.1f} m\n"
                f"Maximum Error with Maneuvering: {summary['max_position_error_with_maneuver']:.1f} m\n"
                f"Maximum Velocity with Maneuvering: {summary['max_velocity_with_maneuver']:.1f} m/s"
            )
            
            self.summary_text = self.fig.text(0.05, 0.92, summary_str, fontsize=12, 
                                            bbox=dict(facecolor='white', alpha=0.8))
            
            # Update status
            self.status_text.set_text("Ready")
            self.status_text.set_bbox(dict(facecolor='lightgreen', alpha=0.8))
            
            # Draw the plot
            self.fig.canvas.draw_idle()
            
        except Exception as e:
            print(f"Error updating plot UI: {e}")
            self.status_text.set_text(f"Error updating plot: {str(e)}")
            self.status_text.set_bbox(dict(facecolor='red', alpha=0.8))
            self.fig.canvas.draw_idle()
        
    def update_plot(self):
        """Initial plot when starting the application"""
        try:
            # Do a fast calculation for the initial plot
            from target_vis import visualize_error_growth
            
            # Create new visualization with reduced sample sizes for responsiveness
            fig, summary = visualize_error_growth(
                self.initial_position, 
                self.initial_velocity, 
                self.t_lock
            )
            
            # Close the figure generated by visualize_error_growth
            plt.close(fig)
            
            # Update our summary text
            if hasattr(self, 'summary_text'):
                self.summary_text.remove()
            
            summary_str = (
                f"Initial Speed: {np.linalg.norm(self.initial_velocity):.1f} m/s\n"
                f"Maximum Error without Maneuvering: {summary['max_position_error_no_maneuver']:.1f} m\n"
                f"Maximum Error with Maneuvering: {summary['max_position_error_with_maneuver']:.1f} m\n"
                f"Maximum Velocity with Maneuvering: {summary['max_velocity_with_maneuver']:.1f} m/s"
            )
            
            self.summary_text = self.fig.text(0.05, 0.92, summary_str, fontsize=12, 
                                            bbox=dict(facecolor='white', alpha=0.8))
            
            # Update status
            self.status_text.set_text("Ready")
            self.status_text.set_bbox(dict(facecolor='lightgreen', alpha=0.8))
            
            # Draw the plot
            self.fig.canvas.draw_idle()
            
        except Exception as e:
            print(f"Error in initial plot: {e}")
            self.status_text.set_text(f"Error: {str(e)}")
            self.status_text.set_bbox(dict(facecolor='red', alpha=0.8))
            self.fig.canvas.draw_idle()
        
    def reset(self, event):
        # Reset all sliders to default values
        self.slider_px.set_val(5000)
        self.slider_py.set_val(2000)
        self.slider_pz.set_val(15000)
        
        self.slider_vx.set_val(-300)
        self.slider_vy.set_val(-250)
        self.slider_vz.set_val(-100)
        
        self.slider_t_lock.set_val(3.0)
        
    def show(self):
        plt.show()

if __name__ == "__main__":
    app = InteractiveTargetApp()
    app.show()