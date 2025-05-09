# VT-Nemesis-AIAA-Code Repository Overview
This repository contains the code and resources for the Virginia Tech Aerospace Engineering Senior Design Nemesis Team (Missile 2), competing in the AIAA Anti-Missile Missile Project.

## Repository Structure

## Actuators/
This folder contains code related to actuator control and implementation for missile control surfaces.

### servo_control.ino:
Arduino code for controlling servos for fin actuation based on control inputs.

### actuator_calibration.py:
Python script for calibrating actuator response and deadband elimination.

## Aerodynamics/
This folder contains code related to aerodynamic analysis of the missile.

### stability.m: 
MATLAB script for analyzing the stability margin of the missile over time, visualizing how stability changes throughout flight until apogee.

### aero_coefficients.m:
MATLAB script for calculating aerodynamic coefficients from simulation data.

## Error Visualization/
This folder contains tools for visualizing and analyzing error propagation for target tracking and interception.

### target_vis.py: 
Python tool for visualizing target error growth over time, including both non-maneuvering and maneuvering scenarios.

### demo_script.py: 
Script to demonstrate various target scenarios, including low-altitude cruise missiles with different lock-on times.

### error_propagation.py:
Python tool for simulating and analyzing error propagation in the guidance system.

## Open Rocket/
This folder contains Open Rocket simulation files and analysis scripts.

### nemesis_missile.ork:
Open Rocket design file for the Nemesis missile.

### flight_simulation_analysis.py:
Python script for post-processing Open Rocket simulation data.

### trajectory_comparison.m:
MATLAB script comparing OpenRocket trajectories with mathematical models.

## Payload/
This folder contains code related to the missile payload system.

### seeker_simulation.py:
Simulation of the seeker head behavior and target acquisition.

### proximity_sensor.ino:
Arduino code for the proximity sensor integration and detonation logic.

### warhead_optimization.m:
MATLAB analysis for optimizing warhead performance based on target parameters.

## Thrust Vectoring Experiment/
This folder contains code and analysis for thrust vectoring experiments using jet vanes.

### VaneExperimentAnalysis.m: 
MATLAB script for analyzing thrust vectoring test data, including calibration, thrust curves, and moment vs. angle of attack relationships.

### multiLCRead.ino: 
Arduino code for reading multiple load cells simultaneously during thrust vectoring experiments.

### vane_design_optimization.py:
Python script for optimizing vane design parameters.

### thermal_analysis.m:
MATLAB script for thermal analysis of jet vanes during operation.

### control_algorithm.py:
Python implementation of the thrust vector control algorithm.

### data_acquisition.ino:
Arduino code for data acquisition during thrust vector testing.

### README.md: 
Brief overview of the thrust vectoring experiment.

## Time to Strike/
This folder contains tools for trajectory analysis and missile intercept calculations.

### trajSets.txt: 
Collection of incoming missile profiles with various parameters including altitude, range, and speeds.

### incomingMissileProfile.m: 
MATLAB function for calculating and visualizing possible impact zones for incoming missiles, including maneuverability analysis.

### interceptor_simulation.m:
MATLAB simulation of interceptor performance against various targets.

### engagement_geometry.py:
Python script for analyzing engagement geometry and optimal launch parameters.

### boost_phase_analysis.m:
Analysis of boost phase intercept opportunities and constraints.

### impact_probability.py:
Calculation of probability of successful intercept based on various parameters.

### old functions/:
#### unpoweredParabolicDescent.m: 
Function for modeling unpowered parabolic descent trajectories.

#### linearPoweredDescent.m: 
Function for modeling linear powered descent trajectories with various parameters.

## Turn Rates/
This folder contains analysis of missile maneuverability and turn performance.

### max_turn_rate.m:
MATLAB script calculating maximum achievable turn rates at different flight conditions.

### g_load_analysis.py:
Python script for analyzing g-loading during various maneuvers.

### structural_limits.m:
Analysis of structural limitations on maximum turn performance.

### pursuit_curves.py:
Simulation of various pursuit curves and their efficiency for interception.