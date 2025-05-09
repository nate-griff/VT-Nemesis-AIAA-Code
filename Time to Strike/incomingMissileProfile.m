function [impact_zone_within_grid, impact_times_within_grid, positions_within_grid, velocities_within_grid] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicsSettings, outputLocation, maxSpeed)
    % az0: initial azimuth angle of the missile (rad)
    % range0: initial range of the missile (m)
    % alt0: initial altitude of the missile (m)
    % vx0, vy0, vz0: initial velocity components of the missile (m/s)
    % graphicsSettings: struct containing horiz_calc, vert_calc, acc_calc, dt, max_iter
    % outputLocation: directory to save output files
    
    % Expand graphicsSettings
    horiz_calc = graphicsSettings.horiz_calc;
    vert_calc = graphicsSettings.vert_calc;
    acc_calc = graphicsSettings.acc_calc;
    dt = graphicsSettings.dt;
    max_iter = graphicsSettings.max_iter;
    numValues = 20;

    % Create output folder name
    [~, a0, ~, ~] = atmoscoesa(alt0); % Speed of sound at the initial altitude
    M = norm([vx0, vy0, vz0])/a0; % Mach number calculation
    folderName = sprintf('Alt_%.1f_m_Mach_%.2f_R0_%.1f_m', alt0, M, range0);
    fullOutputPath = fullfile(outputLocation, folderName);
    if ~exist(fullOutputPath, 'dir')
        mkdir(fullOutputPath);
    end

    % Create data subfolder
    dataOutputPath = fullfile(fullOutputPath, 'data');
    if ~exist(dataOutputPath, 'dir')
        mkdir(dataOutputPath);
    end

    % Open text file for writing
    logFile = fullfile(fullOutputPath, 'output.txt');
    fid = fopen(logFile, 'w');

    % Check for GPU availability
    % useGPU = gpuDeviceCount > 0;
    useGPU = false;
    if useGPU
        gpuDevice(1); % Initialize GPU device
    end
    
    % Constants and initial conditions
    error = 0.05; % Error in the initial position and velocity of the missile
    r_d = 5 * 1609.34; % defense radius (m)

    % Initial position of the missile
    x0 = range0*cos(az0);
    y0 = range0*sin(az0);
    z0 = alt0;

    % Initial position and velocity vectors
    r0 = [x0; y0; z0]; 
    v0 = [vx0; vy0; vz0];

    [hit_pos_nomaneuver, willhit, t_ground] = traj_to_hit(r0, v0); % Check if the missile will hit the defense zone

    % Calculate the maximum g maneuver impact zone
    [impact_zone, impact_times, max_travel_time, positions, velocities] = max_g_maneuver(r0, v0, 3*9.81, maxSpeed*a0, horiz_calc, vert_calc, acc_calc, dt, max_iter, useGPU);

    % Create grid for contour plot
    x = impact_zone(:,1);
    y = impact_zone(:,2);
    z = impact_times;
    
    % Combine x, y, z into data
    data = [x, y, z];
    
    % Sort data ascending on x, y, impact_time
    data_sorted = sortrows(data, [1, 2, 3]);
    
    % Get unique x, y pairs, keeping the one with minimal impact_time
    [~, idx_unique] = unique(data_sorted(:,1:2), 'rows', 'stable');
    data_unique = data_sorted(idx_unique, :);
    
    % Use the unique data
    x = data_unique(:,1);
    y = data_unique(:,2);
    z = data_unique(:,3);
    
    % Ensure x, y, z are double arrays
    if useGPU
        x = gather(x);
        y = gather(y);
        z = gather(z);
    end

    % Define grid over impact zone
    xlin = linspace(min(x), max(x), 100);
    ylin = linspace(min(y), max(y), 100);
    [X, Y] = meshgrid(xlin, ylin);
    
    % Grid the data
    Z = griddata(x, y, z, X, Y, 'natural');
    
    % Find the smallest impact time within the defense radius
    within_defense_zone = sqrt(X.^2 + Y.^2) <= r_d;
    min_impact_time = [];
    min_impact_point = [];
    if any(within_defense_zone, 'all')
        Z_within_defense = Z(within_defense_zone);
        [min_impact_time, idx] = min(Z_within_defense(:));
        [row, col] = find(Z == min_impact_time, 1);
        min_impact_point = [X(row, col), Y(row, col), 0];
    else
        disp('error: no impact points within the defense zone');
        min_impact_time = NaN;
        min_impact_point = [NaN, NaN, NaN];
    end

    if willhit
        r_hit = hit_pos_nomaneuver(1:2); % Hit location from the center of the defense zone
        hit_msg1 = sprintf('The missile is traveling at Mach %.2f and was detected %.1f mi from the center of the defense zone\n', M, norm(r0(1:2)) / 1609.34);
        hit_msg2 = sprintf('If it continues to travel on its current trajectory, it will hit the defense zone in %.1f±%.1f seconds\nat the location (%.1f, %.1f, %.1f) m which is %.1f±%.1f m from the center of the of the defense zone', t_ground,t_ground*error, hit_pos_nomaneuver(1), hit_pos_nomaneuver(2), hit_pos_nomaneuver(3), norm(r_hit), norm(r_hit)*error);
        hit_msg3 = sprintf('If the missile maneuvers at or below 3g, it can hit the ground in %.1f±%.1f seconds\nat the location (%.1f, %.1f, %.1f) m which is %.1f±%.1f m from the center of the defense zone', min_impact_time, min_impact_time*error, min_impact_point(1), min_impact_point(2), 0, norm(min_impact_point), norm(min_impact_point)*error);
        fprintf(fid, '%s\n%s\n%s\n', hit_msg1, hit_msg2, hit_msg3);
    elseif any(within_defense_zone)
        not_hit_msg1 = sprintf('The missile is traveling at Mach %.2f and was detected %.1f mi from the center of the defense zone\n', M, norm(r0(1:2)) / 1609.34);
        not_hit_msg2 = 'The missile is not on a trajectory to hit the defense zone, but it is within maneuvering range\n';
        not_hit_msg3 = sprintf('The missile can hit the ground in %.1f±%.1f seconds at the location (%.1f, %.1f, %.1f) m \nwhich is %.1f±%.1f m from the center of the defense zone', min_impact_time, min_impact_time*error, min_impact_point(1), min_impact_point(2), z0, norm(min_impact_point), norm(min_impact_point)*error);
        fprintf(fid, '%s\n%s\n%s\n', not_hit_msg1, not_hit_msg2, not_hit_msg3);
    else
        not_hit_zone_msg1 = sprintf('The missile is traveling at Mach %.2f and was detected %.1f mi from the center of the defense zone\n', M, norm(r0(1:2)) / 1609.34);
        r_hit_outside = hit_pos_nomaneuver(1:2); % Hit location from the center of the defense zone
        not_hit_zone_msg2 = sprintf('The missile will not hit the defense zone, but will hit the ground in %.1f±%.1f seconds\nat the location (%.1f, %.1f, %.1f) m which is %.1f±%.1f m outside of the defense zone', t_ground, t_ground*error, hit_pos_nomaneuver(1), hit_pos_nomaneuver(2), hit_pos_nomaneuver(3), norm(r_hit_outside), norm(r_hit_outside)*error);
        fprintf(fid, '%s\n%s\n', not_hit_zone_msg1, not_hit_zone_msg2);
    end

    % Plot the full impact zone
    figure1 = figure;
    contourf(X, Y, Z, 128, 'LineStyle', 'none');
    colorbar;
    colormap(flipud(turbo));

    ylabel(colorbar, 'Impact Time (s)','FontName', 'Times New Roman');
    set(gca, 'FontName', 'Times New Roman');
    hold on;
    
    % Plot the hit position without maneuver as an 'x'
    plot(hit_pos_nomaneuver(1), hit_pos_nomaneuver(2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k');
    
    % Plot the hit position with maneuver as a '*'
    if ~isempty(min_impact_point)
        plot(min_impact_point(1), min_impact_point(2),'*', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k');
    end

    % Plot the missile's initial position with an arrow annotation
    quiver(r0(1), r0(2), v0(1)*5, v0(2)*5, 'm', 'LineWidth', 2.5, 'MaxHeadSize', 12);
    

    % Plot the defense zone circle
    theta = linspace(0, 2*pi, 360);
    x_circle = r_d * cos(theta);
    y_circle = r_d * sin(theta);
    plot(x_circle, y_circle, '--', 'LineWidth', 2, 'Color', 'black');

    % Labels and title
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    if acc_calc == 1
        title('Impact Zone under 3G Load');
    else
        title('Impact Zone under 0.1-3G Load');
    end
    sub = {[sprintf('M_0 = %.2f, alt_0 = %.1f m, R_0 = %.1f m, t_{hit, linear} = %.1f sec, t_{hit, shortest} = %.1f sec', M, alt0, range0, t_ground, min_impact_time)], ...
        [sprintf('Time Range = %.2f sec', max_travel_time)]};
    subtitle(sub);
    legend('Impact Zone', 'Hit Position for linear descent', 'Hit Position for Shortest Impact', 'Missile Initial Position','Defense Zone', 'Location', 'best');
    set(gca, 'FontName', 'Times New Roman');

    axis equal;
    grid on;
    hold off;
    savefig(figure1, fullfile(fullOutputPath, 'impact_zone.fig'));

    % Plot the defense zone circle and grid data within the defense zone
    figure2 = figure;
    hold on;

    % Plot the grid data within the defense zone
    within_defense_zone_grid = sqrt(X.^2 + Y.^2) <= r_d * 1.2;
    X_scaled = X .* within_defense_zone_grid;
    Y_scaled = Y .* within_defense_zone_grid;
    Z_scaled = Z .* within_defense_zone_grid;
    contourf(X_scaled, Y_scaled, Z_scaled, 128, 'LineStyle', 'none');
    colorbar;
    colormap(flipud(turbo));
    t_scaled_min = floor(min(Z_scaled(Z_scaled>0)));
    t_scaled_max = max(ceil(max(Z_scaled)));
    clim([t_scaled_min, t_scaled_max]);
    ylabel(colorbar, 'Impact Time (s)','FontName', 'Times New Roman');
    set(gca, 'FontName', 'Times New Roman');
    
    % Plot the defense zone circle
    plot(x_circle, y_circle, '--', 'LineWidth', 2, 'Color', 'black');

    % Plot the hit position without maneuver as an 'x'
    plot(hit_pos_nomaneuver(1), hit_pos_nomaneuver(2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k');
    
    % Plot the hit position with maneuver as a '*'
    if ~isempty(min_impact_point)
        plot(min_impact_point(1), min_impact_point(2),'*', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k');
    end

    % Labels and title
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Impact Zone within Defense Zone');
    sub = {[sprintf('M_0 = %.2f, alt_0 = %.1f m, R_0 = %.1f m, t_{hit, linear} = %.1f sec, t_{hit, shortest} = %.1f sec', M, alt0, range0, t_ground, min_impact_time)], ...
        [sprintf('Time Range = %.2f sec', max_travel_time)]};
    subtitle(sub);
    set(gca, 'FontName', 'Times New Roman');
    axis equal;
    grid on;
    hold off;
    savefig(figure2, fullfile(fullOutputPath, 'impact_zone_within_defense.fig'));

    % Filter impact points within the defense zone grid
    within_defense_zone_grid = sqrt(impact_zone(:,1).^2 + impact_zone(:,2).^2) <= r_d;

    % Create new variables for impact_zone, impact_times, positions, and velocities within the defense zone grid
    impact_zone_within_grid = impact_zone(within_defense_zone_grid, :);
    impact_times_within_grid = impact_times(within_defense_zone_grid);
    positions_within_grid = positions(within_defense_zone_grid);
    velocities_within_grid = velocities(within_defense_zone_grid);

    % Write data to files
    writematrix(impact_zone_within_grid, fullfile(dataOutputPath, 'impact_zone_within_grid.csv'));
    writematrix(impact_times_within_grid, fullfile(dataOutputPath, 'impact_times_within_grid.csv'));
    save(fullfile(dataOutputPath, 'positions_within_grid.mat'), 'positions_within_grid');
    save(fullfile(dataOutputPath, 'velocities_within_grid.mat'), 'velocities_within_grid');

    % Save only a few trajectories
    var_size = length(positions_within_grid);
    idxs = round(linspace(1, var_size, numValues));

    shrunk_pos_data = positions_within_grid(idxs);
    shrunk_vel_data = velocities_within_grid(idxs);

    trajectory_data_set = cell(size(shrunk_pos_data));

    for i = 1:length(idxs)
        trajectory_data_set{i} = [shrunk_pos_data{i}, shrunk_vel_data{i}];
    end
    save(fullfile(dataOutputPath, 'shrunk_PosAndVel_data_export.mat'), "trajectory_data_set")
    disp("Profile complete. Output files saved in: " + fullOutputPath);
end

function [hit_position, willhit, t_ground] = traj_to_hit(r0, v0)
    % willhit: true if the missile is on a trajectory that will hit the defense zone
    % r0: missile position (m) (vector: x, y, z)
    % v0: missile velocity (m/s) (vector: v_x, v_y, v_z)

    % Define the ground hit radius in meters (5 miles)
    hit_radius = 5 * 1609.34;

    % Calculate the time to hit the ground (z = 0)
    if v0(3) == 0
        t_ground = Inf;
    else
        t_ground = -r0(3) / v0(3);
    end
    
    % Calculate the position at the time of hitting the ground
    hit_position = r0 + v0 * t_ground;

    % Check if the hit position is within the hit radius
    distance_to_origin = norm(hit_position(1:2));
    willhit = distance_to_origin <= hit_radius;
end

function [impact_zone, impact_times, max_travel_time, positions, velocities] = max_g_maneuver(r0, v0, max_acc, max_speed, num_theta, num_phi, num_acc, dt, max_iterations, useGPU)
    % max_g_maneuver Simulates missile trajectory under maximum g-force maneuvering.
    %
    % Parameters:
    %   r0            - Initial position of the missile [x, y, z] (m)
    %   v0            - Initial velocity vector of the missile [vx, vy, vz] (m/s)
    %   max_acc       - Maximum acceleration (m/s^2)
    %   max_speed     - Maximum allowable speed (m/s)
    %   num_theta     - Number of theta (azimuth) angles
    %   num_phi       - Number of phi (elevation) angles
    %   num_acc       - Number of acceleration magnitudes
    %   dt            - Time step (s)
    %   max_iterations- Maximum number of iterations for each loop
    %   useGPU        - Boolean flag to use GPU arrays
    %
    % Outputs:
    %   impact_zone   - Impact positions within the defense zone
    %   impact_times  - Impact times corresponding to impact_zone
    %   positions     - Cell array of position trajectories
    %   velocities    - Cell array of velocity trajectories

    max_travel_time = dt * max_iterations;

    % Normalize the initial velocity vector
    v0_norm = v0 / norm(v0);
    epsilon = 1e-6; % Tolerance for floating-point comparisons

    % Edge case check: If v0 is along the x-axis
    if abs(abs(v0_norm(1)) - 1) < epsilon && abs(v0_norm(2)) < epsilon && abs(v0_norm(3)) < epsilon
        % v0 is along the x-axis, no rotation needed
        R = eye(3);
    else
        % Yaw rotation to bring v0 into x-z plane
        psi = atan2(v0_norm(2), v0_norm(1));
        R_yaw = [cos(-psi), -sin(-psi), 0;
                 sin(-psi),  cos(-psi), 0;
                 0,          0,         1];
        v0_yaw = R_yaw * v0_norm;

        % Edge case after yaw rotation: If v0 is now along the x-axis
        if abs(v0_yaw(2)) < epsilon && abs(v0_yaw(3)) < epsilon
            % v0 is along the x-axis after yaw rotation
            R = R_yaw;
        else
            % Pitch rotation to align v0 with x-axis
            theta = atan2(v0_yaw(3), v0_yaw(1));
            R_pitch = [cos(-theta), 0, sin(-theta);
                       0,           1, 0;
                      -sin(-theta), 0, cos(-theta)];
            R = R_pitch * R_yaw;
        end
    end

    % Rotate initial position and velocity
    r0_rot = R * r0;
    v0_rot = R * v0;

    % Initialize rotated position and velocity
    r_init = r0_rot;
    v_init = v0_rot;

    % Define angular resolution for half the hemisphere of accelerations
    num_theta_half = ceil(num_theta / 2);
    theta_vals = linspace(0, pi, num_theta_half);
    phi_vals = linspace(0, pi/2, num_phi);

    % Define acceleration magnitudes
    if num_acc == 1
        acceleration = max_acc;
    else
        acceleration = linspace(1, max_acc, num_acc); % m/s^2, from 1 to max_acc
    end

    % Total number of combinations
    num_combinations = numel(acceleration) * numel(theta_vals) * numel(phi_vals);

    if useGPU
        % Transfer data to GPU
        acceleration = gpuArray(acceleration);
        theta_vals = gpuArray(theta_vals);
        phi_vals = gpuArray(phi_vals);
        r_init = gpuArray(r_init);
        v_init = gpuArray(v_init);
    end

    % -------------------
    % Initialize Progress Tracking
    % -------------------
    total_iterations = num_combinations; % Ensure 'num_combinations' is defined before this point

    % Create a DataQueue for receiving progress updates
    progressQueue = parallel.pool.DataQueue;

    % Initialize a counter for completed iterations
    completed = 0;
    lastPercent = 0;

    % Initialize the waitbar
    hWaitbar = waitbar(0, 'Processing Trajectories...', 'Name', 'Simulation Progress');

    % Define the callback function to update the waitbar
    afterEach(progressQueue, @updateWaitbar);

    % Define the nested callback function
    function updateWaitbar(~)
        completed = completed + 1;
        currentPercent = floor((completed / total_iterations) * 100);
        if currentPercent > lastPercent
            lastPercent = currentPercent;
            waitbar(currentPercent/100, hWaitbar, sprintf('Progress: %d%%', currentPercent));
        end
    end

    % Preallocate arrays to store ground impact points and times
    if useGPU
        impact_zone = gpuArray.zeros(num_combinations, 3, 'double'); % [x, y, z]
        impact_times = gpuArray.zeros(num_combinations, 1, 'double');
    else
        impact_zone = zeros(num_combinations, 3); % [x, y, z]
        impact_times = zeros(num_combinations, 1);
    end

    % Initialize cell arrays to store trajectories
    positions = cell(num_combinations,1);
    velocities = cell(num_combinations,1);

    % Initialize parallel pool once (optional, can be managed outside the function)
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool; %#ok<*NASGU>
    end

    % Parallel loop over all combinations with progress tracking
    parfor idx = 1:num_combinations
        % Extract indices for current combination
        [acc_idx, theta_idx, phi_idx] = ind2sub([numel(acceleration), numel(theta_vals), numel(phi_vals)], idx);

        acc = acceleration(acc_idx);
        theta = theta_vals(theta_idx);
        phi = phi_vals(phi_idx);

        % Compute the acceleration vector for the current angles
        g = 9.81; % gravitational acceleration (m/s^2)
        ax = acc * cos(phi) * cos(theta);
        ay = acc * cos(phi) * sin(theta);
        az = acc * sin(phi) - g;
        a = [ax; ay; az];

        % Initialize position, velocity, and time
        r = r_init;
        v = v_init;
        t = 0;

        % Initialize storage for trajectory
        pos_traj = r';
        vel_traj = v';

        % Simulation loop
        iteration_count = 1;
        while r(3) > 0 && iteration_count < max_iterations
            % Update velocity with total acceleration
            v = v + a * dt;

            % Enforce maximum speed
            speed = norm(v);
            if speed > max_speed
                v = (v / speed) * max_speed;
            end

            % Update position with current velocity
            r = r + v * dt;

            % Update time and iteration count
            t = t + dt;
            iteration_count = iteration_count + 1;

            % Store trajectory
            pos_traj = [pos_traj; r'];
            vel_traj = [vel_traj; v'];
        end

        % Post-simulation: Determine impact
        if r(3) <= 0
            % Hit the ground
            hit_position = r';
            impact_zone(idx, :) = hit_position; % Store [x, y, z]
            impact_times(idx) = t;
        else
            % Did not hit the ground within max iterations
            if useGPU
                impact_zone(idx, :) = gpuArray([NaN, NaN, NaN]);
                impact_times(idx) = gpuArray(NaN);
            else
                impact_zone(idx, :) = [NaN, NaN, NaN];
                impact_times(idx) = NaN;
            end
        end

        % Store trajectories
        if useGPU
            positions{idx} = gather(pos_traj);
            velocities{idx} = gather(vel_traj);
        else
            positions{idx} = pos_traj;
            velocities{idx} = vel_traj;
        end

        % Send progress update
        send(progressQueue, 1);
    end

    % Close the waitbar after completion
    close(hWaitbar);

    % -------------------
    % Mirror the results in the rotated frame (mirror y-coordinate)
    % -------------------
    total_rows = size(impact_zone, 1);
    impact_zone_full = [impact_zone; impact_zone];
    impact_zone_full(total_rows+1:end, 2) = -impact_zone_full(total_rows+1:end, 2); % Mirror y-coordinate
    impact_times_full = [impact_times; impact_times];

    % Transform impact points back to original coordinates
    impact_zone_full_rotated = (R') * impact_zone_full'; % 3xN matrix
    if useGPU
        impact_zone_full_rotated = gather(impact_zone_full_rotated'); % Nx3 matrix
    else
        impact_zone_full_rotated = impact_zone_full_rotated'; % Nx3 matrix
    end

    % -------------------
    % Initialize Progress Tracking for Rotation
    % -------------------
    total_rotations = num_combinations * 2; % Ensure 'num_combinations' is defined before this point

    % Create a DataQueue for receiving progress updates
    rotationQueue = parallel.pool.DataQueue;

    % Initialize a counter for completed rotations
    completed_rotations = 0;
    lastPercent_rotations = 0;

    % Initialize the waitbar for rotations
    hWaitbar_rotations = waitbar(0, 'Rotating Trajectories...', 'Name', 'Rotation Progress');

    % Define the callback function to update the waitbar for rotations
    afterEach(rotationQueue, @updateWaitbarRotations);

    % Define the nested callback function for rotations
    function updateWaitbarRotations(~)
        completed_rotations = completed_rotations + 1;
        currentPercent_rotations = floor((completed_rotations / total_rotations) * 100);
        if currentPercent_rotations > lastPercent_rotations
            lastPercent_rotations = currentPercent_rotations;
            waitbar(currentPercent_rotations/100, hWaitbar_rotations, sprintf('Rotation Progress: %d%%', currentPercent_rotations));
        end
    end

    % Rotate trajectories back to original coordinates
    positions_full = cell(num_combinations * 2,1);
    velocities_full = cell(num_combinations * 2,1);
    parfor idx = 1:num_combinations * 2
        if idx <= num_combinations
            pos_rot = positions{idx};
            vel_rot = velocities{idx};
        else
            original_idx = idx - num_combinations;
            pos_rot = positions{original_idx};
            vel_rot = velocities{original_idx};
            pos_rot(:,2) = -pos_rot(:,2); % Mirror y-coordinate
            vel_rot(:,2) = -vel_rot(:,2); % Mirror y-velocity
        end
        % Rotate back to original frame
        pos_orig = (R') * pos_rot';
        vel_orig = (R') * vel_rot';
        if useGPU
            positions_full{idx} = gather(pos_orig)';
            velocities_full{idx} = gather(vel_orig)';
        else
            positions_full{idx} = pos_orig';
            velocities_full{idx} = vel_orig';
        end
        % Send progress update for rotations
        send(rotationQueue, 1);
    end

    % Close the waitbar for rotations after completion
    close(hWaitbar_rotations);

    % Remove impact points with NaN values
    valid_indices = ~isnan(impact_times_full);
    impact_zone_final = impact_zone_full_rotated(valid_indices, 1:2); % Only [x, y] coordinates
    impact_times_final = impact_times_full(valid_indices);

    % Filter trajectories based on valid indices
    positions = positions_full(valid_indices);
    velocities = velocities_full(valid_indices);
    impact_zone = impact_zone_final;
    impact_times = impact_times_final;
end