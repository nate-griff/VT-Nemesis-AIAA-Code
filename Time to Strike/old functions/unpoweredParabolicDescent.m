clc, clear, close all
%#ok<*NASGU>
%#ok<*GVMIS>

r_d = 5; % defense radius (mi)
r_d = r_d * 1609.34; % defense radius (m)
h_d = 30000; % defense altitude (ft)
h_d = h_d * 0.3048; % defense altitude (m)
g = 9.81; % Acceleration due to gravity (m/s^2)

M = 0.5:0.5:3; % Mach number

h_incoming = 100:100:30000; % incoming missile altitude (ft)
h_incoming = h_incoming * 0.3048; % incoming missile altitude (m)

cutoff_time = 2; % Minimum time needed for intercept (s)

t_out_cent = zeros(length(M), length(h_incoming));
descent_trajectory_x = cell(length(M), length(h_incoming));
descent_trajectory_y = cell(length(M), length(h_incoming));

for j = 1:length(M)
    for i = 1:length(h_incoming)
        t_out_cent(j,i) = descentTime_parabolic_fromEdge(h_incoming(i), r_d, M(j));
        [x, y] = descentTrajCalc(h_incoming(i), t_out_cent(j,i), M(j));
        descent_trajectory_x{j,i} = x;
        descent_trajectory_y{j,i} = y;
    end
end

% Plot descent trajectories for min and max heights for each Mach number
figure;

% Plot descent trajectory for min heights
subplot(2,1,1);
hold on;
for j = 1:length(M)
    % Minimum height
    x_min = descent_trajectory_x{j,1} * 3.28084;
    y_min = descent_trajectory_y{j,1} * 3.28084; % Convert to ft
    t_min = linspace(0, t_out_cent(j,1), length(x_min));
    plot3(t_min, x_min, y_min, '--', 'DisplayName', sprintf('Mach %.1f Min Height', M(j)));
end
ylabel('Range (ft)');
xlabel('Time (s)');
zlabel('Altitude (ft)');
title('Descent Trajectory for Min Heights per Mach Number');
legend('Location', 'best');
grid on;
view(3);
hold off;

% Plot descent trajectory for max heights
subplot(2,1,2);
hold on;
for j = 1:length(M)
    % Maximum height
    x_max = descent_trajectory_x{j,end} * 3.28084;
    y_max = descent_trajectory_y{j,end} * 3.28084; % Convert to ft
    t_max = linspace(0, t_out_cent(j,end), length(x_max));
    plot3(t_max, x_max, y_max, '-', 'DisplayName', sprintf('Mach %.1f Max Height', M(j)));
end
ylabel('Range (ft)');
xlabel('Time (s)');
zlabel('Altitude (ft)');
title('Descent Trajectory for Max Heights per Mach Number');
legend('Location', 'best');
grid on;
view(3);
hold off;

% Plot altutude vs time to center for each Mach number
figure;
plot(t_out_cent, h_incoming * 3.28084, 'LineWidth', 2)
colormap(jet)
ylabel('Incoming missile altitude (ft)')
xlabel('Time to reach the center of the defense zone (s)')
title('Missile Altitude vs. Time to Center for Various Mach Numbers')
legend(arrayfun(@(x) sprintf('Mach %.1f', x), M, 'UniformOutput', false), 'Location', 'southeast')
grid on

% Create a table to show the maximum, minimum, and middle time to reach the target
% from the minimum, middle, and maximum height for each Mach number
max_min_mid_table = array2table(zeros(length(M), 7), ...
    'VariableNames', {'Mach', 'TimeToHitFromMinAlt', 'MinAltInitRange', 'TimeToHitFromMidAlt', 'MidAltInitRange', 'TimeToHitFromMaxAlt', 'MaxAltInitRange'});

for j = 1:length(M)
    % Minimum height
    min_time = t_out_cent(j, 1);
    % Middle height
    mid_idx = ceil(length(h_incoming) / 2);
    mid_time = t_out_cent(j, mid_idx);
    % Maximum height
    max_time = t_out_cent(j, end);

    % Find the index where y first reaches 0
    idx_min = find(descent_trajectory_y{j, 1} == 0, 1);
    idx_mid = find(descent_trajectory_y{j, mid_idx} == 0, 1);
    idx_max = find(descent_trajectory_y{j, end} == 0, 1);

    % Get the corresponding x values
    min_x = descent_trajectory_x{j, 1}(idx_min) * 3.28084 / 5280; % Convert to miles
    mid_x = descent_trajectory_x{j, mid_idx}(idx_mid) * 3.28084 / 5280; % Convert to miles
    max_x = descent_trajectory_x{j, end}(idx_max) * 3.28084 / 5280; % Convert to miles

    % Populate the table
    max_min_mid_table(j, :) = {M(j), min_time, min_x, mid_time, mid_x, max_time, max_x};
end

% Create a figure and uitable
figure;
uitable('Data', max_min_mid_table{:,:}, 'ColumnName', {'Mach', 'Time to Hit from Min Alt (s)', 'Min Alt Initial Range (mi)', 'Time to Hit from Mid Alt (s)', 'Mid Alt Initial Range (mi)', 'Time to Hit from Max Alt (s)', 'Max Alt Initial Range (mi)'}, ...
    'Units', 'Normalized', 'Position', [0, 0, 1, 1]);



function t_out_cent = descentTime_parabolic_fromEdge(h_incoming, r_d, M)
    % t_out_cent: time to reach the center of the defense zone
    % h_incoming: incoming missile altitude (m) (can be vector)
    % r_d: defense radius (m)
    % M: Mach number
    % V: velocity (m/s)
    g = 9.81; % Acceleration due to gravity (m/s^2)
    gamma = 1.4; % Ratio of specific heats
    R = 287; % Specific gas constant (J/kg-K)
    T0 = 288.15; % Standard temperature (K)
    a0 = sqrt(gamma * R * T0); % Speed of sound
    V = M * a0; % Velocity (m/s)

    % Calculate the time to fall from h_incoming to the ground
    t_fall = sqrt(2 * h_incoming / g);

    % Calculate the horizontal distance traveled during the fall
    x_fall = V * t_fall;

    % Calculate the starting horizontal distance
    if x_fall > r_d
        % Missile starts outside the defense radius
        x_start = x_fall - r_d;
    else
        % Missile starts inside the defense radius
        x_start = r_d - x_fall;
    end

    % Calculate the time to travel from the starting point to the center
    t_travel = x_start / V;

    % Total time to reach the center of the defense zone
    t_out_cent = t_fall + t_travel;
end

function [x, y] = descentTrajCalc(h_incoming, t_out_cent, M)
    % h_incoming: incoming missile altitude (m)
    % r_d: defense radius (m)
    % t_out_cent: time to reach the center of the defense zone (s)
    % M: Mach number
    % x: distance from the center of the defense zone (m)
    % y: altitude (m)

    g = 9.81; % Acceleration due to gravity (m/s^2)
    gamma = 1.4; % Ratio of specific heats
    R = 287; % Specific gas constant (J/kg-K)
    T0 = 288.15; % Standard temperature (K)
    a0 = sqrt(gamma * R * T0); % Speed of sound
    V = M * a0; % Initial horizontal velocity (m/s)

    t = linspace(0, t_out_cent, 1000);
    x = V .* t; % Horizontal distance traveled
    y = h_incoming - 0.5 * g * t.^2; % Vertical distance traveled due to gravity

    % Ensure the missile does not go below ground level
    y(y < 0) = 0;
    % % Ensure the missile does not travel beyond the defense radius
    % x(x > r_d) = r_d;
end