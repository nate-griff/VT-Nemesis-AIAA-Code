clc, clear, close all
%#ok<*NASGU>
%#ok<*GVMIS>

r_d = .5; % defense radius (mi) (min before any cutoff issues at 2s is 1.2683)
r_d = r_d * 1609.34; % defense radius (m)
h_d = 30000; % defense altitude (ft)
h_d = h_d * 0.3048; % defense altitude (m)
g = 9.81; % Acceleration due to gravity (m/s^2)

M = 0.5:0.25:3; % Mach number

h_incoming = 100:100:30000; % incoming missile altitude (ft)
h_incoming = h_incoming * 0.3048; % incoming missile altitude (m)

cutoff_time = 2; % Minimum time needed for intercept (s)

t_out_cent = zeros(length(M), length(h_incoming));
descent_trajectory_x = cell(length(M), length(h_incoming));
descent_trajectory_y = cell(length(M), length(h_incoming));

for j = 1:length(M)
    for i = 1:length(h_incoming)
        t_out_cent(j,i) = descentTime_linear_fromEdge(h_incoming(i), r_d, M(j));
        [x, y] = descentTrajCalc(h_incoming(i), r_d, t_out_cent(j,i));
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
yline_value = linspace(min(x_min), max(x_min), 100); % Adjust the range as needed
zline_value = linspace(min(y_min), max(y_min), 100); % Adjust the range as needed
plot3(repmat(cutoff_time, size(yline_value)), yline_value, zline_value, '--r', 'LineWidth', 1.5, 'DisplayName', sprintf('cutoff time (%.1f s)', cutoff_time));
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
view(3)
yline_value = linspace(min(x_max), max(x_max), 100); % Adjust the range as needed
zline_value = linspace(min(y_max), max(y_max), 100); % Adjust the range as needed
plot3(repmat(cutoff_time, size(yline_value)), yline_value, zline_value, '--r', 'LineWidth', 1.5, 'DisplayName', sprintf('cutoff time (%.1f s)', cutoff_time));
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
xline(cutoff_time, '--r', 'LineWidth', 1.5,'Label', sprintf('cutoff time (%.1f s)', cutoff_time));

% Create a table to show the maximum, minimum, and middle time to reach the target
% from the minimum, middle, and maximum height for each Mach number
max_min_mid_table = array2table(zeros(length(M), 4), ...
    'VariableNames', {'Mach', 'TimeToHitFromMinAlt', 'TimeToHitFromMidAlt',  'TimeToHitFromMaxAlt'});

for j = 1:length(M)
    % Minimum height
    min_time = t_out_cent(j, 1);
    % Middle height
    mid_idx = ceil(length(h_incoming) / 2);
    mid_time = t_out_cent(j, mid_idx);
    % Maximum height
    max_time = t_out_cent(j, end);
    % Populate the table
    max_min_mid_table(j, :) = {M(j), min_time, mid_time, max_time};
end

% Create a figure and uitable
figure;
uitable('Data', max_min_mid_table{:,:}, 'ColumnName', {'Mach', 'Time to Hit from Min Alt (s)', 'Time to Hit from Mid Alt (s)', 'Time to Hit from Max Alt (s)'}, ...
    'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

% Preallocate mach_cutoff_table with a reasonable size
mach_cutoff_data = cell(length(M), 2);

for j = 1:length(M)
    max_altitude = -inf;
    for i = 1:length(h_incoming)
        if t_out_cent(j, i) < cutoff_time && h_incoming(i) > max_altitude
            max_altitude = h_incoming(i);
        end
    end
    if max_altitude > -inf
        mach_cutoff_data(j, :) = {M(j), max_altitude * 3.28084}; % Convert altitude to ft
    else
        mach_cutoff_data(j, :) = {M(j), NaN}; % No valid altitude found
    end
end

% Remove rows with NaN values
mach_cutoff_data = mach_cutoff_data(~cellfun(@(x) any(isnan(x)), mach_cutoff_data(:, 2)), :);

% Convert to table
mach_cutoff_table = cell2table(mach_cutoff_data, 'VariableNames', {'Mach', 'MaxInitialAltitude'});
% Display the table
figure;
uitable('Data', mach_cutoff_table{:,:}, 'ColumnName', {'Mach', 'Max Initial Altitude (ft)'}, ...
    'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

function t_out_cent = descentTime_linear_fromEdge(h_incoming, r_d, M)
    % t_out_cent: time to reach the center of the defense zone
    % h_incoming: incoming missile altitude (m) (can be vector)
    % r_d: defense radius (m)
    % M: Mach number
    % V: velocity (m/s)
    gamma = 1.4; % Ratio of specific heats
    R = 287; % Specific gas constant (J/kg-K)
    T0 = 288.15; % Standard temperature (K)
    a0 = sqrt(gamma * R * T0); % Speed of sound
    V = M * a0; % Velocity (m/s)

    theta = atan(h_incoming / r_d);
    distance = sqrt(r_d^2 + h_incoming.^2);
    t_out_cent = distance / V; % time to reach the center of the defense zone (s)
end

function [x, y] = descentTrajCalc(h_incoming, r_d, t_out_cent)
    % h_incoming: incoming missile altitude (m)
    % r_d: defense radius (m)
    % t_out_cent: time to reach the center of the defense zone (s)
    % descent_trajectory: trajectory of the missile during descent
    % x: distance from the center of the defense zone (m)
    % y: altitude (m)

    t = linspace(0, t_out_cent, 1000);
    x = r_d - (r_d / t_out_cent) .* t;
    x(x < 0) = 0;
    y = h_incoming - (h_incoming / t_out_cent) .* t;
    y(y < 0) = 0;
end