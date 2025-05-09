outputLocation = fullfile(pwd, 'Output Files','Report Graphics','SRBM'); % Directory to save output files
az0 = 0;                 % Initial azimuth angle in radians 
range0 = 60000;          % Initial range in meters
alt0 = 24000;            % Initial altitude in meters
mach2 = 2 * 343; % Speed of sound at sea level is approximately 343 m/s
vx0 = -mach2 * (range0 / sqrt(range0^2 + alt0^2)); % Initial velocity in x-direction (m/s)
vz0 = -mach2 * (alt0 / sqrt(range0^2 + alt0^2));   % Initial velocity in z-direction (m/s)
vy0 = 0;                 % Initial velocity in y-direction (m/s)
maxSpeed = 3;            % Maximum speed in mach
graphicSettings = struct(...
    'horiz_calc', 90/2, ...   % Example value for horizontal calculations
    'vert_calc', 22, ...    % Example value for vertical calculations
    'acc_calc', 2, ...      % Example value for acceleration calculations
    'dt', 2.5, ...            % Time step in seconds
    'max_iter', 120 ...     % Maximum number of iterations
);
[~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSettings, outputLocation, maxSpeed);