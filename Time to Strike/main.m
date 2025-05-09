% Define the graphicsSettings struct with required fields
lowResGraphicsSettings = struct(...
    'horiz_calc', 90, ...   % Example value for horizontal calculations
    'vert_calc', 45, ...    % Example value for vertical calculations
    'acc_calc', 4, ...      % Example value for acceleration calculations
    'dt', 1, ...            % Time step in seconds
    'max_iter', 180 ...     % Maximum number of iterations
);

desiredTime = 85; % Desired time in seconds
dtHighRes = 0.1; % Time step in seconds for high resolution calculations
highResGraphicsSettings = struct(...
    'horiz_calc', 360/2, ...   % Example value for horizontal calculations
    'vert_calc', 180/2, ...    % Example value for vertical calculations
    'acc_calc', 15, ...      % Example value for acceleration calculations
    'dt', dtHighRes, ...     % Time step in seconds
    'max_iter', desiredTime/0.5 ...    % Maximum number of iterations
);

graphicSetting = highResGraphicsSettings; % Set the graphic setting to low resolution

% Define the missile parameters

% % YJ-18 missile in configuration of STK low alt cruise missile
% outputLocation = fullfile(pwd, 'Output Files','Report Graphics','YJ-18'); % Directory to save output files
% az0 = 0;                 % Initial azimuth angle in radians 
% range0 = 33000;          % Initial range in meters
% alt0 = 804;              % Initial altitude in meters
% vx0 = -300*3;            % Initial velocity in x-direction (m/s)
% vy0 = 0;                 % Initial velocity in y-direction (m/s)
% vz0 = 0;                 % Initial velocity in z-direction (m/s)
% maxSpeed = 3;            % Maximum speed in mach

% desiredTime = 55; % Desired time in seconds
% dtHighRes = 0.1; % Time step in seconds for high resolution calculations
% highResGraphicsSettings.dt = dtHighRes; % Set the time step for high resolution calculations
% highResGraphicsSettings.max_iter = desiredTime/dtHighRes; % Set the maximum number of iterations for high resolution calculations
% graphicSetting = highResGraphicsSettings; % Set the graphic setting to high resolution

% [~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSetting, outputLocation, maxSpeed);


% % Onyx missile low alt Mach 3
% outputLocation = fullfile(pwd, 'Output Files','Report Graphics','Onyx'); % Directory to save output files
% az0 = 0;                 % Initial azimuth angle in radians 
% range0 = 804;            % Initial range in meters
% alt0 = 15;               % Initial altitude in meters
% vx0 = -3*343;            % Initial velocity in x-direction (m/s)
% vy0 = 0;                 % Initial velocity in y-direction (m/s)
% vz0 = 0;                 % Initial velocity in z-direction (m/s)
% maxSpeed = 3;            % Maximum speed in mach
% 
% 
% desiredTime = 20; % Desired time in seconds
% dtHighRes = .5; % Time step in seconds for high resolution calculations
% highResGraphicsSettings.dt = dtHighRes; % Set the time step for high resolution calculations
% highResGraphicsSettings.max_iter = desiredTime/dtHighRes; % Set the maximum number of iterations for high resolution calculations
% graphicSetting = highResGraphicsSettings; % Set the graphic setting to high resolution
% 
% [~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSetting, outputLocation, maxSpeed);

%KH-102 max speed mach 1 @ max alt
outputLocation = fullfile(pwd, 'DataExportsHighRes','KH-102 Modified'); % Directory to save output files
az0 = 0;                 % Initial azimuth angle in radians 
range0 = 16000;          % Initial range in meters
alt0 = 10000;            % Initial altitude in meters
vx0 = -343*.78;          % Initial velocity in x-direction (m/s)
vy0 = 0;                 % Initial velocity in y-direction (m/s)
vz0 = 0;                 % Initial velocity in z-direction (m/s)
maxSpeed = 1.6;            % Maximum speed in mach

graphicSetting = lowResGraphicsSettings; % Set the graphic setting to low resolution

[~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSetting, outputLocation, maxSpeed);


% SRBM missile in configuration of STK from Lynchburg
% outputLocation = fullfile(pwd, 'Output Files','Report Graphics','SRBM'); % Directory to save output files
% az0 = 0;                 % Initial azimuth angle in radians 
% range0 = 60000;          % Initial range in meters
% alt0 = 24000;            % Initial altitude in meters
% mach2 = 2 * 343; % Speed of sound at sea level is approximately 343 m/s
% vx0 = -mach2 * (range0 / sqrt(range0^2 + alt0^2)); % Initial velocity in x-direction (m/s)
% vz0 = -mach2 * (alt0 / sqrt(range0^2 + alt0^2));   % Initial velocity in z-direction (m/s)
% vy0 = 0;                 % Initial velocity in y-direction (m/s)
% maxSpeed = 2.5;            % Maximum speed in mach

% desiredTime = 120; % Desired time in seconds
% dtHighRes = 1; % Time step in seconds for high resolution calculations
% highResGraphicsSettings.dt = dtHighRes; % Set the time step for high resolution calculations
% highResGraphicsSettings.max_iter = desiredTime/dtHighRes; % Set the maximum number of iterations for high resolution calculationsgraphicSetting = highResGraphicsSettings; % Set the graphic setting to high resolution
% % graphicSetting = highResGraphicsSettings;
% [~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSetting, outputLocation, maxSpeed);

%%
graphicSetting = lowResGraphicsSettings; % Set the graphic setting to low resolution

% Missile variation 1
outputLocation = fullfile(pwd, 'DataExportsHighRes','Variation1'); % Directory to save output files
az0 = pi/6;              % Initial azimuth angle in radians 
range0 = 1609.34;        % Initial range in meters (1 mile)
alt0 = 1524;             % Initial altitude in meters (5000 ft)
maxSpeed = 0.5;          % Maximum speed in mach
speed = 343 * maxSpeed;  % Speed in m/s
vx0 = -speed * (range0 / sqrt(range0^2 + alt0^2 + range0^2)); % Initial velocity in x-direction (m/s)
vy0 = -speed * (range0 / sqrt(range0^2 + alt0^2 + range0^2)); % Initial velocity in y-direction (m/s)
vz0 = -speed * (alt0 / sqrt(range0^2 + alt0^2 + range0^2));   % Initial velocity in z-direction (m/s)

[~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSetting, outputLocation, maxSpeed);

% Missile variation 2
outputLocation = fullfile(pwd, 'DataExportsHighRes','Variation2'); % Directory to save output files
az0 = pi/4;              % Initial azimuth angle in radians 
range0 = 16093.4;        % Initial range in meters (10 miles)
alt0 = 9144;             % Initial altitude in meters (30000 ft)
maxSpeed = 1.5;          % Maximum speed in mach
speed = 343 * maxSpeed;  % Speed in m/s
vx0 = -speed * (range0 / sqrt(range0^2 + alt0^2 + range0^2)); % Initial velocity in x-direction (m/s)
vy0 = -speed * (range0 / sqrt(range0^2 + alt0^2 + range0^2)); % Initial velocity in y-direction (m/s)
vz0 = -speed * (alt0 / sqrt(range0^2 + alt0^2 + range0^2));   % Initial velocity in z-direction (m/s)

[~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSetting, outputLocation, maxSpeed);

% Missile variation 3
outputLocation = fullfile(pwd, 'DataExportsHighRes','Variation3'); % Directory to save output files
az0 = pi/3;              % Initial azimuth angle in radians 
range0 = 32186.9;        % Initial range in meters (20 miles)
alt0 = 6096;             % Initial altitude in meters (20000 ft)
maxSpeed = 2.5;          % Maximum speed in mach
speed = 343 * maxSpeed;  % Speed in m/s
vx0 = -speed * (range0 / sqrt(range0^2 + alt0^2 + range0^2)); % Initial velocity in x-direction (m/s)
vy0 = -speed * (range0 / sqrt(range0^2 + alt0^2 + range0^2)); % Initial velocity in y-direction (m/s)
vz0 = -speed * (alt0 / sqrt(range0^2 + alt0^2 + range0^2));   % Initial velocity in z-direction (m/s)



[~,~,~,~] = incomingMissileProfile(az0, range0, alt0, vx0, vy0, vz0, graphicSetting, outputLocation, maxSpeed);