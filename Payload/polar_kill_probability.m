function [P_f] = polar_kill_probability()
    % Constants and parameters
    C_D = 1.42;            % Drag coefficient
    k = 4.74;              % Shape Factor, g/cm^3
    ro_exp = 1.82;         % Explosive density, g/cm^3
    ro_case = 7.85;        % Case density, g/cm^3
    type = 1;              % Warhead Type, cylinder=1 & sphere=0
    L_cyl = 9.01;
    E = 2926;              % Gurney Constant, m/s
    B =  0.0531;           % Mott Constant, b^1/2in^-7/16
    D = 5.98;               % Case outer diameter, in
    d = 5.48;               % Case inner diameter, in
    V_exp = pi*((d/2)^2)*L_cyl*16.3871+2*pi*((d/2)^2)*16.3871;        % Volume of explosive, cm^3
    V_case = pi*((D/2)^2)*L_cyl*16.3871+2*pi*((D/2)^2)*16.3871-V_exp; % Volume of case, cm^3
    ro_air = 0.001293;     % Air density, g/cm^3
    E_cr = 200;            % Critical impact energy, J
    A_T = 0.292;            % Target area, m^2 3.716122
    r_target = [-10; 0]; % Target position, m -> Asymmetrical position
    v_target = [1029; 0]; % Target velocity, m/s

    % Setup grid (polar-style)
    theta_vals = linspace(0, 2*pi, 360);  % angles (radial directions)
    r_vals = linspace(1, 20, 20);         % radii from detonation center, m

    [Theta, R] = meshgrid(theta_vals, r_vals);
    X = R .* cos(Theta);     % x-coordinates of grid
    Y = R .* sin(Theta);     % y-coordinates of grid

    P_f = zeros(size(R));    % Matrix to hold kill probabilities

    % Constant fragment parameters
    m_exp = ro_exp*V_exp/1000  % Explosive mass, kg
    m_case = ro_case*V_case/1000  % Case mass, kg
    v_i = E*((m_case/m_exp)+0.5)^(-0.5);  % Initial fragment velocity, m/s
    L = (2*(k^(2/3))/C_D/ro_air)*100; %
    t = (D-d)/2;  % Case thickness, in
    Mk = B*(t^(5/16))*(d^(1/3))*(1+t/d);  % Mott Distribution
    N_t = (m_case*2.205)/2/Mk^2;  % Number of fragments
    Q0 = N_t/(4*pi);  % Pieces/m^2

    % Loop through all grid points to calculate kill probability
    for i = 1:numel(R)
        % Position of point relative to detonation center
        r_frag = [X(i); Y(i)];
        r_mag = norm(r_frag);  % Distance from center
        
        % Relative distance between fragment and target
        dist_to_target = norm(r_target-r_frag);  % Euclidean distance
        
        % Kill probability: the closer the fragment to the target, the higher the probability
        % Using an exponential decay based on distance to target
        distance_factor = exp(-dist_to_target / 10);  % 10 is a scaling factor for distance

        % Assuming fragment velocity at this position
        v_frag = v_i*exp(-r_mag/L);  % Fragment velocity (exponentially decaying with distance)
        
        % Calculate relative velocity magnitude between fragment and target
        v_rel = v_frag-v_target;  % Relative velocity vector
        v_rel_mag = norm(v_rel);    % Magnitude of the relative velocity
        
        % Critical mass threshold
        M_cr = 2*E_cr/(v_rel_mag)^2;  % Critical mass from the relative velocity
        
        % Fragment density and kill probability calculation (adjusted for target position)
        q_cr = Q0/(r_mag^2)*exp(-sqrt(2*M_cr/m_case));  % Critical fragment density
        P_f(i) = 1-exp(-q_cr*A_T*distance_factor);  % Kill probability
    end

    % Kill probability heatmap
    figure;
    contourf(X, Y, P_f, 100, 'LineColor', 'none');
    colorbar;
    title('Kill Probability Map');
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    axis equal;
    ax = gca;
    ax.FontSize = 16;
    ax.FontName = 'Times New Roman'
    % Add the target position and velocity to the plot
    hold on;
    % Target velocity vector
    quiver(r_target(1), r_target(2), v_target(1), v_target(2), 0.004, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);  
    % Target position
    plot(r_target(1), r_target(2), 'ro', 'MarkerSize', 6, 'LineWidth', 2);
    hold on;

    % Plot a contour line at kill probability = 0.7 (dotted)
    [C, h] = contour(X, Y, P_f, [0.7 0.7], 'k:');  % 'k:' for black dotted line
    h.LineWidth = 2;
    
     guidance_radius = 0.0375; % meters
    theta_circle = linspace(0, 2*pi, 100);
    x_circle = r_target(1) + guidance_radius * cos(theta_circle);
    y_circle = r_target(2) + guidance_radius * sin(theta_circle);
    plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5); % black dashed circle
    legend('Kill Probability', 'Target Velocity (Direction)', 'Target Position', 'Region With p_{k} of 0.7 or greater', 'Guidance Error');
end
