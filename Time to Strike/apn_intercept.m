function [intercept_trajectories, total_delta_v, closest_distances] = apn_intercept(launcher_pos, launcher_vel, positions_within_grid, velocities_within_grid, N, K1, K2, dt)
    % apn_intercept Implements Augmented Proportional Navigation for Missile Interception
    %
    % Parameters:
    %   launcher_pos          - Initial position of the interceptor [x; y; z] (m)
    %   launcher_vel          - Initial velocity of the interceptor [vx; vy; vz] (m/s)
    %   positions_within_grid - Cell array of [T_i x 3] incoming missile positions
    %   velocities_within_grid- Cell array of [T_i x 3] incoming missile velocities
    %   N                     - Navigation constant
    %   K1, K2                - Augmentation constants for APN
    %   dt                    - Time step (s)
    %
    % Outputs:
    %   intercept_trajectories - Cell array of interceptor trajectories [T x 3]
    %   total_delta_v          - Vector of total delta_v for each maneuver (m/s)
    %   closest_distances      - Vector of closest distances for each maneuver (m)

    num_targets = length(positions_within_grid);
    intercept_trajectories = cell(num_targets,1);
    total_delta_v = zeros(num_targets,1);
    closest_distances = zeros(num_targets,1);
    distance_history = cell(num_targets,1);

    for idx = 1:num_targets
        target_pos = positions_within_grid{idx}.';   % Transpose to [3 x T]
        target_vel = velocities_within_grid{idx}.';  % Transpose to [3 x T]
        T_steps = size(target_pos, 2);

        % Initialize interceptor state
        interceptor_pos = launcher_pos;
        interceptor_vel = launcher_vel;
        delta_v = 0;

        % Initialize trajectory with launcher position
        trajectory = interceptor_pos.'; 

        % Initialize previous LOS and LOS_rate
        prev_LOS = [0; 0; 0];
        prev_LOS_rate = [0; 0; 0];

        % Initialize distance tracking
        prev_distance = inf;  % Initialize with infinity
        closest_distance = inf;  % To store the minimal distance
        dist = zeros(T_steps, 1);

        for t = 1:T_steps
            % Current target state
            tgt_p = target_pos(:, t);
            tgt_v = target_vel(:, t);

            % Relative position and velocity
            rel_pos = tgt_p - interceptor_pos;
            rel_vel = tgt_v - interceptor_vel;

            % Calculate current distance
            current_distance = norm(rel_pos);
            dist(t) = current_distance;

            % Check if the current distance is the closest so far
            if current_distance < closest_distance
                closest_distance = current_distance;
            end

            % Normalize relative position to get LOS unit vector
            if current_distance == 0
                LOS = [0; 0; 0];
            else
                LOS = rel_pos / current_distance;
            end

            % Rate of change of LOS
            if t > 1
                LOS_rate = (LOS - prev_LOS) / dt;
            else
                LOS_rate = [0; 0; 0];
            end

            % Closing velocity
            closing_vel = -dot(rel_vel, LOS);

            % APN Guidance Law
            % Ensure LOS_rate is a column vector
            LOS_rate = LOS_rate(:);

            a_cmd = N * closing_vel * norm(LOS_rate) * LOS + K1 * LOS_rate + K2 * (LOS_rate - prev_LOS_rate) / dt;

            % Update interceptor acceleration
            interceptor_acc = a_cmd;

            % Update velocity and position
            interceptor_vel = interceptor_vel + interceptor_acc * dt;
            interceptor_pos = interceptor_pos + interceptor_vel * dt;

            % Accumulate delta_v
            delta_v = delta_v + norm(interceptor_acc * dt);

            % Store trajectory
            trajectory = [trajectory; interceptor_pos.']; %#ok<AGROW>

            % Update previous LOS and LOS_rate for next iteration
            prev_LOS = LOS;
            prev_LOS_rate = LOS_rate;

            % Update previous distance
            prev_distance = current_distance;
        end

        % Store results for current target
        intercept_trajectories{idx} = trajectory;
        total_delta_v(idx) = delta_v;
        closest_distances(idx) = closest_distance;
        distance_history{idx} = dist;

        % Trim trajectory to the point of closest approach
        [~, closest_distance_idx] = min(distance_history{idx});
        intercept_trajectories{idx} = intercept_trajectories{idx}(1:closest_distance_idx, :);
    end
end