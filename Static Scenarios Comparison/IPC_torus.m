clc; close all; clear;

%% ================= PARAMETERS =================
N_steps = 9000;
dt = 0.01;

% Tuned gains (MATCH SMC / MPC)
K_1 = 0.3;   % linear velocity scale
K_2 = 6;     % angular velocity scale

%% ================= TORUS =================
R0 = 9;
r0 = 3;

rng(1000)
N_pts = 1500;

u_pts = 2*pi*rand(N_pts,1);
v_pts = 2*pi*rand(N_pts,1);

X_obs = (R0 + r0*cos(v_pts)) .* cos(u_pts);
Y_obs = (R0 + r0*cos(v_pts)) .* sin(u_pts);
Z_obs = r0 * sin(v_pts);

points = [X_obs Y_obs Z_obs];

%% ================= START & GOAL =================
start = [R0, 0, 0]';
goal  = [-R0/1.414, R0/1.414, 0]';

pos = start;
r_i = [1;0;0];

traj = pos;
r_hist = r_i;
v_hist = [];

tol = 0.3;

%% ================= VISUALIZATION =================
figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

[u,v] = meshgrid(linspace(0,2*pi,40), linspace(0,2*pi,20));
XT = (R0 + r0*cos(v)).*cos(u);
YT = (R0 + r0*cos(v)).*sin(u);
ZT = r0*sin(v);

surf(XT,YT,ZT,'FaceAlpha',0.1,'EdgeColor','none');
scatter3(points(:,1),points(:,2),points(:,3),3,'r','filled');

plot3(start(1),start(2),start(3),'go','LineWidth',2);
plot3(goal(1),goal(2),goal(3),'bx','LineWidth',2);

h_traj = plot3(pos(1),pos(2),pos(3),'b','LineWidth',3);
h_curr = plot3(pos(1),pos(2),pos(3),'ko','MarkerFaceColor','k');

title('Proposed Control Structure')

sphere_handles = [];

%% ================= MAIN LOOP =================
for k = 1:N_steps


    %% SAMPLE DIRECTIONS
    n_dirs = 30;
    tic
    centers = zeros(3,n_dirs);
    radii = zeros(1,n_dirs);

    for d = 1:n_dirs-1

        randomRow = points(randi(N_pts),:)';
        vec = randomRow - pos;

        if norm(vec) < 1e-6, continue; end

        direction_vec = vec / norm(vec);

        bounds = [];

        for i = 1:N_pts
            p_vec = points(i,:)' - pos;

            if norm(p_vec) < 1e-6, continue; end

            if dot(direction_vec, p_vec / norm(p_vec)) > 0
                val = norm(p_vec) / ...
                    (2 * dot(direction_vec, p_vec / norm(p_vec)));
                bounds = [bounds val];
            end
        end

        if isempty(bounds), continue; end

        r_opt = min(bounds);
        centers(:,d) = r_opt * direction_vec;
        radii(d) = r_opt;
    end

    %% GOAL DIRECTION
    goal_vec = goal - pos;

    if norm(goal_vec) > 1e-6
        dist_goal = norm(goal_vec);
        dir = goal_vec / dist_goal;

        bounds = [];

        for i = 1:N_pts
            p_vec = points(i,:)' - pos;

            if norm(p_vec) < 1e-6, continue; end

            if dot(dir, p_vec / norm(p_vec)) > 0
                val = norm(p_vec) / ...
                    (2 * dot(dir, p_vec / norm(p_vec)));
                bounds = [bounds val];
            end
        end

        if ~isempty(bounds)
            r_opt = min(bounds);

            if r_opt >= dist_goal
                centers(:,n_dirs) = goal_vec;
                radii(n_dirs) = dist_goal;
            else
                centers(:,n_dirs) = r_opt * dir;
                radii(n_dirs) = r_opt;
            end
        end
    end

    %% SELECT BEST
    best_idx = 1; best_dist = inf;

    for d = 1:n_dirs
        if radii(d)==0, continue; end

        candidate_center = pos + centers(:,d);
        dist = norm(goal - candidate_center) - 0.2*radii(d);

        if dist < best_dist
            best_dist = dist;
            best_idx = d;
        end
    end

    center = pos + centers(:,best_idx);
    radius = radii(best_idx);
    ipc_time = toc;
    %% CONTROLLER
    R = pos - center;

    if norm(R) < 1e-6, continue; end

    R_hat = R / norm(R);
    S = dot(R_hat, r_i) + 1;

    % Smooth control (MATCH SMC scale)
    u = -K_1 * tanh(norm(R)) * tanh(5*dot(R, r_i));

    cross_term = cross(r_i, R_hat);
    cross_norm = norm(cross_term);

    if cross_norm < 1e-6
        omega = [0;0;0];
    else
        omega = (-K_2 * tanh(S) ...
            - u * cross_norm) * cross_term / cross_norm;
    end

    %% UPDATE
    pos = pos + u * r_i * dt;
    r_i = r_i + cross(omega, r_i) * dt;
    r_i = r_i / norm(r_i);

    traj(:,end+1) = pos;
    r_hist(:,end+1) = r_i;
    v_hist(end+1) = u;

    %% LIVE VIS
    if mod(k,20)==0

        set(h_traj,'XData',traj(1,:),...
                   'YData',traj(2,:),...
                   'ZData',traj(3,:));

        set(h_curr,'XData',pos(1),...
                   'YData',pos(2),...
                   'ZData',pos(3));

        drawnow limitrate;
    end

    if norm(pos-goal) < tol
        disp('Goal reached!');
        break;
    end
end

%% ================= METRICS =================

diffs = diff(traj,1,2);
path_length = sum(vecnorm(diffs));

time_of_arrival = length(v_hist)*dt;

omega_vec = diff(r_hist,1,2)/dt;
omega_norm = vecnorm(omega_vec);

curvature = omega_norm ./ max(abs(v_hist(1:end)),1e-6);

a_linear = diff(v_hist)/dt;
alpha = diff(omega_norm)/dt;

fprintf('\n===== SPHERE NAV (MATCHED) =====\n');
fprintf('Path length        : %.4f\n', path_length);
fprintf('Time               : %.4f\n', time_of_arrival);
fprintf('Mean curvature     : %.4f\n', mean(curvature));
fprintf('Max curvature      : %.4f\n', max(curvature));

fprintf('\nLinear accel min   : %.4f\n', min(a_linear));
fprintf('Linear accel max   : %.4f\n', max(a_linear));

fprintf('\nAngular accel min  : %.4f\n', min(alpha));
fprintf('Angular accel max  : %.4f\n', max(alpha));

%% ================= PLOTS =================
figure; plot(v_hist); title('Velocity'); grid on;
figure; plot(omega_norm); title('Angular velocity'); grid on;
figure; plot(curvature); title('Curvature'); grid on;
figure; plot(a_linear); title('Linear acceleration'); grid on;
figure; plot(alpha); title('Angular acceleration'); grid on;