clc; close all; clear;

%% ================= PARAMETERS =================
N_steps = 20000;
dt = 0.01;

H = 1;

lambda_w = 0.001;
lambda_orient = 2.0;

u_max = 0.3;
w_max = 6;

d_min = 0.3;
K_near = 1500;

%% ================= TORUS =================
R_major = 9;
r_minor = 3;

rng(1000)
n_obs = 1500;

u_rand = 2*pi*rand(n_obs,1);
v_rand = 2*pi*rand(n_obs,1);

X_obs = (R_major + r_minor*cos(v_rand)) .* cos(u_rand);
Y_obs = (R_major + r_minor*cos(v_rand)) .* sin(u_rand);
Z_obs = r_minor * sin(v_rand);

obstacles = [X_obs Y_obs Z_obs];

%% Torus surface
[u,v] = meshgrid(linspace(0,2*pi,40), linspace(0,2*pi,30));
X_torus = (R_major + r_minor*cos(v)) .* cos(u);
Y_torus = (R_major + r_minor*cos(v)) .* sin(u);
Z_torus = r_minor * sin(v);

%% ================= APF WAYPOINTS =================
start = [R_major, 0, 0]';
goal  = [-R_major/1.414, R_major/1.414, 0];
%goal  = [-3.6, R_major-1, 0]';

k_att = 1.0;
k_rep = 6.0;
d0    = 2.0;
eta   = 0.2;

pos_apf = start';
waypoints = pos_apf;

tic;
for i = 1:400
    F_att = -k_att * (pos_apf - goal);
    
    F_rep = [0 0 0];
    for j = 1:size(obstacles,1)
        d = norm(pos_apf - obstacles(j,:));
        if d < d0 && d > 1e-3
            F_rep = F_rep + k_rep*(1/d - 1/d0)*(1/d^2)*(pos_apf - obstacles(j,:))/d;
        end
    end
    
    F = F_att + F_rep;
    if norm(F)>1e-6
        F = F/norm(F);
    end
    
    pos_apf = pos_apf + eta*F;
    waypoints = [waypoints; pos_apf];
    
    if norm(pos_apf - goal') < 0.3
        break;
    end
end
apfmpc_time = toc;
waypoints = smoothdata(waypoints,1,'movmean',11);

wp_id = 1;
Rc = 0.1;

%% ================= INITIAL STATE =================
pos = start;
%start_pos = waypoints(1,:);
r_i = [1;0;0];

traj = pos;
v_hist = [];
omega_hist = [];

%% ================= LIVE VIS =================
figure; hold on; grid on;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

surf(X_torus,Y_torus,Z_torus,'FaceAlpha',0.1,'EdgeColor','none');
scatter3(obstacles(:,1),obstacles(:,2),obstacles(:,3),6,'k','filled');
scatter3(waypoints(:,1),waypoints(:,2),waypoints(:,3),6,'filled');

plot3(start(1),start(2),start(3),'go','LineWidth',9);
plot3(goal(1),goal(2),goal(3),'bx','LineWidth',9);

h_traj = plot3(pos(1),pos(2),pos(3),'r','LineWidth',6);
h_curr = plot3(pos(1),pos(2),pos(3),'ro','MarkerFaceColor','r');
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('APF MPC')

%% ================= MAIN LOOP =================


for k = 1:N_steps
    
    if wp_id >= size(waypoints,1)
        break;
    end
    
    target = waypoints(wp_id,:)';
    
    if norm(pos - target) < Rc
        wp_id = min(wp_id+1, size(waypoints,1));
    end
    
    %% nearest obstacles
    dists = vecnorm(obstacles' - pos);
    [~,idx] = sort(dists);
    near_pts = obstacles(idx(1:K_near),:);
    
    %% MPC
    U0 = zeros(4,1);
    lb = [-u_max; -w_max; -w_max; -w_max];
    ub = [ u_max;  w_max;  w_max;  w_max];

    tic
    U_opt = fmincon(@(U) cost_fn(U,pos,r_i,target,dt,lambda_w,lambda_orient),...
        U0,[],[],[],[],lb,ub,...
        @(U) constraint_fn(U,pos,r_i,near_pts,d_min,dt),...
        optimoptions('fmincon','Display','off','MaxIterations',30));
    mpc_time = toc;
    u = U_opt(1);
    omega = U_opt(2:4);
    
    %% update
    pos = pos + u*r_i*dt;
    r_i = r_i + cross(omega,r_i)*dt;
    r_i = r_i / norm(r_i);
    
    traj(:,end+1) = pos;
    v_hist(end+1) = u;
    omega_hist(:,end+1) = omega;
    
    %% live plot
    if mod(k,20)==0
        set(h_traj,'XData',traj(1,:),'YData',traj(2,:),'ZData',traj(3,:));
        set(h_curr,'XData',pos(1),'YData',pos(2),'ZData',pos(3));
        drawnow limitrate;
    end
    
    if norm(pos-goal) < 0.3
        disp('Goal reached');
        break;
    end
end

%% ================= METRICS =================

traj = traj;

% Path length
diffs = diff(traj,1,2);
path_length = sum(vecnorm(diffs));

% Time
time_of_arrival = length(v_hist)*dt;

% Angular velocity
omega_norm = vecnorm(omega_hist);

% Curvature
curvature = omega_norm ./ max(abs(v_hist),1e-6);

% Linear acceleration
a_linear = diff(v_hist)/dt;

% Angular acceleration (norm)
alpha = diff(omega_norm)/dt;

fprintf('\n===== METRICS =====\n');
fprintf('Path length        : %.4f\n', path_length);
fprintf('Time of arrival    : %.4f\n', time_of_arrival);
fprintf('Mean curvature     : %.4f\n', mean(curvature));
fprintf('Max curvature      : %.4f\n', max(curvature));

fprintf('\nLinear accel min   : %.4f\n', min(a_linear));
fprintf('Linear accel max   : %.4f\n', max(a_linear));

fprintf('\nAngular accel min  : %.4f\n', min(alpha));
fprintf('Angular accel max  : %.4f\n', max(alpha));

%% ================= PLOTS =================
figure;
plot(v_hist); title('Linear velocity'); grid on;

figure;
plot(omega_norm); title('Angular velocity norm'); grid on;

figure;
plot(curvature); title('Curvature'); grid on;

figure;
plot(a_linear); title('Linear acceleration'); grid on;

figure;
plot(alpha); title('Angular acceleration'); grid on;

%% ================= FUNCTIONS =================
function J = cost_fn(U,pos,r_i,ref,dt,lambda_w,lambda_orient)

    x = pos;
    r = r_i;
    
    u = U(1);
    omega = U(2:4);
    
    x = x + u*r*dt;
    r = r + cross(omega,r)*dt;
    r = r / norm(r);
    
    ref_vec = ref - x;
    if norm(ref_vec)>1e-6
        d = ref_vec / norm(ref_vec);
    else
        d = [0;0;0];
    end
    
    orient_err = 1 - dot(r,d);
    
    J = norm(x-ref)^2 ...
        + lambda_w*norm(omega)^2 ...
        + lambda_orient*orient_err;
end

function [c,ceq] = constraint_fn(U,pos,r_i,near_pts,d_min,dt)

    x = pos;
    r = r_i;
    
    u = U(1);
    omega = U(2:4);
    
    x = x + u*r*dt;
    r = r + cross(omega,r)*dt;
    r = r / norm(r);
    
    c = [];
    for i = 1:size(near_pts,1)
        p = near_pts(i,:)';
        c = [c; d_min^2 - norm(x-p)^2];
    end
    
    ceq = [];
end