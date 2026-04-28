clc; clear; close all;

%% ================= PARAMETERS =================
v0 = -0.3;
Rc = 0.5;
eps_gain = 5;

D_alpha = 0.1;
D_beta  = 0.1;

dt = 0.01;
T  = 200;
N  = T/dt;

sgn = @(x) sign(x + 1e-6);
rng(1000)
%% ================= TORUS =================
R_major = 9;
r_minor = 3;

n_obs = 1500;
u_rand = 2*pi*rand(n_obs,1);
v_rand = 2*pi*rand(n_obs,1);

X_obs = (R_major + r_minor*cos(v_rand)) .* cos(u_rand);
Y_obs = (R_major + r_minor*cos(v_rand)) .* sin(u_rand);
Z_obs = r_minor * sin(v_rand);

obstacles = [X_obs Y_obs Z_obs];

%% smooth torus (visual)
[u,v] = meshgrid(linspace(0,2*pi,40), linspace(0,2*pi,30));
X_torus = (R_major + r_minor*cos(v)) .* cos(u);
Y_torus = (R_major + r_minor*cos(v)) .* sin(u);
Z_torus = r_minor * sin(v);

%% ================= APF =================
start = [R_major, 0, 0];
goal  = [-R_major/1.414, R_major/1.414, 0];

k_att = 1.0;
k_rep = 6.0;
d0    = 2.0;
eta   = 0.2;

pos = start;
waypoints = pos;

for i = 1:400
    
    F_att = -k_att * (pos - goal);
    
    F_rep = [0 0 0];
    for j = 1:size(obstacles,1)
        d = norm(pos - obstacles(j,:));
        if d < d0 && d > 1e-3
            F_rep = F_rep + k_rep*(1/d - 1/d0)*(1/d^2)*(pos - obstacles(j,:))/d;
        end
    end
    
    F = F_att + F_rep;
    
    if norm(F)>1e-6
        F = F/norm(F);
    end
    
    pos = pos + eta*F;
    waypoints = [waypoints; pos];
    
    if norm(pos - goal) < 0.3
        break;
    end
end

waypoints = smoothdata(waypoints,1,'movmean',11);
wp_id = 1;

%% ================= INITIAL STATE =================
start_pos = waypoints(1,:);

x = start_pos(1);
y = start_pos(2);
z = start_pos(3);

dx = waypoints(2,1) - waypoints(1,1);
dy = waypoints(2,2) - waypoints(1,2);
dz = waypoints(2,3) - waypoints(1,3);

alpha = -atan2(dy, dx);
beta  = atan2(dz, sqrt(dx^2 + dy^2));

%% ================= STORAGE =================
traj = [];
u_alpha_hist = [];
u_beta_hist  = [];
v_hist       = [];

K_alpha = 1;
K_beta  = 1;
v = v0;

%% ================= SIMULATION =================
for k = 1:N
    
    if wp_id >= size(waypoints,1)
        break;
    end
    
    target = waypoints(wp_id,:);
    
    dx = target(1) - x;
    dy = target(2) - y;
    dz = target(3) - z;
    
    R = norm([dx dy dz]);
    
    if R < Rc
        wp_id = wp_id + 1;
        if wp_id >= size(waypoints,1)
            break;
        end
        
        target = waypoints(wp_id,:);
        
        dx = target(1) - x;
        dy = target(2) - y;
        dz = target(3) - z;
        
        R0_seg = norm([dx dy dz]);
        Reff = max(R0_seg, Rc);
        
        K_alpha = D_alpha + abs(v)/Reff + eps_gain;
        K_beta  = D_beta  + 2*abs(v)/Reff + eps_gain;
        
        continue;
    end
    
    theta = atan2(dy, dx);
    phi   = atan2(dz, sqrt(dx^2 + dy^2));
    
    s1 = phi + beta;
    s2 = alpha - theta - pi;
    
    u_alpha = -K_alpha * sgn(s2);
    u_beta  = -K_beta  * sgn(s1);
    
    x_dot = v*cos(alpha)*cos(beta);
    y_dot = v*sin(alpha)*cos(beta);
    z_dot = v*sin(beta);
    
    x = x + dt*x_dot;
    y = y + dt*y_dot;
    z = z + dt*z_dot;
    
    alpha = alpha + dt*u_alpha;
    beta  = beta  + dt*u_beta;
    
    traj = [traj; x y z];
    u_alpha_hist = [u_alpha_hist; u_alpha];
    u_beta_hist  = [u_beta_hist; u_beta];
    v_hist       = [v_hist; v];
end

goal_pos = waypoints(end,:);

%% ================= VISUALIZATION =================
figure;
%axis equal;

plot3(traj(:,1),traj(:,2),traj(:,3),'r','LineWidth',2); hold on;
plot3(waypoints(:,1),waypoints(:,2),waypoints(:,3),'m--','LineWidth',3);

surf(X_torus,Y_torus,Z_torus,'FaceAlpha',0.15,'EdgeColor','none');
scatter3(obstacles(:,1),obstacles(:,2),obstacles(:,3),6,'k','filled');

plot3(start_pos(1),start_pos(2),start_pos(3),'go','LineWidth',9);
plot3(goal_pos(1),goal_pos(2),goal_pos(3),'bx','LineWidth',9);

xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; 
title('APF Adaptive SMC');
zlim([-10 10])

%% ================= CONTROL =================
figure;
subplot(2,1,1); plot(u_alpha_hist); title('u_alpha'); grid on;
subplot(2,1,2); plot(u_beta_hist); title('u_beta'); grid on;

figure;
plot(v_hist); title('Velocity'); grid on;

%% ================= METRICS =================

traj_mat = traj';

diffs = diff(traj_mat,1,2);
path_length = sum(vecnorm(diffs));

time_of_arrival = length(v_hist)*dt;

omega_norm = sqrt(u_alpha_hist.^2 + u_beta_hist.^2);

curvature = omega_norm ./ max(abs(v_hist),1e-6);

a_linear = diff(v_hist)/dt;
alpha_acc = diff(omega_norm)/dt;

fprintf('\n===== APF + SMC METRICS =====\n');
fprintf('Path length        : %.4f\n', path_length);
fprintf('Time of arrival    : %.4f\n', time_of_arrival);
fprintf('Mean curvature     : %.4f\n', mean(curvature));
fprintf('Max curvature      : %.4f\n', max(curvature));

fprintf('\nLinear accel min   : %.4f\n', min(a_linear));
fprintf('Linear accel max   : %.4f\n', max(a_linear));

fprintf('\nAngular accel min  : %.4f\n', min(alpha_acc));
fprintf('Angular accel max  : %.4f\n', max(alpha_acc));

%% ================= EXTRA PLOTS =================
figure;
plot(curvature,'LineWidth',1.5); title('Curvature'); grid on;

figure;
plot(a_linear,'LineWidth',1.5); title('Linear Acceleration'); grid on;

figure;
plot(alpha_acc,'LineWidth',1.5); title('Angular Acceleration'); grid on;