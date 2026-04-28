clc; close all; clear;
rng(1);

%% PARAMETERS
N_steps = 1300;
del_T = 0.01;
tol = 0.1;

K_1 = 1;
K_2 = 1;

R0 = 2.0;
r0 = 0.5;

%% EXPERIMENT SETUP
N_mix = 5;
all_traj_cases = cell(N_mix+1,1);
results = struct();

%% START & GOAL
theta = 2*pi*rand;
theta_goal = theta + pi/2;

v_start = (rand-0.5)*pi;
v_goal  = (rand-0.5)*pi;

pos_init = [
    (R0 + r0*cos(v_start))*cos(theta);
    (R0 + r0*cos(v_start))*sin(theta);
    r0*sin(v_start)
];

goal = [
    (R0 + r0*cos(v_goal))*cos(theta_goal);
    (R0 + r0*cos(v_goal))*sin(theta_goal);
    r0*sin(v_goal)
];

%% FIXED TORUS OBSTACLES (SHARED)
N_pts = 270;
u_pts = 2*pi*rand(N_pts,1);
v_pts = 2*pi*rand(N_pts,1);
points_global = torus_from_uv(u_pts,v_pts,R0,r0);

%% BASE NOISE
noise_base.sigma_g = 0.02;
noise_base.sigma_u = 0.02;
noise_base.lambda  = 10;
noise_base.p_imp   = 0.01;
noise_base.imp_mag = 0.5;

%% =========================
%% BASELINE (NO NOISE)
%% =========================
disp('Baseline');
traj = simulate_and_visualize(pos_init, goal, noise_base, true, ...
    'Baseline', 'baseline.gif', points_global);

all_traj_cases{1} = traj;

%% =========================
%% MIXTURE SWEEP
%% =========================
for c = 1:N_mix

    a = rand(4,1);
    a = a / sum(a);

    noise = noise_base;
    noise.a1 = a(1);
    noise.a2 = a(2);
    noise.a3 = a(3);
    noise.a4 = a(4);

    disp(['Mix ', num2str(c), ' -> ', num2str(a','%.2f ')]);

    rng(c); % reproducible noise per run

    traj = simulate_and_visualize(pos_init, goal, noise, false, ...
        sprintf('Mix %d',c), sprintf('mix_%d.gif',c), points_global);

    all_traj_cases{c+1} = traj;
    results(c).a = a;
end

%% =========================
%% FINAL COMPARISON PLOT
%% =========================
figure; hold on; grid on;
title('Trajectory Comparison','FontSize',13);

colors = lines(N_mix+1);
legend_handles = [];

for c = 1:(N_mix+1)

    tr = all_traj_cases{c};

    if c == 1
        lw = 3;
    else
        lw = 1.5;
    end

    h = plot3(tr(1,:), tr(2,:), tr(3,:), ...
        'Color', colors(c,:), 'LineWidth', lw);

    legend_handles = [legend_handles h];
end

plot3(pos_init(1), pos_init(2), pos_init(3),'bo','MarkerFaceColor','b');
plot3(goal(1), goal(2), goal(3),'go','LineWidth',13);

xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; view(45,45);

legend_strings = cell(N_mix+1,1);
legend_strings{1} = 'No Noise';

for c = 1:N_mix
    a = results(c).a;
    legend_strings{c+1} = sprintf('a=[%.2f %.2f %.2f %.2f]', a);
end

legend(legend_handles, legend_strings);

%% =========================
%% SIMULATION FUNCTION
%% =========================
function traj = simulate_and_visualize(pos, goal, noise, no_noise, ...
    fig_title, gif_name, points)

    N_steps = 1300;
    del_T = 0.01;
    tol = 0.1;

    K_1 = 1;
    K_2 = 1;

    r_i = (goal - pos)/norm(goal - pos);
    traj = pos;

    fig = figure('Name',fig_title);
    hold on; grid on; axis equal;
    xlim([-3 3]); ylim([-3 3]); zlim([-2 2]);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    view(45,45);
    title(fig_title,'FontSize',13);

    % show fixed obstacles (context only)
    scatter3(points(:,1),points(:,2),points(:,3),10,[0.7 0.7 0.7],'filled');

    plot3(pos(1),pos(2),pos(3),'bo','MarkerFaceColor','b');
    plot3(goal(1),goal(2),goal(3),'g*');

    h_traj = plot3(pos(1),pos(2),pos(3),'b','LineWidth',2);

    [Xs,Ys,Zs] = sphere(12);
    frame_count = 1;

    camlight headlight
    lighting gouraud

    for k = 1:N_steps

        n_dirs = 25;
        centers = zeros(3,n_dirs);
        radii = zeros(1,n_dirs);

        for d = 1:n_dirs

            vec_true = points(randi(size(points,1)),:)' - pos;

            if no_noise
                vec = vec_true;
            else
                vec = apply_measurement_noise(vec_true, noise);
            end

            if norm(vec) < 1e-6, continue; end

            dir = vec / norm(vec);
            bounds = [];

            for i = 1:size(points,1)
                p = points(i,:)' - pos;

                if no_noise
                    p_vec = p;
                else
                    p_vec = apply_measurement_noise(p, noise);
                end

                if norm(p_vec)<1e-6, continue; end

                if dot(dir, p_vec/norm(p_vec)) > 0
                    val = norm(p_vec)/(2*dot(dir,p_vec/norm(p_vec)));
                    bounds = [bounds val];
                end
            end

            if isempty(bounds), continue; end

            r = min(bounds);
            centers(:,d) = r * dir;
            radii(d) = r;
        end

        [~,idx] = min(vecnorm(goal - (pos + centers)) - 0.2*radii);

        if radii(idx) > norm(goal-pos)
            center = goal;
            radius = norm(goal-pos);
        else
            center = pos + centers(:,idx);
            radius = radii(idx);
        end

        % 🔵 CONTROL SPHERE (persistent)
        surf( ...
            radius*Xs + center(1), ...
            radius*Ys + center(2), ...
            radius*Zs + center(3), ...
            'FaceAlpha',0.01, ...
            'EdgeColor','none', ...
            'FaceColor',[1 1 0], ...
            'AmbientStrength',0.9, ...
            'DiffuseStrength',0.1, ...
            'SpecularStrength',0.1);

        R = pos - center;
        if norm(R)<1e-6, continue; end

        R_hat = R/norm(R);
        S = dot(R_hat,r_i)+1;

        u = -K_1*tanh(norm(R))*sign(dot(R,r_i));

        cross_term = cross(r_i,R_hat);
        cn = norm(cross_term);

        if cn<1e-6
            omega = [0;0;0];
        else
            omega = (-K_2*sign(S)*abs(S)^0.5 - u*cn)*cross_term/cn;
        end

        pos = pos + u*r_i*del_T;
        r_i = r_i + cross(omega,r_i)*del_T;
        r_i = r_i/norm(r_i);

        traj(:,end+1) = pos;

        set(h_traj,'XData',traj(1,:), ...
                   'YData',traj(2,:), ...
                   'ZData',traj(3,:));

        drawnow;

        frame = getframe(fig);
        img = frame2im(frame);
        [imind,cm] = rgb2ind(img,256);

        if frame_count==1
            imwrite(imind,cm,gif_name,'gif','Loopcount',inf,'DelayTime',0.05);
        else
            imwrite(imind,cm,gif_name,'gif','WriteMode','append','DelayTime',0.05);
        end

        frame_count = frame_count + 1;

        if norm(pos-goal)<tol
            break;
        end
    end
end

%% NOISE MODEL
function rho_meas = apply_measurement_noise(rho, p)

    if norm(rho)<1e-9
        rho_meas = rho; return;
    end

    h = rho/norm(rho);

    dg = 10*p.sigma_g*randn;
    du = 10*p.sigma_u*(2*rand-1);
    de = 10*exprnd(1/p.lambda);

    if rand<p.p_imp
        di = 10*p.imp_mag*(2*rand-1);
    else
        di = 0;
    end

    delta = p.a1*dg + p.a2*du + p.a3*de + p.a4*di;

    rho_meas = rho + delta*h;
end

%% TORUS
function pts = torus_from_uv(u,v,R0,r0)
    x = (R0 + r0*cos(v)).*cos(u);
    y = (R0 + r0*cos(v)).*sin(u);
    z = r0*sin(v);
    pts = [x y z];
end