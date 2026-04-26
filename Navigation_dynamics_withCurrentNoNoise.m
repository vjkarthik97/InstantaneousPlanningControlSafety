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

%% CURRENT SWEEP
eps_list = [0, 0.02, 0.05, 0.1];
N_cases = length(eps_list);

all_traj = cell(N_cases,1);
all_v_ctrl = cell(N_cases,1);
all_v_tot  = cell(N_cases,1);
all_omega  = cell(N_cases,1);

% ✅ metrics
all_path_length = zeros(N_cases,1);
all_arrival_time = zeros(N_cases,1);

%% START & GOAL
theta = 2*pi*rand;
theta_goal = theta + pi/2;

pos_init = [R0*cos(theta); R0*sin(theta); 0];
goal     = [R0*cos(theta_goal); R0*sin(theta_goal); 0];

%% =========================
%% RUN ALL CASES
%% =========================
for c = 1:N_cases

    eps_c = eps_list(c);
    disp(['Running eps = ', num2str(eps_c)])

    gif_name = sprintf('eps_%.2f.gif',eps_c);

    [traj, v_ctrl, v_tot, omega_hist, path_len, t_arrival] = ...
        simulate_full_visual(pos_init, goal, eps_c, gif_name);

    all_traj{c} = traj;
    all_v_ctrl{c} = v_ctrl;
    all_v_tot{c}  = v_tot;
    all_omega{c}  = omega_hist;

    all_path_length(c) = path_len;
    all_arrival_time(c) = t_arrival;
end

%% PRINT METRICS
disp('--- Drift Study Metrics ---')
for c = 1:N_cases
    fprintf('eps = %.2f | Path Length = %.4f | Arrival Time = %.4f s\n', ...
        eps_list(c), all_path_length(c), all_arrival_time(c));
end

%% =========================
%% TRAJECTORY PLOT (FIXED)
%% =========================
figure; hold on; grid on;

[U,V] = meshgrid(linspace(0,2*pi,60), linspace(0,2*pi,30));
X = (R0 + r0*cos(V)).*cos(U);
Y = (R0 + r0*cos(V)).*sin(U);
Z = r0*sin(V);

surf(X,Y,Z,'FaceAlpha',0.08,'EdgeColor','none');

colors = lines(N_cases);
h_traj = gobjects(N_cases,1);

for c = 1:N_cases
    tr = all_traj{c};
    h_traj(c) = plot3(tr(1,:),tr(2,:),tr(3,:),...
        'LineWidth',3,'Color',colors(c,:));
end

h_start = plot3(pos_init(1),pos_init(2),pos_init(3),'bo','MarkerFaceColor','b','MarkerSize',9);
h_goal  = plot3(goal(1),goal(2),goal(3),'gp','MarkerFaceColor','g','MarkerSize',9);

axis equal; view(3);

labels = arrayfun(@(x) sprintf('\\epsilon=%.2f',x),eps_list,'UniformOutput',false);

legend([h_traj;h_start;h_goal],[labels(:);{'Start'};{'Goal'}],'Fontsize',13);

title('Trajectories (with Drift)','Fontsize',13);
xlabel('X')
xlabel('Y')
xlabel('Z')

%% =========================
%% LINEAR VELOCITY
%% =========================
figure;
for c = 1:N_cases
    subplot(N_cases,1,c); hold on; grid on;

    t = (1:length(all_v_ctrl{c}))*del_T;

    plot(t, vecnorm(all_v_ctrl{c}),'b','LineWidth',2);
    plot(t, vecnorm(all_v_tot{c}),'r--','LineWidth',2);
    legend('Linear Velocity (control)', 'Linear Velocity (with drift)','Fontsize',13)
    xlabel('Time (in seconds)','Fontsize',13)
    ylabel(['\epsilon=',num2str(eps_list(c))],'Fontsize',13);
end
sgtitle('Linear Velocity');

%% =========================
%% ANGULAR VELOCITY
%% =========================
figure;
for c = 1:N_cases
    subplot(N_cases,1,c); hold on; grid on;

    t = (1:length(all_omega{c}))*del_T;
    plot(t, all_omega{c},'m','LineWidth',2);

    ylabel(['\epsilon=',num2str(eps_list(c))]);
end
sgtitle('Angular Velocity');

%% =========================
%% SIMULATION FUNCTION (UNCHANGED LOGIC)
%% =========================
function [traj, v_ctrl_hist, v_tot_hist, omega_hist, path_len, t_arrival] = ...
    simulate_full_visual(pos, goal, eps_c, gif_name)

    N_steps = 1300; del_T = 0.01; tol = 0.1;
    K_1 = 1; K_2 = 1;

    R0 = 2.0; r0 = 0.5;
    Rotation_M = [cos(-pi/3) -sin(-pi/3) 0; sin(-pi/3) cos(-pi/3) 0;0 0 1];
    r_i = Rotation_M*(goal-pos)/norm(goal-pos);

    traj = pos;
    v_ctrl_hist=[]; v_tot_hist=[]; omega_hist=[];

    % ✅ metrics
    path_len = 0;
    t_arrival = N_steps*del_T;

    N_pts = 200;
    u_pts = 2*pi*rand(N_pts,1);
    v_pts = 2*pi*rand(N_pts,1);

    omega_u = 0.2; omega_v = 0.1;

    for k=1:N_steps

        u_pts = mod(u_pts + omega_u*del_T,2*pi);
        v_pts = mod(v_pts + omega_v*del_T,2*pi);
        points = torus_from_uv(u_pts,v_pts,R0,r0);

        %% ORIGINAL SAMPLING + CONTROL (UNCHANGED)
        n_dirs=25;
        centers=zeros(3,n_dirs+1); radii=zeros(1,n_dirs+1);
        goal_dir = (goal - pos)/norm(goal - pos);

        points(end+1,:) = 100*goal_dir;

        for d=1:n_dirs+1
            if(d==1)
                vec = points(end,:)';
            else
                vec = points(randi(N_pts),:)' - pos;
            end
            if norm(vec)<1e-6, continue; end

            dir = vec/norm(vec);
            bounds=[];

            for i=1:N_pts
                p = points(i,:)' - pos;
                if norm(p)<1e-6, continue; end

                if dot(dir,p/norm(p))>0
                    val = norm(p)/(2*dot(dir,p/norm(p)));
                    bounds=[bounds val];
                end
            end

            if isempty(bounds), continue; end

            r=min(bounds);
            centers(:,d)=r*dir; radii(d)=r;
        end

        [~,idx]=min(vecnorm(goal-(pos+centers)));

        if radii(end) > norm(goal-pos)
            center = goal;
        else
            center = pos + centers(:,idx);
        end

        R = pos-center;
        if norm(R)<1e-6, continue; end

        R_hat = R/norm(R);
        S = dot(R_hat,r_i)+1;

        u = -K_1*tanh(norm(R))*sign(dot(R,r_i));

        cross_term = cross(r_i,R_hat);
        cn = norm(cross_term);

        if cn<1e-6
            omega=[0;0;0];
        else
            omega=(-K_2*sign(S)*abs(S)^0.5 - u*cn)*cross_term/cn;
        end

        %% REVERSED DRIFT
        r_xy=[pos(1);pos(2);0];
        if norm(r_xy)>1e-6
            v_c = eps_c*[r_xy(2);-r_xy(1);0]/norm(r_xy);
        else
            v_c=[0;0;0];
        end

        v_ctrl = u*r_i;
        v_tot  = v_ctrl + v_c;

        v_ctrl_hist(:,end+1)=v_ctrl;
        v_tot_hist(:,end+1)=v_tot;
        omega_hist(end+1)=norm(omega);

        pos = pos + v_tot*del_T;

        % ✅ metrics
        path_len = path_len + norm(v_tot)*del_T;

        r_i = r_i + cross(omega,r_i)*del_T;
        r_i = r_i/norm(r_i);

        traj(:,end+1)=pos;

        if norm(pos-goal)<tol
            t_arrival = k*del_T;
            break;
        end
    end
end

function pts = torus_from_uv(u,v,R0,r0)
    x=(R0+r0*cos(v)).*cos(u);
    y=(R0+r0*cos(v)).*sin(u);
    z=r0*sin(v);
    pts=[x y z];
end