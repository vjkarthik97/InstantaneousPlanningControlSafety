
clc; close all; clear;

%% PARAMETERS
N_steps = 4000;
del_T = 0.01;

K_1 = 1;
K_2 = 1;

%% TORUS
R0 = 2.0;
r0 = 0.5;
tol = 0.1;

%% DENSITY SWEEP
N_pts_list = [50, 150, 300, 350, 400, 500];

%% FIXED START & GOAL (pi/2 apart)
theta0 = 2*pi*rand;
theta_goal = theta0 + 2*pi/3;

rho = R0 + (rand-0.5)*r0;
z = (rand-0.5)*r0;

pos_init  = [rho*cos(theta0); rho*sin(theta0); z];
goal      = [rho*cos(theta_goal); rho*sin(theta_goal); z];

%% LOOP OVER SCENARIOS
for scenario = 1:length(N_pts_list)

    fprintf('\n--- Scenario %d ---\n', scenario);

    %% RESET STATE
    pos = pos_init;
    r_i = [0;1;0];
    traj = pos;

    N_pts = N_pts_list(scenario);

    %% OBSTACLES
    u_pts = 2*pi*rand(N_pts,1);
    v_pts = 2*pi*rand(N_pts,1);

    obs_type = randi([1 3], N_pts, 1);

    omega1 = 0.3 + 0.2*rand(N_pts,1);
    omega2 = 0.2 + 0.2*rand(N_pts,1);

    %% METRICS
    replan_time = [];

    %% FIGURE
    figure('Name',sprintf('Scenario %d',scenario));
    hold on; grid on;
    axis equal;
    xlim([-3 3]); ylim([-3 3]); zlim([-2 2]);
    view(45,45);

    [XT,YT,ZT] = torus_surface(R0,r0);
    surf(XT,YT,ZT,'FaceAlpha',0.05,'EdgeColor','none');

    %% INITIAL POINTS
    points = torus_from_uv(u_pts,v_pts,R0,r0);

    h1 = scatter3(points(obs_type==1,1),points(obs_type==1,2),points(obs_type==1,3),10,'b','filled');
    h2 = scatter3(points(obs_type==2,1),points(obs_type==2,2),points(obs_type==2,3),10,'g','filled');
    h3 = scatter3(points(obs_type==3,1),points(obs_type==3,2),points(obs_type==3,3),10,'r','filled');

    plot3(pos(1),pos(2),pos(3),'bo','MarkerFaceColor','b');
    plot3(goal(1),goal(2),goal(3),'g*','MarkerSize',12);

    %% TRAJECTORY HANDLE
    h_traj = plot3(traj(1,:),traj(2,:),traj(3,:),'b','LineWidth',2);
    h_robot = plot3(pos(1),pos(2),pos(3),'ko','MarkerFaceColor','k');

    %% SPHERES
    sphere_handles = [];
    max_spheres = 20;

    %% MAIN LOOP
    for k = 1:N_steps

        %% UPDATE OBSTACLES
        for i = 1:N_pts
            if obs_type(i)==1
                u_pts(i) = mod(u_pts(i) + omega1(i)*del_T, 2*pi);
            elseif obs_type(i)==2
                u_pts(i) = mod(u_pts(i) + omega1(i)*del_T, 2*pi);
                v_pts(i) = mod(v_pts(i) + 0.3*sin(omega2(i)*k*del_T)*del_T, 2*pi);
            else
                u_pts(i) = mod(u_pts(i) + omega1(i)*cos(v_pts(i))*del_T, 2*pi);
                v_pts(i) = mod(v_pts(i) + omega2(i)*sin(u_pts(i))*del_T, 2*pi);
            end
        end

        points = torus_from_uv(u_pts,v_pts,R0,r0);

        %% TIMER
        tic;

        %% SAMPLING
        n_dirs_random = N_pts;
        n_dirs = n_dirs_random + 1;

        centers = zeros(3,n_dirs);
        radii = zeros(1,n_dirs);

        for d = 1:n_dirs_random
            rand_pt = points(randi(N_pts),:)';
            vec = rand_pt - pos;
            if norm(vec)<1e-6, continue; end

            dir = vec/norm(vec);
            bounds = [];

            for i=1:N_pts
                p_vec = points(i,:)' - pos;
                if norm(p_vec)<1e-6, continue; end

                if dot(dir,p_vec/norm(p_vec))>0
                    val = norm(p_vec)/(2*dot(dir,p_vec/norm(p_vec)));
                    bounds = [bounds val];
                end
            end

            if isempty(bounds), continue; end

            r_opt = min(bounds);
            centers(:,d) = r_opt*dir;
            radii(d) = r_opt;
        end

        %% GOAL DIRECTION
        goal_vec = goal - pos;
        dist_goal = norm(goal_vec);

        if dist_goal > 1e-6
            dir = goal_vec / dist_goal;
            bounds = [];

            for i=1:N_pts
                p_vec = points(i,:)' - pos;
                if norm(p_vec)<1e-6, continue; end

                if dot(dir,p_vec/norm(p_vec))>0
                    val = norm(p_vec)/(2*dot(dir,p_vec/norm(p_vec)));
                    bounds = [bounds val];
                end
            end

            if ~isempty(bounds)
                r_opt = min(bounds);

                if r_opt >= dist_goal
                    centers(:,n_dirs) = goal - pos;
                    radii(n_dirs) = dist_goal;
                else
                    centers(:,n_dirs) = r_opt * dir;
                    radii(n_dirs) = r_opt;
                end
            end
        end

        %% SELECT BEST
        best_idx = 1; best_dist = inf;

        for d=1:n_dirs
            if radii(d)==0, continue; end

            c = pos + centers(:,d);
            dist = norm(goal - c) - 0.2*radii(d);

            if dist < best_dist
                best_dist = dist;
                best_idx = d;
            end
        end

        center = pos + centers(:,best_idx);
        radius = radii(best_idx);

        replan_time(end+1) = toc;

        %% CONTROLLER
        R = pos - center;
        if norm(R)<1e-6, continue; end

        R_hat = R/norm(R);
        S = dot(R_hat,r_i)+1;

        u = -K_1*tanh(norm(R))*sign(dot(R,r_i));

        cross_term = cross(r_i,R_hat);
        cn = norm(cross_term);

        if cn<1e-6
            omega=[0;0;0];
        else
            omega = (-K_2*sign(S)*abs(S)^0.5 - u*cn)*cross_term/cn;
        end

        %% UPDATE STATE
        pos = pos + u*r_i*del_T;
        r_i = r_i + cross(omega,r_i)*del_T;
        r_i = r_i/norm(r_i);

        traj(:,end+1)=pos;

        %% VISUALIZATION
        if mod(k,30)==0

            set(h1,'XData',points(obs_type==1,1),'YData',points(obs_type==1,2),'ZData',points(obs_type==1,3));
            set(h2,'XData',points(obs_type==2,1),'YData',points(obs_type==2,2),'ZData',points(obs_type==2,3));
            set(h3,'XData',points(obs_type==3,1),'YData',points(obs_type==3,2),'ZData',points(obs_type==3,3));

            % Update trajectory
            set(h_traj,'XData',traj(1,:),'YData',traj(2,:),'ZData',traj(3,:));
            set(h_robot,'XData',pos(1),'YData',pos(2),'ZData',pos(3));

            % Fade old spheres
            for i = 1:length(sphere_handles)
                set(sphere_handles(i),'FaceAlpha',0.03);
            end

            % Draw new sphere
            [Xs,Ys,Zs] = sphere(20);
            h = surf(radius*Xs+center(1), ...
                     radius*Ys+center(2), ...
                     radius*Zs+center(3), ...
                     'FaceAlpha',0.25,'EdgeColor','none','FaceColor',[0 0.6 1]);

            sphere_handles = [sphere_handles h];

            if length(sphere_handles) > max_spheres
                delete(sphere_handles(1));
                sphere_handles(1) = [];
            end

            drawnow;
        end

        %% STOP
        if norm(pos-goal)<tol
            disp('Goal reached!');
            break;
        end
    end

    %% SAVE FINAL FIGURE
    saveas(gcf, sprintf('trajectory_scenario_%d.png',scenario));

    %% PRINT METRICS
    fprintf('Avg replanning time: %.6f s\n',mean(replan_time));
    fprintf('Max replanning time: %.6f s\n',max(replan_time));

end

%% FUNCTIONS

function pts = torus_from_uv(u,v,R0,r0)
    x = (R0 + r0*cos(v)).*cos(u);
    y = (R0 + r0*cos(v)).*sin(u);
    z = r0*sin(v);
    pts = [x y z];
end

function [X,Y,Z] = torus_surface(R0,r0)
    [u,v] = meshgrid(linspace(0,2*pi,40), linspace(0,2*pi,20));
    X = (R0 + r0*cos(v)).*cos(u);
    Y = (R0 + r0*cos(v)).*sin(u);
    Z = r0*sin(v);
end

