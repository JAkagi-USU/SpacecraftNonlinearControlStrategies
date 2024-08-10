%% Parameters

addpath utils

mu = 3.9860044188e14;

x0_bounds = [350; 350; 350; 1; 1; 1];
x_max_guide = [425; 425; 425; 5; 5; 5];
u_max_guide = [.01; .01; .01];

x0_mpc = [20; 20; 20; .25; .25; .25];
x_robust_mpc = [60;30;40;4;4;4];
u_robust_mpc = [.006; .006; .006];

x_max_mpc = [75; 75; 75; 5; 5; 5];
u_max_mpc = [.01; .01; .01];

full_error_max = x_max_guide + x_max_mpc;

poles_d = [.85; .85; .95; .95; .9; .9];

Q_guide = zeros(6);
R_guide = eye(3);

Q_mpc = .001*eye(6);
R_mpc = 1*eye(3);

mass = 24;
dt = 10;

N_guide = 800;      % Guidance horizon, in steps
N_mpc = 30;         % MPC horizon, in steps

N_sims = 2;        % Number of sims to run

prop = orbit_prop();

fig_tot_error = figure();
fig_mpc_error = figure();
fig_u = figure();
u_traj_fig = figure();

DV_all = NaN(1,N_sims);

%% Init data collecting

total_errors_all = NaN(6,N_guide+1,N_sims);
mpc_errors_all = NaN(6,N_guide+1,N_sims);
u_all = NaN(3,N_guide,N_sims);
u_all_guide = NaN(3,N_guide,N_sims);



for j = 1:N_sims
    %% Initial Conditions
    
    fprintf('%d\n',j)

    % Reference Orbit
    ref_el = [6878000; 1e-4; rad2deg(25); rad2deg(45); 0; rad2deg(0)];
    ROE = [0; 0.1; .0001454; 0; .0001454; 0];
    [r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
    ref_eci = [r; v];

    % Desired Orbit
    xd_el = ROE_to_El(ROE, ref_el);
    [r, v] = orbel2rv(xd_el(1),xd_el(2),xd_el(3),xd_el(4),xd_el(5),xd_el(6),mu);
    xd_eci = [r;v];
    [r_rel_x, v_rel_x, ~] = rva_relative(ref_eci(1:3)',ref_eci(4:6)',xd_eci(1:3)',xd_eci(4:6)',mu);
    xd_rel = [r_rel_x; v_rel_x];

    % Actual Orbit
    rand_max = randi(2,6,1);
    rand_max(rand_max == 2) = -1;
    x_guide_error =rand_max.*x0_bounds;

    rand_max = randi(2,6,1);
    rand_max(rand_max == 2) = -1;
    x_mpc_error = rand_max.*x0_mpc;

    x0_rel = xd_rel + x_guide_error + x_mpc_error;
    x0_eci = RTN_to_ECI(x0_rel, ref_eci); % Function is correct

    % Initial state for guidance
    x0_guide_rel = xd_rel + x_guide_error;

    % Add mass terms
    ref_eci = [ref_eci; mass];
    xd_eci = [xd_eci; mass];
    x0_eci = [x0_eci; mass];



    %% HCW Dynamics

    n = sqrt(mu/ref_el(1)^3);

    AD_rr = @(n,dt) [4-3*cos(n*dt), 0, 0; 6*(sin(n*dt)-n*dt), 1, 0; 0, 0, cos(n*dt)];
    AD_vr = @(n,dt) [sin(n*dt)/n, 2*(1-cos(n*dt))/n, 0; 2*(cos(n*dt)-1)/n,(4*sin(n*dt)-3*n*dt)/n,0;0,0,sin(n*dt)/n];
    AD_rv = @(n,dt) [3*n*sin(n*dt), 0, 0; 6*n*(cos(n*dt)-1),0,0; 0,0,-n*sin(n*dt)];
    AD_vv = @(n,dt) [cos(n*dt), 2*sin(n*dt),0; -2*sin(n*dt),4*cos(n*dt)-3,0;0,0,cos(n*dt)];

    AD = [AD_rr(n,dt), AD_vr(n,dt);
        AD_rv(n,dt), AD_vv(n,dt)];

    BD_f = @(n,dt) [2*(sin(n*dt/2)/n)^2, -2*(sin(n*dt) - n*dt)/(n^2), 0;
        2*(sin(n*dt) - n*dt)/(n^2), -(8*cos(n*dt) + 3*n^2*dt^2-8)/(2*n^2),0;
        0,0,2*(sin(n*dt/2)/n)^2;
        sin(n*dt)/n, -(2*cos(n*dt)-2)/n,0;
        -4*sin(n*dt/2)^2/n, 4*sin(n*dt)/n-3*dt,0;
        0,0,sin(n*dt)/n];

    BD = BD_f(n, dt);

    %% Controller Setup

    K = -place(AD, BD, poles_d);

    ctrl = sat_control(AD, BD, K);
    ctrl.init_guide(N_guide, Q_guide, R_guide, x_max_guide, u_max_guide)
    ctrl.init_mpc(N_mpc, Q_mpc, R_mpc, x_robust_mpc, u_robust_mpc)


    ctrl.calc_guidance(x0_guide_rel, xd_rel)

    

    


    [total_error, mpc_error] = ctrl.error(1, x0_rel, xd_rel);
    total_errors_all(:,1,j) = total_error;
    mpc_errors_all(:,1,j) = mpc_error;


    N = 800;

    for i = 1:N

        

        [u, u_mpc_combine, u_guide] = ctrl.ctrl(i, x0_rel, xd_rel);
        u_all(:,i,j) = u;
        u_all_guide(:,i,j) = u_guide;

        T = ECI_2_LVLH_T(ref_eci(1:6));
        u_eci = T'*u;

        % Update Dynamics
        ref_eci = prop.prop(dt, ref_eci, zeros(3,1));
        xd_eci = prop.prop(dt, xd_eci, zeros(3,1));
        x0_eci = prop.prop(dt, x0_eci, u_eci);

        [r_rel_x, v_rel_x, ~] = rva_relative(ref_eci(1:3)',ref_eci(4:6)',xd_eci(1:3)',xd_eci(4:6)',mu);
        xd_rel = [r_rel_x; v_rel_x];
        [r_rel_x, v_rel_x, ~] = rva_relative(ref_eci(1:3)',ref_eci(4:6)',x0_eci(1:3)',x0_eci(4:6)',mu);
        x0_rel = [r_rel_x; v_rel_x];

        [total_error, mpc_error] = ctrl.error(i+1, x0_rel, xd_rel);
        total_errors_all(:,i,j) = total_error;
        mpc_errors_all(:,i,j) = mpc_error;


    end    



    figure(fig_tot_error)
    for i = 1:6

        subplot(6,1,i)
        hold on
        plot([0 N_guide]*dt/60,[full_error_max(i) full_error_max(i)],'k-')
        plot([0 N_guide]*dt/60,[-full_error_max(i) -full_error_max(i)],'k-')
        plot((0:N_guide)*dt/60,total_errors_all(i,:,j),'k-')

    end



    figure(fig_mpc_error)
    for i = 1:6

        subplot(6,1,i)
        hold on
        plot([0 N_guide]*dt/60,[x_max_mpc(i) x_max_mpc(i)],'k-')
        plot([0 N_guide]*dt/60,[-x_max_mpc(i) -x_max_mpc(i)],'k-')
        plot((0:N_guide)*dt/60,mpc_errors_all(i,:,j),'k-')

    end

    figure(fig_u)

    for i = 1:3
        subplot(3,1,i)
        hold on
        plot([0 N_guide]*dt/60,[.02 .02],'k-')
        plot([0 N_guide]*dt/60,[-.02 -.02],'k-')
        plot((1:N_guide)*dt/60,u_all(i,:,j),'k-')

    end
    
    figure(u_traj_fig)
    
    for i = 1:3
        subplot(3,1,i)
        hold on
        plot([0 N_guide]*dt/60,[.01 .01],'k-')
        plot([0 N_guide]*dt/60,[-.01 -.01],'k-')
        plot((1:N_guide)*dt/60,ctrl.u_traj(i,:),'k-')
        
    end
    
    drawnow;
end

save("data/TwoBodySimResults.mat","total_errors_all","mpc_errors_all","u_all","u_all_guide")






