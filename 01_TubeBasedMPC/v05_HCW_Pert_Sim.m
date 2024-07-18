addpath utils

x0_bounds = [350; 350; 350; 1; 1; 1];
x_max_guide = [425; 425; 425; 5; 5; 5];
u_max_guide = [.01; .01; .01];

x0_mpc = [20; 20; 20; .25; .25; .25];
x_robust_mpc = [60;30;40;4;4;4];
u_robust_mpc = [.006; .006; .006];

x_max_mpc = [75; 75; 75; 5; 5; 5];
u_max_mpc = [.01; .01; .01];

poles_d = [.85; .85; .95; .95; .9; .9];

W_max = [.1;.1;.1;.02;.02;.02];     % Set max bounds on perturbations

combo = ff2n(6);
combo(combo == 0) = -1;

n = .0011;
dt = 10;

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



K = -place(AD, BD, poles_d);


Q_guide = zeros(6);
R_guide = eye(3);
N_guide = 800;

Q_mpc = .001*eye(6);
R_mpc = eye(3);

N_MPC = 30;     % Set horizon of MPC, in steps

N_sims = 2;     % Set number of sims to run







guide_fig = figure();
mpc_fig = figure();
u_fig = figure();
u_traj_fig = figure();

DV_all = NaN(1,N_sims);
total_errors_all = NaN(6,N_guide+1,N_sims);
mpc_errors_all = NaN(6,N_guide+1,N_sims);
u_all = NaN(3,N_guide,N_sims);
u_all_guide = NaN(3,N_guide,N_sims);



for j = 1:N_sims
    
    fprintf('Sim %d\n',j);
    
    
    
    
    rand_max = randi(2,6,1);
    rand_max(rand_max == 2) = -1;
    
    x0_guide = rand_max.*x0_bounds;
    
    ctrl = sat_control(AD, BD, K);
    ctrl.init_guide(N_guide, Q_guide, R_guide, x_max_guide, u_max_guide)
    ctrl.init_mpc(N_MPC, Q_mpc, R_mpc, x_robust_mpc, u_robust_mpc)
    ctrl.calc_guidance(x0_guide, zeros(6,1))
    
    rand_max = randi(2,6,1);
    rand_max(rand_max == 2) = -1;
    x_act = x0_guide + x0_mpc.*rand_max;
    x_error = x_act - x0_guide;
    
    x_full = NaN(6,N_guide+1);
    x_full(:,1) = x_act;
    
    x_error_all = NaN(6,N_guide+1);
    x_error_all(:,1) = x_error;
    mpc_errors_all(:,1) = x_error;
    
    u_full = NaN(3,N_guide);
    x_guide = NaN(6,N_guide+1);
    
    [total_error, guide_error] = ctrl.error(1, x_act, zeros(6,1));
    x_guide(:,1) = total_error;
    total_errors_all(:,1,j) = total_error;
    
    for i = 1:N_guide
        
        [u, u_mpc_combine, u_guide] = ctrl.ctrl(i, x_act, zeros(6,1));
        
        u_full(:,i) = u;
        u_all(:,i,j) = u;
        u_all_guide(:,i,j) = u_guide;
        
        if any(abs(u_full(:,i)) > .02)
            u_full(:,i)
        end
        
        rand_W = -1 + 2*rand(6,1);
        
        x_act = AD*x_full(:,i) + BD*u_full(:,i) + rand_W.*W_max;
        x_full(:,i+1) = x_act;
        
        [total_error, guide_error] = ctrl.error(i+1, x_act, zeros(6,1));
        x_guide(:,i+1) = total_error;
        total_errors_all(:,i+1,j) = total_error;
        x_error_all(:,i+1) = guide_error;
        mpc_errors_all(:,i+1,j) = guide_error;
        
        
    end
    
    
    figure(guide_fig)
    for i = 1:6
        
        subplot(6,1,i)
        hold on
        plot([0 N_guide*dt/60],[x_max_guide(i) x_max_guide(i)],'k-')
        plot([0 N_guide*dt/60],[-x_max_guide(i) -x_max_guide(i)],'k-')
        plot((0:N_guide)*dt/60,ctrl.x_traj(i,:),'k-')
        
    end
    
    
    
    figure(mpc_fig)
    for i = 1:6
        
        subplot(6,1,i)
        hold on
        plot([0 N_guide*dt/60],[x_max_mpc(i) x_max_mpc(i)],'k-')
        plot([0 N_guide*dt/60],[-x_max_mpc(i) -x_max_mpc(i)],'k-')
        plot((0:N_guide)*dt/60,x_error_all(i,:),'k-')
        
    end
    
    figure(u_fig)
    
    for i = 1:3
        subplot(3,1,i)
        hold on
        plot([0 N_guide]*dt/60,[.02 .02],'k-')
        plot([0 N_guide]*dt/60,[-.02 -.02],'k-')
        plot((1:N_guide)*dt/60,u_full(i,:),'k-')
        
    end
    
    figure(u_traj_fig)
    
    for i = 1:3
        subplot(3,1,i)
        hold on
        plot([0 N_guide]*dt/60,[.01 .01],'k-')
        plot([0 N_guide]*dt/60,[-.01 -.01],'k-')
        plot((1:N_guide)*dt/60,ctrl.u_traj(i,:),'k-')
        
    end
    
    drawnow
    
    
end

save("data/HCWSimResults.mat","total_errors_all","mpc_errors_all","u_all","u_all_guide")

