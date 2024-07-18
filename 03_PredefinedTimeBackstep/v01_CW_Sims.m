addpath utils


M2KM = 1/1000;
S2HR = 1/3600;

x0_bounds = [370; 370; 370; 1.25; 1.25; 1.25].*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];
x0_mpc = [20; 20; 20; .25; .25; .25].*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];

x0_base = x0_bounds + x0_mpc;

n = .0011/S2HR;
dt = 10*S2HR;

N_sim = 800;
T = N_sim*dt;
TF = (N_sim+10)*dt;


Dx = 7e-6*M2KM/S2HR^2;
Dz = 7e-6*M2KM/S2HR^2;




filename = "data/CW_Sims.mat";





%%




N_sims = 5;

sim_idx = randsample(64,N_sims);
init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;


total_error_fig = figure();
total_u_fig = figure();





DV_all = NaN(1,N_sims);


x_sim = cell(1,N_sims);
t_sim = cell(1,N_sims);
u_sim = cell(1,N_sims);



for j = 1:N_sims
    
    fprintf('Sim %d\n',j);
    
    idx = sim_idx(j);
    init_cond = init_cond_idx(idx,:);
    x0_act = x0_base.*init_cond';
    
    
    % Need to manually change the number of additional basis functions
    % within the FB_HCW class itself
    ctrl = FB_HCW(n,T);
    ctrl.init_ctrl(x0_act, Dx, Dz);
    t_span = 0:dt:TF;
    
    
    % Can add drag perturbations in the dyn function 
    [t_all, x_all] = rk4(@(t,x) dyn(t, x, n, @ctrl.calc_u), t_span, x0_act);
    
    n_steps = size(t_all(:),1);
    u_all = NaN(3,n_steps);
    
    for k = 1:n_steps
        t = t_all(k);
        x = x_all(k,:)';
        
        u_all(:,k) = ctrl.calc_u(t,x);
        
    end
    
    
    x_all = x_all.*[1/M2KM; 1/M2KM; 1/M2KM; S2HR/M2KM; S2HR/M2KM; S2HR/M2KM]';
    u_all = u_all*(S2HR^2)/(M2KM);
    
    t_sim{j} = t_all;
    x_sim{j} = x_all;
    u_sim{j} = u_all;
    

    figure(total_error_fig)
    
    for ii = 1:6
       subplot(6,1,ii)
       hold on
       plot(t_all, x_all(:,ii),'k-')
    end
    
    figure(total_u_fig)
    
    for ii = 1:3
       subplot(3,1,ii)
       hold on
       plot(t_all, u_all(ii,:),'k-')
    end

    
    drawnow;
    
end

save(filename,"t_sim","x_sim","u_sim")



figure(total_error_fig)

subplot(6,1,1)
title('State Errors in LVLH Frame (m, m/s)')
plot(t_all(end),0,'kx')
ylabel('X')

subplot(6,1,2)
plot(t_all(end),0,'kx')
ylabel('Y')

subplot(6,1,3)
plot(t_all(end),0,'kx')
ylabel('Z')

subplot(6,1,4)
plot(t_all(end),0,'kx')
ylabel('VX')

subplot(6,1,5)
plot(t_all(end),0,'kx')
ylabel('VY')

subplot(6,1,6)
plot(t_all(end),0,'kx')
ylabel('VZ')
xlabel('Time (hr)')












function [x_dot, u] = dyn(t, x, n, u_func)

M2KM = 1/1000;
S2HR = 1/3600;

drag_a = 0*[0; 5e-6;0]*M2KM/(S2HR^2);


u = u_func(t,x);


A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*n^2 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -n^2 0 0 0];
 
B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];


x_dot = A*x + B*(u + drag_a);

end
