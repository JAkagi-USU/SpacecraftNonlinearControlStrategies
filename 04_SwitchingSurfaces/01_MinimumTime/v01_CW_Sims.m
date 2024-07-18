addpath utils

x0_bounds = [370; 370; 370; 1.25; 1.25; 1.25];
x0_mpc = [20; 20; 20; .25; .25; .25];

x0_base = x0_bounds + x0_mpc;

omega = .0011;
C_xy = [-3*omega^2 0 0 -2*omega;
    0 0 2*omega 0];
C_z = [omega^2 0];

N_step = 500;

switch_f = @(x,v) -sqrt(2*abs(v*x)).*sign(x);
x_lin = linspace(-1000,1000,500);


% Upper and lower bounds on V1 and V2 controls for switching surface
lb = 1e-6;
ub = 100;

dt = .05;

ctrl = switch_surf_ctrl(C_xy, C_z);

phase_fig = figure();




N_sims = 3;

filename = "MinTimeCtrl_CWDyn.mat";

sim_idx = randsample(64,N_sims);
init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;

x_sim = cell(1,N_sims);
t_sim = cell(1,N_sims);
u_sim = cell(1,N_sims);
v_sim = cell(1,N_sims);


for j = 1:N_sims

    
    idx = sim_idx(j);
    init_cond = init_cond_idx(idx,:);
    
    
    x0_all = x0_base.*init_cond';
    
% Randomly initialize V1 and V2 controls for optimization problem
vx = lb + (ub - lb)*rand(1);
wx = lb + (ub - lb)*rand(1);
vy = lb + (ub - lb)*rand(1);
wy = lb + (ub - lb)*rand(1);
vz = lb + (ub - lb)*rand(1);
wz = lb + (ub - lb)*rand(1);




v0 = [vx; wx; vy; wy; vz; wz];



ctrl.optimize_ctrl(x0_all, v0, lb, ub, N_step)
vx_opt = ctrl.vx_all;
vy_opt = ctrl.vy_all;
vz_opt = ctrl.vz_all;

t_max = max([ctrl.tx; ctrl.ty; ctrl.tz]);

v_final = [vx_opt; vy_opt; vz_opt];


t_span = 0:dt:(t_max*1.05);
[t_ode, x_ode] = rk4(@(t,x) dyn(t, x, omega, ctrl), t_span, x0_all);


n_ode = numel(t_ode);
u_ode = NaN(3,n_ode);
for i = 1:n_ode
    x = x_ode(i,:)';
    t = t_ode(i);
    [~, u_ode(:,i)] = dyn(t, x, omega, ctrl);
    
end


J_sol = norm(u_ode(:),inf);




figure(phase_fig)
subplot(1,3,1)
hold on
plot(x_lin, switch_f(x_lin, ctrl.vx_all(2)),'k--')
plot(x_ode(:,1),x_ode(:,4),'b-')
plot(x_ode(1,1), x_ode(1,4), 'b^')

subplot(1,3,2)
hold on
plot(x_lin, switch_f(x_lin, ctrl.vy_all(2)),'k--')
plot(x_ode(:,2),x_ode(:,5),'b-')
plot(x_ode(1,2), x_ode(1,5), 'b^')

subplot(1,3,3)
hold on
plot(x_lin, switch_f(x_lin, ctrl.vz_all(2)),'k--')
plot(x_ode(:,3),x_ode(:,6),'b-')
plot(x_ode(1,3), x_ode(1,6), 'b^')





 t_sim{j} = t_ode;
    x_sim{j} = x_ode;
    u_sim{j} = u_ode;
    v_sim{j} = v_final;




end

filepath = sprintf("data/%s",filename);
save(filepath,"t_sim","x_sim","u_sim","v_sim")





function [x_dot, u] = dyn(t, x, n, ctrl)




A = [ 0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];
B = [ 0 0 0;
    0 0 0;
    0 0 0;
    1 0 0;
    0 1 0;
    0 0 1];



u = ctrl.calc_u(t,x);

x_dot = A*x + B*(u + [-.001; .001; -.001]*0);




end






