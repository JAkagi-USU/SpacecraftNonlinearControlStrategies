addpath utils

x0_bounds = [370; 370; 370; 1.25; 1.25; 1.25];
x0_mpc = [20; 20; 20; .25; .25; .25];

x0_base = x0_bounds + x0_mpc;


filename = "data/CW_Simulations.mat";



T = 300;
TF = T*1.05;
dt = .05;
t_span = 0:dt:TF;


% Select switching function
switch_func = @sign;
switch_func = @(x) tanh(1*x);



rho_1 = T/2; % > 0
rho_2 = .1; % > delta
alpha = 300; % > 0
q = .996; % 0 < q < 1

n = .0011;

param.q = q;
param.rho_1 = rho_1;
param.rho_2 = rho_2;
param.alpha = alpha;
ctrl = CW_Switch_1(param, param, param, n, 'SwitchFunc', switch_func);

N_sims = 5;

sim_idx = randsample(64,N_sims);
init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;

x_sim = cell(1,N_sims);
t_sim = cell(1,N_sims);
u_sim = cell(1,N_sims);

fig_states = figure();
fig_control = figure();


for j = 1:N_sims
    
    fprintf('Sim %d\n',j);
    fprintf('%.2f %.2f %.2f %.2f %.2f %.2f\n\n',x0)

    idx = sim_idx(j);
    init_cond = init_cond_idx(idx,:);
    x0 = x0_base.*init_cond';
    
[t_all, x_all] = rk4(@(t,x) dyn(t, x, ctrl, n), t_span, x0);

n_steps = numel(t_all);
u_all = NaN(3, n_steps);
for i = 1:n_steps
    
    t = t_all(i);
    x = x_all(i,:)';
    

    
    [~, u_all(:,i)] = dyn(t,x,ctrl, n);
    
end


t_sim{j} = t_all;
    x_sim{j} = x_all;
    u_sim{j} = u_all;








figure(fig_states)
for i = 1:6
subplot(6,1,i)

plot(t_all, x_all(:,i), 'k-')
hold on

end

figure(fig_control)
for i = 1:3
subplot(3,1,i)
plot(t_all, u_all(i,:), 'k-')
hold on
% plot(t_all, u_all(i,:), 'k-')
end


end

save(filename,"t_sim","x_sim","u_sim")




function [x_dot, u] = dyn(t, x, ctrl, n)

drag_a = [0; 5e-6;0]*0;

u = ctrl.calc_u(t, x);


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