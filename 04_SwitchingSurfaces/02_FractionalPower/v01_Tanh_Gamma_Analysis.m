addpath utils

x0_bounds = [370; 370; 370; 1.25; 1.25; 1.25];
x0_mpc = [20; 20; 20; .25; .25; .25];

x0_base = x0_bounds + x0_mpc;


filename = "data/TanhGammaAnalysis.mat";


n_Times = 3; % Number of convergence time bounds
n_Tanh = 80; % Number of gamma values

T_range = linspace(600,8000, n_Times);
tanh_range = exp(linspace(-10, 7, n_Tanh));

[T_grid, tanh_grid] = meshgrid(T_range, tanh_range);

N_sims = numel(T_grid);


n = .0011;
dt = .05;





% sim_idx = randsample(64,N_sims);
init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;

x_sim = cell(1,N_sims);
t_sim = cell(1,N_sims);
u_sim = cell(1,N_sims);

fig_states = figure();
fig_control = figure();


for j = 1:N_sims
    
%     idx = sim_idx(j);
    idx = 64;
    init_cond = init_cond_idx(idx,:);
    x0 = x0_base.*init_cond';
    
    fprintf('Sim %d\n',j);
    fprintf('%.2f %.2f %.2f %.2f %.2f %.2f\n\n',x0)
    
    
    T = T_grid(j);
    alpha = tanh_grid(j);
    
    
    
    a1 = 128/T^2;
    a2 = 64/T^2;
    b1 = 128/T^2;
    b2 = 64/T^2;
    k = .1;

    switch_func = @(x) tanh(alpha*x);
%     switch_func = @(x) sign(x);



    param.a1 = a1;
    param.a2 = a2;
    param.b1 = b1;
    param.b2 = b2;
    param.k = k;
    ctrl = CW_Switch_2(param, param, param, n, 'SwitchFunc', switch_func);

    
    t_span = 0:dt:T*1.1;
% [t_all, x_all] = rk4(@(t,x) dyn(t, x, ctrl, n), t_span, x0);
[t_all, x_all] = prop_zero_hold(t_span, x0, ctrl, n);
% [t_all, x_all] = ode45(@(t,x) dyn(t, x, ctrl, n), t_span, x0);

n_steps = numel(t_all);
u_all = NaN(3, n_steps);

for i = 1:n_steps
    
    t = t_all(i);
    x = x_all(i,:)';
    

    
%     [~, u_all(:,i)] = dyn(t,x,ctrl, n);
      [~, u_all(:,i)] = dyn_u(t, x, ctrl.calc_u(t,x), n);  
    
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


save(filename,"t_sim","x_sim","u_sim", "T_grid", "tanh_grid")


function [t_span, x_all] = prop_zero_hold(t_span, x0, ctrl, n)

n_steps = numel(t_span);

x_all = NaN(n_steps, 6);
x_all(1,:) = x0;

for i = 1:n_steps-1
    
    t = t_span(i);
    x = x_all(i,:)';
    u = ctrl.calc_u(t,x);
    
    [~, x_rk4] = rk4(@(t,x) dyn_u(t, x, u, n), [t t_span(i+1)], x);
    
    x_all(i+1,:) = x_rk4(end,:);
    
    
end

end

function [x_dot, u] = dyn_u(t, x, u, n)

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


x_dot = A*x + B*u;

end

% function [x_dot, u] = dyn(t, x, ctrl, n)
% 
% drag_a = [0; 5e-6;0]*0;
% 
% u = ctrl.calc_u(t, x);
% 
% 
% 
% A = [0 0 0 1 0 0;
%      0 0 0 0 1 0;
%      0 0 0 0 0 1;
%      3*n^2 0 0 0 2*n 0;
%      0 0 0 -2*n 0 0;
%      0 0 -n^2 0 0 0];
%  
% B = [0 0 0;
%      0 0 0;
%      0 0 0;
%      1 0 0;
%      0 1 0;
%      0 0 1];
% 
% 
% x_dot = A*x + B*(u + drag_a);
% 
% end
% 
% function x_all = trap_int(x0, t_all, xdot_all)
% 
% n_steps = numel(t_all);
% 
% x_all = NaN(6,n_steps);
% x_all(:,1) = x0;
% 
% for i = 1:n_steps-1
%     a = xdot_all(:,i);
%     b = xdot_all(:,i+1);
%     h = t_all(i+1) - t_all(i);
%     
%     x_all(:,i+1) = x_all(:,i) + h*(a+b)/2;
% 
% end
% 
% 
% 
% end