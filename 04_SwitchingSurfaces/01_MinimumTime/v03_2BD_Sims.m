addpath utils

x0_bounds = [370; 370; 370; 1.25; 1.25; 1.25];
x0_mpc = [20; 20; 20; .25; .25; .25];

x0_base = x0_bounds + x0_mpc;

mu = 3.9860044188e14;

mass = 24;

% Reference Orbit
ref_el = [6878000; 1e-4; rad2deg(25); rad2deg(45); 0; rad2deg(0)];
[r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
ref_eci_0 = [r; v];

n = sqrt(mu/ref_el(1)^3);




% Brunovsky feedback
C_xy = [-3*n^2 0 0 -2*n;
    0 0 2*n 0];
C_z = [n^2 0];

N_step = 500; % Number of points to sample along trajectory for optimziation

switch_f = @(x,v) -sqrt(2*abs(v*x)).*sign(x);
x_lin = linspace(-1000,1000,500);

% 
lb = 1e-6;
ub = 100;

dt = .05;

ctrl = switch_surf_ctrl(C_xy, C_z);

x_fig = figure();
u_fig = figure();
phase_fig = figure();


prop = orbit_prop_2Sat(dt);




N_sims = 3;             % Set number of simulations

filename = "MinTimeCtrl_2BdDyn.mat";


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
    
    x0_eci_0 = RTN_to_ECI(x0_all, ref_eci_0);



    % Add mass terms
    ref_eci = [ref_eci_0; mass];
    x0_eci = [x0_eci_0; mass];



    % Randomly initialize V1 and V2 controls for switching surface
    % optimization
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


x0 = [ref_eci; x0_eci];
    
[t_all, x_all] = prop.prop_2_rk(t_span, x0, ctrl);


n_steps = size(t_all(:),1);
    u_all = NaN(3,n_steps);
    lvlh_error = NaN(6,n_steps);
    
    for k = 1:n_steps
        t = t_all(k);
        x = x_all(k,:)';
        
        x_ref = x_all(k,1:7)';
        x_sat = x_all(k,8:14)';
        
        [r_rel_x, v_rel_x, ~] = rva_relative(x_ref(1:3)',x_ref(4:6)',x_sat(1:3)',x_sat(4:6)',mu);
        lvlh_error(:,k) = [r_rel_x; v_rel_x];        
        
        
        [ ~, u_all(:,k)] = prop.dyn_2(t, x, ctrl);
        
    end






figure(phase_fig)
subplot(1,3,1)
hold on
plot(x_lin, switch_f(x_lin, ctrl.vx_all(2)),'k--')
plot(lvlh_error(1,:),lvlh_error(4,:),'b-')
plot(lvlh_error(1,1), lvlh_error(4,1), 'b^')

subplot(1,3,2)
hold on
plot(x_lin, switch_f(x_lin, ctrl.vy_all(2)),'k--')
plot(lvlh_error(2,:),lvlh_error(5,:),'b-')
plot(lvlh_error(2,1), lvlh_error(5,1), 'b^')

subplot(1,3,3)
hold on
plot(x_lin, switch_f(x_lin, ctrl.vz_all(2)),'k--')
plot(lvlh_error(3,:),lvlh_error(6,:),'b-')
plot(lvlh_error(3,1), lvlh_error(6,1), 'b^')



figure(x_fig)

for i = 1:6
    
    subplot(6,1,i)
    hold on
    plot(t_all, lvlh_error(i,:),'k-')
    
    
    
end


figure(u_fig)

for i = 1:3
    
    subplot(3,1,i)
    hold on
    plot(t_all, u_all(i,:), 'k-')
    
end

 t_sim{j} = t_all;
x_sim{j} = lvlh_error;
u_sim{j} = u_all;
v_sim{j} = v_final;




end

filepath = sprintf("data/%s",filename);
save(filepath,"t_sim","x_sim","u_sim","v_sim")








% 
% function J = switch_obj(v, x, C, N)
% 
% x0 = x([1;3]);
% y0 = x([2;4]);
% 
% v1_all = v([1;2]);
% v2_all = v([3;4]);
% 
% u_all = switch_planning.calc_u(x0, v1_all, y0, v2_all, C, N);
% 
% J = norm(u_all(:),Inf);
% 
% 
% 
% end
% 
% function J = switch_obj_z(v, x0, C, N)
% 
% u_all = switch_planning.calc_u_z(x0, v, C, N);
% 
% J = norm(u_all(:),Inf);
% 
% 
% 
% end












