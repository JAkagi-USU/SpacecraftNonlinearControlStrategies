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


T = 300;
TF = T*1.05;
dt = .05;


% Select switching function
% switch_func = @(x) tanh(x);
switch_func = @(x) sign(x);


rho_1 = T/2; % > 0
rho_2 = .1;


% Number of alpha and q parameters to search over
n_alpha = 15;
n_q = 15;



alpha_range = linspace(1e-3,1000, n_alpha);
q_range = linspace(1e-3, .999, n_q);


[alpha_grid, q_grid] = meshgrid(alpha_range, q_range);

N_sims = numel(alpha_grid);


x_fig = figure();
u_fig = figure();


prop = orbit_prop_2Sat(dt);





filename = "QA_Param_Compare.mat";



x_sim = cell(1,N_sims);
t_sim = cell(1,N_sims);
u_sim = cell(1,N_sims);
error_sim = NaN(N_sims,1);
u_max_sim = NaN(N_sims,1);
conv_sim = NaN(N_sims,1);


for j = 1:N_sims
    
    fprintf("Sim: %d/%d\n", j, N_sims)

    
    q_sim = q_grid(j);
    alpha_sim = alpha_grid(j);
    
    param.q = q_sim;
    param.rho_1 = rho_1;
    param.rho_2 = rho_2;
    param.alpha = alpha_sim;
    ctrl = CW_Switch_1(param, param, param, n, 'SwitchFunc', switch_func);
    
    init_cond = [1;1;1;1;1;1];
    x0_all = x0_base.*init_cond;
    
    x0_eci_0 = RTN_to_ECI(x0_all, ref_eci_0);



    % Add mass terms
    ref_eci = [ref_eci_0; mass];
    x0_eci = [x0_eci_0; mass];




    






t_span = 0:dt:TF;


x0 = [ref_eci; x0_eci];
    
[t_all, x_all] = prop.prop_2_piece(t_span, x0, ctrl);



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


in_bounds = vecnorm(lvlh_error(1:3,:)) <= .1;
if in_bounds(end) == 1
    in_bound_idx = find(diff(in_bounds) ~= 0,1,'last') + 1;
    conv_sim(j) = t_all(in_bound_idx);
end



u_max_sim(j) = max(abs(u_all(:)));

error_sim(j) = norm(lvlh_error(1:3,end));



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






end

filepath = sprintf("data/%s",filename);
save(filepath,"t_sim","x_sim","u_sim","error_sim", "q_grid", "alpha_grid", "conv_sim", "u_max_sim")





% Plot stable
error_sim(isnan(error_sim)) = inf;
stable = error_sim <= 10;
unstable = error_sim > 10;

figure()
hold on
plot(alpha_grid(stable), q_grid(stable), 'kx')
plot(alpha_grid(unstable), q_grid(unstable), 'k.')





















