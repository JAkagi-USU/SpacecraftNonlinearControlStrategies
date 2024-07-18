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


T_conv = 400;
dt = .05;


% Select switching function
% switch_func = @(x) sign(x);
switch_func = @(x) tanh(100*x);



x_fig = figure();
u_fig = figure();


prop = orbit_prop_2Sat(dt);




N_sims = 3;

filename = "data/TwoBody_Simulations.mat";



init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;




sim_idx = randsample(64,N_sims);



x_sim = cell(1,N_sims);
t_sim = cell(1,N_sims);
u_sim = cell(1,N_sims);
conv_sim = NaN(N_sims,1);
T_sim_all = NaN(N_sims, 1);


for j = 1:N_sims
    
    fprintf('Sim: %d/%d\n', j, N_sims)
    

    
    idx = sim_idx(j);
    init_cond = init_cond_idx(idx,:);
    

    
    TF = T_conv*1.05;
    
    
    a1 = 128/T_conv^2;
    a2 = 64/T_conv^2;
    b1 = 128/T_conv^2;
    b2 = 64/T_conv^2;
    k = .1;


    param.a1 = a1;
    param.a2 = a2;
    param.b1 = b1;
    param.b2 = b2;
    param.k = k;
    ctrl = CW_Switch_2(param, param, param, n, 'SwitchFunc', switch_func);
    
    
    x0_all = x0_base.*init_cond'
    
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

save(filename,"t_sim","x_sim","u_sim", "T_sim_all", "conv_sim")




















