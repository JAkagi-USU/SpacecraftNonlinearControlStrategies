addpath utils

file_input = "TwoBody_Sims";
filename = sprintf("data/%s.mat",file_input);

prop = orbit_prop_2Sat();




%% Reference Orbit

mu = 3.9860044188e14;
ref_el = [6878000; 1e-4; rad2deg(25); rad2deg(45); 0; rad2deg(0)];
[r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
ref_eci_0 = [r; v];

n = sqrt(mu/ref_el(1)^3);

mass = 24;


%% Set Up Controller

M2KM = 1/1000;
S2HR = 1/3600;

n_hr = n/S2HR;
a = 2*n_hr;
b = 3*n_hr^2;

InPlane.A0 = [0 1 0 0; b 0 0 a; 0 0 0 1; 0 -a 0 0];
InPlane.B0 = [0;0;0;1];

InPlane.TC = [0, n_hr^2/a, 0, 0];
InPlane.TD = [0 a 0 0; 0 0 a 0; -b 0 1 0; 0 -b 0 1];

CrossTrack.A0 = [0 1; -n_hr^2 0];
CrossTrack.B0 = [0; 1];
CrossTrack.TC = [n_hr^2 0];
CrossTrack.TD = [1 0; 0 1];

ctrl = CommonLyap();
ctrl.setInPlaneDyn(InPlane.A0, InPlane.B0, InPlane.TC, InPlane.TD);
ctrl.setCrossTrackDyn(CrossTrack.A0, CrossTrack.B0, CrossTrack.TC, CrossTrack.TD);

%% Solve Common Lyapunov Problem

R_val = .01;
Q_val = 1;
    
[C_InPlane, ~, ~] = lqr(ctrl.InPlane.AZ, ctrl.InPlane.BZ, Q_val*eye(4), R_val*eye(1));
[C_CrossTrack, ~, ~] = lqr(ctrl.CrossTrack.AZ, ctrl.CrossTrack.BZ,  Q_val*eye(2), R_val*eye(1));

ctrl.setInPlaneGain(-C_InPlane);
ctrl.setCrossTrackGain(-C_CrossTrack);


max_bound = 1e2;
min_bound = 1e-6;

Q_max_bound = 1;
P_max_bound = max_bound;
R_max_bound = max_bound;

Q_min_bound = min_bound;
P_min_bound = min_bound;
R_min_bound = min_bound;


alpha = 1;
Q_weight = .01;
R_weight = 1-alpha;
P_weight = alpha;

P_bounds = [P_min_bound, P_max_bound];
R_bounds = [R_min_bound, R_max_bound];
Q_bounds = [Q_min_bound, Q_max_bound];

ctrl.setInPlaneBounds(P_bounds, Q_bounds, R_bounds);
ctrl.setCrossTrackBounds(P_bounds, Q_bounds, R_bounds);

output = ctrl.solve_push('Q_weight', Q_weight, 'R_weight', R_weight, 'P_weight', P_weight);


warning('on')
if any(eig(ctrl.InPlane.Q) < 0) || any(eig(ctrl.InPlane.P) < 0) || any(eig(ctrl.InPlane.R) < 0) || any(eig(ctrl.CrossTrack.Q) < 0) || any(eig(ctrl.CrossTrack.P) < 0) || any(eig(ctrl.CrossTrack.R) < 0)
    warning('Infeasible Problem')
end




%%

N_sims = 2;


init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;

sim_idx = randsample(64,N_sims);

prop.set_scaling(M2KM, S2HR)


sims_all = cell(1,N_sims);


for i = 1:N_sims
    
    fprintf('Sim: %d/%d\n',i,N_sims);
    
    idx = sim_idx(i);
    init_cond = init_cond_idx(idx,:);
    
    x0_lvlh = init_cond'.*[390; 390; 390; 1.5; 1.5; 1.5];
    x0_eci_0 = RTN_to_ECI(x0_lvlh, ref_eci_0);
    
    % Add mass terms
    ref_eci = [ref_eci_0; mass];
    x0_eci = [x0_eci_0; mass];
    
    x0 = [ref_eci; x0_eci];




    t0 = 0;
    tf = 6*3600;
    dt = .5;

    t_span = t0:dt:tf;

    [t_all, x_all] = prop.prop_2_piece(t_span, x0, @ctrl.ctrl);
    
    
    
    
    n_steps = size(t_all(:),1);
    u_all = NaN(3,n_steps);
    u_comp_all = NaN(4,n_steps);
    lvlh_error = NaN(6,n_steps);
    
    for k = 1:n_steps
        t = t_all(k);
        x = x_all(k,:)';
        
        x_ref = x_all(k,1:7)';
        x_sat = x_all(k,8:14)';
        
        [r_rel_x, v_rel_x, ~] = rva_relative(x_ref(1:3)',x_ref(4:6)',x_sat(1:3)',x_sat(4:6)',mu);
        lvlh_error(:,k) = [r_rel_x; v_rel_x];        
        
        x_error_KM_HR = lvlh_error(:,k).*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];
        
        [u_step, ~, u_comp] = ctrl.ctrl(x_error_KM_HR);
        
        u_comp_all(:,k) = [u_comp.InPlane(:); u_comp.CrossTrack(:)]*(S2HR^2)/M2KM;
        u_all(:,k) = u_step*(S2HR^2)/M2KM;
        
    end
    
    sim.x_all = lvlh_error;
    sim.u_all = u_all;
    sim.t_all = t_all;

    sim.InPlane.Q = ctrl.InPlane.Q;
    sim.InPlane.P = ctrl.InPlane.P;
    sim.InPlane.R = ctrl.InPlane.R;

    sim.CrossTrack.Q = ctrl.CrossTrack.Q;
    sim.CrossTrack.P = ctrl.CrossTrack.P;
    sim.CrossTrack.R = ctrl.CrossTrack.R;
    
    sims_all{i} = sim;
    
    
    
    
end




save(filename, "sims_all")









