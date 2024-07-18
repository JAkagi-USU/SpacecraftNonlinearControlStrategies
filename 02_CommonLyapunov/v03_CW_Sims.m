addpath utils


file_input = "CW_Sims";
filename = sprintf("data/%s.mat",file_input);


M2KM = 1/1000;
S2HR = 1/3600;

n = .0011/S2HR;
a = 2*n;
b = 3*n^2;

InPlane.A0 = [0 1 0 0; b 0 0 a; 0 0 0 1; 0 -a 0 0];
InPlane.B0 = [0;0;0;1];

InPlane.TC = [0, n^2/a, 0, 0];
InPlane.TD = [0 a 0 0; 0 0 a 0; -b 0 1 0; 0 -b 0 1];

CrossTrack.A0 = [0 1; -n^2 0];
CrossTrack.B0 = [0; 1];
CrossTrack.TC = [n^2 0];
CrossTrack.TD = [1 0; 0 1];


ctrl = CommonLyap();
ctrl.setInPlaneDyn(InPlane.A0, InPlane.B0, InPlane.TC, InPlane.TD);
ctrl.setCrossTrackDyn(CrossTrack.A0, CrossTrack.B0, CrossTrack.TC, CrossTrack.TD);




%% Set sim parameters

N_sims = 50;


% Set LQR weights for K feedback gain
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






%% Run sims











state_fig = figure();
ctrl_fig = figure();

sims_all = cell(1,N_sims);


init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;

sim_idx = randsample(64,N_sims);


t0 = 0;
tf = 6;
dt = .5*S2HR;

t_span = t0:dt:tf;
n_steps = numel(t_span);



for i = 1:N_sims
    
    fprintf('\nSim %d/%d\n',i,N_sims);
    
    idx = sim_idx(i);
    init_cond = init_cond_idx(idx,:);
    
    x0 = [390; 390; 390; 1.5; 1.5; 1.5];
    x0 = init_cond'.*x0.*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];
    

x_all = NaN(n_steps,6);
u_all = NaN(n_steps-1,3);
t_all = t_span;

x_all(1,:) = x0;
for j = 1:n_steps-1
    
    x = x_all(j,:)';
    t1 = t_span(j);
    t2 = t_span(j+1);
    
    
    u = ctrl.ctrl(x);
    u_all(j,:) = u;
    [t_ode, x_ode] = rk4(@(t,x) dyn_u(t,x,n,u), [t1 t2], x);
    
    x_all(j+1,:) = x_ode(end,:);
    
end

sim.x_all = x_all.*[1/M2KM 1/M2KM 1/M2KM S2HR/M2KM S2HR/M2KM S2HR/M2KM];
sim.u_all = u_all.*[S2HR^2/M2KM S2HR^2/M2KM S2HR^2/M2KM];
sim.t_all = t_span/S2HR;
sim.dt = dt/S2HR;

Q_IP_L = min(eig(ctrl.InPlane.Q));
u_bound_IP = norm(C_InPlane)/sqrt(Q_IP_L);
u_bound_m_s_IP = u_bound_IP*(S2HR^2)/M2KM;

Q_CT_L = min(eig(ctrl.CrossTrack.Q));
u_bound_CT = norm(C_CrossTrack)/sqrt(Q_CT_L);
u_bound_m_s_CT = u_bound_CT*(S2HR^2)/M2KM;


sim.InPlane.Q = ctrl.InPlane.Q;
sim.InPlane.P = ctrl.InPlane.P;
sim.InPlane.R = ctrl.InPlane.R;
sim.InPlane.obj = ctrl.InPlane.obj;

sim.InPlane.uBound = u_bound_m_s_IP;

sim.CrossTrack.Q = ctrl.CrossTrack.Q;
sim.CrossTrack.P = ctrl.CrossTrack.P;
sim.CrossTrack.R = ctrl.CrossTrack.R;
sim.CrossTrack.obj = ctrl.CrossTrack.obj;

sim.CrossTrack.uBound = u_bound_m_s_CT;

sim.Q_val = Q_val;
sim.R_val = R_val;

sim.Q_max_bound = Q_max_bound;
sim.P_max_bound = P_max_bound;
sim.R_max_bound = R_max_bound;

sim.Q_min_bound = Q_min_bound;
sim.P_min_bound = P_min_bound;
sim.R_min_bound = R_min_bound;

sim.Q_weight = Q_weight;
sim.R_weight = R_weight;
sim.P_weight = P_weight;


sims_all{i} = sim;

end


    save(filename, "sims_all")










function xdot = dyn_u(t, x, n, u)

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
      
      
      xdot = A*x + B*u;



end






























