addpath utils


file_input = "PQR_Survey";
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


% Initialize the lower bounds on eigenvalues of Q
N_Q = 10;
N_Q = 2;
Q_max_all = 10.^(linspace(-1,2,N_Q));

% Initialize the balance between the P and R weights 
N_alpha = 10;
N_alpha = 2;
alpha_all = linspace(0,1,N_alpha);

% Set the LQR ratio for the feedback term
QR_ratio = 100;
Q_val = 1;

% Initialize the overall lower bounds
N_bounds = 10;
N_bounds = 2;
min_bound_all = linspace(1e-3, 1e1, N_bounds);

[Q_grid, alpha_grid, bound_grid] = ndgrid( Q_max_all, alpha_all, min_bound_all);

N_sims = numel(Q_grid);

max_bound = 1e5;











%% Run sims

bounds = calc_bound();

sims_fig = figure();
error_fig = figure();

sims_all = cell(1,N_sims);

for i = 1:N_sims
    
    fprintf('\nSim %d/%d\n',i,N_sims);
    
%     QR_ratio = QR_grid(i);
    R_val = Q_val/QR_ratio;
    
    [C_InPlane, ~, ~] = lqr(ctrl.InPlane.AZ, ctrl.InPlane.BZ, Q_val*eye(4), R_val*eye(1));
    [C_CrossTrack, ~, ~] = lqr(ctrl.CrossTrack.AZ, ctrl.CrossTrack.BZ,  Q_val*eye(2), R_val*eye(1));
    
    ctrl.setInPlaneGain(-C_InPlane);
    ctrl.setCrossTrackGain(-C_CrossTrack);


Q_max_bound = Q_grid(i);
P_max_bound = max_bound;
R_max_bound = max_bound;

Q_min_bound = bound_grid(i);
P_min_bound = bound_grid(i);
R_min_bound = bound_grid(i);



Q_weight = 1;
R_weight = 1-alpha_grid(i);
P_weight = alpha_grid(i);

P_bounds = [P_min_bound, P_max_bound];
R_bounds = [R_min_bound, R_max_bound];
Q_bounds = [Q_min_bound, Q_max_bound];

ctrl.setInPlaneBounds(P_bounds, Q_bounds, R_bounds);
ctrl.setCrossTrackBounds(P_bounds, Q_bounds, R_bounds);

output = ctrl.solve_push('Q_weight', Q_weight, 'R_weight', R_weight, 'P_weight', P_weight);

[Q_u_bound_IP, Q_l_bound_IP] = bounds.solve_eig(ctrl.InPlane.Q);
[Q_u_bound_CT, Q_l_bound_CT] = bounds.solve_eig(ctrl.CrossTrack.Q);

[R_u_bound_IP, R_l_bound_IP] = bounds.solve_eig(ctrl.InPlane.R);
[R_u_bound_CT, R_l_bound_CT] = bounds.solve_eig(ctrl.CrossTrack.R);

[P_u_bound_IP, P_l_bound_IP] = bounds.solve_eig(ctrl.InPlane.P);
[P_u_bound_CT, P_l_bound_CT] = bounds.solve_eig(ctrl.CrossTrack.P);

fprintf('Q: %.4f, R/P: %.4f\n',Q_u_bound_IP, R_u_bound_IP/P_l_bound_IP)
fprintf('Q: %.4f, R/P: %.4f\n',Q_u_bound_CT, R_u_bound_CT/P_l_bound_CT)

fprintf("%s\n",output.InPlane.out.info)
fprintf("%s\n",output.CrossTrack.out.info)

warning('on')
if any(eig(ctrl.InPlane.Q) < 0) || any(eig(ctrl.InPlane.P) < 0) || any(eig(ctrl.InPlane.R) < 0) || any(eig(ctrl.CrossTrack.Q) < 0) || any(eig(ctrl.CrossTrack.P) < 0) || any(eig(ctrl.CrossTrack.R) < 0)
    warning('Infeasible Problem')
    continue
end

% fprintf("Q: Det: %.3f Eig: %.3f\n", det(ctrl.InPlane.Q), max(abs(eig(ctrl.InPlane.Q))))
% fprintf("R: Det: %.3f Eig: %.3f\n", det(ctrl.InPlane.R), max(abs(eig(ctrl.InPlane.R))))
% fprintf("P: Det: %.3f Eig: %.3f\n\n", det(ctrl.InPlane.P), max(abs(eig(ctrl.InPlane.P))))

warning('on')
switch output.InPlane.out.problem
    case 1
        warning('Infeasible Problem')
        continue
        
    case 4
%         warning('In Plane Numerical Problem')
        
    case 5
%         warning('Cross Track Lack of Progress')
        
end

switch output.CrossTrack.out.problem
    case 1
        warning('Infeasible Problem')
        continue
        
    case 4
%         warning('Cross Track Numerical Problem')
        
    case 5
%         warning('Cross Track Lack of Progress')
        
end

% ctrl.check()





t0 = 0;
tf = 8;
dt = 5*S2HR;

t_span = t0:dt:tf;
n_steps = numel(t_span);


x0 = [390; 390; 390; 1.5; 1.5; 1.5];
x0 = x0.*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];


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

sim.Q_IP = [Q_l_bound_IP Q_u_bound_IP];
sim.Q_CT = [Q_l_bound_CT Q_u_bound_CT];

sim.P_IP = [P_l_bound_IP P_u_bound_IP];
sim.P_CT = [P_l_bound_CT P_u_bound_CT];

sim.R_IP = [R_l_bound_IP R_u_bound_IP];
sim.R_CT = [R_l_bound_CT R_u_bound_CT];

sim.InPlane.Q = ctrl.InPlane.Q;
sim.InPlane.P = ctrl.InPlane.P;
sim.InPlane.R = ctrl.InPlane.R;
sim.InPlane.obj = ctrl.InPlane.obj;

sim.CrossTrack.Q = ctrl.CrossTrack.Q;
sim.CrossTrack.P = ctrl.CrossTrack.P;
sim.CrossTrack.R = ctrl.CrossTrack.R;
sim.CrossTrack.obj = ctrl.CrossTrack.obj;

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









function [xdot, u] = dyn(t, x, n, ctrl)

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
      
      u = ctrl.ctrl(x);
      
      xdot = A*x + B*u;



end

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






























