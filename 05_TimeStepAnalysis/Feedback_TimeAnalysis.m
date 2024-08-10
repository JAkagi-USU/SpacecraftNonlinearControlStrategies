% Generates the plots for the timestep analysis comparison.
% The three different approaches (uncontrolled, discrete, and continuous)
% can be changed by choosing which K feedback gain is applied to the
% control law

clear all
addpath utils


mu = 3.9860044188e14;

x0_base = [390; 390; 390; 1.5; 1.5; 1.5];
mass = 24;

% Reference Orbit
ref_el = [6878000; 1e-4; deg2rad(25); deg2rad(45); 0; deg2rad(0)];
[r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
ref_eci_0 = [r; v];

n = sqrt(mu/ref_el(1)^3);


% Define the continuous CW dynamic matrices
A_cont = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];

B_cont = [0 0 0;
    0 0 0;
    0 0 0;
    1 0 0;
    0 1 0;
    0 0 1];

% Define the discrete CW matrices as a function of n and dt
AD_rr = @(n,dt) [4-3*cos(n*dt), 0, 0; 6*(sin(n*dt)-n*dt), 1, 0; 0, 0, cos(n*dt)];
AD_vr = @(n,dt) [sin(n*dt)/n, 2*(1-cos(n*dt))/n, 0; 2*(cos(n*dt)-1)/n,(4*sin(n*dt)-3*n*dt)/n,0;0,0,sin(n*dt)/n];
AD_rv = @(n,dt) [3*n*sin(n*dt), 0, 0; 6*n*(cos(n*dt)-1),0,0; 0,0,-n*sin(n*dt)];
AD_vv = @(n,dt) [cos(n*dt), 2*sin(n*dt),0; -2*sin(n*dt),4*cos(n*dt)-3,0;0,0,cos(n*dt)];

BD_f = @(n,dt) [2*(sin(n*dt/2)/n)^2, -2*(sin(n*dt) - n*dt)/(n^2), 0;
    2*(sin(n*dt) - n*dt)/(n^2), -(8*cos(n*dt) + 3*n^2*dt^2-8)/(2*n^2),0;
    0,0,2*(sin(n*dt/2)/n)^2;
    sin(n*dt)/n, -(2*cos(n*dt)-2)/n,0;
    -4*sin(n*dt/2)^2/n, 4*sin(n*dt)/n-3*dt,0;
    0,0,sin(n*dt)/n];


% Set poles for both the closed loop feedback gain
poles_d = -[.06 .05 .04 .03 .02 .01];

% Calculate continuous control K
K_cont = -place(A_cont, B_cont, poles_d);

% Set simulation time
T = 5400;

% Initialize propagator and control
prop = orbit_prop_2Sat();
ctrl = FeedbackCtrl();




timesteps_all = [.001 .005 .01 .05 .1 .5 1 5 10];
t_save = 0:max(timesteps_all):T;

N_sims = numel(timesteps_all);


eci_sim = NaN(6,numel(t_save),N_sims);






for j = 1:N_sims

    fprintf('Sim: %d/%d\n',j,N_sims);

    % Initial conditions
    x0_eci_0 = RTN_to_ECI(x0_base, ref_eci_0);
    ref_eci = [ref_eci_0; mass];
    x0_eci = [x0_eci_0; mass];
    x0 = [ref_eci; x0_eci];

    dt = timesteps_all(j);
    t_span = 0:dt:T;

    % Initialize discrete AD and BD matrices
    AD = [AD_rr(n,dt), AD_vr(n,dt);
        AD_rv(n,dt), AD_vv(n,dt)];
    BD = BD_f(n, dt);

    % Calculate discrete control K
    K_disc = -place(AD, BD, poles_d);


    % Choose which approach to use for the analysis
%     ctrl.setK(K_disc); % Discrete control law
%     ctrl.setK(K_cont); % Continuous control law
    ctrl.setK(K_cont*0); % Uncontrolled simulation

    % Propagate simulation
    [t_all, x_all] = prop.prop_2_piece(t_span, x0, ctrl);

    n_steps = size(t_all(:),1);

    % Save off trajectory points corresponding to the save timestep
    for k = 1:n_steps
        t = t_all(k);
        x = x_all(k,:)';

        x_ref = x_all(k,1:7)';
        x_sat = x_all(k,8:14)';

        save_idx = find(t == t_save);
        if save_idx
            eci_sim(:,save_idx,j) = x_sat(1:6);
        end

    end

end


%% Plot results



colorlist = [31,120,180;
    51,160,44;
    227,26,28;
    255,127,0;
    171,217,233;
    106,61,154;
    177,89,40;
    0,0,0]/255;


figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);
ax = gca;
set(gca,'ColorOrder',colorlist,'nextplot','replacechildren')

for j = 2:N_sims

    cmp_error = vecnorm(eci_sim(1:3,:,1) - eci_sim(1:3,:,j));

    line_name = sprintf("%.3f s",timesteps_all(j));
    semilogy(t_save/60, cmp_error, 'DisplayName', line_name, 'LineWidth',1.5)
    hold on

end

ylim([10^-10 10^3]) % Discrete
% ylim([10^-6 10^3]) % Continuous

grid on
legend('Location','best','NumColumns',2)

xlabel('Time (min)')
ylabel('Position Error (m)')






















