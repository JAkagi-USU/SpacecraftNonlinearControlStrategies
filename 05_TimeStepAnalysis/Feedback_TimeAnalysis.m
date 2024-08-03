clear all

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



% poles_d = -[.9 .85 .8 .75 .7 .65];
poles_d = -[.06 .05 .04 .03 .02 .01];
K_cont = -place(A_cont, B_cont, poles_d);


T = 3600;
T = 5400;


dt = 0;
prop = orbit_prop_2Sat(dt);




timesteps_all = [.001 .005 .01 .05 .1 .5 1 5 10];
% timesteps_all = [.5 10];
t_save = 0:max(timesteps_all):T;

N_sims = numel(timesteps_all);


eci_error = NaN(6,numel(t_save),N_sims);




% x_sim = cell(1,N_sims);
% t_sim = cell(1,N_sims);
% u_sim = cell(1,N_sims);
% v_sim = cell(1,N_sims);

ctrl = FeedbackCtrl();



for j = 1:N_sims

    fprintf('Sim: %d/%d\n',j,N_sims);

    
    
    % Initial conditions
    x0_all = x0_base;
    x0_eci_0 = RTN_to_ECI(x0_all, ref_eci_0);
    
    ref_eci = [ref_eci_0; mass];
    x0_eci = [x0_eci_0; mass];
    x0 = [ref_eci; x0_eci];

    dt = timesteps_all(j);
    t_span = 0:dt:T;

    AD = [AD_rr(n,dt), AD_vr(n,dt);
        AD_rv(n,dt), AD_vv(n,dt)];
    BD = BD_f(n, dt);

    K_disc = -place(AD, BD, poles_d);

%     ctrl.setK(K_disc);
    ctrl.setK(K_cont);
%     ctrl.setK(K_cont*0);



    
    [t_all, x_all] = prop.prop_2_piece(t_span, x0, ctrl);








n_steps = size(t_all(:),1);


for k = 1:n_steps
    t = t_all(k);
    x = x_all(k,:)';

    x_ref = x_all(k,1:7)';
    x_sat = x_all(k,8:14)';

    save_idx = find(t == t_save);
    if save_idx
%         eci_error(:,save_idx,j) = x_ref(1:6) - x_sat(1:6);
        eci_error(:,save_idx,j) = x_sat(1:6);
    end

end













end


%%


% error_val = NaN(N_sims-1,numel(t_save));

colorlist = [165,0,38;
%              215,48,39;
             244,109,67;
%              253,174,97;
             254,224,144;
%              224,243,248;
             171,217,233;
%              116,173,209;
             69,117,180;
             49,54,149]/255;

% colorlist = [204,236,230;
%              153,216,201;
%              102,194,164;
%              65,174,118;
%              35,139,69;
%              0,109,44;
%              0,68,27]/255;

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
% ax.ColorOrder = colorlist;
set(gca,'ColorOrder',colorlist,'nextplot','replacechildren')

for j = 2:N_sims

    cmp_error = vecnorm(eci_error(1:3,:,1) - eci_error(1:3,:,j));

%     cmp_error(cmp_error == 0) = 1e-100;
%     cmp_error(1) = 0;
line_name = sprintf("%.3f s",timesteps_all(j));
semilogy(t_save/60, cmp_error, 'DisplayName', line_name, 'LineWidth',1.5)
hold on

end

ylim([10^-10 10^3]) % Discrete
% ylim([10^-6 10^3]) % Continuous

grid on
legend('Location','best','NumColumns',2)
% legend('Location','southwest')

xlabel('Time (min)')
ylabel('Position Error (m)')

% plot(t_save, error_val)





















