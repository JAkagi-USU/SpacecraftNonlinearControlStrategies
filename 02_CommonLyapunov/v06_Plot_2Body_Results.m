
addpath utils

% Select data file
filename = "data/Dissertation_TwoBody_CommonLyap.mat";
% filename = "data/TwoBody_Sims.mat";
load(filename)

savefigs = true;

M2KM = 1/1000;
S2HR = 1/3600;

mu = 3.9860044188e14;
ref_el = [6878000; 1e-4; rad2deg(25); rad2deg(45); 0; rad2deg(0)];
[r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
ref_eci_0 = [r; v];

n = sqrt(mu/ref_el(1)^3);

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

R_val = .01;
Q_val = 1;
    
[C_InPlane, ~, ~] = lqr(ctrl.InPlane.AZ, ctrl.InPlane.BZ, Q_val*eye(4), R_val*eye(1));
[C_CrossTrack, ~, ~] = lqr(ctrl.CrossTrack.AZ, ctrl.CrossTrack.BZ,  Q_val*eye(2), R_val*eye(1));

%%
sim_num = 1;


Q_IP_L = min(eig(sims_all{sim_num}.InPlane.Q));
P_IP_L = min(eig(sims_all{sim_num}.InPlane.P));
R_IP_U = max(eig(sims_all{sim_num}.InPlane.R));

Q_CT_L = min(eig(sims_all{sim_num}.CrossTrack.Q));
P_CT_L = min(eig(sims_all{sim_num}.CrossTrack.P));
R_CT_U = max(eig(sims_all{sim_num}.CrossTrack.R));

RP_IP = R_IP_U/P_IP_L;
RP_CT = R_CT_U/P_CT_L;

u_bound_IP = norm(C_InPlane)/sqrt(Q_IP_L);
u_bound_m_s_IP = u_bound_IP*(S2HR^2)/M2KM;

u_bound_CT = norm(C_CrossTrack)/sqrt(Q_CT_L);
u_bound_m_s_CT = u_bound_CT*(S2HR^2)/M2KM;


%%

n_sims = size(sims_all,2);

conv_time_IP = NaN(1,n_sims);
conv_time_CT = NaN(1,n_sims);
conv_idx_IP = NaN(1,n_sims);
conv_idx_CT = NaN(1,n_sims);
max_ctrl_IP = NaN(1,n_sims);
max_ctrl_CT = NaN(1,n_sims);
T0_all_IP = NaN(1,n_sims);
T0_all_CT = NaN(1,n_sims);
TF_IP = NaN(1,n_sims);
TF_CT = NaN(1,n_sims);

conv_bound = .1;



n_total_sims = 0;

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    x_all = sims_all{i}.x_all;
    t_all = sims_all{i}.t_all;
    
    x0 = x_all(:,1).*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];
    
    
    in_bounds_IP = vecnorm(x_all(1:2,:)) <= conv_bound;
    in_bounds_CT = abs(x_all(3,:)) <= conv_bound;
    
    if in_bounds_IP(end) == 1
        in_bound_idx = find(diff(in_bounds_IP) ~= 0,1,'last') + 1;
        conv_idx_IP(i) = in_bound_idx;
        conv_time_IP(i) = t_all(in_bound_idx);
    else
        conv_time_IP(i) = Inf;
        conv_idx_IP(i) = numel(t_all);
    
    end
    max_ctrl_IP(i) = max(abs(sims_all{i}.u_all(:,2)));
    
    x = x0([1;4;2;5]);
    Bx = InPlane.TD\x;
    T0_all_IP(i) = CommonLyap.calc_T(Bx, sims_all{i}.InPlane.Q);
    TF_IP(i) = T0_all_IP(i)*RP_IP;
    
    
     if in_bounds_CT(end) == 1
        in_bound_idx = find(diff(in_bounds_CT) ~= 0,1,'last') + 1;
        conv_idx_CT(i) = in_bound_idx;
        conv_time_CT(i) = t_all(in_bound_idx);
        
     else
         conv_time_CT(i) = Inf;
         conv_idx_CT(i) = numel(t_all);
     end
     max_ctrl_CT(i) = max(abs(sims_all{i}.u_all(:,3)));
     
     
     z = x0([3;6]);
     Bz = CrossTrack.TD\z;
     T0_all_CT(i) = CommonLyap.calc_T(Bz, sims_all{i}.CrossTrack.Q);
     TF_CT(i) = T0_all_CT(i)*RP_CT;
    
     n_total_sims = n_total_sims + 1;
end

fprintf('IP Convergence Bound: %.3f +/- %.3f\n', mean(TF_IP), std(TF_IP));
fprintf('CT Convergence Bound: %.3f +/- %.3f\n', mean(TF_CT), std(TF_CT));

%% All error and control plots

IP_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);
CT_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);

IP_u_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);
CT_u_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);

IP_max = 0;
CT_max = 0;

for i = 1:n_sims
    
    x_all = sims_all{i}.x_all;
    t_all = sims_all{i}.t_all;
    u_all = sims_all{i}.u_all;
    
    figure(IP_fig)
    hold on
    IP_error = vecnorm(x_all(1:2,:));
%     plot(t_all(1:conv_idx_IP(i))/3600, IP_error(1:conv_idx_IP(i)), 'k')
    plot(t_all/3600, IP_error, 'k')
    
    figure(CT_fig)
    hold on
    CT_error = abs(x_all(3,:));
%     plot(t_all(1:conv_idx_IP(i))/3600, CT_error(1:conv_idx_IP(i)), 'k')
    plot(t_all/3600, CT_error, 'k')
    

    figure(IP_u_fig)
    hold on
%     plot(t_all(1:conv_idx_IP(i))/3600, abs(u_all(2,1:conv_idx_IP(i))),'k')
    plot(t_all/3600, abs(u_all(2,:)),'k')
    
    if max(abs(u_all(2,:))) > IP_max
        IP_max = max(abs(u_all(2,:)));
    end

    figure(CT_u_fig)
    hold on
%     plot(t_all(1:conv_idx_IP(i))/3600, abs(u_all(3,1:conv_idx_IP(i))),'k')
    plot(t_all/3600, abs(u_all(3,:)),'k')
    
    if max(abs(u_all(3,:))) > CT_max
        CT_max = max(abs(u_all(3,:)));
    end
    
    
end

figure(IP_fig)
set(gca, 'YScale', 'log')
ylabel('Error (m)')
xlabel('Time (hr)')
title('In-Plane Error')
grid on

if savefigs
    saveas(gcf, 'figs/IP_Error', 'epsc')
    saveas(gcf, 'figs/IP_Error', 'png')
end

figure(CT_fig)
set(gca, 'YScale', 'log')
ylabel('Error (m)')
xlabel('Time (hr)')
title('Cross-Track Error')
grid on

if savefigs
    saveas(gcf, 'figs/CT_Error', 'epsc')
    saveas(gcf, 'figs/CT_Error', 'png')
end


figure(IP_u_fig)
plot([0 t_all(end)/3600], [u_bound_m_s_IP u_bound_m_s_IP], 'r--')

set(gca, 'YScale', 'log')
ylabel('Control (m/s^2)')
xlabel('Time (hr)')
title('In-Plane Control Magnitude')
grid on

if savefigs
    saveas(gcf, 'figs/IP_Ctrl', 'epsc')
    saveas(gcf, 'figs/IP_Ctrl', 'png')
end

figure(CT_u_fig)

plot([0 t_all(end)/3600], [u_bound_m_s_CT u_bound_m_s_CT], 'r--')

set(gca, 'YScale', 'log')
ylabel('Control (m/s^2)')
xlabel('Time (hr)')
title('Cross-Track Control Magnitude')
grid on

if savefigs
    saveas(gcf, 'figs/CT_Ctrl', 'epsc')
    saveas(gcf, 'figs/CT_Ctrl', 'png')
end






    
    