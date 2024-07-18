
filename = "data/PQR_Survey.mat";
load(filename)

savefigs = false;

%% Calc K
M2KM = 1/1000;
S2HR = 1/3600;

x0 = [390; 390; 390; 1.5; 1.5; 1.5];
x0 = x0.*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];

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


% NEEDS TO MATCH THE ORIGINAL SIM DATA TO CORRECTLY CALCULATE U BOUND
QR_ratio = 100;
Q_val = 1;
R_val = Q_val/QR_ratio;

ctrl = CommonLyap();
ctrl.setInPlaneDyn(InPlane.A0, InPlane.B0, InPlane.TC, InPlane.TD);
ctrl.setCrossTrackDyn(CrossTrack.A0, CrossTrack.B0, CrossTrack.TC, CrossTrack.TD);

[C_InPlane, ~, ~] = lqr(ctrl.InPlane.AZ, ctrl.InPlane.BZ, Q_val*eye(4), R_val*eye(1));
[C_CrossTrack, ~, ~] = lqr(ctrl.CrossTrack.AZ, ctrl.CrossTrack.BZ,  Q_val*eye(2), R_val*eye(1));
    

%%

n_sims = size(sims_all,2);

conv_time_IP = NaN(1,n_sims);
conv_time_CT = NaN(1,n_sims);
max_ctrl_IP = NaN(1,n_sims);
max_ctrl_CT = NaN(1,n_sims);
T0_all_IP = NaN(1,n_sims);
T0_all_CT = NaN(1,n_sims);

conv_bound = .1;

n_total_sims = 0;

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    x_all = sims_all{i}.x_all;
    t_all = sims_all{i}.t_all;
    
    
    in_bounds_IP = vecnorm(x_all(:,1:2),2,2) <= conv_bound;
    in_bounds_CT = vecnorm(x_all(:,3),2,2) <= conv_bound;
    
    if in_bounds_IP(end) == 1
        in_bound_idx = find(diff(in_bounds_IP) ~= 0,1,'last') + 1;
        conv_time_IP(i) = t_all(in_bound_idx);
    else
        conv_time_IP(i) = Inf;
    
    end
    max_ctrl_IP(i) = max(abs(sims_all{i}.u_all(:,2)));
    
    x = x0([1;4;2;5]);
    Bx = InPlane.TD\x;
    T0_all_IP(i) = CommonLyap.calc_T(Bx, sims_all{i}.InPlane.Q);
    
    
     if in_bounds_CT(end) == 1
        in_bound_idx = find(diff(in_bounds_CT) ~= 0,1,'last') + 1;
        conv_time_CT(i) = t_all(in_bound_idx);
        
     else
         conv_time_CT(i) = Inf;
     end
     max_ctrl_CT(i) = max(abs(sims_all{i}.u_all(:,3)));
     
     
     z = x0([3;6]);
     Bz = CrossTrack.TD\z;
     T0_all_CT(i) = CommonLyap.calc_T(Bz, sims_all{i}.CrossTrack.Q);
    
     n_total_sims = n_total_sims + 1;
end

cmap = jet;







%% Conv Time Against Q, R/P, IP

figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    RP_IP = sims_all{i}.R_IP(2)/sims_all{i}.P_IP(1);
    Q_IP = sims_all{i}.Q_IP(1);
    
    if conv_time_IP(i) < Inf
        scatter(RP_IP, Q_IP, [], log10(conv_time_IP(i)/3600), 'filled')
    else
        scatter(RP_IP, Q_IP, 'k.')
    end

end


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlim([30 10^5])
xticks([10^2, 10^3, 10^4, 10^5])

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')

title('In-Plane Convergence Time (hr, log_{10})')

grid on

colorbar

if savefigs
    saveas(gcf, '03_figs/12_figs/IP_ConvTime', 'epsc')
    saveas(gcf, '03_figs/12_figs/IP_ConvTime', 'png')
end

%% Conv Time Max Against Q, R/P, IP

figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

TF_max = 20;

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    RP_IP = sims_all{i}.R_IP(2)/sims_all{i}.P_IP(1);
%     PR_IP = sims_all{i}.R_IP(1)/sims_all{i}.P_IP(2);
    Q_IP = sims_all{i}.Q_IP(1);
    
    TF = T0_all_IP(i)*RP_IP;
    

    scatter(RP_IP, Q_IP, [], log10(TF), 'filled')


end


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlim([30 10^5])
xticks([10^2, 10^3, 10^4, 10^5])

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')

title('In-Plane Max Convergence Time (hr, log_{10})')

grid on

colorbar

if savefigs
    saveas(gcf, '03_figs/12_figs/IP_ConvBound', 'epsc')
    saveas(gcf, '03_figs/12_figs/IP_ConvBound', 'png')
end


%% Conv Time Against Q, R/P, CT


figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
%     RP_IP = sims_all{i}.R_IP(1)/sims_all{i}.P_IP(2);
    RP_CT = sims_all{i}.R_CT(2)/sims_all{i}.P_CT(1);
%     Q_IP = sims_all{i}.Q_IP(1);
    Q_CT = sims_all{i}.Q_CT(1);
    
    if conv_time_CT(i) < Inf
%     scatter(RP_IP, Q_IP, [], conv_time_IP(i)/3600, 'filled')
    scatter(RP_CT, Q_CT, [], log10(conv_time_CT(i)/3600), 'filled')
    else
        scatter(RP_CT, Q_CT, 'k.')
    end

end

xlim([4,20])

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')
title('Cross-Track Convergence Time (hr, log_{10})')

grid on

colorbar

if savefigs
    saveas(gcf, '03_figs/12_figs/CT_ConvTime', 'epsc')
    saveas(gcf, '03_figs/12_figs/CT_ConvTime', 'png')
end

%% Max Conv Time Against Q, R/P, CT


figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    RP_CT = sims_all{i}.R_CT(2)/sims_all{i}.P_CT(1);
    Q_CT = sims_all{i}.Q_CT(1);
    
    TF = T0_all_CT(i)*RP_CT;
    
    scatter(RP_CT, Q_CT, [], log10(TF), 'filled')


end

xlim([4,20])

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')
title('Cross-Track Max Convergence Time (hr, log_{10})')

grid on

colorbar

if savefigs
    saveas(gcf, '03_figs/12_figs/CT_ConvBound', 'epsc')
    saveas(gcf, '03_figs/12_figs/CT_ConvBound', 'png')
end


%% Max Ctrl Against Q, R/P, IP

figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

count_IP = 0;

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    RP_IP = sims_all{i}.R_IP(2)/sims_all{i}.P_IP(1);
    Q_IP = sims_all{i}.Q_IP(1);

    u_bound = norm(C_InPlane)/sqrt(Q_IP);
    u_bound_m_s = u_bound*(S2HR^2)/M2KM;
    
        scatter(RP_IP, Q_IP, [], max_ctrl_IP(i), 'filled')
        
        if max_ctrl_IP(i) > u_bound_m_s
            fprintf('Sim: %d Max Bound: %.3e Max Control: %.3e\n', i, u_bound_m_s, max_ctrl_IP(i))
            count_IP = count_IP+1;
            scatter(RP_IP, Q_IP, [],'ro')
        end


end


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlim([30 10^5])
xticks([10^2, 10^3, 10^4, 10^5])

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')

title('In-Plane Max Control (m/s^2)')

grid on

colorbar

if savefigs
    saveas(gcf, '03_figs/12_figs/IP_CtrlMax', 'epsc')
    saveas(gcf, '03_figs/12_figs/IP_CtrlMax', 'png')
end

%% Ctrl Bound Against Q, R/P, IP

figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    RP_IP = sims_all{i}.R_IP(2)/sims_all{i}.P_IP(1);
    Q_IP = sims_all{i}.Q_IP(1);
    
    u_bound = norm(C_InPlane)/sqrt(Q_IP);
    u_bound_m_s = u_bound*(S2HR^2)/M2KM;
    
        scatter(RP_IP, Q_IP, [], u_bound_m_s, 'filled')


end


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlim([30 10^5])
xticks([10^2, 10^3, 10^4, 10^5])

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')

title('In-Plane Control Bound (m/s^2)')

grid on

colorbar

if savefigs
    saveas(gcf, '03_figs/12_figs/IP_CtrlBound', 'epsc')
    saveas(gcf, '03_figs/12_figs/IP_CtrlBound', 'png')
end

%% Max Ctrl Against Q, R/P, CT


figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
%     RP_IP = sims_all{i}.R_IP(1)/sims_all{i}.P_IP(2);
    RP_CT = sims_all{i}.R_CT(2)/sims_all{i}.P_CT(1);
%     Q_IP = sims_all{i}.Q_IP(1);
    Q_CT = sims_all{i}.Q_CT(1);
    
    u_bound = norm(C_CrossTrack)/sqrt(Q_CT);
    u_bound_m_s = u_bound*(S2HR^2)/M2KM;
    
%     scatter(RP_IP, Q_IP, [], conv_time_IP(i)/3600, 'filled')
    scatter(RP_CT, Q_CT, [], max_ctrl_CT(i), 'filled')
    
    if max_ctrl_CT(i) > u_bound_m_s
        scatter(RP_CT, Q_CT, 'ro')
    end


end

xlim([4, 20])

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')

title('Cross-Track Max Control (m/s^2)')

grid on

colorbar


if savefigs
    saveas(gcf, '03_figs/12_figs/CT_CtrlMax', 'epsc')
    saveas(gcf, '03_figs/12_figs/CT_CtrlMax', 'png')
end


%% Ctrl Bound Against Q, R/P, CT


figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4])
hold on

for i = 1:n_sims
    
    if isempty(sims_all{i})
        continue
    end
    
    RP_CT = sims_all{i}.R_CT(2)/sims_all{i}.P_CT(1);
    Q_CT = sims_all{i}.Q_CT(1);
    
    u_bound = norm(C_CrossTrack)/sqrt(Q_CT);
    u_bound_m_s = u_bound*(S2HR^2)/M2KM;
    
    scatter(RP_CT, Q_CT, [], u_bound_m_s, 'filled')


end

xlim([4, 20])

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

ylabel('\lambda_{min}(Q)')
xlabel('\lambda_{max}(R)/\lambda_{min}(P)')

title('Cross-Track Control Bound (m/s^2)')

grid on

colorbar

if savefigs
    saveas(gcf, '03_figs/12_figs/CT_CtrlBound', 'epsc')
    saveas(gcf, '03_figs/12_figs/CT_CtrlBound', 'png')
end





