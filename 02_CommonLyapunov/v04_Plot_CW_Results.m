
% Select data file
load data/Dissertation_CW_CommonLyap.mat
% load data/CW_Sims.mat

save_file = true;

N_sims = size(sims_all,2);

inplane_log = figure('DefaultAxesFontSize',12);
crosstrack_log = figure('DefaultAxesFontSize',12);
inplane_u_log = figure('DefaultAxesFontSize',12);
crosstrack_u_log = figure('DefaultAxesFontSize',12);


for j = 1:N_sims
    
    x_single = sims_all{j}.x_all;
    u_single = sims_all{j}.u_all;
    t_single = sims_all{j}.t_all;

    
    figure(inplane_log);
    hold on
    plot(t_single/3600,vecnorm(x_single(:,1:2),2,2),'k-')
    
    figure(crosstrack_log);
    hold on
    plot(t_single/3600,abs(x_single(:,3)),'k-')
    
    figure(inplane_u_log);
    hold on
    plot(t_single(1:end-1)/3600,vecnorm(u_single(:,1:2),2,2),'k-')
    
    figure(crosstrack_u_log);
    hold on
    plot(t_single(1:end-1)/3600,abs(u_single(:,3)),'k-')
    

    
end

IP_bound = sims_all{1}.InPlane.uBound;
CT_bound = sims_all{1}.InPlane.uBound;


figure(inplane_log)
title('In-Plane Position Error')
ylabel('Error (m)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(inplane_log, 'figs/HCW_NoPert_IP_Log.png')
    saveas(inplane_log, 'figs/HCW_NoPert_IP_Log.eps')
end


figure(crosstrack_log)
title('Cross-Track Position Error')
ylabel('Error (m)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(crosstrack_log, 'figs/HCW_NoPert_CT_Log.png')
    saveas(crosstrack_log, 'figs/HCW_NoPert_CT_Log.eps')
end


figure(inplane_u_log)
plot([t_single(1) t_single(end)]/3600, [IP_bound IP_bound], 'r--')
title('In-Plane Control Norm')
ylabel('|u| (m/s^2)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on

if save_file
    saveas(inplane_u_log, 'figs/HCW_NoPert_IP_u_Log.png')
    saveas(inplane_u_log, 'figs/HCW_NoPert_IP_u_Log.eps')
end

figure(crosstrack_u_log)
plot([t_single(1) t_single(end)]/3600, [CT_bound CT_bound], 'r--')
title('Cross-Track Control Norm')
ylabel('|u| (m/s^2)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(crosstrack_u_log, 'figs/HCW_NoPert_CT_u_Log.png')
    saveas(crosstrack_u_log, 'figs/HCW_NoPert_CT_u_Log.eps')
end

