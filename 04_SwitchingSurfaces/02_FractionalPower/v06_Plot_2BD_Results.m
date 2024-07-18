
% Select data file
load data/Dissertation_TwoBody_Results.mat
% load data/TwoBody_Simulations.mat

save_file = true;

N_sims = size(x_sim,2);

inplane_log = figure('DefaultAxesFontSize',12);
crosstrack_log = figure('DefaultAxesFontSize',12);
inplane_u_log = figure('DefaultAxesFontSize',12);
crosstrack_u_log = figure('DefaultAxesFontSize',12);


for j = 1:N_sims

    x_single = x_sim{j};
    u_single = u_sim{j};
    t_single = t_sim{j};

    
    figure(inplane_log);
    hold on
    plot(t_single,vecnorm(x_single(1:2,:)),'k-')
    
    figure(crosstrack_log);
    hold on
    plot(t_single,abs(x_single(3,:)),'k-')
    
    figure(inplane_u_log);
    hold on
    plot(t_single,vecnorm(u_single(1:2,:)),'k-')
    
    figure(crosstrack_u_log);
    hold on
    plot(t_single,abs(u_single(3,:)),'k-')
    

    
end


figure(inplane_log)
title('In-Plane Position Error')
ylabel('Error (m)')
xlabel('Time (s)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(inplane_log, 'figs/2BD_IP_Log.png')
    saveas(inplane_log, 'figs/2BD_IP_Log.eps')
end


figure(crosstrack_log)
title('Cross-Track Position Error')
ylabel('Error (m)')
xlabel('Time (s)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(crosstrack_log, 'figs/2BD_CT_Log.png')
    saveas(crosstrack_log, 'figs/2BD_CT_Log.eps')
end


figure(inplane_u_log)
title('In-Plane Control Norm')
ylabel('|u|_2 (m/s^2)')
xlabel('Time (s)')
set(gca, 'YScale', 'log')
grid on

if save_file
    saveas(inplane_u_log, 'figs/2BD_IP_u_Log.png')
    saveas(inplane_u_log, 'figs/2BD_IP_u_Log.eps')
end

figure(crosstrack_u_log)
title('Cross-Track Control Norm')
ylabel('|u| (m/s^2)')
xlabel('Time (s)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(crosstrack_u_log, 'figs/2BD_CT_u_Log.png')
    saveas(crosstrack_u_log, 'figs/2BD_CT_u_Log.eps')
end

