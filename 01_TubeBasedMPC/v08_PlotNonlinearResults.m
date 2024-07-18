
% Select data file
load data/Dissertation_TwoBody_MPC.mat
% load data/TwoBodySimResults.mat

save_file = true;

N = size(u_all_guide,3);
N_guide = size(u_all_guide,2);
dt = 10;

inplane_log = figure('DefaultAxesFontSize',12);
crosstrack_log = figure('DefaultAxesFontSize',12);
inplane_u_log = figure('DefaultAxesFontSize',12);
crosstrack_u_log = figure('DefaultAxesFontSize',12);
fig_tot_error = figure('DefaultAxesFontSize',12);


for j = 1:N


    
    figure(inplane_log);
    hold on
    plot((0:N_guide)*dt/60,vecnorm(total_errors_all(1:2,:,j)),'k-')
    
    figure(crosstrack_log);
    hold on
    plot((0:N_guide)*dt/60,abs(total_errors_all(3,:,j)),'k-')
    
    figure(inplane_u_log);
    hold on
    plot((1:N_guide)*dt/60,vecnorm(u_all(1:2,:,j)),'k-')
    
    figure(crosstrack_u_log);
    hold on
    plot((1:N_guide)*dt/60,abs(u_all(3,:,j)),'k-')
    
    
    figure(fig_tot_error)
    for i = 1:6

        subplot(6,1,i)
        hold on
        plot((0:N_guide)*dt/60,total_errors_all(i,:,j),'k-')

    end
    

    
end


figure(inplane_log)
title('In-Plane Position Error')
ylabel('Error (m)')
xlabel('Time (min)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(inplane_log, 'figs/2BD_IP_Log.png')
    saveas(inplane_log, 'figs/2BD_IP_Log.eps')
end


figure(crosstrack_log)
title('Cross-Track Position Error')
ylabel('Error (m)')
xlabel('Time (min)')
set(gca, 'YScale', 'log')
yticks([10^-6 10^-4 10^-2 10^0 10^2 10^4])
grid on


if save_file
    saveas(crosstrack_log, 'figs/2BD_CT_Log.png')
    saveas(crosstrack_log, 'figs/2BD_CT_Log.eps')
end



figure(inplane_u_log)
title('In-Plane Control Norm')
ylabel('|u|_2 (m/s^2)')
xlabel('Time (min)')
set(gca, 'YScale', 'log')
grid on

if save_file
    saveas(inplane_u_log, 'figs/2BD_IP_u_Log.png')
    saveas(inplane_u_log, 'figs/2BD_IP_u_Log.eps')
end

figure(crosstrack_u_log)
title('Cross-Track Control Norm')
ylabel('|u| (m/s^2)')
xlabel('Time (min)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(crosstrack_u_log, 'figs/2BD_CT_u_Log.png')
    saveas(crosstrack_u_log, 'figs/2BD_CT_u_Log.eps')
end


figure(fig_tot_error)
subplot(6,1,1)
title('Total Error')
ylabel('X (m)')
subplot(6,1,2)
ylabel('Y (m)')
subplot(6,1,3)
ylabel('Z (m)')
subplot(6,1,4)
ylabel('VX (m/s)')
subplot(6,1,5)
ylabel('VY (m/s)')
subplot(6,1,6)
ylabel('VZ (m/s)')
xlabel('Time (min)')

if save_file
    saveas(fig_tot_error, 'figs/2BD_full_error.png')
    saveas(fig_tot_error, 'figs/2BD_full_error.eps')
end
