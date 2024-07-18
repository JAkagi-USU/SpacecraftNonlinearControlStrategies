
load data/Dissertation_TwoBody_Results.mat
% load data/MinTimeCtrl_2BdDyn.mat

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
    plot(t_single/3600,vecnorm(x_single(1:2,:)),'k-')
    
    figure(crosstrack_log);
    hold on
    plot(t_single/3600,abs(x_single(3,:)),'k-')
    
    figure(inplane_u_log);
    hold on
    plot(t_single/3600,vecnorm(u_single(1:2,:)),'k-')
    
    figure(crosstrack_u_log);
    hold on
    plot(t_single/3600,abs(u_single(3,:)),'k-')
    

    
end


%% Phase Plot

phase_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,12,4]);

switch_f = @(x,v) -sqrt(2*abs(v*x)).*sign(x);
x_lin = linspace(-1500,1500,1000);

for j = 1:N_sims
    
    x_ode = x_sim{j}';
    v_ode = v_sim{j};
    
    figure(phase_fig)
    subplot(1,3,1)
    
    hold on
    plot(x_lin, switch_f(x_lin, v_ode(2)),'k:')
    plot(x_ode(:,1),x_ode(:,4),'k-')
    plot(x_ode(1,1), x_ode(1,4), 'k^')

    subplot(1,3,2)
    hold on

    plot(x_lin, switch_f(x_lin, v_ode(4)),'k:')
    plot(x_ode(:,2),x_ode(:,5),'k-')
    plot(x_ode(1,2), x_ode(1,5), 'k^')

    subplot(1,3,3)
    hold on

    plot(x_lin, switch_f(x_lin, v_ode(6)),'k:')
    plot(x_ode(:,3),x_ode(:,6),'k-')
    plot(x_ode(1,3), x_ode(1,6), 'k^')
end



subplot(1,3,1)
xlim([-500, 500])
ylim([-15, 15])

grid on

xlabel('X (m)')
ylabel('VX (m/s)')

subplot(1,3,2)
xlim([-500, 500])
ylim([-15, 15])

xlabel('Y (m)')
ylabel('VY (m/s)')

grid on


subplot(1,3,3)
xlim([-1100, 1100])
ylim([-4.5, 4.5])

xlabel('Z (m)')
ylabel('VZ (m/s)')

grid on

if save_file
    saveas(gcf,'figs/2BD_Phase','epsc');
    saveas(gcf,'figs/2BD_Phase','png');
end




figure(inplane_log)
title('In-Plane Position Error')
ylabel('Error (m)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(inplane_log, 'figs/2BD_IP_Log.png')
    saveas(inplane_log, 'figs/2BD_IP_Log.eps')
end


figure(crosstrack_log)
title('Cross-Track Position Error')
ylabel('Error (m)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(crosstrack_log, 'figs/2BD_CT_Log.png')
    saveas(crosstrack_log, 'figs/2BD_CT_Log.eps')
end


figure(inplane_u_log)
title('In-Plane Control Norm')
ylabel('|u|_2 (m/s^2)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on

if save_file
    saveas(inplane_u_log, 'figs/2BD_IP_u_Log.png')
    saveas(inplane_u_log, 'figs/2BD_IP_u_Log.eps')
end

figure(crosstrack_u_log)
title('Cross-Track Control Norm')
ylabel('|u| (m/s^2)')
xlabel('Time (hr)')
set(gca, 'YScale', 'log')
grid on


if save_file
    saveas(crosstrack_u_log, 'figs/2BD_CT_u_Log.png')
    saveas(crosstrack_u_log, 'figs/2BD_CT_u_Log.eps')
end

