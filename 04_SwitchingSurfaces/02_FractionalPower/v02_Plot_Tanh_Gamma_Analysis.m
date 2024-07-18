

save_plots = true;

filename = "data/TanhGammaAnalysis.mat";
load(filename);



%% Final error as function of gamma



final_error_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);

N_sims = numel(T_grid);

for i = 1:N_sims
    
    x_all = x_sim{i};
    t_all = t_sim{i};
    T = T_grid(i);
    gamma = tanh_grid(i);
    
    [~, t_idx] = min(abs(t_all - T));
    T_error = norm(x_all(t_idx,1:3));
    
    
    loglog(gamma, T_error, 'ko')
    hold on
    
    
    
end

xlim([10^-5, 10^5])
xticks([10^-4, 10^-2, 10^0, 10^2, 10^4])

grid on
ylabel('Error (m)')
xlabel('\gamma')


if save_plots
    saveas(gcf,'figs/Tanh_Gamma_Error','epsc');
    saveas(gcf,'figs/Tanh_Gamma_Error','png');
end


%% Final ctrl as function of gamma

final_error_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);

N_sims = numel(T_grid);

for i = 1:N_sims
    
    u_all = u_sim{i};
    t_all = t_sim{i};
    T = T_grid(i);
    gamma = tanh_grid(i);
    
    [~, t_idx] = min(abs(t_all - T));
    T_error = norm(u_all(1:3, t_idx));
    
    
    loglog(gamma, T_error, 'ko')
    hold on
    
    
    
end

xlim([10^-5, 10^5])
xticks([10^-4, 10^-2, 10^0, 10^2, 10^4])

grid on
ylabel('|u|_2 (m/s)')
xlabel('\gamma')

if save_plots
    saveas(gcf,'figs/Tanh_Gamma_Ctrl','epsc');
    saveas(gcf,'figs/Tanh_Gamma_Ctrl','png');
end
