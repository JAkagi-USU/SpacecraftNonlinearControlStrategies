addpath data

save_plots = true;

filename = "QA_Param_Compare.mat";
load(filename);


%% Converged Sims - Sign



error_sim(isnan(error_sim)) = inf;
stable = error_sim <= .01;
unstable = error_sim > .01;

stable_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,6,4]);
hold on
plot(alpha_grid(stable), q_grid(stable), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'DisplayName','Converged')
plot(alpha_grid(unstable), q_grid(unstable), 'k.', 'DisplayName', 'Unconverged')




ylim([0, 1.05])

legend('Location','eastoutside')

xlabel('\alpha')
ylabel('q')
title('System Convergence')


if save_plots
    saveas(gcf,'figs/Convergence_Plot','epsc');
    saveas(gcf,'figs/Convergence_Plot','png');
end

%% Convergence Time

error_sim(isnan(error_sim)) = inf;
conv_sim(isnan(conv_sim)) = inf;
error_stable = error_sim <= .01;
conv_stable = conv_sim < inf;

stable = conv_stable & error_stable;

stable_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);
hold on
scatter(alpha_grid(:), q_grid(:), 'k.')
scatter(alpha_grid(stable), q_grid(stable), [], conv_sim(stable), 'filled')

colorbar

xlabel('\alpha')
ylabel('q')


title('Convergence Time (s)')


if save_plots
    saveas(gcf,'figs/Convergence_Time','epsc');
    saveas(gcf,'figs/Convergence_Time','png');
end



%% Max Control

error_sim(isnan(error_sim)) = inf;
conv_sim(isnan(conv_sim)) = inf;
error_stable = error_sim <= .01;
conv_stable = conv_sim < inf;

stable = conv_stable & error_stable;% & (u_max_sim < 100);

u_max_under_idx = (u_max_sim < 100) & stable;
u_max_over_idx = (u_max_sim >= 100) & stable;



stable_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);
hold on
scatter(alpha_grid(:), q_grid(:), 'k.')
scatter(alpha_grid(u_max_under_idx), q_grid(u_max_under_idx), [], u_max_sim(u_max_under_idx), 'filled')
scatter(alpha_grid(u_max_over_idx), q_grid(u_max_over_idx), 'ro', 'filled')


ylim([0, 1.05])

colorbar

xlabel('\alpha')
ylabel('q')

title('Max Control (m/s^2)')


if save_plots
    saveas(gcf,'figs/Control_Plot','epsc');
    saveas(gcf,'figs/Control_Plot','png');
end


%% Terminal Error and Ctrl Plots

n_sims = numel(t_sim);
terminal_error = NaN(1,n_sims);
terminal_ctrl = NaN(1,n_sims);

for i = 1:n_sims
    [~, t_idx] = min(abs(t_sim{i} - 300));
    terminal_error(i) = norm(x_sim{i}(1:3,t_idx));
    terminal_ctrl(i) = norm(u_sim{i}(:,t_idx));
end

error_stable = terminal_error <= .01;

stable_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);
hold on
scatter(alpha_grid(:), q_grid(:), 'k.')
scatter(alpha_grid(error_stable), q_grid(error_stable), [], log10(terminal_error(error_stable)), 'filled')

colorbar

xlabel('\alpha')
ylabel('q')

title('Terminal Error (m, log_{10})')

if save_plots
    saveas(gcf,'figs/Terminal_Error','epsc');
    saveas(gcf,'figs/Terminal_Error','png');
end

stable_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);
hold on
scatter(alpha_grid(:), q_grid(:), 'k.')
scatter(alpha_grid(error_stable), q_grid(error_stable), [], log10(terminal_ctrl(error_stable)), 'filled')

colorbar

xlabel('\alpha')
ylabel('q')

title('Terminal Control (m/s^2, log_{10})')


if save_plots
    saveas(gcf,'figs/Terminal_Control','epsc');
    saveas(gcf,'figs/Terminal_Control','png');
end








