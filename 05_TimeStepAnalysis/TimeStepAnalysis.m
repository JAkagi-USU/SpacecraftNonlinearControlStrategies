%% Parameters
mu = 3.9860044188e14;
mass = 24;
N_dt = 12;
N_inc = 4;

N_sims = N_dt*N_inc;

ref_el = [6878000; 1e-4; deg2rad(25); deg2rad(45); 0; deg2rad(0)];

orbit_period = 2*pi*sqrt(ref_el(1)^3/mu);

T_full = orbit_period;

dt_all = 10.^linspace(-2,2,N_dt);
inc_all = linspace(0,deg2rad(90),N_inc);

[dt_grid, inc_grid] = ndgrid(dt_all, inc_all);

N_steps_max = floor(T_full/min(dt_all))+2;

prop = orbit_prop();

traj_plot = figure();
hold on




%% Init data collecting




states_all = NaN(7,N_steps_max,N_sims);
t_all   = NaN(1,N_steps_max,N_sims);



for j = 1:N_sims
    %% Initial Conditions
    dt = dt_grid(j);

    fprintf('%d: %.3e s\n',j, dt)

    % Reference Orbit
    ref_el(3) = inc_grid(j);
    [r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
    states_all(:,1,j) = [r; v; mass];
    t_all(:,1,j) = 0;

    

    

    N_steps = floor(T_full/dt);
    mod_dt = mod(T_full,dt);

    counter = 10;

    for i = 1:N_steps

        perc = i*100/N_steps;
        if perc >= counter
            fprintf('%.0f%% Done\n', perc)
            counter = counter + 10;
        end

        % Update Dynamics
        x0_eci = prop.prop(dt, states_all(:,i,j), zeros(3,1));
        states_all(:,i+1,j) = x0_eci;

        t_all(:,i+1,j) = t_all(:,i,j) + dt;

    end    

    if mod_dt > 0
        x0_eci = prop.prop(mod_dt, states_all(:,i+1,j), zeros(3,1));
        states_all(:,i+2,j) = x0_eci;
        t_all(:,i+2,j) = t_all(:,i+1,j) + mod_dt;
    end


    figure(traj_plot)
    plot3(states_all(1,:,j), states_all(2,:,j), states_all(3,:,j))


    [row, col] = find(~isnan(states_all(:,:,j)),1,'last');
    states_all(:,1,j) - states_all(:,col,j)


    drawnow;
end

%%
error_fig = figure('DefaultAxesFontSize',12,'units','inch','position',[1,1,5,4]);



ref_idx = 1;
for j = 2:N_sims

    

    [~, ref_col] = find(~isnan(states_all(:,:,ref_idx)),1,'last');

[~, col] = find(~isnan(states_all(:,:,j)),1,'last');

state_error = norm(states_all(1:3,ref_col,ref_idx) - states_all(1:3,col,j));

loglog(dt_grid(j), state_error, 'ko')
    hold on

    if mod(j,N_dt) == 0
        ref_idx = j+1;
    end

end

grid on

xlabel('Simulation dt (s)')
ylabel('Position Error (m)')






