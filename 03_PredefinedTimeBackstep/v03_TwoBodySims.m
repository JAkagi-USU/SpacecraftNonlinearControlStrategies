addpath utils

M2KM = 1/1000;
S2HR = 1/3600;

x0_bounds = [370; 370; 370; 1.25; 1.25; 1.25];
x0_mpc = [20; 20; 20; .25; .25; .25];

x0_base = x0_bounds + x0_mpc;


mu = 3.9860044188e14;

mass = 24;


dt = 10;
guide_dt = dt*S2HR;

N_sim = 800;
N_end = N_sim + 100;
T = N_sim*dt;
TF = (N_end)*dt;




% Reference Orbit
ref_el = [6878000; 1e-4; rad2deg(25); rad2deg(45); 0; rad2deg(0)];
[r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
ref_eci_0 = [r; v];%.*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];

n = sqrt(mu/ref_el(1)^3);




prop = orbit_prop_2Sat(dt);
prop.set_scaling(M2KM, S2HR)



filename = "data/TwoBodySims.mat";


%%



% Set the number of simulations
N_sims = 3;

sim_idx = randsample(64,N_sims);
init_cond_idx = ff2n(6);
init_cond_idx(init_cond_idx == 0) = -1;



total_error_fig = figure();
total_u_fig = figure();





DV_all = NaN(1,N_sims);


x_sim = cell(1,N_sims);
t_sim = cell(1,N_sims);
u_sim = cell(1,N_sims);
lvlh_sim = cell(1,N_sims);



for j = 1:N_sims
    
    fprintf('Sim %d\n',j);
    
    idx = sim_idx(j);
    init_cond = init_cond_idx(idx,:);
    
    
    x0_lvlh = x0_base.*init_cond';
    
    x0_eci_0 = RTN_to_ECI(x0_lvlh, ref_eci_0);



    % Add mass terms
    ref_eci = [ref_eci_0; mass];
    x0_eci = [x0_eci_0; mass];
    


    

    %% FB

    
    Dx = .0001*M2KM/S2HR^2;
    Dy = .0001*M2KM/S2HR^2;
    Dz = .0001*M2KM/S2HR^2;
    
    x0_lvlh_scaled = x0_lvlh.*[M2KM; M2KM; M2KM; M2KM/S2HR; M2KM/S2HR; M2KM/S2HR];
    
    
    % Manually change the number of additial basis functions in the
    % FB_HCW_XYZ class (flexStates variable)
    fb_ctrl = FB_HCW_XYZ(n/S2HR,T*S2HR);
    fb_ctrl.init_ctrl(x0_lvlh_scaled, Dx, Dy, Dz);
    fprintf('Control Initiated\n');
    
    x0 = [ref_eci; x0_eci];
    
    [t_all, x_all] = prop.prop_2(TF, dt, x0, fb_ctrl);
    

    n_steps = size(t_all(:),1);
    u_all = NaN(3,n_steps);
    lvlh_error = NaN(6,n_steps);
    
    for k = 1:n_steps
        t = t_all(k);
        x = x_all(k,:)';
        
        x_ref = x_all(k,1:7)';
        x_sat = x_all(k,8:14)';
        
        [r_rel_x, v_rel_x, ~] = rva_relative(x_ref(1:3)',x_ref(4:6)',x_sat(1:3)',x_sat(4:6)',mu);
        lvlh_error(:,k) = [r_rel_x; v_rel_x];        
        
        
        [ ~, u_all(:,k)] = prop.dyn_2(t, x, fb_ctrl);
        
    end
    
    
    x_ref_all = x_all(:,1:7);
    x_sat_all = x_all(:,8:14);
    
    final_error = x_ref_all(1:6,end) - x_sat_all(1:6,end);
    final_mass = x_sat_all(7,end);
    
    
    lvlh_error(:,end)
    
    
    t_sim{j} = t_all;
    x_sim{j} = x_all;
    u_sim{j} = u_all;
    lvlh_sim{j} = lvlh_error;
    DV_all(j) = prop.params.Isp*prop.params.g0*log(mass/final_mass);
    
    

    figure(total_error_fig)
    
    for ii = 1:6
       subplot(6,1,ii)
       hold on
       plot(t_all, lvlh_error(ii,:),'k-')
    end
    
    figure(total_u_fig)
    
    for ii = 1:3
       subplot(3,1,ii)
       hold on
       plot(t_all, u_all(ii,:),'k-')
    end

    
    drawnow;
    
end

save(filename,"t_sim","x_sim","u_sim","lvlh_sim","DV_all")



figure(total_error_fig)

subplot(6,1,1)
title('State Errors in LVLH Frame (m, m/s)')
hold on; grid on
ylabel('X')

subplot(6,1,2)
hold on; grid on
ylabel('Y')

subplot(6,1,3)
hold on; grid on
ylabel('Z')

subplot(6,1,4)
hold on; grid on
ylabel('VX')

subplot(6,1,5)
hold on; grid on
ylabel('VY')

subplot(6,1,6)
hold on; grid on
ylabel('VZ')
xlabel('Time (hr)')






