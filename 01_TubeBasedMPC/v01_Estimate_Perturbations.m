%% Parameters

addpath utils

mu = 3.9860044188e14;
mass = 24;
dt = 10; % 10 second step

AD_rr = @(n,dt) [4-3*cos(n*dt), 0, 0; 6*(sin(n*dt)-n*dt), 1, 0; 0, 0, cos(n*dt)];
AD_vr = @(n,dt) [sin(n*dt)/n, 2*(1-cos(n*dt))/n, 0; 2*(cos(n*dt)-1)/n,(4*sin(n*dt)-3*n*dt)/n,0;0,0,sin(n*dt)/n];
AD_rv = @(n,dt) [3*n*sin(n*dt), 0, 0; 6*n*(cos(n*dt)-1),0,0; 0,0,-n*sin(n*dt)];
AD_vv = @(n,dt) [cos(n*dt), 2*sin(n*dt),0; -2*sin(n*dt),4*cos(n*dt)-3,0;0,0,cos(n*dt)];

BD_f = @(n,dt) [2*(sin(n*dt/2)/n)^2, -2*(sin(n*dt) - n*dt)/(n^2), 0;
    2*(sin(n*dt) - n*dt)/(n^2), -(8*cos(n*dt) + 3*n^2*dt^2-8)/(2*n^2),0;
    0,0,2*(sin(n*dt/2)/n)^2;
    sin(n*dt)/n, -(2*cos(n*dt)-2)/n,0;
    -4*sin(n*dt/2)^2/n, 4*sin(n*dt)/n-3*dt,0;
    0,0,sin(n*dt)/n];


%% 2 Body Dynamics

prop = orbit_prop();

%% Orbit Parameters
min_el = [6728000; 0; deg2rad(0); deg2rad(0); deg2rad(0); deg2rad(0)];
max_el = [7378000; .01; deg2rad(90); deg2rad(360); deg2rad(360); deg2rad(360)];

xd_max_rel = [1000; 1000; 1000; 50; 50; 50];
xd_min_rel = -xd_max_rel;

x0_max_rel = [1000; 1000; 1000; 10; 10; 10];
x0_min_rel = -x0_max_rel;

max_u = [.25; .25; .25];
min_u = -max_u;



N_samples = 50000;
error_all = NaN(6,N_samples);

for i = 1:N_samples
    
    if mod(i,1000) == 0
        fprintf('%d/%d Complete\n',i, N_samples)
    end


%% Generate Sample Points
ref_el = min_el + (max_el - min_el).*rand(6,1);
[r, v] = orbel2rv(ref_el(1),ref_el(2),ref_el(3),ref_el(4),ref_el(5),ref_el(6),mu);
ref_eci_rv = [r; v];

xd_rel = xd_min_rel + (xd_max_rel - xd_min_rel).*rand(6,1);
x0_rel = xd_rel + x0_min_rel + (x0_max_rel - x0_min_rel).*rand(6,1);

u_rel = min_u + (max_u - min_u).*rand(3,1);

xd_eci_rv = RTN_to_ECI(xd_rel, ref_eci_rv); % Function is correct
x0_eci_rv = RTN_to_ECI(x0_rel, ref_eci_rv); % Function is correct

ref_eci = [ref_eci_rv; mass];
x0_eci = [x0_eci_rv; mass];
xd_eci = [xd_eci_rv; mass];



%% HCW Dynamics

n = sqrt(mu/ref_el(1)^3);

AD = [AD_rr(n,dt), AD_vr(n,dt);
    AD_rv(n,dt), AD_vv(n,dt)];
BD = BD_f(n, dt);

%% Propagate

T = ECI_2_LVLH_T(ref_eci(1:6));
u_eci = T'*u_rel;

ref_eci = prop.prop(dt, ref_eci, zeros(3,1));
x0_eci = prop.prop(dt, x0_eci, u_eci);
xd_eci = prop.prop(dt, xd_eci, zeros(3,1));

[r_rel_x, v_rel_x, ~] = rva_relative(ref_eci(1:3)',ref_eci(4:6)',x0_eci(1:3)',x0_eci(4:6)',mu);
x0_rel_ECI = [r_rel_x; v_rel_x];

[r_rel_x, v_rel_x, ~] = rva_relative(ref_eci(1:3)',ref_eci(4:6)',xd_eci(1:3)',xd_eci(4:6)',mu);
xd_rel_ECI = [r_rel_x; v_rel_x];


x0_rel_LVLH = AD*x0_rel + BD*u_rel;
xd_rel_LVLH = AD*xd_rel;

ECI_diff = x0_rel_ECI - xd_rel_ECI;
LVLH_diff = x0_rel_LVLH - xd_rel_LVLH;

error_all(:,i) = ECI_diff - LVLH_diff;



end

max_pert = max(abs(error_all),[],2)

% 50000 Samples
% max_pert =
% 
%    0.071128799696453
%    0.072760138893273
%    0.017424039137723
%    0.019086639797326
%    0.017631978292535
%    0.000514871031797

