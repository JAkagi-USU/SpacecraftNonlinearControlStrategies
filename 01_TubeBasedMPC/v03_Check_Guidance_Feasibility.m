clear all

addpath utils

x0_bounds = [350; 350; 350; 1; 1; 1];
x_max = [425; 425; 425; 5; 5; 5];
u_max = [.01; .01; .01];

combo = ff2n(6);
combo(combo == 0) = -1;

n = .0011;
dt = 10;

AD_rr = @(n,dt) [4-3*cos(n*dt), 0, 0; 6*(sin(n*dt)-n*dt), 1, 0; 0, 0, cos(n*dt)];
AD_vr = @(n,dt) [sin(n*dt)/n, 2*(1-cos(n*dt))/n, 0; 2*(cos(n*dt)-1)/n,(4*sin(n*dt)-3*n*dt)/n,0;0,0,sin(n*dt)/n];
AD_rv = @(n,dt) [3*n*sin(n*dt), 0, 0; 6*n*(cos(n*dt)-1),0,0; 0,0,-n*sin(n*dt)];
AD_vv = @(n,dt) [cos(n*dt), 2*sin(n*dt),0; -2*sin(n*dt),4*cos(n*dt)-3,0;0,0,cos(n*dt)];

AD = [AD_rr(n,dt), AD_vr(n,dt);
      AD_rv(n,dt), AD_vv(n,dt)];
  
BD_f = @(n,dt) [2*(sin(n*dt/2)/n)^2, -2*(sin(n*dt) - n*dt)/(n^2), 0;
        2*(sin(n*dt) - n*dt)/(n^2), -(8*cos(n*dt) + 3*n^2*dt^2-8)/(2*n^2),0;
        0,0,2*(sin(n*dt/2)/n)^2;
        sin(n*dt)/n, -(2*cos(n*dt)-2)/n,0;
        -4*sin(n*dt/2)^2/n, 4*sin(n*dt)/n-3*dt,0;
        0,0,sin(n*dt)/n];
    
BD = BD_f(n, dt);


Q = zeros(6);
R = eye(3);
N = 800;
guidance_lin = LINPROG_solver(AD, BD, N, Q, R, x_max, u_max);

for i = 1:64
    
    x0 = combo(i,:)'.*x0_bounds;
    [x_all, u_all,status] = guidance_lin.solve(x0);
    
    if all(x_all == 0)
        warning('Failure: %d\n',i)
    else
        fprintf('Success: %d\n',i)
    end
    
end