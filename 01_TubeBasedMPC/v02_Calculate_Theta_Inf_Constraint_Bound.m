%% Initialize Dynamics

addpath utils

n = .0011;
dt = 10;

u_max = [.01; .01; .01];
x_max = [75; 75; 75; 5; 5; 5];

W_max = [.1;.1;.1;.02;.02;.02];



%% Set up dynamics and constraint matrices

n_states = 6;
n_control = 3;

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


W_vertices = ff2n(6);
W_vertices(W_vertices == 0) = -1;
W_vertices = W_vertices.*W_max';


C = [eye(6) zeros(6,3);
     -eye(6) zeros(6,3);
     zeros(3,6) eye(3);
     zeros(3,6) -eye(3)];
 
D = [x_max; x_max; u_max; u_max];
 
cx = C(:,1:n_states);
cu = C(:,n_states+1:end);

%% Calculate Gain K

poles_d = [.85; .85; .95; .95; .9; .9];
K = -place(AD, BD, poles_d);

% Check pole placement

AK = AD + BD*K;
assert(all(abs(eig(AK)) < 1))


%% Calculate specific alpha for N

N = 250;

A_n = AK^N;

max_AW = max(vecnorm(A_n*W_vertices',inf,1));
max_KAW = max(vecnorm(K*A_n*W_vertices',inf,1));
max_w = max(vecnorm(W_vertices',inf,1));
max_Kw = max(vecnorm(K*W_vertices',inf,1));

alpha = max(max_AW/max_w, max_KAW/max_Kw)



%% Calculate theta for given alpha and N



n_rows = size(cx,1);
n_states = size(cx,2);

theta = zeros(n_rows, 1);

options = optimoptions('linprog','Display','none');

for i = 1:n_rows
    A_eq = -eye(n_rows);
    b_eq = zeros(n_rows,1);

    for j = 0:N-1

        A_eq = [A_eq (cx*AK^j+cu*K*AK^j)];

    end

    ub = [inf(n_rows,1); repmat(W_max,N,1)];
    lb = [-inf(n_rows,1); repmat(-W_max,N,1)];

    f = [zeros(n_rows,1); zeros(n_states*N,1)];
    f(i) = 1;

    [sol, J, flag] = linprog(-f, [], [], A_eq, b_eq, lb, ub,[],options);

    theta(i) = sol(i);

end


D_prime = D - (1-alpha)^-1*theta;
D_prime([1;2;3;4;5;6;13;14;15])

