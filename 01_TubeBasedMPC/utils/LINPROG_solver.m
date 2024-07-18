classdef LINPROG_solver < handle
   
    % L1 Solver
    
    properties
        AD;
        BD;
        N;
        
        A_eq;
        b_eq;
        
        A_ineq;
        b_ineq;
        
        ub;
        lb;
        
        objective;
        opts;
        
        n_states;
        n_ctrl;
        
        
    end
    
    methods
       
        function obj = LINPROG_solver(AD, BD, N, Q, R, UB_X, UB_U)
            
            obj.AD = AD;
            obj.BD = BD;
            obj.N = N;
            
            obj.n_states = size(AD,1);
            obj.n_ctrl = size(BD,2);
            
            A_eye = kron(eye(N),eye(obj.n_states));
            A_AD  = kron(diag(ones(N-1,1),-1),-AD);
            A_BD  = kron(eye(N),-BD);
            A_eps_0 = zeros(N*obj.n_states,N*(obj.n_states+obj.n_ctrl));
            A_final = [zeros(obj.n_states,(N-1)*obj.n_states) eye(obj.n_states) zeros(obj.n_states,N*obj.n_ctrl) zeros(obj.n_states,N*(obj.n_states+obj.n_ctrl))];
            
            obj.A_eq = [A_eye+A_AD, A_BD, A_eps_0; A_final];
            obj.b_eq = [zeros(obj.n_states*N,1); zeros(obj.n_states,1)];
            
            A_Q = kron(eye(N),Q);
            A_Q(end-(obj.n_states-1):end,end-(obj.n_states-1):end) = 0; % Zero weight for terminal error
            A_R = kron(eye(N),R);
            
            A_eps = [A_Q zeros(N*obj.n_states,N*obj.n_ctrl); zeros(N*obj.n_ctrl, N*obj.n_states) A_R];
            
            obj.A_ineq = [A_eps -eye(N*(obj.n_states+obj.n_ctrl)); -A_eps -eye(N*(obj.n_states+obj.n_ctrl))];
            obj.b_ineq = zeros(2*N*(obj.n_states+obj.n_ctrl),1);
            
            obj.ub = [repmat(UB_X,N,1); repmat(UB_U,N,1); inf(N*(obj.n_states+obj.n_ctrl),1)];
            obj.lb = [repmat(-UB_X,N,1); repmat(-UB_U,N,1); zeros(N*(obj.n_states+obj.n_ctrl),1)];
            
            
            obj.objective = [zeros(N*(obj.n_states+obj.n_ctrl),1); ones(N*(obj.n_states+obj.n_ctrl),1)];
            obj.opts = optimoptions('linprog','Display','none');
             
 
            
        end
        
        function [x_sol, u_sol, status] = solve(obj, x0)
            

            obj.b_eq(1:obj.n_states) = obj.AD*x0;
            
            [sol, J, flag] = linprog(obj.objective,obj.A_ineq,obj.b_ineq,obj.A_eq,obj.b_eq,obj.lb,obj.ub, obj.opts);
            
            if flag == 1
                status = 1;
                x_sol = reshape(sol(1:obj.N*obj.n_states),obj.n_states,obj.N);
                u_sol = reshape(sol(obj.N*obj.n_states+1:obj.N*(obj.n_states+obj.n_ctrl)),obj.n_ctrl,obj.N);
            else
                status = -1;
                x_sol = zeros(obj.n_states,obj.N);
                u_sol = zeros(obj.n_ctrl,obj.N);
            end
                        
            
        end
        
    end
    
    
    
    
    
end