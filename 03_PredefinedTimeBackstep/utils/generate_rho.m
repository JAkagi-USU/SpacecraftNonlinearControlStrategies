classdef generate_rho < handle
    
    % Generates tuning function from t = 0 to t = 1
    properties
       rho; 
    end
    
    methods
        function obj = generate_rho(x0, df_dt, yd, T, basis, flexStates)
            
            % Generate initial conditions for rho
            n_states = size(x0,1);
            P0 = NaN(n_states,1);
            
            for j = 0:n_states-1
                P0(j+1) = x0(j+1) - yd{j+1}(0);
                
                sum = 0;
                for k = 1:j
                    sum = sum + df_dt(x0,k,j-k);
                end
                
                P0(j+1) = P0(j+1) + sum;
                
            end
            
            syms x
            
            
            n_der = n_states+1;
            n_coef = 2*n_states+1+flexStates;
            leg_f = sym(zeros(n_der,n_coef));
            
            for i = 1:n_coef
                
                for j = 1:n_der
                    L = basis(i-1,x);
                    for k = 1:j-1
                        L = diff(L,x);
                    end
                    leg_f(j,i) = L;
                    
                end
            end
            
            rho_func = matlabFunction(leg_f);
            
            P_all = [P0; zeros(n_states+1,1)];
            
            P_0 = rho_func(0);
            P_T = rho_func(T);
            
            Q = [P_0(1:n_states,:); P_T(1:n_states+1,:)];
            
            if flexStates == 0
                
                coef = Q\P_all;
                
            else
                
                tmp = min(abs(Q(:,1:2*n_states+1)\P_all));
                
                coef0 = [Q(:,1:2*n_states+1)\P_all; tmp*rand(flexStates,1)];
                
                rho_x = rho_func(x);
                rho_n = rho_x(end,:);
                
                options = optimoptions('fmincon','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-20,'StepTolerance',1e-20);
%                 options = optimoptions('fmincon','Display','off');
                
                [coef, J, flag, output] = fmincon(@(x)generate_rho.opt_rho(x,rho_n,T),coef0,[],[],Q,P_all,[],[],[],options);
                
                flag
                fprintf("J: %.3f\n",J);
                
            end
            
            assert(not(any(isnan(coef))))       % Check for NaNs in coefficients
            assert(not(any(abs(coef) == inf)))  % Check for Inf in coefficients
            
            obj.rho = cell(n_der,1);
            for i = 1:n_der
                obj.rho{i} = matlabFunction(leg_f(i,:)*coef);
            end
            
            
            
        end
        
    end
    
    methods(Static)
        
        function J = opt_rho(coef, rho_n,T)
            
            t_all = linspace(0,T,1000);
            n1_f = matlabFunction(rho_n*coef);
            
            try
                J = max(abs(n1_f(t_all)));
            catch
                J = 1e10;
            end
            
        end
        
    end
    
    
    
end

