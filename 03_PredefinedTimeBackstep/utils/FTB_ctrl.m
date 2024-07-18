% Finite Time Backstepping

classdef FTB_ctrl < handle
    properties
        D;          % Disturbance gain
        yd;          % Desired trajectory (t)
        T;          % Convergence Time
        T_y;        % Final time desired trajectory is valid
                
        error_states; % Error States
        alpha_states; % Alpha States
        
        df_dt;      % Dynamics
        g;          % Control input function
        rho;        % Tuning function

    end
    
    methods
        function obj = FTB_ctrl(df_dt, g, rho, yd, D, T, T_y)
            
            obj.df_dt = df_dt;
            obj.g = g;
            obj.rho = rho;
            obj.yd = yd;
            obj.D = D;
            obj.T = T;
            obj.T_y = T_y;
            
        end
        
        function u = calc_u(obj, t, x)
            
            n_states = size(x,1);
            
            obj.error_states = NaN(1,n_states);
            obj.alpha_states = NaN(1,n_states);
            
            for i = 1:n_states
                obj.alpha_states(i) = obj.calc_alpha(t,x,i);
                obj.error_states(i) = obj.calc_error(t,x,i);
            end
            
            u = (obj.g(x)^-1)*(-obj.df_dt(x,n_states,0) + obj.alpha_dot(t, x, n_states) - obj.D*sign(obj.error_states(n_states)) - obj.error_states(n_states-1));
            
        end

        function alpha = calc_alpha(obj, t, x, j)
            % Eq. 19
            
            if j == 1
                alpha = NaN;
            elseif j == 2
                alpha = -obj.df_dt(x,1,0) + obj.calc_yd(1+1,t) + obj.calc_rho(1+1,t);
            else
                alpha = -obj.df_dt(x,j-1,0) + obj.alpha_dot(t, x, j-1) - obj.error_states(j-2);
            end
            
        end
        
        function val = calc_rho(obj, n_f,  t)
           
            if t >= obj.T
                val = 0;
            else
                val = obj.rho{n_f}(t);
            end
            
            
        end
        
        function val = calc_yd(obj, n_f, t)
            
            if t >= obj.T_y
                val = 0;
            else
                val = obj.yd{n_f}(t);
            end
           
            
            
        end
        
        function alpha_d = alpha_dot(obj, t, x, n)
            
            alpha_d = 0;
            
            j = n - 1;
            % Dynamics portion
            for i = 1:j
                alpha_d = alpha_d - obj.df_dt(x,i,j+1-i);
            end
            
            % Desired Trajectory
            alpha_d = alpha_d + obj.calc_yd(j+1+1,t);
            alpha_d = alpha_d + obj.calc_rho(j+1+1,t);
            
            % Error Derivatives
            for k = 1:j-1
                
                alpha_d = alpha_d - obj.error_dot(k,j-k);
                
            end

        end
        
        function error = calc_error(obj, t, x, j)
            
            if j == 1
                error = x(1) - obj.calc_yd(0+1,t) - obj.calc_rho(0+1,t);
            else
                error = x(j) - obj.alpha_states(j);
            end
            
        end
        
        function e_dot = error_dot(obj, sub, sup)
            
%             assert(sup > 0)
%             assert(sub < obj.n_states)
%             assert(sub > 0)
            
            if sub == 1
                if sup == 1
                    e_dot = obj.error_states(2);
                else
                    e_dot = obj.error_dot(2,sup-1);
                end
            else
                if sup == 1
                    e_dot = obj.error_states(sub+1)-obj.error_states(sub-1);
                else
                    e_dot = obj.error_dot(sub+1,sup-1)...
                        - obj.error_dot(sub-1,sup-1);
                end
            end

            
        end

    end
    

end