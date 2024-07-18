classdef sat_control < handle
   
    properties
        
        AD;
        BD;
        K;
        
        x_traj;
        u_traj;
        mpc_x_pred = NaN(6,1);
        
        guide_law;
        mpc_law;
        

        N_guide;
        N_current = 0;
        
    end
    
    methods
       
        function obj = sat_control(AD, BD, K)
            
            obj.AD = AD;
            obj.BD = BD;
            obj.K = K;
            
        end
        
        function init_guide(obj, N, Q, R, x_max, u_max)
            
            obj.N_guide = N;
            obj.guide_law = LINPROG_solver(obj.AD, obj.BD, N, Q, R, x_max, u_max);  
            
        end
        
        function init_mpc(obj, N, Q, R, x_max, u_max)
            
            obj.mpc_law = LINPROG_solver(obj.AD, obj.BD, N, Q, R, x_max, u_max);
            
        end
        
        function calc_guidance(obj, x0_LVLH, xd_LVLH)
            
            x0_total_error = x0_LVLH - xd_LVLH;
            
            [obj.x_traj, obj.u_traj, status] = obj.guide_law.solve(x0_total_error);
            
            obj.x_traj = [x0_total_error obj.x_traj];
                
            if status == 1
                fprintf('Guidance Success\n')
            else
                warning('Guidance Failure')
            end
            
        end
        
        function [total_error, guide_error] = error(obj, N, x0_LVLH, xd_LVLH)
            
            total_error = x0_LVLH - xd_LVLH;
            
            xd_guide = obj.x_traj(:,N);
            guide_error = total_error - xd_guide;
            
        end
        
        function [u, u_mpc_combine, u_guide] = ctrl(obj, N, x0_LVLH, xd_LVLH)
            
            x0_total_error = x0_LVLH - xd_LVLH;
            
            xd_guide = obj.x_traj(:,N);
            guide_error = x0_total_error - xd_guide;
            
            [x_mpc, u_mpc, status] = obj.mpc_law.solve(guide_error);
            
            if status ~= 1
                warning('MPC Failed')
            end
            
            mpc_feedback_error = guide_error - obj.mpc_x_pred;
            
            if any(isnan(mpc_feedback_error))
                mpc_feedback_error = zeros(6,1);
            end
            
            u_feedback = obj.K*mpc_feedback_error;
            u = obj.u_traj(:,N) + u_feedback + u_mpc(:,1);
            u_mpc_combine = u_feedback + u_mpc(:,1);
            u_guide = obj.u_traj(:,N);
        
            obj.mpc_x_pred = x_mpc(:,1);
                
            
            
        end
        
        
        
    end
    
    
    
end