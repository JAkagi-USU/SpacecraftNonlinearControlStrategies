classdef switch_planning
    
    methods(Static)
        
        function [t_all, x_all, v_all, ts, tf] = prop_states(x0, v1, v2, N_samples)
            
            v1 = switch_planning.calc_v1_mag(x0, v1, v2);
            [ts, v2] = switch_planning.calc_surf_time(x0, v1, v2);
            tf = switch_planning.calc_tf(x0, v1, v2, ts);
            
            N_1 = round(ts*N_samples/(tf+ts));
            N_2 = N_samples - N_1;
            
            if N_1 == N_samples
                N_1 = N_samples - 2;
                N_2 = 2;
            elseif N_2 == N_samples
                N_1 = 2;
                N_2 = N_samples - 2;
            end
            
            t_all_s = linspace(0,ts,N_1+1);
            t_all_f = linspace(0,tf,N_2);
            
            t_all = [t_all_s(1:end-1) t_all_f+ts];
            v1_all = repmat(v1,1,numel(t_all_s(1:end-1)));
            v2_all = repmat(v2,1,numel(t_all_f));
            v_all = [v1_all v2_all];
            
            x_1_all = switch_planning.calc_x(t_all_s, x0, v1);
            x_2_all = switch_planning.calc_x(t_all_f, x_1_all(:,end), v2);
            
            x_all = [x_1_all(:,1:end-1) x_2_all];
            
            
        end
        
        function tf = calc_tf(x0, v1, v2, ts)
            
            xs = switch_planning.calc_x(ts, x0, v1);
            
            if xs(2) == 0
                tf = 0;
            else
                tf = -xs(2)/v2;
            end
            
        end
        
        function x_all = calc_x(t, x0, v)
           
            x1 = x0(1) + x0(2).*t + .5*v.*t.^2;
            x2 = x0(2) + v.*t;
            
            x_all = [reshape(x1,1,numel(x1));
                     reshape(x2,1,numel(x2))];
            
        end
        
        function v = calc_v1_mag(x0, v1, v2)
            
            s = switch_planning.surf_x2(v2,x0(1));
            
            if abs(x0(2) - s) < 1e-6
                v = -sign(x0(2))*abs(v1);
            elseif x0(2) > s
                v = -abs(v1);
            elseif x0(2) < s
                v = abs(v1);
            end
            
        end
        
        function [t, v2] = calc_surf_time(x, v1, v2)
            
            if v1 == 0
                t = 0;
                v2 = -sign(x(2))*abs(v2);
                return
            end

            a = -2*x(2)*(v1-v2);
            b = sqrt((4*v2*x(2)^2 - 8*v1*v2*x(1))*(v2 - v1));
            c = 2*(v1^2 - v1*v2);
            
            t_opt = abs([(a + b)/c; (a - b)/c]);
            
            
            z1 = x(1) + x(2)*t_opt(1) + .5*v1*t_opt(1)^2;
            z2 = x(1) + x(2)*t_opt(2) + .5*v1*t_opt(2)^2;
            
            x2_t1 = switch_planning.surf_x2(v2,z1);
            x2_t2 = switch_planning.surf_x2(v2,z2);
            
            if abs(x(2) + v1*t_opt(1) - x2_t1) < 1e-6
                t = t_opt(1);
            elseif abs(x(2) + v1*t_opt(2) - x2_t2) < 1e-6
                t = t_opt(2);
            else
                v2 = -v2;
                
                a = -2*x(2)*(v1-v2);
                b = sqrt((4*v2*x(2)^2 - 8*v1*v2*x(1))*(v2 - v1));
                c = 2*(v1^2 - v1*v2);
                
                t_opt = abs([(a + b)/c; (a - b)/c]);
                
                
                z3 = x(1) + x(2)*t_opt(1) + .5*v1*t_opt(1)^2;
                z4 = x(1) + x(2)*t_opt(2) + .5*v1*t_opt(2)^2;
                
                x2_t3 = switch_planning.surf_x2(v2,z3);
                x2_t4 = switch_planning.surf_x2(v2,z4);
                
                if abs(x(2) + v1*t_opt(1) - x2_t3) < 1e-6
                    t = t_opt(1);
                elseif abs(x(2) + v1*t_opt(2) - x2_t4) < 1e-6
                    t = t_opt(2);
                else
                    t = NaN;
                end
                
            end
            
            
            
        end
        
        function x2 = surf_x2(v,x1)
            x2 = -sqrt(2*abs(v*x1)).*sign(x1);
        end
        
        
        function [t_all, x_all, u_all] = calc_full_u(x0, v1_all, y0, v2_all, C, N_step)
            
            t1_switch = [NaN;NaN];
            v1_all(1) = switch_planning.calc_v1_mag(x0, v1_all(1), v1_all(2));
            [t1_switch(1), v1_all(2)] = switch_planning.calc_surf_time(x0, v1_all(1), v1_all(2));
            t1_switch(2) = switch_planning.calc_tf(x0, v1_all(1), v1_all(2), t1_switch(1)) + t1_switch(1);
            
            t2_switch = [NaN;NaN];
            v2_all(1) = switch_planning.calc_v1_mag(y0, v2_all(1), v2_all(2));
            [t2_switch(1), v2_all(2)] = switch_planning.calc_surf_time(y0, v2_all(1), v2_all(2));
            t2_switch(2) = switch_planning.calc_tf(y0, v2_all(1), v2_all(2), t2_switch(1)) + t2_switch(1);
            

            t_switch = sort(unique([t1_switch; t2_switch]));
            n_switch = numel(t_switch);
            
            t0 = 0;
            x_step_1 = x0;
            x_step_2 = y0;
            t_all = [];
            x_all = [];
            u_all = [];
            
            for i = 1:n_switch
                
                
               
                tf = t_switch(i);
                
                if tf <= t1_switch(1)
                    v1 = v1_all(1);
                elseif tf <= t1_switch(2)
                    v1 = v1_all(2);
                elseif tf > t1_switch(2)
                    v1 = 0;
                end
                
                if tf <= t2_switch(1)
                    v2 = v2_all(1);
                elseif tf <= t2_switch(2)
                    v2 = v2_all(2);
                elseif tf > t2_switch(2)
                    v2 = 0;
                end
                
                dt_all = linspace(0, tf - t0, N_step);
                
                if all(dt_all == 0)
                    continue
                end
                
                t_all = [t_all dt_all+t0];
                t0 = tf;
                
                
                
                x1_all = switch_planning.calc_x(dt_all, x_step_1, v1);
                x2_all = switch_planning.calc_x(dt_all, x_step_2, v2);
                
                x_all_step = [x1_all(1,:); x2_all(1,:); x1_all(2,:); x2_all(2,:)];

                u_all = [u_all [v1; v2] + C*x_all_step];
                x_all = [x_all x_all_step];
                
                x_step_1 = x_all([1;3],end);
                x_step_2 = x_all([2;4],end);
                
                
                
            end
        end
        
        function u_all = calc_u_z(z0, v_all, C, N_step)
            
            t1_switch = [NaN;NaN];
            v_all(1) = switch_planning.calc_v1_mag(z0, v_all(1), v_all(2));
            [t1_switch(1), v_all(2)] = switch_planning.calc_surf_time(z0, v_all(1), v_all(2));
            t1_switch(2) = switch_planning.calc_tf(z0, v_all(1), v_all(2), t1_switch(1)) + t1_switch(1);
            
            % First segment
            dt_all_1 = linspace(0, t1_switch(1), N_step);
            x1_all = switch_planning.calc_x(dt_all_1, z0, v_all(1));
            
            % Second Segment
            dt_all_2 = linspace(t1_switch(1), t1_switch(2), N_step);
            x2_all = switch_planning.calc_x(dt_all_2 - t1_switch(1), x1_all(:,end), v_all(2));
            
            u_all = [v_all(1)+C*x1_all, v_all(2)+C*x2_all];
            
        end
        
        function [t_all, x_all, u_all] = calc_all_z(z0, v_all, C, N_step)
            
            t1_switch = [NaN;NaN];
            v_all(1) = switch_planning.calc_v1_mag(x0, v_all(1), v_all(2));
            [t1_switch(1), v_all(2)] = switch_planning.calc_surf_time(x0, v_all(1), v_all(2));
            t1_switch(2) = switch_planning.calc_tf(x0, v_all(1), v_all(2), t1_switch(1)) + t1_switch(1);
            
            % First segment
            dt_all_1 = linspace(0, t1_switch(1), N_step);
            x1_all = switch_planning.calc_x(dt_all_1, z0, v_all(1));
            
            % Second Segment
            dt_all_2 = linspace(t1_switch(1), t1_switch(2), N_step);
            x2_all = switch_planning.calc_x(dt_all_2 - t1_switch(1), x1_all(:,end), v_all(2));
            
            u_all = [v_all(1)+C*x1_all, v_all(2)+C*x2_all];
            x_all = [x1_all x2_all];
            t_all = [dt_all_1 dt_all_2];
            
        end
        
        
        function [tx, ty, v1_all, v2_all] = calc_switch(x0, v1_all, y0, v2_all)
            
             tx = [NaN;NaN];
            v1_all(1) = switch_planning.calc_v1_mag(x0, v1_all(1), v1_all(2));
            [tx(1), v1_all(2)] = switch_planning.calc_surf_time(x0, v1_all(1), v1_all(2));
            tx(2) = switch_planning.calc_tf(x0, v1_all(1), v1_all(2), tx(1)) + tx(1);
            
            ty = [NaN;NaN];
            v2_all(1) = switch_planning.calc_v1_mag(y0, v2_all(1), v2_all(2));
            [ty(1), v2_all(2)] = switch_planning.calc_surf_time(y0, v2_all(1), v2_all(2));
            ty(2) = switch_planning.calc_tf(y0, v2_all(1), v2_all(2), ty(1)) + ty(1);
            
        end
        
        function [tx, ty, tz, v1_all, v2_all, v3_all] = calc_switch_3(x0, v1_all, y0, v2_all, z0, v3_all)
            
             tx = [NaN;NaN];
            v1_all(1) = switch_planning.calc_v1_mag(x0, v1_all(1), v1_all(2));
            [tx(1), v1_all(2)] = switch_planning.calc_surf_time(x0, v1_all(1), v1_all(2));
            tx(2) = switch_planning.calc_tf(x0, v1_all(1), v1_all(2), tx(1)) + tx(1);
            
            ty = [NaN;NaN];
            v2_all(1) = switch_planning.calc_v1_mag(y0, v2_all(1), v2_all(2));
            [ty(1), v2_all(2)] = switch_planning.calc_surf_time(y0, v2_all(1), v2_all(2));
            ty(2) = switch_planning.calc_tf(y0, v2_all(1), v2_all(2), ty(1)) + ty(1);
            
            tz = [NaN;NaN];
            v3_all(1) = switch_planning.calc_v1_mag(z0, v3_all(1), v3_all(2));
            [tz(1), v3_all(2)] = switch_planning.calc_surf_time(z0, v3_all(1), v3_all(2));
            tz(2) = switch_planning.calc_tf(z0, v3_all(1), v3_all(2), tz(1)) + tz(1);
            
        end
            
            function [u_all] = calc_u(x0, v1_all, y0, v2_all, C, N_step)
            
            t1_switch = [NaN;NaN];
            v1_all(1) = switch_planning.calc_v1_mag(x0, v1_all(1), v1_all(2));
            [t1_switch(1), v1_all(2)] = switch_planning.calc_surf_time(x0, v1_all(1), v1_all(2));
            t1_switch(2) = switch_planning.calc_tf(x0, v1_all(1), v1_all(2), t1_switch(1)) + t1_switch(1);
            
            t2_switch = [NaN;NaN];
            v2_all(1) = switch_planning.calc_v1_mag(y0, v2_all(1), v2_all(2));
            [t2_switch(1), v2_all(2)] = switch_planning.calc_surf_time(y0, v2_all(1), v2_all(2));
            t2_switch(2) = switch_planning.calc_tf(y0, v2_all(1), v2_all(2), t2_switch(1)) + t2_switch(1);
            
            
            
            
            
            
            
            t_switch = sort(unique([t1_switch; t2_switch]));
            n_switch = numel(t_switch);
            
            t0 = 0;
            x_step_1 = x0;
            x_step_2 = y0;
            u_all = [];
            
            for i = 1:n_switch
                
                
               
                tf = t_switch(i);
                
                if tf <= t1_switch(1)
                    v1 = v1_all(1);
                elseif tf <= t1_switch(2)
                    v1 = v1_all(2);
                elseif tf > t1_switch(2)
                    v1 = 0;
                end
                
                if tf <= t2_switch(1)
                    v2 = v2_all(1);
                elseif tf <= t2_switch(2)
                    v2 = v2_all(2);
                elseif tf > t2_switch(2)
                    v2 = 0;
                end
                
                dt_all = linspace(0, tf - t0, N_step);
                
                if all(dt_all == 0)
                    continue
                end
                
                t0 = tf;
                
                
                
                x1_all = switch_planning.calc_x(dt_all, x_step_1, v1);
                x2_all = switch_planning.calc_x(dt_all, x_step_2, v2);
                
                x_all_step = [x1_all(1,:); x2_all(1,:); x1_all(2,:); x2_all(2,:)];

                u_all = [u_all [v1; v2] + C*x_all_step];
                
                x_step_1 = x_all_step([1;3],end);
                x_step_2 = x_all_step([2;4],end);
                
                
                
            end
            
            
            
            
        end
        
        
        
    end
    
    
    
    
end