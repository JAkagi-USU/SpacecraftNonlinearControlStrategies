classdef orbit_prop_2Sat < handle
    
    properties
        mu;
        params; 
        
        M2KM;
        S2HR;
        
        alt_density = [0 1.225;
                       1000 1.112;
                       2000 1.007;
                       3000 .9093;
                       4000 .8194;
                       5000 .7364;
                       6000 .6601;
                       7000 .5900;
                       8000 .8258;
                       9000 .4671;
                       10000 .4135;
                       15000 .1948;
                       20000 .08891;
                       25000 .04008;
                       30000 .01841;
                       40000 .003996;
                       50000 .001027;
                       60000 .0003097;
                       70000 .00008283;
                       80000 .00001846;
                       90000 3.416e-6;
                       100000 5.604e-7;
                       110000 9.708e-8;
                       120000 2.222e-8;
                       130000 8.152e-9;
                       140000 3.831e-9;
                       150000 2.076e-9;
                       160000 1.233e-9;
                       170000 7.815e-10;
                       180000 5.194e-10;
                       190000 3.581e-10;
                       200000 2.541e-10;
                       250000 6.073e-11;
                       300000 1.916e-11;
                       350000 7.014e-12;
                       400000 2.803e-12;
                       450000 1.814e-12;
                       500000 5.215e-13;
                       550000 2.384e-13;
                       600000 1.137e-13;
                       650000 5.712e-14;
                       700000 3.070e-14;
                       750000 1.788e-14;
                       800000 1.136e-14;
                       850000 7.824e-15;
                       900000 5.759e-15;
                       950000 4.453e-15;
                       1000000 3.561e-15];
                   
                   
           A_f = @(n) [0 0 0 1 0 0;
                        0 0 0 0 1 0;
                        0 0 0 0 0 1;
                        3*n^2 0 0 0 2*n 0;
                        0 0 0 -2*n 0 0;
                        0 0 -n^2 0 0 0];
            B = [0 0 0;
                      0 0 0;
                      0 0 0
                      1 0 0
                      0 1 0;
                      0 0 1];
        
    end
    
    methods
       
        function obj = orbit_prop_2Sat()
            
            params.mu = 3.9860044188e14;
            params.CD = 2.2;
            params.A = .048;
            params.J2 = 0.1082635854e-2;
            params.omega_E =  [0;0;7.2921159e-5];
            params.Isp = 160;
            params.g0 = 9.81;
            params.R_E = 6.3781363000e6;
            
            
            obj.params = params;
            
            obj.M2KM = 1;
            obj.S2HR = 1;
            
        end
        
        function set_scaling(obj, M2KM, S2HR)
            
            obj.M2KM = M2KM;
            obj.S2HR = S2HR;
            
        end
        
        function rho = atm_density(obj, h)
            
            rho = interp1(obj.alt_density(:,1), obj.alt_density(:,2), h);
            if isnan(rho)
                rho = 0;
            end
        end
        
        function [t_all, x_all] = prop_2_piece(obj, t_span, x0, finite_ctrl)
            
            h = norm(x0(1:3)) - obj.params.R_E;
            obj.params.rho = obj.atm_density(h);
            
            n_steps = numel(t_span);
            
            t_all = t_span;
            x_all = NaN(n_steps, numel(x0));
            x_all(1,:) = x0;
            
            for i = 1:n_steps-1
                
                x = x_all(i,:)';
                t = t_all(i);
                tf = t_all(i+1);
                
                x_ref = x(1:7);
                x_sat = x(8:14);
                
                [r_rel_x, v_rel_x, ~] = rva_relative(x_ref(1:3)',x_ref(4:6)',x_sat(1:3)',x_sat(4:6)',obj.params.mu);
                x0_rel = [r_rel_x; v_rel_x];
                
                u_lvlh = finite_ctrl.calc_u(t, x0_rel);
            
                [~, x_step] = rk4(@(t,x) obj.dyn_2_const(t, x, u_lvlh), [t tf], x);
%                 [~, x_step] = euler(@(t,x) obj.dyn_2_const(t, x, u_lvlh), [t tf], x);
                
                x_all(i+1,:) = x_step(end,:);
                
            end
            
            
            
        end
        
        
        
        function [t_all, x_all] = prop_2_rk(obj, t_span, x0, finite_ctrl)
            
            h = norm(x0(1:3)) - obj.params.R_E;
            obj.params.rho = obj.atm_density(h);
            
            
%             [t_all, x_all] = ode45(@(t,x) obj.dyn_2(t, x, finite_ctrl), [t0 max(t_span)], x0);
            [t_all, x_all] = rk4(@(t,x) obj.dyn_2(t, x, finite_ctrl), t_span, x0);
            
            
            
        end
        
        function [t_all, x_all] = prop_2_ode(obj, t_span, x0, finite_ctrl)
            
            h = norm(x0(1:3)) - obj.params.R_E;
            obj.params.rho = obj.atm_density(h);
            
            
            [t_all, x_all] = ode45(@(t,x) obj.dyn_2(t, x, finite_ctrl), [t0 max(t_span)], x0);
            
            
            
        end
        
        function [t_all, x_all] = prop_2_euler(obj, t_span, x0, finite_ctrl)
            
            h = norm(x0(1:3)) - obj.params.R_E;
            obj.params.rho = obj.atm_density(h);
            
            
            [t_all, x_all] = euler(@(t,x) obj.dyn_2(t, x, finite_ctrl), t_span, x0);
            
            
            
        end
        
        
        
        function accel = calc_accel(obj, x0, u)
            
            h = norm(x0(1:3)) - obj.params.R_E;
            obj.params.rho = obj.atm_density(h);
            
            accel = orbit_prop.dyn(0, x0, u, obj.params);
            
        end
        
        function [x_dot, u_lvlh] = dyn_2_const(obj, t, x, u_lvlh)
           
           x_ref = x(1:7);
           x_sat = x(8:14);
           
           
           

            [r_rel_x, v_rel_x, ~] = rva_relative(x_ref(1:3)',x_ref(4:6)',x_sat(1:3)',x_sat(4:6)',obj.params.mu);
            x0_rel = [r_rel_x; v_rel_x];
            

           
           x_ref_dot = obj.dyn_1(t, x_ref, zeros(3,1));
           
           r = x_sat(1:3);
           v = x_sat(4:6);
           m = x_sat(7);
           
           r_norm = norm(r);
           
           v_rel = v - cross(obj.params.omega_E,r);
           a_d = -.5*obj.params.rho*norm(v_rel)*obj.params.CD*obj.params.A*v_rel/m;
           
           zr_5 = (5*r(3)^2)/(r_norm^2);
           a_J2_1 = 3*obj.params.J2 * obj.params.mu * obj.params.R_E^2/(2*r_norm^5);
           a_J2_2 = [r(1)*(zr_5 - 1);
                     r(2)*(zr_5 - 1);
                     r(3)*(zr_5 - 3)];
            a_J2 = a_J2_1*a_J2_2;
           
            

            u_lvlh = u_lvlh.*repmat(obj.S2HR^2/obj.M2KM,3,1);
            
            
            T = ECI_2_LVLH_T(x_ref(1:6));
            u_eci = T'*u_lvlh;
            

           thrust = u_eci*m;
           
           r_dot = v;
           v_dot = -obj.params.mu*r/r_norm^3 + u_eci + a_J2 + a_d;
           m_dot = -norm(thrust,1)/(obj.params.Isp*obj.params.g0);
           x_sat_dot = [r_dot; v_dot; m_dot];
           
           x_dot = [x_ref_dot; x_sat_dot];
           
           
            
       end 

        
       function [x_dot, u_lvlh] = dyn_2(obj, t, x, finite_ctrl)
           
           x_ref = x(1:7);
           x_sat = x(8:14);
           
           
           

            [r_rel_x, v_rel_x, ~] = rva_relative(x_ref(1:3)',x_ref(4:6)',x_sat(1:3)',x_sat(4:6)',obj.params.mu);
            x0_rel = [r_rel_x; v_rel_x];
            

           
           x_ref_dot = obj.dyn_1(t, x_ref, zeros(3,1));
           
           r = x_sat(1:3);
           v = x_sat(4:6);
           m = x_sat(7);
           
           r_norm = norm(r);
           
           v_rel = v - cross(obj.params.omega_E,r);
           a_d = -.5*obj.params.rho*norm(v_rel)*obj.params.CD*obj.params.A*v_rel/m;
           
           zr_5 = (5*r(3)^2)/(r_norm^2);
           a_J2_1 = 3*obj.params.J2 * obj.params.mu * obj.params.R_E^2/(2*r_norm^5);
           a_J2_2 = [r(1)*(zr_5 - 1);
                     r(2)*(zr_5 - 1);
                     r(3)*(zr_5 - 3)];
            a_J2 = a_J2_1*a_J2_2;
           
            
            t = t*obj.S2HR;
            x0_rel = x0_rel.*[repmat(obj.M2KM,3,1); repmat(obj.M2KM/obj.S2HR, 3,1)];
            u_lvlh = finite_ctrl.calc_u(t, x0_rel);
            u_lvlh = u_lvlh.*repmat(obj.S2HR^2/obj.M2KM,3,1);
            
            
            T = ECI_2_LVLH_T(x_ref(1:6));
            u_eci = T'*u_lvlh;
            

           thrust = u_eci*m;
           
           r_dot = v;
           v_dot = -obj.params.mu*r/r_norm^3 + u_eci + a_J2 + a_d*100;
           m_dot = -norm(thrust,1)/(obj.params.Isp*obj.params.g0);
           x_sat_dot = [r_dot; v_dot; m_dot];
           
           x_dot = [x_ref_dot; x_sat_dot];
           
           
            
       end 

       function x_dot = dyn_1(obj, t, x, u_eci)
           
           r = x(1:3);
           v = x(4:6);
           m = x(7);
           
           r_norm = norm(r);
           
           v_rel = v - cross(obj.params.omega_E,r);
           a_d = -.5*obj.params.rho*norm(v_rel)*obj.params.CD*obj.params.A*v_rel/m;
           
           zr_5 = (5*r(3)^2)/(r_norm^2);
           a_J2_1 = 3*obj.params.J2 * obj.params.mu * obj.params.R_E^2/(2*r_norm^5);
           a_J2_2 = [r(1)*(zr_5 - 1);
                     r(2)*(zr_5 - 1);
                     r(3)*(zr_5 - 3)];
            a_J2 = a_J2_1*a_J2_2;
           
            
            
           
           T = u_eci*m;
           
           r_dot = v;
           v_dot = -obj.params.mu*r/r_norm^3 + u_eci + a_J2 + a_d;
           m_dot = -norm(T,1)/(obj.params.Isp*obj.params.g0);
           x_dot = [r_dot; v_dot; m_dot];
           
           
            
       end 
       
       function [position,isterminal,direction] = switch_surface(t,x,v)
           
          s_x2 = -sqrt(2*abs(v(1)*x(1))).*sign(x(1));
          s_y2 = -sqrt(2*abs(v(2)*x(2))).*sign(x(2));
          s_z2 = -sqrt(2*abs(v(3)*x(3))).*sign(x(3));
           
          position = [x(4) - s_x2; x(5) - s_y2; x(6) - s_z2]; % The value that we want to be zero
          isterminal = [1; 1; 1];  % Halt integration 
          direction = [0; 0; 0];   % The zero can be approached from either direction
        end
        
       
    end
    
    
    
    
end