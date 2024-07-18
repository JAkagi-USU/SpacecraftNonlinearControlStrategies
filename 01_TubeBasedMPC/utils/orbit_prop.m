classdef orbit_prop < handle
    
    properties
        mu;
        params; 
        
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
        
    end
    
    methods
       
        function obj = orbit_prop()
            
            params.mu = 3.9860044188e14;
            params.CD = 2.2;
            params.A = .048;
            params.J2 = 0.1082635854e-2;
            params.omega_E =  [0;0;7.2921159e-5];
            params.Isp = 160;
            params.g0 = 9.81;
            params.R_E = 6.3781363000e6;
            
            
            obj.params = params;
            
        end
        
        function rho = atm_density(obj, h)
            
            rho = interp1(obj.alt_density(:,1), obj.alt_density(:,2), h);
            if isnan(rho)
                rho = 0;
            end
        end
        
        function x_new = prop(obj, dt, x0, u)
            
            h = norm(x0(1:3)) - obj.params.R_E;
            obj.params.rho = obj.atm_density(h);
            
            [~, x_all] = ode45(@(t,x) orbit_prop.dyn(t, x, u, obj.params), [0 dt], x0);
            
            x_new = x_all(end,:)';
            
            
        end
        
        
    end
    
    methods(Static)
        
        
       function x_dot = dyn(t, x, u_eci, params)
           
           r = x(1:3);
           v = x(4:6);
           m = x(7);
           
           r_norm = norm(r);
           
           v_rel = v - cross(params.omega_E,r);
           a_d = -.5*params.rho*norm(v_rel)*params.CD*params.A*v_rel/m;
           
           zr_5 = (5*r(3)^2)/(r_norm^2);
           a_J2_1 = 3*params.J2 * params.mu * params.R_E^2/(2*r_norm^5);
           a_J2_2 = [r(1)*(zr_5 - 1);
                     r(2)*(zr_5 - 1);
                     r(3)*(zr_5 - 3)];
            a_J2 = a_J2_1*a_J2_2;
           
            
            
           
           T = u_eci*m;
           
           r_dot = v;
           v_dot = -params.mu*r/r_norm^3 + u_eci + a_J2 + a_d;
           m_dot = -norm(T,1)/(params.Isp*params.g0);
           x_dot = [r_dot; v_dot; m_dot];
           
           
            
        end 
    end
    
    
    
    
end