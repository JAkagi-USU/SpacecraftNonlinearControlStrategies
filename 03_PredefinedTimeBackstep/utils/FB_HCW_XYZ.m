classdef FB_HCW_XYZ < handle
    
    properties
        T_x;
        T_y;
        
        dy_x;
        dy_y;
        dy_z;
        
        T;
        n;
        
        rho_x;
        rho_y;
        rho_z;
        
        ctrl_x;
        ctrl_y;
        ctrl_z;
        
        g;
        
        
    end
    
    methods
        
        function obj = FB_HCW_XYZ(n, T)
            
            obj.T = T;
            obj.n = n;
            
            obj.T_x = [3*n^2 0 0 2*n];
            obj.T_y = [0 0 -2*n 0];
            
            obj.dy_x = cell(2+1,1);
            for i = 0:2
                obj.dy_x{i+1} = @(t) 0;
            end
            
            obj.dy_y = cell(2+1,1);
            for i = 0:2
                obj.dy_y{i+1} = @(t) 0;
            end
            
            obj.dy_z = cell(2+1,1);
            for i = 0:2
                obj.dy_z{i+1} = @(t) 0;
            end
            
            obj.g = @(x) 1;
            
            
            
            
            
        end
        
        function init_ctrl(obj, x0, Dx, Dy, Dz)
            
            x = x0([1;4]);
            y = x0([2;5]);
            z = x0([3;6]);
            
            flexStates = 3;
            basis = @(x,y) legendreP(x,y);
            
            
            obj.rho_x = generate_rho(x, @obj.df_x, obj.dy_x, obj.T, basis, flexStates).rho;
            obj.rho_y = generate_rho(y, @obj.df_y, obj.dy_x, obj.T, basis, flexStates).rho;
            obj.rho_z = generate_rho(z, @obj.df_z, obj.dy_z, obj.T, basis, flexStates).rho;
            
            obj.ctrl_x = FTB_ctrl(@obj.df_x, obj.g, obj.rho_x, obj.dy_x, Dx, obj.T, obj.T);
            obj.ctrl_y = FTB_ctrl(@obj.df_y, obj.g, obj.rho_y, obj.dy_y, Dy, obj.T, obj.T);
            obj.ctrl_z = FTB_ctrl(@obj.df_z, obj.g, obj.rho_z, obj.dy_z, Dz, obj.T, obj.T);
            
        end
        
        function df = df_x(obj, x, state, derivative)
            
            df = 0;
            
        end
        
        function df = df_y(obj, x, state, derivative)
            
            df = 0;
            
        end
        
        function df = df_z(obj, x, state, derivative)
            
            df = 0;
            if state == 2
                if derivative == 0
                    df = -obj.n^2*x(1);
                elseif derivative == 1
                    df = -obj.n^2*x(2); % Double check this is correct although I don't believe it is ever actually used
                end
            end
            
        end
        
        function u = calc_u(obj, t, x0)
            
            x = x0([1;4]);
            y = x0([2;5]);
            z = x0([3;6]);
            
            xy = x0([1;2;4;5]);
            
            v_x = obj.ctrl_x.calc_u(t,x);
            v_y = obj.ctrl_y.calc_u(t,y);
            u_z = obj.ctrl_z.calc_u(t,z);
            
            u_x = v_x - obj.T_x*xy;
            u_y = v_y - obj.T_y*xy;
            
            
            
            u = [u_x;u_y;u_z];
            
            
        end
        
        
    end
    
    
end