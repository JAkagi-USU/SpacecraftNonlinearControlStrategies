classdef FB_HCW < handle
    
    properties
        TD;
        TC;
        
        dy_xy;
        dy_z;
        
        T;
        n;
        
        rho_xy;
        rho_z;
        
        ctrl_xy;
        ctrl_z;
        
        g;
        
        
    end
    
    methods
        
        function obj = FB_HCW(n, T)
            
            obj.T = T;
            obj.n = n;
            
            a = 2*n;
            b = 3*n^2;
            obj.TD = [0 a 0 0; 0 0 a 0; -b 0 1 0; 0 -b 0 1];
            obj.TC = [0, n^2/a, 0, 0];
            
            obj.dy_xy = cell(4+1,1);
            for i = 0:4
                obj.dy_xy{i+1} = @(t) 0;
            end
            
            obj.dy_z = cell(2+1,1);
            for i = 0:2
                obj.dy_z{i+1} = @(t) 0;
            end
            
            obj.g = @(x) 1;
            
            
            
            
            
        end
        
        function init_ctrl(obj, x0, Dx, Dz)
            
            xy = x0([1;4;2;5]);
            z = x0([3;6]);
            
            flexStates = 0;
            basis = @(x,y) legendreP(x,y);
            
            xy_z = obj.TD\xy;
            
            obj.rho_xy = generate_rho(xy_z, @obj.df_xy, obj.dy_xy, obj.T, basis, flexStates).rho;
            obj.rho_z = generate_rho(z, @obj.df_z, obj.dy_z, obj.T, basis, flexStates).rho;
            
            obj.ctrl_xy = FTB_ctrl(@obj.df_xy, obj.g, obj.rho_xy, obj.dy_xy, Dx, obj.T, obj.T);
            obj.ctrl_z = FTB_ctrl(@obj.df_z, obj.g, obj.rho_z, obj.dy_z, Dz, obj.T, obj.T);
            
        end
        
        function df = df_xy(obj, x, state, derivative)
            
            df = 0;
            
        end
        
        function df = df_z(obj, x, state, derivative)
            
            df = 0;
            if state == 2
                if derivative == 0
                    df = -obj.n^2*x(1);
                elseif derivative == 1
                    df = -obj.n^2;
                end
            end
            
        end
        
        function u = calc_u(obj, t, x)
            
            xy = x([1;4;2;5]);
            z = x([3;6]);
            
            xy_z = obj.TD\xy;
            
            v_y = obj.ctrl_xy.calc_u(t,xy_z);
            u_y = v_y + obj.TC*xy;
            
            u_z = obj.ctrl_z.calc_u(t,z);
            
            u = [0;u_y;u_z];
            
            
        end
        
        
    end
    
    
end