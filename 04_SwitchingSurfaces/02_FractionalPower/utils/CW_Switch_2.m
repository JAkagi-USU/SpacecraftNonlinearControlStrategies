classdef CW_Switch_2
   
    properties
        ctrl_x;
        ctrl_y;
        ctrl_z;
        
        T_x;
        T_y;
        T_z;
        
    end
    
    methods
        function obj = CW_Switch_2(x_param, y_param, z_param, n, varargin)
            
            switch_func = @sign;
            
            nvarargin = nargin - 4;
            
            for setting = 1:2:nvarargin
                switch varargin{setting}
                    case 'SwitchFunc'
                        switch_func = varargin{setting + 1};
                    
                    otherwise
                        warning('Option not found');
                end
            end
            
            obj.ctrl_x = ctrl_switch_2(x_param.a1, x_param.a2, x_param.b1, x_param.b2, x_param.k, 'SwitchFunc', switch_func);
            obj.ctrl_y = ctrl_switch_2(y_param.a1, y_param.a2, y_param.b1, y_param.b2, y_param.k, 'SwitchFunc', switch_func);
            obj.ctrl_z = ctrl_switch_2(z_param.a1, z_param.a2, z_param.b1, z_param.b2, z_param.k, 'SwitchFunc', switch_func);
            
            obj.T_x = [3*n^2 0 0 2*n];
            obj.T_y = [0 0 -2*n 0];
            obj.T_z = [n^2 0];
            
        end
        
        function u = calc_u(obj, t, x0)
            
            x = x0([1;4]);
            y = x0([2;5]);
            z = x0([3;6]);
            
            vx = obj.ctrl_x.calc_u(x);
            vy = obj.ctrl_y.calc_u(y);
            vz = obj.ctrl_z.calc_u(z);
            
            u_xy = [vx; vy] + [obj.T_x; obj.T_y]*x0([1;2;4;5]);
            u_z = vz + obj.T_z*z;
            
            u = [u_xy; u_z];
            
            
            
            
        end
       
       
    end
        
        
    
    
    
    
    
    
    
end