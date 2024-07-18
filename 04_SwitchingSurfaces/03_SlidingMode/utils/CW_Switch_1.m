classdef CW_Switch_1
    
    properties
        ctrl_x;
        ctrl_y;
        ctrl_z;
        
        T_x;
        T_y;
        T_z;
        
    end
    
    methods
        function obj = CW_Switch_1(x_param, y_param, z_param, n, varargin)
            
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
            
            obj.ctrl_x = ctrl_switch(x_param.q, x_param.rho_1, x_param.rho_2, x_param.alpha, 'SwitchFunc', switch_func);
            obj.ctrl_y = ctrl_switch(y_param.q, y_param.rho_1, y_param.rho_2, y_param.alpha, 'SwitchFunc', switch_func);
            obj.ctrl_z = ctrl_switch(z_param.q, z_param.rho_1, z_param.rho_2, z_param.alpha, 'SwitchFunc', switch_func);
            
            obj.T_x = [3*n^2 0 0 2*n];
            obj.T_y = [0 0 -2*n 0];
            obj.T_z = [n^2 0];
            
        end
        
        function u = calc_u(obj, t, x0)
            
            %             M2KM = 1/1000;
            %             S2HR = 1/3600;
            
            x0 = x0;% - [0; 150; 0; 0; 0; 0];
            
            x = x0([1;4]);
            y = x0([2;5]);
            z = x0([3;6]);
            
            vx = obj.ctrl_x.calc_u(x);
            vy = obj.ctrl_y.calc_u(y);
            vz = obj.ctrl_z.calc_u(z);
            
            %             vx = obj.ctrl_x.calc_u(x.*[M2KM; M2KM/S2HR])*(S2HR^2/M2KM);
            %             vy = obj.ctrl_y.calc_u(y.*[M2KM; M2KM/S2HR])*(S2HR^2/M2KM);
            %             vz = obj.ctrl_z.calc_u(z.*[M2KM; M2KM/S2HR])*(S2HR^2/M2KM);
            
            u_xy = [vx; vy] + [obj.T_x; obj.T_y]*x0([1;2;4;5]);
            u_z = vz + obj.T_z*z;
            
            u = [u_xy; u_z];
            
            
            
            
        end
        
        
    end
    
    
    
    
    
    
    
    
    
end