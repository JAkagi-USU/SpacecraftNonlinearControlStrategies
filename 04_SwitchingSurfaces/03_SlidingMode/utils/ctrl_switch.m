% Eq 14 From Jimenez-Rodriguez 2019

classdef ctrl_switch
    
    properties
        q;
        rho1;
        rho2;
        alpha;
        
        switch_func = @sign;
        
    end
    
    methods
        
        function obj = ctrl_switch(q, rho1, rho2, alpha, varargin)
            
            nvarargin = nargin - 4;
            
            for setting = 1:2:nvarargin
                switch varargin{setting}
                    case 'SwitchFunc'
                        obj.switch_func = varargin{setting + 1};
                    otherwise
                        warning('Option not found');
                end
            end
            
            obj.q = q;
            obj.rho1 = rho1;
            obj.rho2 = rho2;
            obj.alpha = alpha;
            
            
        end
        
        function [y, yd, ydd] = w(obj, r)
            
            y = (r^obj.q)/(r^obj.q + obj.alpha);
            yd = (obj.alpha*obj.q*r^(obj.q-1))/((obj.alpha + r^obj.q)^2);
            ydd = -((-obj.alpha*obj.q + obj.alpha + (obj.q + 1)*r^obj.q)*obj.alpha*obj.q*r^(obj.q-2))/((obj.alpha + r^obj.q)^3);
            
        end
        
        function [y, yd] = v(obj, x)
            [~, wd, wdd] = obj.w(abs(x));
            y = -obj.switch_func(x)/(obj.rho1*wd);
            
            yd = wdd/(obj.rho1*wd^2);
            
        end
        
        function s = sd(obj, x)
            
            [v1, ~] = obj.v(x(1));
            s = x(2) + obj.sign_norm(obj.sign_norm(x(2),2) - 2*obj.sign_norm(v1,2), 1/2);
        end
        
        function u = calc_u(obj, x)
            
            sd = obj.sd(x);
            [v2, ~] = obj.v(sd);
            [v1, v1_d] = obj.v(x(1));
            u = v2 - obj.rho2*obj.switch_func(sd) - obj.s(x) + 2*abs(v1)*v1_d*obj.switch_func(sd);
        end
        
        function y = sign_norm(obj, x, h)
            
            y = obj.switch_func(x)*abs(x)^h;
            
        end
        
        function y = s(obj, x)
            
            [v1, ~] = obj.v(x(1));
            y = x(2) - v1;
            
            
            
            
        end
        
        
        
        
        
    end
    
end



