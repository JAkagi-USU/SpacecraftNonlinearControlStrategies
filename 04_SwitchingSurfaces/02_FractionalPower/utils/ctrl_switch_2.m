% Eq 17 From Jimenez-Rodriguez 2019

classdef ctrl_switch_2
    
    properties
        a1;
        a2;
        b1;
        b2;
        k;
        
        switch_func = @sign;
        
    end
    
    methods
        
        function obj = ctrl_switch_2(a1, a2, b1, b2, k, varargin)
            
            nvarargin = nargin - 5;
            
            for setting = 1:2:nvarargin
                switch varargin{setting}
                    case 'SwitchFunc'
                        obj.switch_func = varargin{setting + 1};
                    otherwise
                        warning('Option not found');
                end
            end
            
            obj.a1 = a1;
            obj.a2 = a2;
            obj.b1 = b1;
            obj.b2 = b2;
            obj.k = k;
            
            
        end
        
        function s = sigma(obj, x)
            
            x1 = x(1);
            x2 = x(2);
            
%             s = x2 + obj.norm(obj.norm(x2,2) + obj.a1*x1 + obj.b1*obj.norm(x1,3),1/2); % Control law found in original paper but seems to have a typo
            s = x2 + obj.norm(obj.norm(x1,2) + obj.a1*x1 + obj.b1*obj.norm(x1,3),1/2);
            
        end
        
        function u = calc_u(obj, x)
           x1 = x(1);
           s = obj.sigma(x);
           
            u = -(obj.a1 + 3*obj.b1*x1^2 + 2*obj.k)*obj.switch_func(s)/2 - obj.norm(obj.a2*s + obj.b2*obj.norm(s,3),1/2);
            
        end
        
        
        
        function y = norm(obj, x, h)
            
            y = obj.switch_func(x)*abs(x)^h;
            
        end
        
        
        
        
        
    end
    
end



