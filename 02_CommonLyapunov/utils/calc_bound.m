classdef calc_bound < handle
    
    properties
       options = optimoptions('fmincon','Display','off'); 
    end

    methods
       
        function obj = calc_bound()
            
        end
        
        function [sol_upper, sol_lower] = solve(obj, x0, Q)
            
            sol_upper = fmincon(@(x) x, x0, [],[],[],[],0,[],@(x)obj.con_upper(x,Q), obj.options);
            sol_lower = fmincon(@(x) -x, x0, [],[],[],[],0,[],@(x)obj.con_lower(x,Q), obj.options);
            
        end
        
        function [sol_upper, sol_lower] = solve_eig(obj, Q)
           
            eigenvalues = eig(Q);
            sol_upper = eigenvalues(end);
            sol_lower = eigenvalues(1);
            
        end
        
    end



methods(Static)
    
   
    function [c, ceq] = con_lower(x, Q)

        eig_val = -eig(Q-x*eye(size(Q)));
        ceq = 0;
        c = eig_val;

    end
    
    function [c, ceq] = con_upper(x, Q)

        eig_val = -eig(-Q+x*eye(size(Q)));
        ceq = 0;
        c = eig_val;

    end
    
end



end