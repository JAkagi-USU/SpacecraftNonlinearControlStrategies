classdef switch_surf_ctrl < handle
    
    properties
       
        C_xy;
        C_z;
        C_full;
        tx; 
        ty;
        tz;
        vx_all;
        vy_all;
        vz_all;
        options;
        switch_state;
        
    end
    
    methods
       
        function obj = switch_surf_ctrl(C_xy, C_z)
            
            obj.C_xy = C_xy;
            obj.C_z = C_z;
            obj.C_full = [C_xy zeros(2,2); zeros(1,4) C_z];
            obj.options = optimoptions('fmincon','Display','off');
            
        end
        
        function optimize_ctrl(obj, xyz0, v0, lb, ub, N_step)
            
            xy0 = xyz0([1;2;4;5]);
            x0 = xyz0([1;4]);
            y0 = xyz0([2;5]);
            z0 = xyz0([3;6]);
            
            vxy = v0([1;2;3;4]);
            vz = v0([5;6]);
            
            LB = repmat(lb,4,1);
            UB = repmat(ub,4,1);
            
            [v_sol, J_sol, flag] = fmincon(@(v) switch_surf_ctrl.switch_obj(v, xy0, obj.C_xy, N_step), vxy,[],[],[],[],LB, UB, [], obj.options);
            
            LB = repmat(lb,2,1);
            UB = repmat(ub,2,1);
            
            [v_sol_z, J_sol_z, flag_z] = fmincon(@(v) switch_surf_ctrl.switch_obj_z(v, z0, obj.C_z, N_step), vz,[],[],[],[],LB, UB, [], obj.options);
            
            [obj.tx, obj.ty, obj.tz, obj.vx_all, obj.vy_all, obj.vz_all] = switch_planning.calc_switch_3(x0, v_sol(1:2), y0, v_sol(3:4), z0, v_sol_z);
            
            
        end
        
        function set_v(obj, x, v0)
            
            x0 = x([1;4]);
            y0 = x([2;5]);
            z0 = x([3;6]);
            
            vx = v0([1;2]);
            vy = v0([3;4]);
            vz = v0([5;6]);
           
             [obj.tx, obj.ty, obj.tz, obj.vx_all, obj.vy_all, obj.vz_all] = switch_planning.calc_switch_3(x0, vx, y0, vy, z0, vz);
            
        end
        
        function u = calc_u(obj, t, x)
            
            x0 = x([1;4]);
            y0 = x([2;5]);
            z0 = x([3;6]);
            
            if t == 0
               obj.switch_state = ones(3,1); 
               s = -sign(switch_planning.surf_x2(obj.vx_all(2),x0(1)));
               if x0(2) > s
                    obj.switch_state(1) = 1;
                elseif x0(2) <= s
                    obj.switch_state(1) = -1;
               end
                
               s = -sign(switch_planning.surf_x2(obj.vy_all(2),y0(1)));
               if y0(2) > s
                    obj.switch_state(2) = 1;
                elseif y0(2) <= s
                    obj.switch_state(2) = -1;
               end
                
               s = -sign(switch_planning.surf_x2(obj.vz_all(2),z0(1)));
               if z0(2) > s
                    obj.switch_state(3) = 1;
                elseif z0(2) <= s
                    obj.switch_state(3) = -1;
               end
                
            end
            
            
            
            eps = 1e-3;
            
            if abs(obj.switch_state(1)) == 1
                s = switch_planning.surf_x2(obj.vx_all(2),x0(1));
                
                if abs(x0(2) - s) < eps
                    obj.switch_state(1) = 2;
                end

%                 if x0(2) > s
%                     cur_state = 1;
%                 elseif x0(2) <= s
%                     cur_state = -1;
%                 end
%                 
%                 
%                 if cur_state ~= obj.switch_state(1)
%                     obj.switch_state(1) = 2;
%                 end
                
            end
            
            if abs(obj.switch_state(2)) == 1
                s = switch_planning.surf_x2(obj.vy_all(2),y0(1));
                
                if abs(y0(2) - s) < eps
                    obj.switch_state(2) = 2;
                end

%                 if y0(2) > s
%                     cur_state = 1;
%                 elseif y0(2) <= s
%                     cur_state = -1;
%                 end
%                 
%                 if cur_state ~= obj.switch_state(2)
%                     obj.switch_state(2) = 2;
%                 end
            end
            
            if abs(obj.switch_state(3)) == 1
                s = switch_planning.surf_x2(obj.vz_all(2),z0(1));
                
                if abs(z0(2) - s) < eps
                    obj.switch_state(3) = 2;
                end

%                 if z0(2) > s
%                     cur_state = 1;
%                 elseif z0(2) <= s
%                     cur_state = -1;
%                 end
%                 
%                 if cur_state ~= obj.switch_state(3)
%                     obj.switch_state(3) = 2;
%                 end
            end
            
            
            
            
            if abs(obj.switch_state(1)) == 1
                v1 = switch_planning.calc_v1_mag(x0, obj.vx_all(1), obj.vx_all(2));
            else
                v1 = switch_planning.calc_v1_mag(x0, obj.vx_all(2), obj.vx_all(2));
            end
            
            if abs(obj.switch_state(2)) == 1
                v2 = switch_planning.calc_v1_mag(y0, obj.vy_all(1), obj.vy_all(2));
            else
                v2 = switch_planning.calc_v1_mag(y0, obj.vy_all(2), obj.vy_all(2));
            end
            
            if abs(obj.switch_state(3)) == 1
                v3 = switch_planning.calc_v1_mag(z0, obj.vz_all(1), obj.vz_all(2));
            else
                v3 = switch_planning.calc_v1_mag(z0, obj.vz_all(2), obj.vz_all(2));
            end
            
            v_full = [v1; v2; v3];
            x_shift = x([1;2;4;5;3;6]);
            u = v_full + obj.C_full*x_shift;
            
            
            
            
            
        end
        
    end
    
    methods(Static)
    
        function J = switch_obj(v, x, C, N)

x0 = x([1;3]);
y0 = x([2;4]);

v1_all = v([1;2]);
v2_all = v([3;4]);

u_all = switch_planning.calc_u(x0, v1_all, y0, v2_all, C, N);

J = norm(u_all(:),Inf);



end

function J = switch_obj_z(v, x0, C, N)

u_all = switch_planning.calc_u_z(x0, v, C, N);

J = norm(u_all(:),Inf);



end
        
    end
    
    
    
    
    
    
    
    
    
end