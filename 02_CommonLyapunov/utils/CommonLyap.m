classdef CommonLyap < handle
    
    properties
        
        InPlane;
        CrossTrack;
        
        opts;
        
        
        
    end
    
    methods
        function obj = CommonLyap()
            
            
            obj.opts        = sdpsettings;
            obj.opts.solver = 'sdpt3';
            obj.opts.verbose = 0;
            obj.opts.sdpt3.maxit = 5000;
            obj.opts.sdpt3.steptol = 1e-8;
            
            
        end
        
        function setInPlaneDyn(obj, A0, B0, TC, TD)
            
            obj.InPlane.A0 = A0;
            obj.InPlane.B0 = B0;
            
            obj.InPlane.TC = TC;
            obj.InPlane.TD = TD;
            
            obj.InPlane.AZ = obj.InPlane.TD\(obj.InPlane.A0+obj.InPlane.B0*obj.InPlane.TC)*obj.InPlane.TD;
            obj.InPlane.BZ = obj.InPlane.TD\obj.InPlane.B0;
            
        end
        
        function setCrossTrackDyn(obj, A0, B0, TC, TD)
            
            obj.CrossTrack.A0 = A0;
            obj.CrossTrack.B0 = B0;
            
            obj.CrossTrack.TC = TC;
            obj.CrossTrack.TD = TD;
            
            obj.CrossTrack.AZ = obj.CrossTrack.TD\(obj.CrossTrack.A0+obj.CrossTrack.B0*obj.CrossTrack.TC)*obj.CrossTrack.TD;
            obj.CrossTrack.BZ = obj.CrossTrack.TD\obj.CrossTrack.B0;
            
        end
        
        function setInPlaneGain(obj, K)
            
            obj.InPlane.K = K;
            obj.InPlane.AT = obj.InPlane.AZ + obj.InPlane.BZ*K;
            
        end
        
        function setCrossTrackGain(obj, K)
            
            obj.CrossTrack.K = K;
            obj.CrossTrack.AT = obj.CrossTrack.AZ + obj.CrossTrack.BZ*K;
            
        end
        
        function setInPlaneBounds(obj, P, Q, R)
            
            obj.InPlane.P_bound = P;
            obj.InPlane.Q_bound = Q;
            obj.InPlane.R_bound = R;
            
        end
        
        function setCrossTrackBounds(obj, P, Q, R)
            
            obj.CrossTrack.P_bound = P;
            obj.CrossTrack.Q_bound = Q;
            obj.CrossTrack.R_bound = R;
            
        end
        
        function sol = GeneralSolve(obj, AT, P_bound, Q_bound, R_bound)
            
            n = max(size(AT));
            M = diag(-n:1:-1);
            
            Q_LMI           = sdpvar(n,n);
            R_LMI           = sdpvar(n,n);
            P_LMI           = sdpvar(n,n);
            LMI         = [Q_LMI*AT+AT.'*Q_LMI == -P_LMI, ...
                Q_LMI*M + M*Q_LMI == -R_LMI, ...
                eye(n)*R_bound(2)>=R_LMI>=eye(n)*R_bound(1), ...
                eye(n)*P_bound(2)>=P_LMI>=eye(n)*P_bound(1), ...
                eye(n)*Q_bound(2)>=Q_LMI>=eye(n)*Q_bound(1)];
            sol.out = optimize(LMI,0,obj.opts);
            
            sol.Q = double(Q_LMI);
            sol.P = double(P_LMI);
            sol.R = double(R_LMI);
            
            
            
        end
        
        function sol = GeneralSolve_opt(obj, AT, P_bound, Q_bound, R_bound, param)
            
            n = max(size(AT));
            M = diag(-n:1:-1);
            
            Q_min           = sdpvar(1);
            P_min           = sdpvar(1);
            R_min           = sdpvar(1);
            
            
            switch param.OptVar
                case "Q"
                    opt_J = Q_min;
                case "P"
                    opt_J = P_min;
                case "R"
                    opt_J = R_min;
                otherwise
                    opt_J = 0;
            end
            
            
            Q_LMI           = sdpvar(n,n);
            R_LMI           = sdpvar(n,n);
            P_LMI           = sdpvar(n,n);
            
            LMI         = [Q_LMI*AT+AT.'*Q_LMI == -P_LMI, ...
                Q_LMI*M + M*Q_LMI == -R_LMI, ...
                eye(n)*R_bound(2)>=R_LMI>=eye(n)*R_min, ...
                eye(n)*P_bound(2)>=P_LMI>=eye(n)*P_min, ...
                eye(n)*Q_bound(2)>=Q_LMI>=eye(n)*Q_min, ...
                Q_min >= Q_bound(1), ...
                R_min >= R_bound(1), ...
                P_min >= P_bound(1)];
            sol.out = optimize(LMI,-opt_J,obj.opts);
            
            sol.Q = double(Q_LMI);
            sol.P = double(P_LMI);
            sol.R = double(R_LMI);
            sol.obj = double(opt_J);
            
            
            
        end
        
        function sol = GeneralSolve_target(obj, AT, P_bound, Q_bound, R_bound, param)
            
            n = max(size(AT));
            M = diag(-n:1:-1);
            
            Q_min           = sdpvar(1);
            P_min           = sdpvar(1);
            R_min           = sdpvar(1);
            Q_max           = sdpvar(1);
            P_max           = sdpvar(1);
            R_max           = sdpvar(1);
            
            opt_J = 0;
            
            LMI = [];
            
            Q_weight = 1;
            R_weight = 1;
            P_weight = 1;
            
            fields = fieldnames(param);
            
            for i = 1:numel(fields)
                
               switch fields{i}
                   case "Q_weight"
                       Q_weight = param.Q_weight;
                   case "R_weight"
                       R_weight = param.R_weight;
                   case "P_weight"
                       P_weight = param.P_weight;
               end
                
            end
            
            for i = 1:numel(fields)
               
                switch fields{i}
                    case "Q_target"
                        target = param.Q_target;
                        Qe1 = sdpvar(1);
                        Qe2 = sdpvar(1);
                        LMI = [LMI, ...
                               -Qe1 <= Q_max - target, ...
                               -Qe1 <= -Q_max + target, ...
                               -Qe2 <= Q_min - target, ...
                               -Qe2 <= -Q_min + target, ...
                               Qe1 >= 0, ...
                               Qe2 >= 0];
                        opt_J = opt_J + Q_weight*Qe1 + Q_weight*Qe2;
                        
                    case "P_target"
                        
                        target = param.P_target;
                        Pe1 = sdpvar(1);
                        Pe2 = sdpvar(1);
                        LMI = [LMI, ...
                               -Pe1 <= P_max - target, ...
                               -Pe1 <= -P_max + target, ...
                               -Pe2 <= P_min - target, ...
                               -Pe2 <= -P_min + target, ...
                               Pe1 >= 0, ...
                               Pe2 >= 0];
                        opt_J = opt_J + P_weight*Pe1 + P_weight*Pe2;
                        
                    case "R_target"
                        
                        target = param.R_target;
                        Re1 = sdpvar(1);
                        Re2 = sdpvar(1);
                        LMI = [LMI, ...
                               -Re1 <= R_max - target, ...
                               -Re1 <= -R_max + target, ...
                               -Re2 <= R_min - target, ...
                               -Re2 <= -R_min + target, ...
                               Re1 >= 0, ...
                               Re2 >= 0];
                        opt_J = opt_J + R_weight*Re1 + R_weight*Re2;
                        
                end
                
            end
            
            
            
            
            Q_LMI           = sdpvar(n,n);
            R_LMI           = sdpvar(n,n);
            P_LMI           = sdpvar(n,n);
            
            LMI         = [LMI, Q_LMI*AT+AT.'*Q_LMI == -P_LMI, ...
                Q_LMI*M + M*Q_LMI == -R_LMI, ...
                eye(n)*R_max>=R_LMI>=eye(n)*R_min, ...
                eye(n)*P_max>=P_LMI>=eye(n)*P_min, ...
                eye(n)*Q_max>=Q_LMI>=eye(n)*Q_min, ...
                Q_min >= Q_bound(1), ...
                R_min >= R_bound(1), ...
                P_min >= P_bound(1), ...
                Q_max <= Q_bound(2), ...
                R_max <= R_bound(2), ...
                P_max <= P_bound(2)];
            sol.out = optimize(LMI,opt_J,obj.opts);
            
            sol.Q = double(Q_LMI);
            sol.P = double(P_LMI);
            sol.R = double(R_LMI);
            sol.obj = double(opt_J);
            
            
            
        end
        
        function sol = GeneralSolve_target_ind(obj, AT, P_bound, Q_bound, R_bound, param)
            
            n = max(size(AT));
            M = diag(-n:1:-1);
            
            Q_min           = sdpvar(1);
            P_min           = sdpvar(1);
            R_min           = sdpvar(1);
            Q_max           = sdpvar(1);
            P_max           = sdpvar(1);
            R_max           = sdpvar(1);
            
            opt_J = 0;
            
            LMI = [];
            
            Q_weight_max = 1;
            R_weight_max = 1;
            P_weight_max = 1;
            Q_weight_min = 1;
            R_weight_min = 1;
            P_weight_min = 1;
            
            fields = fieldnames(param);
            
            for i = 1:numel(fields)
                
               switch fields{i}
                   case "Q_weight_max"
                       Q_weight_max = param.Q_weight_max;
                   case "Q_weight_min"
                       Q_weight_min = param.Q_weight_min;
                   case "P_weight_max"
                       P_weight_max = param.P_weight_max;
                   case "P_weight_min"
                       P_weight_min = param.P_weight_min;
                   case "R_weight_max"
                       R_weight_max = param.R_weight_max;
                   case "R_weight_min"
                       R_weight_min = param.R_weight_min;
                   
               end
                
            end
            
            for i = 1:numel(fields)
               
                switch fields{i}
                    case "Q_target_max"
                        Q_target_max = param.Q_target_max;
                        Qe_max = sdpvar(1);
                        LMI = [LMI, ...
                               -Qe_max <= Q_max - Q_target_max, ...
                               -Qe_max <= -Q_max + Q_target_max, ...
                               Qe_max >= 0];
                        opt_J = opt_J + Q_weight_max*Qe_max;
                        
                    case "Q_target_min"
                        Q_target_min = param.Q_target_min;
                        Qe_min = sdpvar(1);
                        LMI = [LMI, ...
                               -Qe_min <= Q_min - Q_target_min, ...
                               -Qe_min <= -Q_min + Q_target_min, ...
                               Qe_min >= 0];
                        opt_J = opt_J + Q_weight_min*Qe_min;
                        
                    case "P_target_max"
                        P_target_max = param.P_target_max;
                        Pe_max = sdpvar(1);
                        LMI = [LMI, ...
                               -Pe_max <= P_max - P_target_max, ...
                               -Pe_max <= -P_max + P_target_max, ...
                               Pe_max >= 1];
                        opt_J = opt_J + P_weight_max*Pe_max;
                        
                    case "P_target_min"
                        P_target_min = param.P_target_min;
                        Pe_min = sdpvar(1);
                        LMI = [LMI, ...
                               -Pe_min <= P_min - P_target_min, ...
                               -Pe_min <= -P_min + P_target_min, ...
                               Pe_min >= 1];
                        opt_J = opt_J + P_weight_min*Pe_min;
                        
                    case "R_target_max"
                        R_target_max = param.R_target_max;
                        Re_max = sdpvar(1);
                        LMI = [LMI, ...
                               -Re_max <= R_max - R_target_max, ...
                               -Re_max <= -R_max + R_target_max, ...
                               Re_max >= 0];
                        opt_J = opt_J + R_weight_max*Re_max;
                        
                    case "R_target_min"
                        R_target_min = param.R_target_min;
                        Re_min = sdpvar(1);
                        LMI = [LMI, ...
                               -Re_min <= R_min - R_target_min, ...
                               -Re_min <= -R_min + R_target_min, ...
                               Re_min >= 0];
                        opt_J = opt_J + R_weight_min*Re_min;
                        
                        
                end
                
            end
            
            
            
            
            Q_LMI           = sdpvar(n,n);
            R_LMI           = sdpvar(n,n);
            P_LMI           = sdpvar(n,n);
            
            LMI         = [LMI, Q_LMI*AT+AT.'*Q_LMI == -P_LMI, ...
                Q_LMI*M + M*Q_LMI == -R_LMI, ...
                eye(n)*R_max>=R_LMI>=eye(n)*R_min, ...
                eye(n)*P_max>=P_LMI>=eye(n)*P_min, ...
                eye(n)*Q_max>=Q_LMI>=eye(n)*Q_min, ...
                Q_min >= Q_bound(1), ...
                R_min >= R_bound(1), ...
                P_min >= P_bound(1), ...
                Q_max <= Q_bound(2), ...
                R_max <= R_bound(2), ...
                P_max <= P_bound(2)];
            sol.out = optimize(LMI,opt_J,obj.opts);
            
            fprintf("Max: Target: %.4f\n\tPe_max: %.4f\n\tWeight: %.4f\n\tObj: %.4f\n", P_target_max, double(Pe_max), P_weight_max, P_weight_max*double(Pe_max));
            fprintf("Min: Target: %.4f\n\tPe_min: %.4f\n\tWeight: %.4f\n\tObj: %.4f\n", P_target_min, double(Pe_min), P_weight_min, P_weight_min*double(Pe_min));
            
            sol.Q = double(Q_LMI);
            sol.P = double(P_LMI);
            sol.R = double(R_LMI);
            sol.obj = double(opt_J);
            
            
            
        end
        
         function sol = GeneralSolve_push(obj, AT, P_bound, Q_bound, R_bound, param)
            
            n = max(size(AT));
            M = diag(-n:1:-1);
            
            LMI = [];
            
            Q_weight = 1;
            R_weight = 1;
            P_weight = 1;
            
            fields = fieldnames(param);
            
            for i = 1:numel(fields)
                
               switch fields{i}
                   case "Q_weight"
                       Q_weight = param.Q_weight;
                   case "P_weight"
                       P_weight = param.P_weight;
                   case "R_weight"
                       R_weight = param.R_weight;
               end
                
            end
            
            
            
            
            Q_LMI           = sdpvar(n,n);
            R_LMI           = sdpvar(n,n);
            P_LMI           = sdpvar(n,n);
            
            Q_min           = sdpvar(1);
            P_min           = sdpvar(1);
            R_max           = sdpvar(1);
            
            LMI         = [LMI, Q_LMI*AT+AT.'*Q_LMI == -P_LMI, ...
                            Q_LMI*M + M*Q_LMI == -R_LMI, ...
                            eye(n)*R_max>=R_LMI>=eye(n)*R_bound(1), ...
                            eye(n)*P_bound(2)>=P_LMI>=eye(n)*P_min, ...
                            eye(n)*Q_bound(2)>=Q_LMI>=eye(n)*Q_min, ...
                            Q_min >= Q_bound(1), ...
                            Q_min <= Q_bound(2), ...
                            P_min >= P_bound(1), ...
                            P_min <= P_bound(2), ...
                            R_max <= R_bound(2),...
                            R_max >= R_bound(1)];
            
            opt_J = -Q_weight*Q_min -P_weight*P_min + R_weight*R_max;
            
            sol.out = optimize(LMI,opt_J,obj.opts);
            
            sol.Q = double(Q_LMI);
            sol.P = double(P_LMI);
            sol.R = double(R_LMI);
            sol.obj = double(opt_J);
            
            
            
        end
        
        function dt_max = finite_time(obj, x0)
            
            n_IP = 4;
            n_CT = 2;
           
            x = x0([1;4;2;5]);
            z = x0([3;6]);
            Bx = obj.InPlane.TD\x;
            Bz = obj.CrossTrack.TD\z;
            
            [T_IP,~,ef] = fzero(@obj.Tsolve,[1e-10, 20000],[],obj.InPlane.Q,Bx);
            [T_CT,~,ef] = fzero(@obj.Tsolve,[1e-10, 20000],[],obj.CrossTrack.Q,Bz);
            
            del_IP = diag(T_IP.^(-n_IP:1:-1));
            y_IP   = del_IP*Bx;
            
            del_CT = diag(T_CT.^(-n_CT:1:-1));
            y_CT   = del_CT*Bz;
            
            
            
            M_IP = diag(-n_IP:1:-1);
            IP_num = (y_IP'*(obj.InPlane.Q*obj.InPlane.AT + obj.InPlane.AT'*obj.InPlane.Q)*y_IP);
            IP_den = (y_IP'*(M_IP*obj.InPlane.AT + M_IP*obj.InPlane.Q)*y_IP);
            T_dot_IP = -IP_num/IP_den;
            dt_IP = T_IP/T_dot_IP;
            
            M_CT = diag(-n_CT:1:-1);
            CT_num = (y_CT'*(obj.CrossTrack.Q*obj.CrossTrack.AT + obj.CrossTrack.AT'*obj.CrossTrack.Q)*y_CT);
            CT_den = (y_CT'*(M_CT*obj.CrossTrack.AT + M_CT*obj.CrossTrack.Q)*y_CT);
            T_dot_CT = -CT_num/CT_den;
            dt_CT = T_CT/T_dot_CT;
            
            dt_max = max(dt_IP, dt_CT);
            
            
            
        end
        
        function result = solve_push(obj, varargin)
            
            param.OptVar = "None";
            if ~isempty(varargin)
                for i = 1:2:length(varargin) % work for a list of name-value pairs
                    param.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
                
            end
            
            InPlaneOut = obj.GeneralSolve_push(obj.InPlane.AT, obj.InPlane.P_bound, obj.InPlane.Q_bound, obj.InPlane.R_bound, param);
            
            CrossTrackOut = obj.GeneralSolve_push(obj.CrossTrack.AT, obj.CrossTrack.P_bound, obj.CrossTrack.Q_bound, obj.CrossTrack.R_bound, param);
            
            obj.InPlane.Q = InPlaneOut.Q;
            obj.InPlane.P = InPlaneOut.P;
            obj.InPlane.R = InPlaneOut.R;
            obj.InPlane.obj = InPlaneOut.obj;
            
            obj.CrossTrack.Q = CrossTrackOut.Q;
            obj.CrossTrack.P = CrossTrackOut.P;
            obj.CrossTrack.R = CrossTrackOut.R;
            obj.CrossTrack.obj = CrossTrackOut.obj;
            
            result.InPlane.out = InPlaneOut.out;
            result.CrossTrack.out = CrossTrackOut.out;
            
            
        end
        
        function result = solve_target_ind(obj, varargin)
            
            param.OptVar = "None";
            if ~isempty(varargin)
                for i = 1:2:length(varargin) % work for a list of name-value pairs
                    param.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
                
            end
            
            InPlaneOut = obj.GeneralSolve_target_ind(obj.InPlane.AT, obj.InPlane.P_bound, obj.InPlane.Q_bound, obj.InPlane.R_bound, param);
            
            CrossTrackOut = obj.GeneralSolve_target_ind(obj.CrossTrack.AT, obj.CrossTrack.P_bound, obj.CrossTrack.Q_bound, obj.CrossTrack.R_bound, param);
            
            obj.InPlane.Q = InPlaneOut.Q;
            obj.InPlane.P = InPlaneOut.P;
            obj.InPlane.R = InPlaneOut.R;
            obj.InPlane.obj = InPlaneOut.obj;
            
            obj.CrossTrack.Q = CrossTrackOut.Q;
            obj.CrossTrack.P = CrossTrackOut.P;
            obj.CrossTrack.R = CrossTrackOut.R;
            obj.CrossTrack.obj = CrossTrackOut.obj;
            
            result.InPlane.out = InPlaneOut.out;
            result.CrossTrack.out = CrossTrackOut.out;
            
            
        end
        
        function result = solve_target(obj, varargin)
            
            param.OptVar = "None";
            if ~isempty(varargin)
                for i = 1:2:length(varargin) % work for a list of name-value pairs
                    param.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
                
            end
            
            InPlaneOut = obj.GeneralSolve_target(obj.InPlane.AT, obj.InPlane.P_bound, obj.InPlane.Q_bound, obj.InPlane.R_bound, param);
            
            CrossTrackOut = obj.GeneralSolve_target(obj.CrossTrack.AT, obj.CrossTrack.P_bound, obj.CrossTrack.Q_bound, obj.CrossTrack.R_bound, param);
            
            obj.InPlane.Q = InPlaneOut.Q;
            obj.InPlane.P = InPlaneOut.P;
            obj.InPlane.R = InPlaneOut.R;
            obj.InPlane.obj = InPlaneOut.obj;
            
            obj.CrossTrack.Q = CrossTrackOut.Q;
            obj.CrossTrack.P = CrossTrackOut.P;
            obj.CrossTrack.R = CrossTrackOut.R;
            obj.CrossTrack.obj = CrossTrackOut.obj;
            
            result.InPlane.out = InPlaneOut.out;
            result.CrossTrack.out = CrossTrackOut.out;
            
            
        end
        
        function result = solve_opt(obj, varargin)
            
            param.OptVar = "None";
            if ~isempty(varargin)
                for i = 1:2:length(varargin) % work for a list of name-value pairs
                    param.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
                
            end
            
            InPlaneOut = obj.GeneralSolve_opt(obj.InPlane.AT, obj.InPlane.P_bound, obj.InPlane.Q_bound, obj.InPlane.R_bound, param);
            
            CrossTrackOut = obj.GeneralSolve_opt(obj.CrossTrack.AT, obj.CrossTrack.P_bound, obj.CrossTrack.Q_bound, obj.CrossTrack.R_bound, param);
            
            obj.InPlane.Q = InPlaneOut.Q;
            obj.InPlane.P = InPlaneOut.P;
            obj.InPlane.R = InPlaneOut.R;
            obj.InPlane.obj = InPlaneOut.obj;
            
            obj.CrossTrack.Q = CrossTrackOut.Q;
            obj.CrossTrack.P = CrossTrackOut.P;
            obj.CrossTrack.R = CrossTrackOut.R;
            obj.CrossTrack.obj = CrossTrackOut.obj;
            
            result.InPlane.out = InPlaneOut.out;
            result.CrossTrack.out = CrossTrackOut.out;
            
            
        end
        
        function result = solve(obj)
            
            InPlaneOut = obj.GeneralSolve(obj.InPlane.AT, obj.InPlane.P_bound, obj.InPlane.Q_bound, obj.InPlane.R_bound);
            
            CrossTrackOut = obj.GeneralSolve(obj.CrossTrack.AT, obj.CrossTrack.P_bound, obj.CrossTrack.Q_bound, obj.CrossTrack.R_bound);
            
            obj.InPlane.Q = InPlaneOut.Q;
            obj.InPlane.P = InPlaneOut.P;
            obj.InPlane.R = InPlaneOut.R;
            
            obj.CrossTrack.Q = CrossTrackOut.Q;
            obj.CrossTrack.P = CrossTrackOut.P;
            obj.CrossTrack.R = CrossTrackOut.R;
            
            result.InPlane.out = InPlaneOut.out;
            result.CrossTrack.out = CrossTrackOut.out;
            
            
        end
        
        
        function [u, T_all, u_all] = ctrl(obj, x0)
            
            % Assume incoming control as x y z xd yd zd
            
            x = x0([1;4;2;5]);
            z = x0([3;6]);
            
            % In Plane Control
            %             if norm(x) < 1e-6
            %                 uy = 0;
            %             else
            
            n = size(obj.InPlane.A0,1);
            
            try
                Bx = obj.InPlane.TD\x;
                [T,~,ef] = fzero(@obj.Tsolve,[1e-10, 20000],[],obj.InPlane.Q,Bx);
                T_all.InPlane = T;
            catch
                T_all.InPlane = NaN;
                ef = 0;
            end
            
            
            if ef == 1
                del = diag(T.^(-n:1:-1));
                y   = del*Bx;
                v   = obj.InPlane.K*y;
            else
                v = 0;
            end
            
            uy = v+obj.InPlane.TC*x;
            
            u_all.InPlane = [v; obj.InPlane.TC*x];
            %             end
            
            % Cross Track Control
            %             if norm(z) < 1e-6
            %                 uz = 0;
            %             else
            n = size(obj.CrossTrack.A0,1);
            
            try
                Bz = obj.CrossTrack.TD\z;
                [T,~,ef] = fzero(@obj.Tsolve,[1e-10, 20000],[],obj.CrossTrack.Q,Bz);
                T_all.CrossTrack = T;
            catch
                T_all.CrossTrack = NaN;
                ef = 0;
            end
            
            
            if ef == 1
                del = diag(T.^(-n:1:-1));
                y   = del*Bz;
                v   = obj.CrossTrack.K*y;
            else
                v = 0;
            end
            
            uz = v+obj.CrossTrack.TC*z;
            
            u_all.CrossTrack = [v; obj.CrossTrack.TC*z];
            %             end
            
            u = [0;uy;uz];
            
            
        end
        
        function [u, u_IP, u_CT] = ctrl_K(obj, x0)
            
            % Assume incoming control as x y z xd yd zd
            
            x = x0([1;4;2;5]);
            z = x0([3;6]);            
            
            Bx = obj.InPlane.TD\x;
            u_IP_K = obj.InPlane.K*Bx; 
            u_IP_B = obj.InPlane.TC*x;
            u_IP = [u_IP_K, u_IP_B];
            
            Bz = obj.CrossTrack.TD\z;
            u_CT_K   = obj.CrossTrack.K*Bz;            
            u_CT_B = obj.CrossTrack.TC*z;
            u_CT = [u_CT_K, u_CT_B];
            
            u = [0;u_IP_K+u_IP_B;u_CT_K+u_CT_B];
            
            
        end
        
        function status = check(obj)
            
            n = max(size(obj.InPlane.A0));
            M = diag(-n:1:-1);
            
            % All checks should be zero...
            ck1 = obj.InPlane.Q*obj.InPlane.AT + obj.InPlane.AT.'*obj.InPlane.Q + obj.InPlane.P;
            ck2 = obj.InPlane.Q*M + M*obj.InPlane.Q + obj.InPlane.R;
            ck3 = sum( eig(obj.InPlane.Q)>0 ) - n;
            ck4 = sum( eig(obj.InPlane.R)>0 ) - n;
            
            x = ones(n,1)*.001;
            
            [T,~,~] = fzero(@obj.Tsolve,norm(x),[],obj.InPlane.Q,x);
            
            del = diag(T.^(-n:1:-1));
            ck5 = x.'*del*obj.InPlane.Q*del*x-1;
            
            status.InPlane = [det(ck1) det(ck2) ck3 ck4 det(ck5)];
            
            
            %% Cross Track
            n = max(size(obj.CrossTrack.A0));
            M = diag(-n:1:-1);
            
            % All checks should be zero...
            ck1 = obj.CrossTrack.Q*obj.CrossTrack.AT + obj.CrossTrack.AT.'*obj.CrossTrack.Q + obj.CrossTrack.P;
            ck2 = obj.CrossTrack.Q*M + M*obj.CrossTrack.Q + obj.CrossTrack.R;
            ck3 = sum( eig(obj.CrossTrack.Q)>0 ) - n;
            ck4 = sum( eig(obj.CrossTrack.R)>0 ) - n;
            
            x = ones(n,1)*.001;
            
            [T,~,~] = fzero(@obj.Tsolve,norm(x),[],obj.CrossTrack.Q,x);
            
            del = diag(T.^(-n:1:-1));
            ck5 = x.'*del*obj.CrossTrack.Q*del*x-1;
            
            status.CrossTrack = [det(ck1) det(ck2) ck3 ck4 det(ck5)];
            
        end
        
        
    end
    
    
    
    
    
    
    methods(Static)
        function F = Tsolve(T,Q,x)
            
            n = max(size(x));
            del = diag(T.^(-n:1:-1));
            F = x.'*del*Q*del*x-1;
            
        end
        
        function T = calc_T(x0, Q)
            
            [T,~,~] = fzero(@CommonLyap.Tsolve,[1e-10, 20000],[],Q,x0);
            
        end
        
        
        
        
    end
    
    
    
    
end