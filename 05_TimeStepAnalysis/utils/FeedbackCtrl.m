classdef FeedbackCtrl < handle

    properties
    K
    end

    methods

        function obj = FeedbackCtrl()

        end

        function obj = setK(obj,K)
            obj.K = K;
        end

        function u = calc_u(obj, t, x)
            u = obj.K*x;
        end
    end


end