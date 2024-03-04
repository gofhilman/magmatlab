classdef Filter
    %FILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
    end
    
    methods
        function obj = Filter(type)
            %FILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.type = type;
        end

    end

    enumeration
        none                (0);
        normalization       (1);
        noiseCancelling     (2);
    end
end

