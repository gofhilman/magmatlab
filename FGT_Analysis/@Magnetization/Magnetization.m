classdef Magnetization
    %MAGNETIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type;
    end
    
    methods
        function obj = Magnetization(type)
            %MAGNETIZATION Construct an instance of this class
            %   Detailed explanation goes here
            obj.type = type;
        end
    end

    enumeration
        none            (0);
        inPlane         (1);
        outOfPlane      (2);
        paramagnetic    (3);
    end
end

