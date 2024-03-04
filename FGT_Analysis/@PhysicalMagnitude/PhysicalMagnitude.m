classdef PhysicalMagnitude
    %PHYSICALMAGNITUDE - Class to store value, error and unit of magnitude
    %   Detailed explanation goes here
    
    properties
        value = 0;
        error = 0;
        unit;
    end
    
    methods
        function obj = PhysicalMagnitude(value, error, unit)
            %PHYSICALMAGNITUDE Construct an instance of this class
            %   Detailed explanation goes here

            if nargin == 0, return; end
            
            obj.value = value;
            obj.error = error;
            obj.unit = unit;
        end
    end
end

