classdef DampingParameters
    %DAMPINGPARAMETERS - Class to store the damping fit parameters
    %   This class stores the damping parameters after the fit
    %
    %   Properties
    %       damping - Gilbert damping coefficient
    %           PhysicalMagnitude
    %       inhomogeneousDamping - inhomogeneousDamping
    %           PhysicalMagnitude
    
    properties
        damping;
        inhomogeneousDamping;
    end
    
    methods
        function obj = DampingParameters(damping, inhDamping)
            %DAMPINGPARAMETERS - Construct an instance of this class
            %
            %   Syntax
            %       obj = DAMPINGPARAMETERS(damping, inhDamping)
            %
            %   Input Arguments
            %       damping - Gilbert damping coefficient
            %           PhysicalMagnitude
            %       inhDamping - inhomogeneousDamping
            %           PhysicalMagnitude
            %
            %   Output Arguments
            %       obj - Instance of KittelParameters
            %           KittelParameters

            obj.damping = damping;
            obj.inhomogeneousDamping = inhDamping;
        end
    end
end