classdef KittelParameters
    %KITTELPARAMETERS - Class to store the Kittel fit parameters
    %   This class stores the Kittel parameters after the fit
    %
    %   Properties
    %       gyromagneticRatio - Gyromagnetic ratio
    %           PhysicalMagnitude
    %       anisotropyField - Anisotropy field
    %           PhysicalMagnitude
    %       effectiveMagnetization - Effective magnetization
    %           PhysicalMagnitude
    
    properties
        gyromagneticRatio;
        anisotropyField;
        effectiveMagnetization;
        magneticSusceptibility;
    end
    
    methods
        function obj = KittelParameters(gyroPhysM, aniPhysM, magnPhysM, ...
                                        magnSuscM)
            %KITTELPARAMETERS - Construct an instance of this class
            %
            %   Syntax
            %       obj = KittelParameters(gyroPhysM, aniPhysM, magnPhysM)
            %
            %   Input Arguments
            %       gyroPhysM - Gyromagnetic ratio physical magnitude
            %           PhysicalMagnitude
            %       aniPhysM - Anisotropy field physical magnitude
            %           PhysicalMagnitude
            %       magnPhysM - Effective magnetization physical magnitude
            %           PhysicalMagnitude
            %       magnSuscM - Magnetic susceptibility physical magnitude
            %           PhysicalMagnitude
            %
            %   Output Arguments
            %       obj - Instance of KittelParameters
            %           KittelParameters

            obj.gyromagneticRatio = gyroPhysM;
            obj.anisotropyField = aniPhysM;
            obj.effectiveMagnetization = magnPhysM;
            obj.magneticSusceptibility = magnSuscM;
        end
    end
end