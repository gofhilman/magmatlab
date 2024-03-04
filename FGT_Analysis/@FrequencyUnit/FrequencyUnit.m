classdef FrequencyUnit
    properties
        value;
        tag;
    end
    
    methods
        function obj = FrequencyUnit(value, tag)
            obj.value = value;
            obj.tag = tag;
        end

        function uStruct = mtimes(scale, unit)
            %MTIMES - Overload * operator to multiply units and scales

            uStruct = UnitStruct(scale.value * unit.value, ...
                                 append(scale.tag, unit.tag));
        end

        function uStruct = mrdivide(unit1, unit2)
            %MTIMES - Overload / operator to divide units

            uStruct = UnitStruct(unit1.value / unit2.value, ...
                                 append(unit1.tag, '/', unit2.tag));
        end
    end
    
    enumeration
        hertz   (1, 'Hz');
    end
end