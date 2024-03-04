classdef MagneticUnit
    properties
        value;
        tag;
    end
    
    methods
        function obj = MagneticUnit(value, tag)
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
    
    enumeration % Number is factor conversion tesla -> unit
        tesla           (1, 'T');
        ampere_metre    (2500000 / pi, 'A/m');
        gauss           (1E4, 'G');
        oersted         (1E4, 'Oe');
    end
end