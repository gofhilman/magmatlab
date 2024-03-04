classdef Scale
    %SCALE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        value;
        tag;
    end

    methods
        function obj = Scale(value, tag)
            obj.value = value;
            obj.tag = tag;
        end

        function uStruct = mtimes(scale, unit)
            %MTIMES - Overload * operator to multiply units and scales

            uStruct = UnitStruct(scale.value * unit.value, ...
                                 append(scale.tag, unit.tag));
        end
    end

    enumeration 
        peta    (1E+15, 'P');
        tera    (1E+12, 'T');
        giga    (1E+09, 'G');
        mega    (1E+06, 'M');
        kilo    (1E+03, 'k');
        hecto   (1E+02, 'h');
        deca    (1E+01, 'da');
        one     (1E+00, '');
        deci    (1E-01, 'd');
        centi   (1E-02, 'c');
        milli   (1E-03, 'm');
        micro   (1E-06, 'μ');
        nano    (1E-09, 'n');
        pico    (1E-12, 'p');
        femto   (1E-15, 'f');
    end
end

