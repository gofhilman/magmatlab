classdef FileParameters < handle
    %FILEPARAMETERS - Parameters of the input file   
    %   The parameters of the input file are numbers of the columns
    %   containing the field, frequency, gain and indexing values
    %
    %   Creation
    %       P = FILEPARAMETERS(fieldCol, freqCol, gainCol, indexingCol)
    %
    %   Input Arguments
    %       fieldCol - Number of the field column
    %           scalar
    %       freqCol - Number of the frequency column
    %           scalar
    %       gainCol - Number of the gain (S12) column
    %           scalar
    %       indexingCol - Number of the indexing variable column
    %           scalar
    %
    %   Output Arguments
    %       P - Output parameters
    %           File parameters
    properties (Access = private)
        frequency;
        field;
        gain;
        indexing;
    end
    
    methods
        function obj = FileParameters(freqCol, fieldCol, ...
                                      gainCol, indexingCol)
            %FILEPARAMETERS - Construct an instance of this class
            obj.frequency = freqCol;
            obj.field = fieldCol;
            obj.gain = gainCol;
            obj.indexing = indexingCol;
        end

% ------------------------------------------------------------------------

        function out = getfrequency(obj)
            %GETFREQUENCY - Return the value of the frequency column
            %   Output the number of frequency column stored in variable
            %
            %   Syntax
            %       c = GETFREQUENCY(obj)
            %
            %   Input parameters
            %       obj - The file parameters
            %           File parameters
            %
            %   Output parameters
            %       c - The value of the frequency parameter
            %           scalar
            out = obj.frequency;
        end

        function out = getfield(obj)
            %GETFIELD - Return the value of the field column
            %   Output the number of field column stored in variable
            %
            %   Syntax
            %       c = GETFIELD(obj)
            %
            %   Input parameters
            %       obj - The file parameters
            %           File parameters
            %
            %   Output parameters
            %       c - The value of the field parameter
            %           scalar
            out = obj.field;
        end

        function out = getgain(obj)
            %GETGAIN - Return the value of the gain column
            %   Output the number of gain column stored in variable
            %
            %   Syntax
            %       c = GETGAIN(obj)
            %
            %   Input parameters
            %       obj - The file parameters
            %           File parameters
            %
            %   Output parameters
            %       c - The value of the gain parameter
            %           scalar
            out = obj.gain;
        end

        function out = getindexing(obj)
            %GETINDEXING - Return the value of the indexing value column
            %   Output the number of indexing column stored in variable
            %
            %   Syntax
            %       c = GETINDEXING(obj)
            %
            %   Input parameters
            %       obj - The file parameters
            %           File parameters
            %
            %   Output parameters
            %       c - The value of the indexing parameter
            %           scalar
            out = obj.indexing;
        end

        function out = getparameters(obj)
            %GETPARAMETERS - Return an array with the file parameters
            %   Output an array containing the file parameters
            %
            %   Syntax
            %       [freqC, fieldC, gainC, indexingC] = GETPARAMETERS(obj)
            %
            %   Input parameters
            %       obj - The file parameters
            %           File parameters
            %
            %   Output parameters
            %       freqC - Number of the frequency column
            %           scalar
            %       fieldC - Number of the field column
            %           scalar
            %       gainC - Number of the gain (S12) column
            %           scalar
            %       indexingC - Number of the indexing variable column
            %           scalar
            out = [obj.frequency, obj.field, obj.gain, obj.indexing];
        end

% ------------------------------------------------------------------------

        function setfrequencycolumn(obj, freqCol)
            %SETFREQUENCYCOLUMN - Set the value of frequency parameter
            %   Changes the value of the frequency parameter
            %
            %   Syntax
            %       SETFREQUENCYCOLUMN(obj, freqCol)
            %
            %   Input Arguments
            %       freqCol - New value for frequency parameter
            %           scalar
            obj.frequency = freqCol;
        end

        function setfieldcolumn(obj, fieldCol)
            %SETFIELDCOLUMN - Set the value of field parameter
            %   Changes the value of the field parameter
            %
            %   Syntax
            %       SETFIELDCOLUMN(obj, fieldCol)
            %
            %   Input Arguments
            %       fieldCol - New value for field parameter
            %           scalar
            obj.field = fieldCol;
        end

        function setgaincolumn(obj, gainCol)
            %SETGAINCOLUMN - Set the value of gain parameter
            %   Changes the value of the gain parameter
            %
            %   Syntax
            %       SETGAINCOLUMN(obj, gainCol)
            %
            %   Input Arguments
            %       gainCol - New value for gain parameter
            %           scalar
            obj.gain = gainCol;
        end

        function setindexingcolumn(obj, indexingCol)
            %SETINDEXINGCOLUMN - Set the value of indexing value parameter
            %   Changes the value of the indexing value parameter
            %
            %   Syntax
            %       SETINDEXINGCOLUMN(obj, indexingCol)
            %
            %   Input Arguments
            %       indexingCol - New value for indexing parameter
            %           scalar
            obj.indexing = indexingCol;
        end
    end
end