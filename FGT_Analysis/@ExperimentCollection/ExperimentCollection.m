classdef ExperimentCollection < handle
    %EXPERIMENTCOLLECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        experimentArray = Experiment;
        nExperiments = 0;
        indexingValues = [];
    end
    
    methods
        function obj = ExperimentCollection(fileNames, fileParams)
            %EXPERIMENTCOLLECTION Construct an instance of this class
            %   Detailed explanation goes here

            % Check that the number of file parameters is equal
            % to the number of file names. Or if there is only one
            % file parameter, make nFiles copies of it.

            if nargin == 0, return; end

            nFiles = length(fileNames);
            if length(fileParams) ~= nFiles
                if length(fileParams) ~= 1
                    disp(['ERROR The number of file parameters ', ...
                          'must be equal to the number of files or 1.']);
                    return
                end
                fileParams = repelem(fileParams, nFiles);
            end

            % Create Experiment array and write number of experiments
            obj.experimentArray(nFiles) = Experiment();
            for iFile = 1:nFiles
                obj.experimentArray(iFile) = Experiment(fileNames(iFile), ...
                                                        fileParams(iFile));
            end

            obj.nExperiments = nFiles;

            % Create indexing array
            obj.indexingValues = zeros(nFiles,1);
            for iFile = 1:nFiles
                obj.indexingValues = getindex(obj.experimentArray(iFile));
            end

        end
        
        function addexperiment(obj, experiment)
            %ADDEXPERIMENT - Add an experiment to the experiment collection
            %   Add an experiment to the experiment collection 
            %   and update the pertinent variables.
            %
            %   Syntax
            %       ADDEXPERIMENT(obj, experiment)
            %
            %   Input Arguments
            %       obj - Experiment collection to add the experiment to
            %           ExperimentCollection
            %       experiment - Experiment to be added
            %           Experiment
        
            obj.nExperiments = obj.nExperiments + 1;
            obj.experimentArray(obj.nExperiments) = experiment;
            obj.indexingValues(obj.nExperiments) = experiment.getindex();
        end

        function removeexperiment(obj, n)
            %REMOVEEXPERIMENT - Remove an experiment from the collection
            %   Remove an existing experiment from the collection
            %   by specifying the position it occupies in the array
            %
            %   Syntax
            %       REMOVEEXPERIMENT(obj, n)
            %
            %   Input Arguments
            %       obj - Experiment collection to remove experiment from
            %           ExperimentCollection
            %       n - Experiment number(s) in the experiment array
            %           scalar | vector

            % Check if n is a valid number
            if ~isnumeric(n) || n > obj.nExperiments || n < 1
                disp(['ERROR Experiment number must be a valid ', ...
                     'number in the array']);
                return;
            end

            % Remove experiment from array
            obj.experimentArray(n) = [];
            obj.indexingValues(n) = [];
            obj.nExperiments = obj.nExperiments - length(n);
        end

        function convertfieldunits(obj, newUnits)
            %CONVERTFIELDUNITS - Convert the units of magnetic field
            %
            %   Syntax
            %       CONVERTFIELDUNITS(obj, newUnits)
            %
            %   Input Arguments
            %       obj - The experiment
            %           ExperimentCollection
            %       newUnits - Units to convert to
            %           MagneticUnit | UnitStruct         

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).convertfieldunits(newUnits);
            end
        end

        function convertfrequencyunits(obj, newUnits)
            %CONVERTFREQUENCYUNITS - Convert the frequency scale of units
            %
            %   Syntax
            %       CONVERTFREQUENCYUNITS(obj, newUnits)
            %
            %   Input Arguments
            %       obj - The experiment
            %           ExperimentCollection
            %       newUnits - Scale to convert to
            %           FrequencyUnit | UnitStruct         

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).convertfrequencyunits(newUnits);
            end
        end

        function resonanceanalysis(obj, iExp)
            %RESONANCEANALYSIS - Makes the resonance analysis
            %   This function saves frequency, ressonant field,
            %   FWHM and the errors involved in the resonanceParameters 
            %   property.
            %
            %   Syntaxis
            %       RESONANCEANALYSIS(obj, iExp)
            %
            %   Input Arguments
            %       obj - Experiment collection
            %           ExperimentCollection
            %       iExp - Experiment number in experimentArray
            %           scalar

            if nargin == 2 && isnumeric(iExp)
                disp(['Resonance analysis of experiment #', ...
                     num2str(iExp), ' with indexing parameter,', ...
                     num2str(obj.experimentArray(iExp).getindex), ...
                     ' ', obj.experimentArray(iExp).indexingUnits]);
                obj.experimentArray(iExp).resonanceanalysis();
                return
            end

            for iExp = 1:obj.nExperiments
                disp(['Resonance analysis of experiment #', ...
                     num2str(iExp), ' with indexing parameter,', ...
                     num2str(obj.experimentArray(iExp).getindex), ...
                     ' ', obj.experimentArray(iExp).indexingUnits]);
                obj.experimentArray(iExp).resonanceanalysis();
            end
        end

        function makefieldpositive(obj, expNumber)
            %MAKEFIELDPOSITIVE - Make the field values positive (abs(.))
            %   Takes the absolute value of the field
            %
            %   Syntax
            %       MAKEFIELDPOSITIVE(obj, iExp)
            %
            %   Input Arguments
            %       obj - Experiment collection
            %           ExperimentCollection
            %       iExp - Experiment number in experimentArray
            %           scalar

            if nargin == 2 && isnumeric(expNumber)
                obj.experimentArray(expNumber).makefieldpositive();
                return
            end

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).makefieldpositive();
            end
        end

        function makefieldnegative(obj, expNumber)
            %MAKEFIELDNEGATIVE - Make the field values negative (-abs(.))
            %   Takes the absolute value of the field and changes sign.
            %
            %   Syntax
            %       MAKEFIELDNEGATIVE(obj, iExp)
            %
            %   Input Arguments
            %       obj - Experiment collection
            %           ExperimentCollection
            %       iExp - Experiment number in experimentArray
            %           scalar

            if nargin == 2 && isnumeric(expNumber)
                obj.experimentArray(expNumber).makefieldnegative();
                return
            end

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).makefieldnegative();
            end
        end

        function makeinplanekittelfit(obj, iExp)
            %MAKEINPLANEKITTELFIT - Make Kittel's law fit w/ resonance data 
            %   Makes the best Kittel equation parameter fit with the
            %   experiment resonance data.
            %   Kittel in plane equation:
            %       frequency = ɣ/2π μ sqrt((H + Hani) * (H + Hani + Meff))
            %
            %       Meff = Ms - Hu
            %   where Ms is the saturation magnetization and Hu the 
            %   uniaxial out-of-plane anisotropy field.
            %
            %   Syntax
            %       MAKEINPLANEKITTELFIT(obj, iExp)
            %
            %   Input Arguments
            %       obj - Experiment collection
            %           ExperimentCollection
            %       iExp - Experiment number in experimentArray
            %           scalar

            if nargin == 2 && isnumeric(iExp)
                obj.experimentArray(iExp).makeinplanekittelfit();
                return
            end

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).makeinplanekittelfit();
            end
        end

        function makeoutofplanekittelfit(obj, iExp)
            %MAKEOUTOFPLANEKITTELFIT - Make Kittel's law fit with data
            %   Makes the best Kittel equation parameter fit with the
            %   experiment resonance data.
            %   Kittel out-of-plane equation:
            %       frequency = ɣ/2π * μ0 (H - Meff)
            %
            %       Meff = Ms + Hani - Hu
            %   where Ms is the saturation magnetization and Hani the
            %   magnetocrystalline anisotropy field and Hu the uniaxial
            %   out-of-plane anisotropy field.
            %
            %   Syntax
            %       MAKEOUTOFPLANEKITTELFIT(obj, iExp)
            %
            %   Input Arguments
            %       obj - Experiment collection
            %           ExperimentCollection
            %       iExp - Experiment number in experimentArray
            %           scalar

            if nargin == 2 && isnumeric(iExp)
                obj.experimentArray(iExp).makeoutofplanekittelfit();
                return
            end

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).makeoutofplanekittelfit();
            end
        end

        function makedampingfit(obj, iExp)
            %MAKEDAMPINGFIT - Make damping fit and get damping parameters
            %
            %   Syntax
            %       MAKEDAMPINGFIT(obj, iExp)
            %
            %   Input Arguments
            %       obj - Experiment
            %           Experiment
            %       iExp - Experiment number in experimentArray
            %           scalar

            if nargin == 2 && isnumeric(iExp)
                obj.experimentArray(iExp).makedampingfit();
                return
            end

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).makedampingfit();
            end
        end

% ------------------------------------------------------------------------

        function filter(obj, isequalindex)
            %FILTER - Filter the experiments by index value
            %   Filter the experiments by index value. Equality
            %   of indices is determined by the input anonymous
            %   functions.
            %
            %   Syntax
            %       FILTER(obj, isequalindex)
            %
            %   Input Arguments
            %       obj - The experiment collection to filter
            %           ExperimentCollection
            %       isequalindex - Function with two inputs that states
            %                      equality between two values
            %           function_handle
            %
            %   Example
            %       
            %       expColl.filter(@(x,y) abs(x,y) < 1)
            %       Filters values in different experiments by their
            %       indexing values. The indexing values are considered
            %       to be equal if they differ by less than 1.

            % Filter all experiments
            for iExp = 1:obj.nExperiments
                expVar = obj.experimentArray(1);
                % Get data from experiment
                expData = expVar.getdata();
                indexingColumn = expVar.getindexingcolumn();
                fileParam = FileParameters(1,2,3,indexingColumn);

                % Filter by indexing parameter all the rows of
                % the data matrix
                while ~isempty(expData)
                    % Take the first index and look for other "equal"
                    % indices in the matrix
                    refIdx = expData(1, indexingColumn);
                    [idxData, expData] = popequalidx(expData, ...
                                                     indexingColumn, ...
                                                     refIdx, ...
                                                     isequalindex);
                    
                    % Create and add experiment with data
                    if ~isempty(idxData)
                        newExp = Experiment(idxData, fileParam);
                        newExp.setfieldunits(expVar.getfieldunits());
                        newExp.setfrequencyunits(expVar.getfrequencyunits());
                        newExp.setindexingunits(expVar.getindexingunits());

                        obj.addexperiment(newExp);

                    end
                end

                % Remove the filtered experiment
                obj.removeexperiment(1);
            end
        end
        
        function filterbytemperature(obj, delta)
            %FILTERBYTEMPERATURE - Filter experiments by temperature
            %   Filter experiments by temperature. Two temperatures
            %   are considered equal if |T1 - T2| < delta.
            %
            %   Syntax
            %       FILTERBYTEMPERATURE(obj,delta)
            %
            %   Input Arguments
            %       obj - Experiment collection to filter
            %           ExperimentCollection
            %       delta - Interval of temperature to be considered equal
            %           scalar

            if nargin == 1, delta = 1; end

            % Set indexing column to fourth column
            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setindexingcolumn(4);
            end

            obj.filter(@(x,y) abs(x - y) < delta);
        end
    
        function filterbyfieldsign(obj)
            %FILTERBYFIELDSIGN - Filter experiments by field sign (+/-)
            %   Filter experiments by field sign: positive and negative.
            %
            %   Syntax
            %       FILTERBYFIELDSIGN(obj)
            %
            %   Input Arguments
            %       obj - Experiment collection to filter
            %           ExperimentCollection

            % Set indexing column to field column
            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setindexingcolumn(2);
            end

            % Filter by sign
            obj.filter(@(x,y) 0.5 * abs(sign(x) + sign(y)));

            % Set indexing parameter to +1 or -1 matching the sign
            for iExp = 1:obj.nExperiments
                expVar = obj.experimentArray(iExp);
                indexValue = sign(expVar.getindex);

                expVar.setindex(indexValue);
                obj.indexingValues(iExp) = indexValue;
            end
        end

        function filterpositivefield(obj)
            %FILTERPOSITIVEFIELD - Filter experiments by positive field
            %   Keep only the part of the experiments with positive
            %   magnetic field values.
            %
            %   Syntax
            %       FILTERPOSITIVEFIELD(obj)
            %
            %   Input Arguments
            %       obj - Experiment collection to filter
            %           ExperimentCollection

            % Set indexing column to field column
            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setindexingcolumn(2);
            end

            % Filter by positive values
            obj.filter(@(x,y) 0.5 * (1 + sign(y)) );

            % Set indexing parameter to +1
            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setindex(1);
                obj.indexingValues(iExp) = 1;
            end
        end

        function filternegativefield(obj)
            %FILTERNEGATIVEFIELD - Filter experiments by negative field
            %   Keep only the part of the experiments with negative
            %   magnetic field values.
            %
            %   Syntax
            %       FILTERNEGATIVEFIELD(obj)
            %
            %   Input Arguments
            %       obj - Experiment collection to filter
            %           ExperimentCollection

            % Set indexing column to field column
            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setindexingcolumn(2);
            end

            % Filter by negative values
            obj.filter(@(x,y) 0.5 * (1 - sign(y)) );

            % Set indexing parameter to -1
            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setindex(-1);
                obj.indexingValues(iExp) = -1;
            end
        end

% ------------------------------------------------------------------------

        function setindexingvalues(obj, newIndexingValues)
            %SETINDEXINGVALUES - Set the value to indexing values
            %   Set the value of indexing values of experiments.
            %   The value will also be changed in the corresponding
            %   experiment variable.
            %
            %   Syntax
            %       SETINDEXINGVALUES(obj, newIndexingValues)
            %
            %   Input Arguments
            %       obj - The experiment collection
            %           ExperimentCollection
            %       newIndexingValues - New indexing values
            %           vector

            if length(newIndexingValues) ~= obj.nExperiments
                disp(['ERROR Input indexing values must be the ' ...
                     'same length as existing indexing values'])
                return;
            end

            obj.indexingValues = newIndexingValues;
            for idx = 1:obj.nExperiments
                obj.experimentArray(idx).setindex(newIndexingValues(idx));
            end
        end
    
        function setfieldunits(obj, unit)
            %SETFIELDUNITS - Set field unit of measurement for experiments
            %
            %   Syntax
            %       SETFIELDUNITS(obj, unit)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           ExperimentCollection
            %       unit - Magnetic unit
            %           MagneticUnit | UnitStruct
            %
            %   Example
            %       SETFIELDUNITS(obj, MagneticUnit.gauss, Scale.milli)
            %       Set the field units to mG (milligauss).

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setfieldunits(unit);
            end
        end

        function setfrequencyunits(obj, unit)
            %SETFREQUENCYUNITS - Set frequency scale of measurements
            %
            %   Syntax
            %       SETFREQUENCYUNITS(obj, unit)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           Experiment
            %       unit - Units of frequency
            %           FrequencyUnit | UnitStruct

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setfrequencyunits(unit);
            end
        end

        function setindexingunits(obj, unitName)
            %SETINDEXINGUNITS - Set indexing unit of measurements
            %
            %   Syntax
            %       SETINDEXINGUNITS(obj, unitName)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           ExperimentCollection
            %       unitName - New unit of measurement name
            %           string

            for iExp = 1:obj.nExperiments
                obj.experimentArray(iExp).setindexingunits(unitName);
            end
        end

    end
end

function [idxData, data] = popequalidx(data, idxCol, refIdx, isequalindex)
    %POPEQUALIDX - Pops the rows of data with the specified index
    %   Returns the rows of data with the same index as the reference
    %   index (to the eyes of the isequalidx function) and deletes
    %   the rows from the data
    %
    %   Syntax
    %       [idxData, data] = POPEQUALIDX(data, idxCol, refIdx, isequalidx)
    %
    %   Input Arguments
    %       data - Data to pop the data from
    %           matrix
    %       idxCol - Indexing column number
    %           scalar
    %       refIdx - Reference index. Rows with this index will be popped
    %           scalar
    %       isequalidx - Function of equality between indices
    %           function_handle
    %
    %   Output Arguments
    %       idxData - Data with the same index parameter as refIdx
    %           matrix
    %       data - Data with different index parameter as refIdx
    %           matrix

    indexingColumn = data(:,idxCol);
    isEqualIndex = logical(isequalindex(refIdx, indexingColumn));
    
    idxData = data(isEqualIndex, :);
    data(isEqualIndex, :) = [];

    if isempty(idxData)
        data = [];
    end
end