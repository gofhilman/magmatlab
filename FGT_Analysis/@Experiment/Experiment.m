classdef Experiment < handle
    %EXPERIMENT - Experiment class variable   
    %   This class stores the data from an experiment, the
    %   number of frequencies and an indexing parameter like
    %   temperature, angle...
    %
    %   If no indexing variable is proportioned, a mean will be
    %   calculated from a specified column in the file parameters.
    %   If no such column was specified, the value will be 0.
    %
    %   Creation
    %       E = EXPERIMENT(fileName, fileParams)
    %       E = EXPERIMENT(fileName, fileParams, index)
    %
    %       E = EXPERIMENT(dataArray, fileParams)
    %       E = EXPERIMENT(dataArray, fileParams, index)
    %
    %   Input Arguments
    %       fileName - Name of the file containing the data
    %           string scalar
    %       fileParams - File column parameters
    %           FileParameters
    %       index - Value of the indexing variable (optional)
    %           scalar
    %       dataArray - Array of data ordered as in fileParams
    %           matrix
    %
    %   Output Arguments
    %       E - Experiment class
    %           Experiment
    
    properties 
        data = [];
        resonanceParameters = []; % [freq, hr, fwhm, dhr, dfwhm]
        crossingPoints = [];
        magnetization = Magnetization.none;

        nFrequencies = 0;
        indexingColumn = 0;
        indexingParameter = 0;
        
        fieldUnits;
        frequencyUnits;
        indexingUnits = '';

        kittelParameters;
        dampingParameters;
    end
    
    methods
        function E = Experiment(inputData, fileParams)
            %EXPERIMENT - Construct an instance of this class
            %
            %   Syntax
            %       E = EXPERIMENT(fileName, fileParams)
            %       E = EXPERIMENT(fileName, fileParams, index)
            %
            %       E = EXPERIMENT(dataArray, fileParams)
            %       E = EXPERIMENT(dataArray, fileParams, index)
            %
            %   Input Arguments
            %       fileName - Name of the file containing the data
            %           string scalar
            %       fileParams - File column parameters
            %           FileParameters
            %       index - Value of the indexing variable (optional)
            %           scalar
            %       dataArray - Array of data ordered as in fileParams
            %           matrix
            %
            %   Output Arguments
            %       E - Experiment class
            %           Experiment

            % If not enough input arguments, return empty experiment
            if nargin < 2, return; end

            % Import data, sort columns and rows.
            % Calculate the number of frequency values
            % and magnetic field values.
            if isstring(inputData)
                E.data = sortcolumns(load(inputData), fileParams);
            elseif isnumeric(inputData)
                E.data = sortcolumns(inputData, fileParams);
            else
                return;
            end
            E.data = sortrows(E.data);

            E.nFrequencies = countfrequencies(E.data);

            % Write indexing parameter
            indexingCol = fileParams.getindexing();
            if indexingCol > 0 % Valid indexing column
                E.indexingColumn = indexingCol;
                E.indexingParameter = mean(E.data(:,indexingCol));
            else
                E.indexingParameter = 0;
            end

        end
       
% ------------------------------------------------------------------------

        function convertfieldunits(obj, newUnits)
            %CONVERTFIELDUNITS - Convert the units of magnetic field
            %
            %   Syntax
            %       CONVERTFIELDUNITS(obj, newUnits)
            %
            %   Input Arguments
            %       obj - The experiment
            %           Experiment
            %       newUnits - Units to convert to
            %           MagneticUnit | UnitStruct         

            if isempty(obj.fieldUnits)
                disp('ERROR Field has no units. Specify units first');
                return
            end

            conversionFactor = newUnits.value / obj.fieldUnits.value;
            obj.data(:,2) =  conversionFactor * obj.data(:,2);
            obj.fieldUnits = newUnits;
        end

        function convertfrequencyunits(obj, newUnits)
            %CONVERTFREQUENCYUNITS - Convert the frequency scale of units
            %
            %   Syntax
            %       CONVERTFREQUENCYUNITS(obj, newScale)
            %
            %   Input Arguments
            %       obj - The experiment
            %           Experiment
            %       newUnits - Scale to convert to
            %           FrequencyUnit | UnitStruct        

            if isempty(obj.frequencyUnits)
                disp('ERROR Frequency has no scale set yet');
                return
            end

            conversionFactor = obj.frequencyUnits.value / newUnits.value;
            obj.data(:,1) =  conversionFactor * obj.data(:,1);
            obj.frequencyUnits = newUnits;
        end

        function resonanceanalysis(obj)
            %RESONANCEANALYSIS - Makes the resonance analysis
            %   This function saves frequency, ressonant field,
            %   FWHM and the errors involved in the resonanceParameters 
            %   property.
            %
            %   Syntaxis
            %       RESONANCEANALYSIS(obj)
            %
            %   Input Arguments
            %       obj - Experiment to make resonance analysis
            %           Experiment

            experimentData = obj.data;

            % Get frequency values
            freqValues = unique(experimentData(:,1));

            % Find resonance parameters for each frequency and
            % store them in the resonanceParameters property
            for iFreq = 1:obj.nFrequencies
                freq = freqValues(iFreq);

                % Extract data at given frequency
                monoFreqData = filterdatabyfreq(experimentData, freq);

                % Find resonance parameters
                resonanceParams = findresonance(monoFreqData, ...
                                                obj.fieldUnits, ...
                                                obj.frequencyUnits);

                % If found, append to existing resonance parameters
                if ~isempty(resonanceParams)
                    obj.resonanceParameters = [obj.resonanceParameters; ...
                                               freq, resonanceParams];
                end
            end
        end

        function crossP = getcrossingpoints(obj, filter)
            
            % Plot SC Image
            gain = obj.data(:, Column.gain);
            field = unique(obj.data(:, Column.field));
            freq = flip(unique(obj.data(:, Column.frequency)), 1);

            nFreq = obj.nFrequencies;
            nField = floor(size(gain, 1) / nFreq);

            scData = flip(reshape(gain, [], nFreq)', 1);

            switch filter.type
                case 1
                    for i = 1:nFreq
                        scData(i,:) = scData(i,:) - mean(scData(i,:));
                    end
                case 2
                    for i = 1:nFreq
                        smData = smoothdata(scData(i,:), "gaussian", floor(nField / 5));
                        scData(i,:) = scData(i,:) - smData;
                    end
            end

            size(scData)
            figure
            imagesc(field, freq, scData);
            axis xy;
            set(gca,'CLim',[-0.5 0]);

            % Get points
            crossP = [];
            while true
                x = ginput(1);
                if isempty(x)
                    break;
                end
                [~, fieldValIdx] = min(abs(field - x(1)));
                fieldVal = field(fieldValIdx);
                [~, freqValIdx] = min(abs(freq - x(2)));
                freqVal = freq(freqValIdx);

                crossP = [crossP; fieldVal, freqVal];
                hold on
                plot(fieldVal, freqVal, '.', 'Color', 'white');
                hold off
            end
            obj.crossingPoints = crossP;
        end

        function makecrossingfit(obj, dfR)

            damping = @(H) 2 * obj.dampingParameters.damping.value / obj.kittelParameters.gyromagneticRatio.value ...
                           * H + obj.dampingParameters.inhomogeneousDamping.value;
            switch obj.magnetization.type 
                case 0
                    disp('ERROR No magnetization type found.')
                    return
                case 1
                    kittelEq = @(H) obj.kittelParameters.gyromagneticRatio.value ...
                               .* sqrt((H + obj.kittelParameters.anisotropyField.value) .* ...
                               (H + obj.kittelParameters.anisotropyField.value + ...
                                obj.kittelParameters.effectiveMagnetization.value));
                    deltaFreqFMR = @(H) obj.kittelParameters.gyromagneticRatio.value .* ...
                                   sqrt(1 + (obj.kittelParameters.gyromagneticRatio.value .* obj.kittelParameters.effectiveMagnetization.value ./ ...
                                        (4 .* pi .* kittelEq(H))).^2) .* damping(H);
                case 2
                    kittelEq = @(H) obj.kittelParameters.gyromagneticRatio.value ...
                               .* abs(H - obj.kittelParameters.effectiveMagnetizatio.valuen);
                    deltaFreqFMR = @(H) obj.kittelParameters.gyromagneticRatio.value .* ...
                                    damping(H);
                case 3
                    kittelEq = @(H) obj.kittelParameters.gyromagneticRatio.value ...
                               .* sqrt(1 + obj.kittelParameters.magneticSusceptibility.value) ...
                               .* abs(H);
                    deltaFreqFMR = @(H) obj.kittelParameters.gyromagneticRatio.value .* ...
                                    sqrt(1 + obj.kittelParameters.magneticSusceptibility.value) .* ...
                                    damping(H);
            end
            
            complexFreq = @(H) kittelEq(H) - 1i .* deltaFreqFMR(H);
            complexRess = @(fR, dfR) fR - 1i .* dfR;

            % couplingEq1 = @(param, H) real((complexFreq(H) + complexRess(param(1), param(2))) ./ 2 - ...
            %                           1i .* param(3) + ...
            %                           sqrt(((complexFreq(H) - complexRess(param(1), param(2))) ./ 2).^2 - ...
            %                           (param(3) .* exp(1i * param(4))).^2));
            % couplingEq2 = @(param, H) real((complexFreq(H) + complexRess(param(1), param(2))) ./ 2 - ...
            %                           1i .* param(3) - ...
            %                           sqrt(((complexFreq(H) - complexRess(param(1), param(2))) ./ 2).^2 - ...
            %                           (param(3) .* exp(1i .* param(4))).^2));
            % 
            couplingEq = @(fR, dfR, g, phi, H, f) (real(sqrt(((complexFreq(H) - complexRess(fR, dfR)) ./ 2).^2 - ...
                                                   (g .* exp(1i .* phi)).^2))).^2 - ...
                                                  (f - real((complexFreq(H) + complexRess(fR, dfR)) ./ 2 - 1i .* g)).^2;

            fun = @(param, H) [couplingEq1(param, H), couplingEq2(param, H)];
            fields = obj.crossingPoints(:, 1);
            % disp(size(fields))
            frequencies = obj.crossingPoints(:,2);
            dataToFit = obj.crossingPoints;
            z = zeros(size(fields));

            ft = fittype(@(fR, g, phi, H, f) couplingEq(fR, dfR, g, phi, H, f), 'independent', {'H', 'f'}, 'dependent', 'z');
            opts = fitoptions('Method', 'NonlinearLeastSquares');
            opts.Display = 'Off';
            opts.MaxFunEvals = 5000;
            opts.MaxIter = 5000;
            opts.StartPoint = [20.5 2 1.57];
            opts.Lower = [0, 0, 0];
            opts.Upper = [30, 10, 3.14];
            opts.TolFun = 1e-14;

            [fitresult, gof] = fit(dataToFit, z, ft, opts)
            % options = optimoptions("lsqcurvefit", "MaxIterations", 5000, "MaxFunctionEvaluations", 5000);
            % [x, res] = lsqcurvefit(fun, [20.5 0.6 1 1.57], fields, frequencies, [0, 0.6, 0.7658, 1.57], [30, 0.6, 2, 1.57], options)
            % % [x, res] = lsqcurvefit(couplingEq1, [20 0.5 1 1.57], fields, frequencies(:,1), [0, 0, 0, 0], [30, 5, 10, 3.14])
            % % fR, dfR, g, phi

            fR = fitresult.fR;
            % dfR = 0.6;
            g = fitresult.g;
            phi = fitresult.phi;
            
            hold on;
            fimplicit(@(x, y) couplingEq(fR, dfR, g, phi, x, y), [-1, 0, 15, 30]);
            plot(fields, frequencies, 'ko');
            hold off;
            % plot(fields, frequencies, 'ko', sort(fields), fun([fR, dfR, g, phi], sort(fields)), 'b-')
        end

        function makekittel(obj)
            
            switch obj.magnetization.type 
                case 0
                    disp('ERROR No magnetization type found.')
                    return
                case 1
                    makeinplanekittelfit(obj);
                case 2
                    makeoutofplanekittelfit(obj);
                case 3
                    makeparamagnetickittelfit(obj);
            end            
        end

        function kittelParams = makeinplanekittelfit(obj, startEffMag, ...
                                                     startAniField, ...
                                                     startGyroR)
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
            %       kittelParams = MAKEINPLANEKITTELFIT(obj, startEffMag,
            %                                           startAniField,
            %                                           startGyroR)
            %
            %   Input Arguments
            %       obj - Experiment
            %           Experiment
            %       startEffMag - Starting value of effective magnetization
            %           scalar
            %       startAniField - Starting value of anisotropy field
            %           scalar
            %       startGyroR - Starting value for gyromagnetic ratio
            %           scalar
            %
            %   Output Arguments
            %       kittelParams - Kittel fit parameters
            %           KittelParameters

            if isempty(obj.resonanceParameters)
                disp(['ERROR Resonance parameters not found. ' ...
                     'Use "obj.resonanceanalysis()".']);
                return;
            end

            % Calculate start point value
            field = obj.resonanceParameters(:, Column.resonanceField);
            frequency = obj.resonanceParameters(:,Column.frequency);
            if nargin < 3
                startPoint = [-0.2*sign(field(1)), 0, 28];
            else
                startPoint = [startEffMag, startAniField, startGyroR];
            end
            

            % Calculate fit parameters for Kittel equation
            [fitresult, ~] = fitinplanekittel(field, ...
                                                frequency, ...
                                                startPoint, ...
                                                obj.fieldUnits, ...
                                                obj.frequencyUnits);

            % Calculate confidence intervals
            confInt = confint(fitresult);
            magConfInt = 0.5 * abs(confInt(1,1) - confInt(2,1));
            aniConfInt = 0.5 * abs(confInt(1,2) - confInt(2,2));
            gyroConfInt = 0.5 * abs(confInt(1,3) - confInt(2,3));
            

            % Create physical magnitude instances to store data
            gyromagneticRatio = PhysicalMagnitude(fitresult.gamma, ...
                                                  gyroConfInt, ...
                                                  obj.frequencyUnits / ...
                                                  obj.fieldUnits);
            anisotropyField = PhysicalMagnitude(fitresult.aF, ...
                                                aniConfInt, ...
                                                obj.fieldUnits);
            effectiveMagnetization = PhysicalMagnitude(fitresult.Meff, ...
                                                       magConfInt, ...
                                                       obj.fieldUnits);
            magneticSusceptibility = PhysicalMagnitude;

            % Create Kittel parameters instance
            kittelParams = KittelParameters(gyromagneticRatio, ...
                                            anisotropyField, ...
                                            effectiveMagnetization, ...
                                            magneticSusceptibility);

            % Set this class Kittel parameters to instance
            obj.kittelParameters = kittelParams;
        end

        function kittelParams = makeoutofplanekittelfit(obj, ...
                                                        startEffMag, ...
                                                        startGyroR)
            %MAKEOUTOFPLANEKITTELFIT - Make Kittel's law fit with data
            %   Makes the best Kittel equation parameter fit with the
            %   experiment resonance data.
            %   Kittel out of plane equation:
            %       frequency = ɣ/2π * μ0 |H - Meff|
            %
            %       Meff = Ms + Hani - Hu
            %   where Ms is the saturation magnetization and Hani the
            %   magnetocrystalline anisotropy field and Hu the uniaxial
            %   out-of-plane anisotropy field.
            %
            %   Syntax
            %       kittelParams = MAKEOUTOFPLANEKITTELFIT(obj, 
            %                                              startEffMag,
            %                                              startGyroR)
            %
            %   Input Arguments
            %       obj - Experiment
            %           Experiment
            %       startEffMag - Starting value of effective magnetization
            %           scalar
            %       startGyroR - Starting value for gyromagnetic ratio
            %           scalar
            %
            %   Output Arguments
            %       kittelParams - Kittel fit parameters
            %           KittelParameters

            if isempty(obj.resonanceParameters)
                disp(['ERROR Resonance parameters not found. ' ...
                     'Use "obj.resonanceanalysis()".']);
                return;
            end

            % Calculate start point value
            field = obj.resonanceParameters(:, Column.resonanceField);
            frequency = obj.resonanceParameters(:,Column.frequency);
            if nargin < 2
                startPoint = [0, 28]; 
            else
                startPoint = [startEffMag, startGyroR];
            end

            % Calculate fit parameters for Kittel equation
            [fitresult, ~] = fitoutofplanekittel(field, ...
                                                 frequency, ...
                                                 startPoint, ...
                                                 obj.fieldUnits, ...
                                                 obj.frequencyUnits);

            % Calculate confidence intervals
            confInt = confint(fitresult);
            magConfInt = 0.5 * abs(confInt(1,1) - confInt(2,1));
            gyroConfInt = 0.5 * abs(confInt(1,2) - confInt(2,2));
            

            % Create physical magnitude instances to store data
            gyromagneticRatio = PhysicalMagnitude(fitresult.gamma, ...
                                                  gyroConfInt, ...
                                                  obj.frequencyUnits / ...
                                                  obj.fieldUnits);
            
            effectiveMagnetization = PhysicalMagnitude(fitresult.Meff, ...
                                                       magConfInt, ...
                                                       obj.fieldUnits);

            anisotropyField = PhysicalMagnitude;
            magneticSusceptibility = PhysicalMagnitude;

            % Create Kittel parameters instance
            kittelParams = KittelParameters(gyromagneticRatio, ...
                                            anisotropyField, ...
                                            effectiveMagnetization, ...
                                            magneticSusceptibility);

            % Set this class Kittel parameters to instance
            obj.kittelParameters = kittelParams;
        end

        function kittelParams = makeparamagnetickittelfit(obj, ...
                                                          startSuscept, ...
                                                          startGyroR)
            %MAKEPARAMAGNETICKITTELFIT - Make Kittel's law fit with data
            %   Makes the best Kittel equation parameter fit with the
            %   experiment resonance data. The material is above Curie
            %   temperature and thus, the magnetization is proportional
            %   to the applied field by M = χ * H.
            %
            %   Paramegnetic Kittel equation:
            %       frequency = ɣ/2π * μ0 sqrt(1 + χ) * H
            %
            %   where χ may have a temperature dependence.
            %
            %   Syntax
            %       kittelParams = MAKEPARAMAGNETICKITTELFIT(obj, 
            %                                                startSuscept,
            %                                                startGyroR)
            %
            %   Input Arguments
            %       obj - Experiment
            %           Experiment
            %       startSuscept - Starting value of susceptibility
            %           scalar
            %       startGyroR - Starting value for gyromagnetic ratio
            %           scalar
            %
            %   Output Arguments
            %       kittelParams - Kittel fit parameters
            %           KittelParameters

            if isempty(obj.resonanceParameters)
                disp(['ERROR Resonance parameters not found. ' ...
                     'Use "resonanceanalysis()" method first.']);
                return;
            end

            % Calculate start point value
            field = obj.resonanceParameters(:, Column.resonanceField);
            frequency = obj.resonanceParameters(:,Column.frequency);
            if nargin < 2
                startPoint = [0, 28]; 
            else
                startPoint = [startSuscept, startGyroR];
            end

            % Calculate fit parameters for Kittel equation
            [fitresult, ~] = fitparamagnetickittel(field, ...
                                                   frequency, ...
                                                   startPoint, ...
                                                   obj.fieldUnits, ...
                                                   obj.frequencyUnits);

            % Calculate confidence intervals
            confInt = confint(fitresult);
            susConfInt = 0.5 * abs(confInt(1,1) - confInt(2,1));
            gyroConfInt = 0;
            

            % Create physical magnitude instances to store data
            gyromagneticRatio = PhysicalMagnitude(fitresult.gamma, ...
                                                  gyroConfInt, ...
                                                  obj.frequencyUnits / ...
                                                  obj.fieldUnits);
            
            effectiveMagnetization = PhysicalMagnitude;

            anisotropyField = PhysicalMagnitude;

            magneticSusceptibility = PhysicalMagnitude(fitresult.X, ...
                                                       susConfInt, ...
                                                       Scale.one);

            % Create Kittel parameters instance
            kittelParams = KittelParameters(gyromagneticRatio, ...
                                            anisotropyField, ...
                                            effectiveMagnetization, ...
                                            magneticSusceptibility);

            % Set this class Kittel parameters to instance
            obj.kittelParameters = kittelParams;
        end

        function dampingParams = makedampingfit(obj)
            %MAKEDAMPINGFIT - Make damping fit and get damping parameters
            %
            %   Syntax
            %       dampingParms = MAKEDAMPINGFIT(obj)
            %
            %   Input Arguments
            %       obj - Experiment
            %           Experiment
            %
            %   Output Arguments
            %       dampingParams - Damping parameters
            %           DampingParameters

            % Get fwhm, frequency and gyromagnetic ratio
            fwhm = obj.resonanceParameters(:, Column.fwhm);
            frequency = obj.resonanceParameters(:, Column.frequency);
            if isempty(obj.kittelParameters)
                gamma = 28;
            else
                gamma = obj.kittelParameters.gyromagneticRatio.value;
            end

            % Calculate fit parameters for Kittel equation
            [fitresult, ~] = dampingfit(frequency, ...
                                          fwhm, ...
                                          obj.frequencyUnits, ...
                                          obj.fieldUnits);

            % Calculate confidence intervals
            confInt = confint(fitresult);
            dampingConfInt = 0.5 * abs(confInt(1,1) - confInt(2,1));
            inhDampingConfInt = 0.5 * abs(confInt(1,2) - confInt(2,2));
            

            % Create physical magnitude instances to store data
            damping = PhysicalMagnitude(2 * gamma * fitresult.p1, ...
                                        2 * gamma * dampingConfInt, ...
                                        Scale.one);
            
            inhomogeneousDamping = PhysicalMagnitude(fitresult.p2, ...
                                                     inhDampingConfInt, ...
                                                     obj.fieldUnits);

            % Create damping parameters instance
            dampingParams = DampingParameters(damping, ...
                                              inhomogeneousDamping);

            % Set class damping parameters value to instance
            obj.dampingParameters = dampingParams;
        end

        function makefieldpositive(obj)
            %MAKEFIELDPOSITIVE - Make the field values positive (abs(.))
            %   Takes the absolute value of the field
            %
            %   Syntax
            %       MAKEFIELDPOSITIVE(obj)
            %
            %   Input Arguments
            %       obj - Experiment
            %           Experiment

            obj.data(:, Column.field) = abs(obj.data(:, Column.field));
        end

        function makefieldnegative(obj)
            %MAKEFIELDNEGATIVE - Make the field values negative (-abs(.))
            %   Takes the absolute value of the field and changes sign
            %
            %   Syntax
            %       MAKEFIELDNEGATIVE(obj)
            %
            %   Input Arguments
            %       obj - Experiment
            %           Experiment

            obj.data(:, Column.field) = - abs(obj.data(:, Column.field));
        end

% -----------------------------------------------------------------------%
%                               GETTERS                                  %
% -----------------------------------------------------------------------%

        function D = getdata(obj)
            %GETDATA - Return the data array of the experiment
            %   Output the data array
            %
            %   Syntax
            %       D = GETDATA(obj)
            %
            %   Input parameters
            %       obj - The experiment variable
            %           Experiment
            %
            %   Output parameters
            %       D - Data array
            %           matrix

            D = obj.data;
        end

        function i = getindex(obj)
            %GETINDEX - Return the indexing parameter of the experiment
            %   Output the indexing parameter.
            %
            %   Syntax
            %       i = GETINDEX(obj)
            %
            %   Input parameters
            %       obj - The experiment variable
            %           Experiment
            %
            %   Output parameters
            %       i - Indexing parameter
            %           scalar
            
            i = obj.indexingParameter;
        end

        function n = getnfrequencies(obj)
            %GETNFREQUENCIES - Return the number of frequencies
            %   Return the number of frequencies of the experiment.
            %
            %   Syntax
            %       n = GETNFREQUENCIES(obj)
            %
            %   Input parameters
            %       obj - The experiment variable
            %           Experiment
            %
            %   Output parameters
            %       n - Number of frequencies
            %           scalar
            
            n = obj.nFrequencies;
        end

        function c = getindexingcolumn(obj)
            %GETINDEXINGCOLUMN - Return the indexing column number
            %   Output the indexing column number of the experiment.
            %
            %   Syntax
            %       c = GETINDEXINGCOLUMN(obj)
            %
            %   Input parameters
            %       obj - The experiment variable
            %           Experiment
            %
            %   Output parameters
            %       c - Indexing column number
            %           scalar
            
            c = obj.indexingColumn;
        end

        function u = getfieldunits(obj)
            %GETFIELDUNITS - Get the field unit of measurement
            %
            %   Syntax
            %       u = GETFIELDUNITS(obj)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           Experiment
            %
            %   Output Arguments
            %       u - Unit of measurement name
            %           string

            u = obj.fieldUnits;
        end

        function u = getfrequencyunits(obj)
            %GETFREQUENCYUNITS - Get the frequency unit of measurement
            %
            %   Syntax
            %       u = GETFREQUENCYUNITS(obj)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           Experiment
            %
            %   Output Arguments
            %       u - Unit of measurement
            %           FrequencyUnit | UnitStruct

            u = obj.frequencyUnits;
        end

        function u = getindexingunits(obj)
            %GETINDEXINGUNITS - Get the indexing parameter units
            %
            %   Syntax
            %       u = GETINDEXINGUNITS(obj)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           Experiment
            %
            %   Output Arguments
            %       u - Unit of measurement name
            %           string

            u = obj.indexingUnits;
        end

% -----------------------------------------------------------------------%
%                               SETTERS                                  %
% -----------------------------------------------------------------------%

        function setindex(obj, newIndexingParameter)
            %SETINDEX - Set the indexing parameter of the experiment
            %   Set the indexing parameter of the experiment.
            %
            %   Syntax
            %       SETINDEX(obj, newIndexingParameter)
            %
            %   Input parameters
            %       obj - The experiment variable
            %           Experiment
            %       newIndexingParameter - New indexing parameter
            %           scalar
            
            obj.indexingParameter = newIndexingParameter;
        end
    
        function setindexingcolumn(obj, newIndexingColumn)
            %SETINDEXINGCOLUMN - Set the indexing column of the experiment
            %   Set the indexing column of the experiment.
            %
            %   Syntax
            %       SETINDEXINGCOLUMN(obj, newIndexingColumn)
            %
            %   Input parameters
            %       obj - The experiment variable
            %           Experiment
            %       newIndexingColumn - New indexing column number
            %           scalar
            
            obj.indexingColumn = newIndexingColumn;
        end
    
        function setfieldunits(obj, unit)
            %SETFIELDUNITS - Set field unit of measurement
            %
            %   Syntax
            %       SETFIELDUNITS(obj, unit)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           Experiment
            %       unit - Magnetic unit
            %           MagneticUnit | UnitStruct
            %
            %   Example
            %       SETFIELDUNITS(obj, MagneticUnit.gauss, Scale.milli)
            %       Set the field units to mG (milligauss).

            obj.fieldUnits = unit;
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

            obj.frequencyUnits = unit;
        end

        function setindexingunits(obj, unitName)
            %SETINDEXINGUNITS - Set indexing parameter unit of measurement
            %
            %   Syntax
            %       SETINDEXINGUNITS(obj, unitName)
            %
            %   Input Arguments
            %       obj - Experiment class
            %           Experiment
            %       unitName - New unit of measurement name
            %           string

            if ~ischar(unitName)
                disp('ERROR Input is not a character array')
                return
            end

            obj.indexingUnits = unitName;
        end
    
        function setmagnetization(obj, magType)
            obj.magnetization = magType;
        end
    end
end


% -----------------------------------------------------------------------%
%                         AUXILIARY FUNCTIONS                            %
% -----------------------------------------------------------------------%

function outData = sortcolumns(data, fileParams)
    %SORTCOLUMNS - Sort the columns of the input data
    %   Sort the columns of the input data so that they go
    %   in the following order: frequency, field, gain and
    %   (if exists) indexing parameter.
    %
    %   Syntax
    %       outData = SORTCOLUMNS(data, fileParams)
    %
    %   Input Arguments
    %       data - Data file containing the data to be sorted
    %           matrix
    %       fileParams - File parameters
    %           FileParameters
    %
    %   Output Arguments
    %       outData - Sorted data
    %           matrix

    outData = zeros(size(data,1),4);

    outData(:,1) = data(:, fileParams.getfrequency());
    outData(:,2) = data(:, fileParams.getfield());
    outData(:,3) = data(:, fileParams.getgain());

    % If indexing parameter column exists, append it
    if fileParams.getindexing() > 0
        outData(:,4) = data(:,4);
    end
end

function nFrequencies = countfrequencies(data)
    %COUNTFREQUENCIES - Count the number of different frequency values
    %
    %   Syntax
    %       nFrequencies = COUNTFREQUENCIES(data)
    %
    %   Input Arguments
    %       data - Completely sorted data array
    %           matrix
    %
    %   Output Arguments
    %       nFrequencies - Number of differenct frequency values
    %           scalar
    nFrequencies = length(unique(data(:,1)));
end

function croppedData = cropdata(data, fieldUnits, frequencyUnits)
    %CROPDATA - Crops the data to the selected interval
    %   Through the input of the user, this function crops the
    %   data to yield the selected portion of the data.
    %
    %   Syntax
    %       croppedData = CROPDATA(data)
    %
    %   Input Arguments
    %       data - Data to crop
    %           matrix
    %       fieldUnits - Units of magnetic field
    %           MagneticUnit | UnitStruct
    %       frequencyUnits - Units of frequency
    %           FrequencyUnit | UnitStruct
    % 
    %   Output Arguments
    %       croppedData - Cropped data
    %           matrix

    croppedData = []; % Default return value

    frequencyValue = data(1,1);
    fieldData = data(:,2);
    gainData = data(:,3);

    smoothData = smoothdata(gainData, 'rloess', 500);
    % Plot data to ask input from user
    figure('Name', 'Select a field interval to crop data (2 clicks)');
    plot(fieldData, gainData - smoothData, '.');%, fieldData, smoothData, 'r-');
    title(['Frequency = ',num2str(frequencyValue),' ',frequencyUnits.tag]);
    xlabel('Magnetic field (', fieldUnits.tag,')', 'interpreter', 'latex');
    ylabel('Gain (dB)', 'interpreter', 'latex')
    box on;
    grid on;

    % Get input
    try [xlims, ~] = ginput(2); catch, return; end
    close(gcf);

    if length(xlims) < 2, return; end

    % Find closest data indices to limits
    [~, idx1] = min(abs(fieldData - xlims(1)));
    [~, idx2] = min(abs(fieldData - xlims(2)));
    
    minIdx = min(idx1, idx2);
    maxIdx = max(idx1, idx2);

    % Crop data
    croppedData = data(minIdx:maxIdx, :);
end

function startPoint = getstartpoint(field, gain)
    %GETSTARTPOINT - Calculate the optimal starting values for fit
    %   Calculate the optimal starting point for the fit.
    %
    %   Syntax
    %       startPoint = GETSTARTPOINT(field, gain)
    %
    %   Input Arguments
    %       field - Magnetic field data (x variable)
    %           vector
    %       gain - Gain data (y variable)
    %           vector

    % Find the minimum gain value and its index
    [minGain, minGainIdx] = min(gain);
    gainMean = 0.5 * (gain(1) + gain(end));

    % Width of peak FWHM (Ansatz: approx. 1/2 of the interval width)
    fwhm = 0.5 * abs(field(end) - field(1));

    % Resonant H field (peak minimum)
    hr = field(minGainIdx);

    % Symmetric part (approx. area of peak A = 1/2 * width * height)
    A = 0.5 * fwhm * ((gain(1) + gain(end))/2 - minGain);

    % Antisymmertric part (approx. twice area of peak)
    B = 2 * A;

    % Slope of the surrounding of the peak
    k = (gain(1) - gain(end)) / (field(1) - field(end));

    % Line equation: k * hr + offset = mean(gain)
    offset = gainMean - k * hr;

    startPoint = [A, B, fwhm, hr, k, offset];
end

function [fitresult, gof] = fitlorentziantodata(xData, yData, startPoint)
    %FITLORENTIZANTODATA - Fit Lorentzian curve to data
    %   Fit Lorentzian curve to data with given start point.
    %
    %   Syntax
    %   [fitresult, gof] = FITLORENTIZANTODATA(xData, yData, startPoint)
    %
    %   Input Arguments
    %       xData - x axis data
    %           vector
    %       yData - y axis data
    %           vector
    %       startPoint - Start point parameters for fit
    %           vector
    %   Output Arguments
    %       fitresult - Fit result
    %           cfit
    %       gof - Goodness-of-fit statistics
    %           gof structure   


    % Set up fittype and options.
    ft = fittype(['A / (4 * (x-hr)^2 + fwhm^2) - ' ...
                  '(fwhm * B * (x-hr)) / (4 * (x-hr)^2 + fwhm^2) + ' ...
                  'offset + k * x'], 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.StartPoint = startPoint;
    opts.TolFun = 1e-14;

    % Fit model to data
    [fitresult, gof] = fit(xData, yData, ft, opts);

    % Plot fit with data.
    figure( 'Name', 'Lorentzian fit over data' );
    plot(fitresult, xData, yData);
    title(['r^2 = ', num2str(gof.rsquare)]);
    close(gcf);
end

function [fitresult, gof] = fitinplanekittel(xData, yData, startPoint, ...
                                             fieldU, freqU)
    %FITINPLANEKITTEL - Fit in plane Kittel equation to data
    %   Fit in plane Kittel equation to data with given start point.
    %
    %   Syntax
    %   [fitresult, gof] = FITINPLANEKITTEL(xData, yData, startPoint,
    %                                       fieldU, freqU)
    %
    %   Input Arguments
    %       xData - x axis data
    %           vector
    %       yData - y axis data
    %           vector
    %       startPoint - Start point parameters for fit
    %           vector
    %       fieldU - Magnetic field units
    %           MagneticUnit | UnitStruct
    %       freqU - Frequency units
    %           FrequencyUnit | UnitStruct
    %
    %   Output Arguments
    %       fitresult - Fit result
    %           cfit
    %       gof - Goodness-of-fit statistics
    %           gof structure   

    % Set up fittype and options.
    ft = fittype( 'gamma*sqrt((x+aF)*(x+aF+Meff))', 'independent', 'x', ...
                  'dependent', 'y' ); % in plane
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = startPoint;
    opts.Display = 'Off';
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.TolFun = 1e-14;

    % Fit model to data
    [fitresult, gof] = fit(xData, yData, ft, opts)

    % Plot fit with data.
    figure( 'Name', 'Kittel fit' );
    plot(fitresult, xData, yData);
    xlabel(['$\mu_0 H$ (', fieldU.tag, ')'], 'Interpreter', 'latex');
    ylabel(['Frequency (', freqU.tag, ')'], 'Interpreter', 'latex');
    title(['r^2 = ', num2str(gof.rsquare)]);
end

function [fitresult, gof] = fitoutofplanekittel(xData, yData, ...
                                                startPoint, ...
                                                fieldU, freqU)
    %FITOUTOFPLANEKITTEL - Fit out of plane Kittel equation to data
    %   Fit out of plane Kittel equation to data with given start point.
    %
    %   Syntax
    %   [fitresult, gof] = FITOUTOFPLANEKITTEL(xData, yData, startPoint,
    %                                          fieldU, freqU)
    %
    %   Input Arguments
    %       xData - x axis data
    %           vector
    %       yData - y axis data
    %           vector
    %       startPoint - Start point parameters for fit
    %           vector
    %       fieldU - Magnetic field units
    %           MagneticUnit | UnitStruct
    %       freqU - Frequency units
    %           FrequencyUnit | UnitStruct
    %
    %   Output Arguments
    %       fitresult - Fit result
    %           cfit
    %       gof - Goodness-of-fit statistics
    %           gof structure   

    % Set up fittype and options.
    ft = fittype( 'gamma * abs(x - Meff)', 'independent', 'x', ...
                  'dependent', 'y' ); % in plane
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = startPoint;
    opts.Display = 'Off';
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.TolFun = 1e-14;

    % Fit model to data
    [fitresult, gof] = fit(xData, yData, ft, opts);

    % Plot fit with data.
    figure( 'Name', 'Kittel fit' );
    plot(fitresult, xData, yData);
    xlabel(['$\mu_0 H$ (', fieldU.tag, ')'], 'Interpreter', 'latex');
    ylabel(['Frequency (', freqU.tag, ')'], 'Interpreter', 'latex');
    title(['r^2 = ', num2str(gof.rsquare)]);
end

function [fitresult, gof] = fitparamagnetickittel(xData, yData, ...
                                                  startPoint, ...
                                                  fieldU, freqU)
    %FITPARAMAGNETICKITTEL - Fit out of plane Kittel equation to data
    %   Fit out of plane Kittel equation to data with given start point.
    %
    %   Syntax
    %   [fitresult, gof] = FITPARAMAGNETICKITTEL(xData, yData, startPoint,
    %                                            fieldU, freqU)
    %
    %   Input Arguments
    %       xData - x axis data
    %           vector
    %       yData - y axis data
    %           vector
    %       startPoint - Start point parameters for fit
    %           vector
    %       fieldU - Magnetic field units
    %           MagneticUnit | UnitStruct
    %       freqU - Frequency units
    %           FrequencyUnit | UnitStruct
    %
    %   Output Arguments
    %       fitresult - Fit result
    %           cfit
    %       gof - Goodness-of-fit statistics
    %           gof structure   

    % Set up fittype and options.
    ft = fittype( 'gamma * sqrt(1 + X) * abs(x)', 'independent', 'x', ...
                  'dependent', 'y' ); % in plane
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = startPoint;
    opts.Lower = [-Inf, startPoint(2)];
    opts.Upper = [ Inf, startPoint(2)];
    opts.Display = 'Off';
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.TolFun = 1e-14;

    % Fit model to data
    [fitresult, gof] = fit(xData, yData, ft, opts);

    % Plot fit with data.
    figure( 'Name', 'Kittel fit' );
    plot(fitresult, xData, yData);
    xlabel(['$\mu_0 H$ (', fieldU.tag, ')'], 'Interpreter', 'latex');
    ylabel(['Frequency (', freqU.tag, ')'], 'Interpreter', 'latex');
    title(['r^2 = ', num2str(gof.rsquare)]);
end

function [fitresult, gof] = dampingfit(xData, yData, freqU, fieldU)
    %DAMPINGFIT - Fit damping equation to data
    %   Fit line equation to data to find damping parameters.
    %
    %   Syntax
    %   [fitresult, gof] = DAMPINGFIT(xData, yData, freqU, fieldU)
    %
    %   Input Arguments
    %       xData - x axis data
    %           vector
    %       yData - y axis data
    %           vector
    %       freqU - Frequency units
    %           FrequencyUnit | UnitStruct
    %       fieldU - Magnetic field units
    %           MagneticUnit | UnitStruct
    %
    %   Output Arguments
    %       fitresult - Fit result
    %           cfit
    %       gof - Goodness-of-fit statistics
    %           gof structure

    % Set up fittype and options.
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [0 0];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'Kittel fit' );
    plot(fitresult, xData, yData);
    xlabel(['Frequency (', freqU.tag, ')'], 'Interpreter', 'latex');
    ylabel(['$\mu_0\Delta H$ (', fieldU.tag, ')'], 'Interpreter', 'latex');
    title(['r^2 = ', num2str(gof.rsquare)]);
end

function [hrError, fwhmError] = getconfint(fitResult)
    %GETCONFINT - Calculate the confidence of the ressonant field and FWHM
    %   Calculate half the confidence intervals length of
    %   ressonant field and FWHM
    %
    %   Syntax
    %       [hrError, fwhmError] = GETCONFINT(gof)
    %
    %   Input Arguments
    %       gof - Goodness-of-fit statistics
    %           gof structure
    %
    %   Output Arguments
    %       hrError - Half confidence interval length for ressonant field
    %           scalar
    %       fwhmError - Half confidence interval length for FWHM
    %           scalar

    confInt = confint(fitResult);

    % Confidence interval limits for fwhm are located in the
    % 3rd column. The 4th column is for the resonance field
    fwhmError = 0.5 * abs(confInt(1, 3) - confInt(2, 3));
    hrError = 0.5 * abs(confInt(1, 4) - confInt(2, 4));
end

function resonanceParams = findresonance(data, fieldUnits, frequencyUnits)
    %FINDRESONANCE - Find the resonance parameters of data
    %   Finds the resonance parameters (ressonant field,
    %   field FWHM, symmetry and antisymmetry parameters
    %   among others) asking the user input.
    %
    %   The resonance parameter output vector is structured
    %   as follows:
    %       [ressonant field, fwhm, ...
    %        error in ressonant field, error in fwhm]
    %
    %   Syntax
    %       resonanceParams = FINDRESONANCE(data)
    %
    %   Input Arguments
    %       data - Data of A SINGLE FREQUENCY value
    %           matrix
    %       fieldUnits - Units of magnetic field
    %           MagneticUnit | UnitStruct
    %       frequencyUnits - Scale of frequency
    %           FrequencyUnit | UnitStruct
    %
    %   Output Arguments
    %       resonanceParams - Resonance parameters
    %           vector

    % Default return value
    resonanceParams = [];

    isBadFit = 1;
    while isBadFit
        % Crop data to make a better fit. The crop must be close
        % to the peak region
        cropData = cropdata(data, fieldUnits, frequencyUnits);
        if isempty(cropData), return; end
    
        % Copy data to variables
        fieldData = cropData(:, 2);
        gainData = cropData(:, 3);
    
        % Calculate starting point values for fit
        startPoint = getstartpoint(fieldData, gainData);
    
        % Fit model to data
        [fitResult, gof] = fitlorentziantodata(fieldData, ...
                                               gainData, ...
                                               startPoint);

        % Check if the fit is bad
        isBadFit = (gof.rsquare < 0.85);
    end
    
    % If the fit is good enough, save the results
    [hrError, fwhmError] = getconfint(fitResult);
    resonanceParams = [fitResult.hr, fitResult.fwhm, ...
                        hrError, fwhmError];   
end

function monoFreqData = filterdatabyfreq(data, freq)
    %FILTERDATABYFREQ - Filter data by frequency value
    %   Return the data rows with the given frequency value
    %
    %   Syntax
    %       monoFreqData = FILTERDATABYFREQ(data, freq)
    %
    %   Input Arguments
    %       data - Data to be filtered
    %           matrix
    %       freq - Filter frequency
    %           scalar
    %
    %   Output Arguments
    %       monoFreqData - Data filtered by frequency
    %           matrix

    freqData = data(:, 1);
    monoFreqData = data(freqData == freq,:);
end