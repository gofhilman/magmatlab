% Set the column number: Frequency, Field, S12, Index
fp = FileParameters(3, 2, 4, 0);

% Read experiment from file
FGT_exp = ExperimentCollection("FGT300K.txt", fp);

% Set units of the data in file
FGT_exp.setfieldunits(MagneticUnit.gauss);
FGT_exp.setfrequencyunits(FrequencyUnit.hertz);

% Convert units
FGT_exp.convertfieldunits(MagneticUnit.tesla);
FGT_exp.convertfrequencyunits(Scale.giga * FrequencyUnit.hertz);

% Filter by positive or negative field
%   If you want the positive field, comment the 
%   first line and uncomment the second line
FGT_exp.filternegativefield();
% FGT_exp.filterpositivefield();

%% Plot S12 vs Field for each frequency and select peaks
% To go to the next frequency select a range for the peak,
% press ENTER or close the window
FGT_exp.resonanceanalysis(1);

%% Make kittel plot
% Magnetization can be set to inPlane, outOfPlane, paramagnetic
FGT_exp.experimentArray(1).setmagnetization(Magnetization.inPlane);
FGT_exp.experimentArray(1).makekittel();

%% Make 2D plot with normalization filter (Marc's approach)
% When the plot pops up, press ENTER to be able to
% zoom in, and modify the plot options
FGT_exp.experimentArray(1).getcrossingpoints(Filter.normalization)

%% Make 2D plot with noise cancelling filter (Pablo's approach)
% When the plot pops up, press ENTER to be able to
% zoom in, and modify the plot options
% You may need to change the Colormap Limits to see properly.
% To do that, go Edit -> Colormap... -> Set Colormap Limits
% And set CLim Minimum to -0.01
FGT_exp.experimentArray(1).getcrossingpoints(Filter.noiseCancelling)