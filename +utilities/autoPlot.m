function autoPlot(setup, s, f, t)
%AUTOPLOT plots a spectrum or spectra computed by the autofft function.
%
% Copyright (c) 2023, Luboš Smolík
% v1.0.0beta (build 21. 11. 2023)
%
% This code is published under BSD-3-Clause License.
%
% autoPlot(setup, s, f)
% autoPlot(setup, s, f, t)
%
% autoPlot(setup, s, f) plots a spectrum or spectra stored in s assuming
%   frequencies stored in a vector f. Visualisation parameters are
%   determined by PlotLayout and EngineeringUnit fields in a structure 
%   array setup.
%   
% autoPlot(setup, s, f, t) plots a spectrogram or spectrograms stored in s
%   assuming frequencies stored in a vector f and time stamps stored in a
%   vector t. Visualisation parameters are determined by PlotLayout and
%   EngineeringUnit fields in a structure array setup.
%
% --- TO DO: Detailed explanation of field values ---

% Validate number of input arguments
narginchk(3,4)

%% Validate inputs
% Check if a "EngineeringUnit" field exists and if it is nonempty
if ~isfield(setup, "EngineeringUnit") || strtrim(setup.EngineeringUnit) == ""
    setup.EngineeringUnit = "EU";
end

% Check if a "Averaging" field exists and if it is nonempty
if ~isfield(setup, "Averaging") || strtrim(setup.Averaging) == ""
    setup.Averaging = "linear";
end

% Check if a "PlotLayout" field exists and if it is valid
if ~isfield(setup, "PlotLayout") || ...
   (setup.PlotLayout == "stacked" && setup.Averaging == "none")
    setup.PlotLayout = "tiled";
elseif setup.PlotLayout == "none"
    % Return control to the invoking script or function
    return;
end

% Check if a "SpectralUnit" field exists and if it is nonempty
if ~isfield(setup, "SpectralUnit") || strtrim(setup.SpectralUnit) == ""
    setup.SpectralUnit = "power";
end

% Check if a "dbReference" field exists and if it is nonempty
if ~isfield(setup, "dbReference") || isempty(setup.dbReference)
    setup.dbReference = NaN;
end

%% Determine spectral unit and set y-axis label
if isnan(setup.dbReference)
    dblab = "";
else
    dblab = " dB/" + num2str(setup.dbReference, "%.2e");
end

switch lower(setup.SpectralUnit)
    case "rms"                          % RMS magnitude
        unitlab = "Autospectrum" + dblab + " (" + setup.EngineeringUnit ...
                   + ", RMS)";
    case {"pk", "0-pk", "peak"}         % 0-peak magnitude
        unitlab = "Autospectrum" + dblab + " (" + setup.EngineeringUnit ...
                   + ", 0-Pk)";
    case {"pp", "pk-pk", "peak2peak"}   % Peak-peak magnitude
        unitlab = "Autospectrum" + dblab + " (" + setup.EngineeringUnit ...
                   + ", Pk-Pk)";
    case {"asd", "psd"}                 % Power spectral density
        unitlab = "Power Spectral Density" + dblab + " (" ...
                   + setup.EngineeringUnit + ")^2 / Hz";
    case {"rsd", "rmssd"}               % Root mean square spectral density
        unitlab = "Power Spectral Density" + dblab + " (" ...
                   + setup.EngineeringUnit + ") / (Hz^1/2)";
    otherwise                           % Autospectrum / power spectrum
        unitlab = "Autospectrum" + dblab + " (" + setup.EngineeringUnit ...
                   + ")^2";
end

%% Plot data
switch setup.PlotLayout
    % Plots stacked in one axes
    case "stacked"
        % Initialise figure
        fig = figure("Color", "w");
        ax  = axes(fig);
	    setAxes(ax, unitlab);
        legend(ax);
	    
        % Cycle through spectra
        for col = 1:size(s, ndims(s))
            plot(ax, f, s(:, col), "DisplayName", "signal " + col);
        end

    % Graphs tiled in one figure
    case "tiled"
        % Initialise figure
        fig = figure("color", "w");

        % Try to use a tiled layout (R2019a or newer)
        try
            tiledlayout(fig, "flow")

            % Cycle through results
            for col = 1:size(s, ndims(s))
                % Initialise new tile
                 ax = nexttile;

                % Plot spectrograms as surfaces or spectra as plots
		        if setup.Averaging == "none"
                    setAxes(ax, "Time (s)");
                    surface(f, t, squeeze(s(:, :, col)).', ...
                            "DisplayName", "signal " + col, ...
                            "EdgeColor","none", "FaceColor", "interp");
                     cbar = colorbar; 
                    title(cbar, unitlab);
                else
                    setAxes(ax, unitlab);
	                plot(f, s(:, col), "DisplayName", "signal " + col);
                end

                legend(ax);
            end
            
        % Use subplot if tiledlayout does not exist (R2018b or older)
        catch
            % Determine the number of tile rows and columns
            rows = 1;
            cols = size(s, ndims(s));

            while rows < cols
                rows = rows + 1;
                cols = ceil(size(s, ndims(s)) / rows);
            end

            % Cycle through results
            for col = 1:size(s, ndims(s))
                % Initialise new subplot
                ax = subplot(rows, cols, col);

                % Plot spectrograms as surfaces or spectra as plots
		        if setup.Averaging == "none"
                    setAxes(ax, "Time (s)");
                    surface(f, t, squeeze(s(:, :, col)).', ...
                            "DisplayName", "signal " + col, ...
                            "EdgeColor","none", "FaceColor", "interp");
                    cbar = colorbar; 
                    title(cbar, unitlab);
                else
                    setAxes(ax, unitlab);
	                plot(f, s(:, col), "DisplayName", "signal " + col);
                end

                legend(ax);
            end
        end

    % Graphs displayed in separate figures
    case "separated"
        % Cycle through results
        for col = 1:size(s, ndims(s))
            % Initialise new figure
            figure("Color", "w");

            % Plot spectrograms as surfaces or spectra as plots
		    if setup.Averaging == "none"
		        surface(f, t, squeeze(s(:, :, col)).', ...
                        "DisplayName", "signal " + col, ...
                        "EdgeColor","none", "FaceColor", "interp");
                setAxes(gca, "Time (s)");
                cbar = colorbar; 
                title(cbar, unitlab);
            else
	            plot(f, s(:, col), "DisplayName", "signal " + col);
                setAxes(gca, unitlab);
		    end
		    
            % Set axes and show legend
	        legend(gca);
        end	
end

%% Local functions
function setAxes(target, str)
% setAxes(target, label) sets properties of axes stored in target and uses
%   str to label y-axis.

    target.Box = "on";
    target.FontName = "Times New Roman";
    target.Layer = "top";
    target.NextPlot = "add";
    target.XLabel.String = "Frequency (Hz)";
    target.YLabel.String = str;
end

% End of main function
end