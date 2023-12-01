function autoPlot(setup, s, f, t)
%AUTOPLOT plots a spectrum or spectra computed by the autofft function.
%
% Copyright (c) 2023, Luboš Smolík
% v1.0.0 (build 1. 12. 2023)
%
% This code is published under BSD-3-Clause License.
%
% autoPlot(setup, s, f)
% autoPlot(setup, s, f, t)
%
% autoPlot(setup, s, f) plots a spectrum or spectra stored in s assuming
%   frequencies stored in a vector f. Visualisation parameters are
%   determined by 'PlotLayout' and 'EngineeringUnit' fields in a structure 
%   array setup.
%   
% autoPlot(setup, s, f, t) plots a spectrogram or spectrograms stored in s
%   assuming frequencies stored in a vector f and time stamps stored in a
%   vector t. Visualisation parameters are determined by PlotLayout and
%   EngineeringUnit fields in a structure array setup.
%
%
% Construction of setup:
%	setup = struct('param', 'value', ...);
%
% For a detailed list of valid name-value pair arguments for setup, see
% documentation of autofft. Below are two name-value pair arguments which
% are important for this function:
%
%   - 'EngineeringUnit' - [ character {'EU'} | string ]
%     Specifies a unit of measure of xs. This parameter is used to generate
%     labels for data visualisation.
%
%   - 'PlotLayout' - [ 'none' | 'separated' | 'stacked' | 'tiled' ]
%      Specifies, which layout is used to visualise results.
%      - 'none'      - Does not visualise any results. It is a {default}
%                      option if auttofft is called with no output.
%      - 'separated' - Visualise each spectrum or spectrogram in a separate
%                      figure. It is a {default} option if auttofft is
%                      called with at least one output variable.
%      - 'stacked'   - Plots spectra into single axes. This value cannot be
%                      used to visualise spectrograms.
%      - 'tiled'     - Creates a panel for each spectrum or spectrogram in
%                      a single figure.
%

%% Validate input variables
%  Validate number of input arguments
narginchk(3,4)

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
	    setAxes(ax, [f(1), f(end)], unitlab);
	    
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
                    setAxes(ax, [f(1), f(end)], "Time (s)");
                    surface(f, t, squeeze(s(:, :, col)).', ...
                            "DisplayName", "signal " + col, ...
                            "EdgeColor","none", "FaceColor", "interp");
                     cbar = colorbar; 
                    title(cbar, unitlab);
                else
                    setAxes(ax, [f(1), f(end)], unitlab);
	                plot(f, s(:, col), "DisplayName", "signal " + col);
                end
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
                    setAxes(ax, [f(1), f(end)], "Time (s)");
                    surface(f, t, squeeze(s(:, :, col)).', ...
                            "DisplayName", "signal " + col, ...
                            "EdgeColor","none", "FaceColor", "interp");
                    cbar = colorbar; 
                    title(cbar, unitlab);
                else
                    setAxes(ax, [f(1), f(end)], unitlab);
	                plot(f, s(:, col), "DisplayName", "signal " + col);
                end
            end
        end

    % Graphs displayed in separate figures
    case "separated"
        % Determine dimension for the loop indexing
        if setup.Averaging == "none"
            dim = 3;
        else
            dim = 2;
        end

        % Cycle through results
        for col = 1:size(s, dim)
            % Initialise new figure
            figure("Color", "w");

            % Plot spectrograms as surfaces or spectra as plots
		    if setup.Averaging == "none"
		        surface(f, t, squeeze(s(:, :, col)).', ...
                        "DisplayName", "signal " + col, ...
                        "EdgeColor","none", "FaceColor", "interp");
                setAxes(gca, [f(1), f(end)], "Time (s)");
                cbar = colorbar; 
                title(cbar, unitlab);
            else
	            plot(f, s(:, col), "DisplayName", "signal " + col);
                setAxes(gca, [f(1), f(end)], unitlab);
		    end
        end	
end

%% Local functions
function setAxes(target, xlims, str)
% setAxes(target, label) sets properties of axes stored in target and uses
%   str to label y-axis.
    target.Box = "on";
    target.FontName = "Times New Roman";
    target.Layer = "top";
    target.NextPlot = "add";
    target.XLabel.String = "Frequency (Hz)";
    target.YLabel.String = str;
    target.XLim = xlims;
    legend(target);
end

% End of main function
end