function autoPlot(setup, s, f, t)
%AUTOPLOT plots a spectrum or spectra computed by the autofft function.
%
% Copyright (c) 2023-2025, Luboš Smolík
% v1.0.1 (build 7. 8. 2025)
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

% CHANGELOG
% v1.0.1 - Fixed bug when the tiled layout was selected by the user for
%          time-frequency analysis results from only one channel.
%        - Fixed figure background color in Matlab R2025a dark mode.
%        - Minor code refactoring.

%% Validate input arguments
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
    labdb = "";
else
    labdb = " dB/" + num2str(setup.dbReference, "%.2e");
end

switch lower(setup.SpectralUnit)
    case "rms"                          % RMS magnitude
        labunit = "Autospectrum" + labdb + " (" + setup.EngineeringUnit ...
                   + ", RMS)";
    case {"pk", "0-pk", "peak"}         % 0-peak magnitude
        labunit = "Autospectrum" + labdb + " (" + setup.EngineeringUnit ...
                   + ", 0-Pk)";
    case {"pp", "pk-pk", "peak2peak"}   % Peak-peak magnitude
        labunit = "Autospectrum" + labdb + " (" + setup.EngineeringUnit ...
                   + ", Pk-Pk)";
    case {"asd", "psd"}                 % Power spectral density
        labunit = "Power Spectral Density" + labdb + " (" ...
                   + setup.EngineeringUnit + ")^2 / Hz";
    case {"rsd", "rmssd"}               % Root mean square spectral density
        labunit = "Power Spectral Density" + labdb + " (" ...
                   + setup.EngineeringUnit + ") / (Hz^1/2)";
    otherwise                           % Autospectrum / power spectrum
        labunit = "Autospectrum" + labdb + " (" + setup.EngineeringUnit ...
                   + ")^2";
end

%% Plot data
% Determine which dimension correspond to the channel
if setup.Averaging == "none"
    nchannel = 3;
else
    nchannel = 2;
end

% Prepare layout and plot data
switch setup.PlotLayout
    % Plots stacked in one axes
    case "stacked"
        % Initialise figure
        fig = figure;
        ax  = axes(fig);
	    setAxes(ax, [f(1), f(end)], labunit);
	    
        % Cycle through channels
        for ch = 1:size(s, nchannel)
            plot(ax, f, s(:, ch), "DisplayName", "signal " + ch);
        end

    % Graphs tiled in one figure
    case "tiled"
        % Initialise figure
        fig = figure;

        % Try to use a tiled layout 
        try
            % Use tiledlayout (R2019a or newer)
            tiledlayout(fig, "flow")
            tiledlayoutexists = true;
        catch
            % Use subplot if tiledlayout does not exist (R2018b or older)
            tiledlayoutexists = false;

            % Determine the number of tile rows and columns
            rows = 1;
            cols = size(s, nchannel);

            while rows < cols
                rows = rows + 1;
                cols = ceil(size(s, nchannel) / rows);
            end
        end

        % Cycle through channels
        for ch = 1:size(s, nchannel)
            % Initialise new tile or subplot
            if tiledlayoutexists == true
                ax = nexttile;
            else
                ax = subplot(rows, cols, ch);
            end
            
            % Plot spectrograms as surfaces or spectra as plots
            if setup.Averaging == "none"
                setAxes(ax, [f(1), f(end)], "Time (s)");
                surface(f, t, squeeze(s(:, :, ch)).', ...
                        "DisplayName", "signal " + ch, ...
                        "EdgeColor","none", "FaceColor", "interp");
                cbar = colorbar; 
                title(cbar, labunit);
            else
                setAxes(ax, [f(1), f(end)], labunit);
                plot(f, s(:, ch), "DisplayName", "signal " + ch);
            end
        end

    % Graphs displayed in separate figures
    case "separated"
        % Cycle through channels
        for ch = 1:size(s, nchannel)
            % Initialise new figure
            fig = figure;
            ax = axes(fig);

            % Plot spectrograms as surfaces or spectra as plots
		    if setup.Averaging == "none"
		        surface(ax, f, t, squeeze(s(:, :, ch)).', ...
                        "DisplayName", "signal " + ch, ...
                        "EdgeColor","none", "FaceColor", "interp");
                setAxes(ax, [f(1), f(end)], "Time (s)");
                cbar = colorbar; 
                title(cbar, labunit);
            else
	            plot(ax, f, s(:, ch), "DisplayName", "signal " + ch);
                setAxes(ax, [f(1), f(end)], labunit);
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