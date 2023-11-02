function autoPlot(setup, s, f, t)

    % Validate number of input arguments
    narginchk(3,4)
    
    % Change plot layout from stacked to tiled for spectrograms
    if setup.PlotLayout == "stacked" && setup.Averaging == "none"
       setup.PlotLayout = "tiled";
    end
    
    % Determine spectral unit and set y-axis label

    % Plot data
    switch setup.PlotLayout
        case "stacked"
            % Stacked plots
	        fig = figure;
		    ax  = initAxes(fig, sLabel);
		    legend(ax);
		    
            for col = 1:size(s,2)
                plot(ax, f, s(:, col), DisplayName = "Signal " + col);
            end
		    
        case "tiled"
            % Tiled plots
            fig = figure;
    
            % Use a tiled layout (R2019a or newer)
            try
                tiledlayout(fig, "flow")
    
            % Use subplot if tiledlayout does not exist (R2018b or older)
            catch
    
            end
       
        otherwise
            for col = 1:size(s,2)
                % Initialize figure
                fig = figure;
                
                % Plot spectra as plots and spectrograms as surfaces
			    if setup.Averaging == "none"
			        ax = initAxes(fig, "Time (s)");
			        surface(ax, f, t, squeeze(s(:, :, col)), ...
                            DisplayName = "Signal " + col, ...
                            EdgeColor = "none", FaceColor = "interp");
                    cbar = colorbar; 
                    title(cbar, slabel);
                else
                    ax = initAxes(fig, slabel);
		            plot(ax, f, s(:, col), DisplayName = "Signal " + col);
			    end
			    
		        legend(ax);
            end	
    end

    function ax = initAxes(target, label)
	    ax = axes(target);
	    ax.Box = "on";
	    ax.FontName = "Times New Roman";
        ax.Layer = "top";
	    ax.NextPlot = "add";
	    ax.XLabel.String = "Frequency (Hz)";
        ax.YLabel.String = label;
    end

% End of main function
end
