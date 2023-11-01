function autoPlot(setup, s, f, t)

% Validate number of input arguments
narginchk(3,4)

% 
if setup.PlotLayout == "stacked" && setup.Averaging == "none"
   setup.PlotLayout = "tiled"
end

% Determine spectral unit and set y-axis label

% Plot data
switch setup.PlotLayout
    case "stacked"
	    fig = figure;
		ax  = initAxes(fig, sLabel);
		legend(ax);
		
	    for col = 1:size(s,2)
		    plot(ax, f, s(:, col), DisplayName = "Signal " + col);
		end
		
	  case "tiled"
        % Try to use tiledlayout
        try
            tiledlayout("flow")
        % Use subplot if tiledlayout does not exist (R2018b or older)
        catch
	
	   end
   
    otherwise
	    for col = 1:size(s,2)
		    fig = figure;
		    			
			if setup.Averaging == "none"
			    ax  = initAxes(fid, "Time (s)");
			    surface(ax, f, t, transpose(s));
			else
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



% Initialize figure
f = figure;

for 
