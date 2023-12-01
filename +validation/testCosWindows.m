% Test functions
tFun = {@utilities.autoBlackmanHarris, @utilities.autoFlatTop, ...
      @utilities.autoHann, @utilities.autoHamming, @utilities.autoUniform};

% Reference functions - requires the Signal Processin Toolbox
rFun = {@blackmanharris, @flattopwin, @hann, @hamming, @rectwin};

% Number of test samples
n = [1:10, 1000*rand(1,10) + 24, round(logspace(3,5,10))];

% Numeric types
nType = {@double, @single, @int64, @uint64};

% Array to store maximum errors
errSym = zeros(length(n), length(tFun));
errPer = zeros(length(n), length(tFun));

% Validate test functions
for i = 1:length(tFun)
    % Empty input should return 0×1 empty double column vector
    if ~isempty(tFun{i}([]))
        warning(func2str(tFun{i}) + "([]) should return empty vector.");
        warning("Instead, it returned " + string(tFun{i}([])));
    end

    % 0 should return 0×1 empty double column vector
    if ~isempty(tFun{i}(0))
        warning(func2str(tFun{i}) + "(0) should return empty vector.");
        warning("Instead, it returned " + string(tFun{i}(0)));
    end

    % Verify that custom functions return the same values as built-ins
    for j = 1:length(n)
        % Generate the numeric type
        k = randi([1 length(nType)]);

        % Symmetric windows
        errSym(j, i) = max(abs(tFun{i}(nType{k}(n(j))) - rFun{i}(nType{k}(n(j)))));
        
        % Periodic windows (skip a uniform window)
        if ~isequal(rFun{5}, @rectwin)
            errPer(j, i) = max(abs(tFun{i}(nType{k}(n(j)), 'per') - ...
                                   rFun{i}(nType{k}(n(j)), 'per')));
        end
    end

    % Plot data
    figure;
    tiledlayout(2,1);

    ax = nexttile;
    bar(errSym(:, i));
    ax.XTick = 1:30;
    ax.XTickLabel = string(round(n,2));
    ax.XTickLabelRotation = 90;
    xlabel("Input value");
    ylabel("Maximum error");
    title(func2str(tFun{i}) + " - symmetric window");

    nexttile;
    bar(errSym(:, i));
    ax.XTick = 1:30;
    ax.XTickLabel = string(round(n,2));
    ax.XTickLabelRotation = 90;
    xlabel("Input value");
    ylabel("Maximum error");
    title(func2str(tFun{i}) + " - periodic window");
end