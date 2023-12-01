% Test functions
tFun = @utilities.autoKaiser;

% Reference functions - requires the Signal Processin Toolbox
rFun = @kaiser;

% Number of test samples
n = [1:10, 1000*rand(1,10) + 24, round(logspace(3,5,10))];

% Numeric types
nType = {@double, @single, @int64, @uint64};

% Values for beta
beta = [0, rand(1,3), 10*rand(1,3) + 1, 100*rand(1,3)+10];

% Array to store maximum errors
err = zeros(length(n), length(beta));

%% Validate test functions
% Empty input should return 0×1 empty double column vector
if ~isempty(tFun([]))
    warning(func2str(tFun) + "([]) should return empty vector.");
    warning("Instead, it returned " + string(tFun([])));
end

% 0 should return 0×1 empty double column vector
if ~isempty(tFun(0))
    warning(func2str(tFun) + "(0) should return empty vector.");
    warning("Instead, it returned " + string(tFun(0)));
end

% Initialize figure
figure;
tiledlayout(5,2);

% Verify that custom functions return the same values as built-ins
for i = 1:length(beta)
    for j = 1:length(n) 
    % Generate the numeric type
    k = randi([1 length(nType)]);

    % Symmetric windows
    err(j, i) = max(abs(tFun(nType{k}(n(j)), beta(i)) - ...
                        rFun(nType{k}(n(j)), beta(i))));
    end
    % Plot data
    ax = nexttile;
    bar(err(:, i));
    ax.XTick = 1:30;
    ax.XTickLabel = string(round(n,2));
    ax.XTickLabelRotation = 90;
    xlabel("Input value");
    ylabel("Maximum error");
    title("Kaiser window with beta = " + string(beta(i)));
end