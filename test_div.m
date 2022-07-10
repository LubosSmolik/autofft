n = round(logspace(0, 7, 45));
tdiv = zeros(length(n), 1);
trdiv = zeros(length(n), 1);
tmult = zeros(length(n), 1);

for i = 1:length(n)
    r = rand(1, n(i));

    tic;
        r = 2 * r / (1024);
    tdiv(i) = toc;

    r = rand(1, n(i));

    tic
        r = 2 * r ./ (1024);
    trdiv(i) = toc;

    r = rand(1, n(i));

    tic
        r = r .* (2 / 1024);
    tmult(i) = toc;
end

figure;
axes(Box = "on", NextPlot = "add", XScale = "log", YScale = "log");
plot(n, tdiv, DisplayName = "/", LineWidth = 1);
plot(n, trdiv, DisplayName = "./", LineWidth = 1);
plot(n, tmult, DisplayName = ".*", LineWidth = 1);
