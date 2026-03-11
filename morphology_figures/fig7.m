%% Load data
T = readtable('fascicle_data.csv');

%% Extract only summary rows
S = T(strcmp(T.mask_type,'summary'), :);

%% Define nerve groups
Sciatic = ["39-2","41-2","44-1","46-2","48-2"];
Tibial  = ["39-1","41-1","44-2","46-1","48-1"];
cVN     = ["40-1","43-1","45-1","47-1"];

nerveGroups = {Sciatic, Tibial, cVN};
nerveNames  = {'Sciatic','Tibial','cVN'};

%% Metrics to plot
wholeA = S.whole_area_um2;     % Whole nerve area
numF   = S.num_fascicles;      % Number of fascicles
outerA = S.outer_area_um2;     % Surrounding tissue area
fascA  = S.total_area_um2;     % All fascicles area
axonA  = S.axon_area_um2;      % Axon area

% Percentages
rawMat = [outerA, fascA, axonA];
pctMat = rawMat ./ sum(rawMat,2) * 100;

% Combine into 8-column matrix
allMetrics = [wholeA, numF, outerA, fascA, axonA, pctMat];

%% Column axis labels (only bottom row)
xLabels = ["Whole Nerve Area (µm^2)", ...
           "# Fascicles", ...
           "Outer Area (µm^2)", ...
           "All Fascicles Area (µm^2)", ...
           "Axon Area (µm^2)", ...
           "Outer Area (%)", ...
           "All Fascicles Area (%)", ...
           "Axon Area (%)"];

nCols = numel(xLabels);

%% Create figure
figure('Color','w','Position',[100 100 1800 900]);

%% Loop through nerve groups (rows)
for r = 1:3
    group = nerveGroups{r};
    idx = ismember(string(S.subject), group);

    metrics_r = allMetrics(idx, :);

    % Collapse to a single y level (true 1D scatter)
    y_r = ones(sum(idx),1);

    for c = 1:nCols
        subplot(3, nCols, (r-1)*nCols + c); hold on; box on;

        scatter(metrics_r(:,c), y_r, 50, 'filled');

        % Y-axis: no ticks or labels, just row title on leftmost column
        set(gca,'YTick',[]);
        if c == 1
            ylabel(nerveNames{r}, 'FontWeight','bold');
        end

        % X-axis: numbers on all rows, labels only on bottom row
        if r == 3
            xlabel(xLabels(c));
        end

        grid on;
        ax = gca;
        ax.XAxis.Exponent = 0;

    end
end
