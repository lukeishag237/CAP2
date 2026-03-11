%% Load data
T = readtable('fascicle_data.csv');

%% Define nerve groups in the order you want them on the x-axis
Sciatic = ["39-2","41-2","44-1","46-2","48-2"];
Tibial  = ["39-1","41-1","44-2","46-1","48-1"];
cVN     = ["40-1","43-1","45-1","47-1"];

orderedSubjects = [Sciatic, Tibial, cVN];

%% -----------------------------
%  Extract fascicle-level rows
% -----------------------------
F = T(strcmp(T.mask_type,'i'), :);
F = F(ismember(string(F.subject), orderedSubjects), :);

subjectsF = string(F.subject);
majorF = F.major_axis_um;
minorF = F.minor_axis_um;
areaDF = F.area_diam_um;

%% Convert subject labels to numeric x positions
% Map each subject to its index in orderedSubjects
xF = arrayfun(@(s) find(orderedSubjects == s), subjectsF);

%% -----------------------------
%  Extract whole-nerve rows from SUMMARY
% -----------------------------
S = T(strcmp(T.mask_type,'summary'), :);
S = S(ismember(string(S.subject), orderedSubjects), :);

subjectsS = string(S.subject);
majorS = S.whole_major_axis_um;
minorS = S.whole_minor_axis_um;
areaDS = S.whole_area_diam_um;
xS = arrayfun(@(s) find(orderedSubjects == s), subjectsS);


%% -----------------------------
%  Create figure
% -----------------------------
figure('Color','w','Position',[100 100 1500 900]);

%% ============================
%   TOP PLOT: Fascicle diameters
% ============================
subplot(4,1,1); hold on; box on;

jit = 0.12;

h1 = scatter(xF - jit, minorF, 35, 'filled', 'MarkerFaceColor',[0.2 0.4 1]);
h2 = scatter(xF,       majorF, 35, 'filled', 'MarkerFaceColor',[1 0.4 0.2]);
h3 = scatter(xF + jit, areaDF, 35, 'filled', 'MarkerFaceColor',[0.2 0.7 0.2]);

set(gca,'XTick',1:numel(orderedSubjects),'XTickLabel',orderedSubjects);
ylabel('Fascicle Diameter (µm)');
grid on;

legend([h1 h2 h3], {'Minor Axis','Major Axis','Area-Derived Diameter'}, 'Location','best');

%% Add vertical dividers (not in legend)
nSci = numel(Sciatic);
nTib = numel(Tibial);

xline(nSci + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
xline(nSci + nTib + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');

%% Add section labels above the plot
yl = ylim;
text(nSci/2, yl(2)*1.02, 'Sciatic', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib/2 + 0.5, yl(2)*1.02, 'Tibial', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib + numel(cVN)/2 + 0.5, yl(2)*1.02, 'cVN', 'HorizontalAlignment','center', 'FontWeight','bold');

%% ============================
%   BOTTOM PLOT: Whole nerve
% ============================
subplot(4,1,2); hold on; box on;

jit = 0.12;

scatter(xS - jit, minorS, 60, 'filled', 'MarkerFaceColor',[0.2 0.4 1]);
scatter(xS,       majorS, 60, 'filled', 'MarkerFaceColor',[1 0.4 0.2]);
scatter(xS + jit, areaDS, 60, 'filled', 'MarkerFaceColor',[0.2 0.7 0.2]);

set(gca,'XTick',1:numel(orderedSubjects),'XTickLabel',orderedSubjects);
ylabel('Whole-Nerve Diameter (µm)');
grid on;

%% Add vertical dividers (not in legend)
xline(nSci + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
xline(nSci + nTib + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');

%% Add section labels
yl = ylim;
text(nSci/2, yl(2)*1.02, 'Sciatic', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib/2 + 0.5, yl(2)*1.02, 'Tibial', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib + numel(cVN)/2 + 0.5, yl(2)*1.02, 'cVN', 'HorizontalAlignment','center', 'FontWeight','bold');

%% ============================================================
%   NEW FIGURE: Fascicle Count + Stacked Area Composition
%% ============================================================

% figure('Color','w','Position',[100 100 1500 900]);

%% Extract summary rows again (already filtered earlier)
subjectsS = string(S.subject);
numF = S.num_fascicles;
axonA = S.axon_area_um2;
wholeA = S.whole_area_um2;
outerA = S.outer_area_um2;

% Convert subject labels to numeric x positions
xS = arrayfun(@(s) find(orderedSubjects == s), subjectsS);

%% ============================
%   TOP SUBPLOT: num_fascicles
% ============================
subplot(4,1,3); hold on; box on;

scatter(xS, numF, 80, 'filled', 'MarkerFaceColor',[0.1 0.1 0.9]);

set(gca,'XTick',1:numel(orderedSubjects),'XTickLabel',orderedSubjects);
ylabel('Number of Fascicles');
grid on;

% Add vertical dividers (not in legend)
nSci = numel(Sciatic);
nTib = numel(Tibial);

xline(nSci + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
xline(nSci + nTib + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');

% Add section labels
yl = ylim;
text(nSci/2, yl(2)*1.02, 'Sciatic', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib/2 + 0.5, yl(2)*1.02, 'Tibial', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib + numel(cVN)/2 + 0.5, yl(2)*1.02, 'cVN', 'HorizontalAlignment','center', 'FontWeight','bold');

%% ============================================
%   BOTTOM SUBPLOT: Stacked Area Composition (%)
%% ============================================

subplot(4,1,4); hold on; box on;

% Build matrix for stacked bar chart
rawAreas = [axonA, wholeA, outerA];

% Convert to percentages
rowSums = sum(rawAreas, 2);
pctAreas = rawAreas ./ rowSums * 100;   % each row sums to 100%

% Reorder rows to match orderedSubjects (via xS)
[~, order] = sort(xS);

pctAreas = pctAreas(order, :);
axonA    = axonA(order);
wholeA   = wholeA(order);
outerA   = outerA(order);
numF     = numF(order);   % optional, to keep subplot 3 perfectly aligned too
subjectsS = subjectsS(order);
xS       = xS(order);

hb = bar(pctAreas, 'stacked', 'BarWidth', 0.8);

% Color each component
hb(1).FaceColor = [0.2 0.4 1];   % Axon Area
hb(2).FaceColor = [0.2 0.7 0.2]; % All Fascicles Area
hb(3).FaceColor = [1 0.4 0.2];   % Surrounding Tissue Area

set(gca,'XTick',1:numel(orderedSubjects),'XTickLabel',orderedSubjects);
ylabel('Percent of Total Area (%)');
grid on;

legend({'Axon Area','All Fascicles Area','All Surrounding Tissue Area'}, ...
       'Location','bestoutside');

% Add vertical dividers (not in legend)
xline(nSci + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
xline(nSci + nTib + 0.5, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');

% Add section labels
yl = ylim;
text(nSci/2, yl(2)*1.02, 'Sciatic', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib/2 + 0.5, yl(2)*1.02, 'Tibial', 'HorizontalAlignment','center', 'FontWeight','bold');
text(nSci + nTib + numel(cVN)/2 + 0.5, yl(2)*1.02, 'cVN', 'HorizontalAlignment','center', 'FontWeight','bold');
