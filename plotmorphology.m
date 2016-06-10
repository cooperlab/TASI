function [F, H1legend] = plotmorphology...
    (Summary_directory, Group, FieldName, FeatureCol, YLabel, FigName,...
    TimeInterval, LegendLabel, Linestyle, Marker, fontsize, Width)
% This function is to plot 3 figures for dynamic morphology features.
% Figure1 is the individual curve for each spheroid in each group.
% Figure2 is the mean feature curve with shaded STD for each group.
% Figure3 is the mean feature curve with shaded 95% CI for each group.
% Figure1 is named as [FigName]; Figure2 is named as [FigName-STD];
% Figure3 is named as [FigName-95CI]. Each figure will saved in 3 formats,
% '.fig', '.tiff', and '.eps' formats under Summary_directory folder.
% The output data structure F is also saved under Summary_directory folder.
%
% Inputs:
% Summary_directory: Input and Output folder for loading Grouped morphology
% data and saving output figures and data.
%
% Group: N*1 structure with 3 fields: 'fill', 'unfill' and 'Ng'.
% row1/2/3 is the filled or unfilled morphology data for each group.
% eg: for group1, I have 3 spheroids, then Group(1).fill is a 3*1 cell,
% each cell is the filled morphology matrix for each spheroid in group1.
% The morphology matrix is a matrix with size [timepoint, feature number].
% Each column is a feature, each row is a timepoint.
% Similarly format for 'unfill' field.
% 'Ng': number of spheroids in each group.
% If not specified, a pop-up window will let you select the Group.mat.
%
% FieldName: a string, must be either 'fill' or 'unfill', indicating which
% types of spheroid feature you want to analyze.
%
% FeatureCol: an integer, the column number of the feature you want to
% analyze. If not specified, an error will occur, but the feature name for
% each column will be displayed in the error for you to select the correct
% feature column.
%
% YLabel: a string, the label for y-axis.
%
% FigName: a string, the name you want to save the figure.
%
% TimeInterval: a number, the time interval (in minute) for image sequence.
%
% LegendLabel: 1*N or N*1 cell. Each cell contains a string, which is the
% legend label for each group. N is the number of groups.
% e.g. LegendLabel = {'Leader','Follower','Parental'};
%
% Linestyle: 1*N or N*1 cell. Each cell contains a string, which is the
% line style for plotting group1/2/3... curves.
% e.g. Linestyle = {'-','--',':'};
%
% Marker: 1*N or N*1 cell. Each cell contains a string, which is the
% marker for plotting group1/2/3... curves.
% e.g. Marker = {'x','none','o'};
%
% fontsize: 1*N or N*1 vector, the fontsize for
% [3 figures, figure1 legend, figure2/3 legend]. If not specified, the
% default = [14 10 14].
%
% Width: Linewidth for all curves. If not specified, default = 2.
%
% Outputs:
% F: 1*N structure with 15 fields: 'Field', 'FeatureCol', 'feature',
% 'NormFeature', 'Mean', 'STD', 'CI', 'X', 'Ng', 'YLabel', 'FigName',
% 'LegendLabel', 'Linestyle', 'Marker', 'Width'.
% F(1).field is group1 data, F(2).field is group2 data, F(3).Field is
% group3 data, etc.
%
% 'Field': 'fill' or 'unfill' for the feature plotted.
%
% 'FeatureCol': column of feature from morphology data.
%
% 'feature': a Nt*Ng matrix for the feature you want to plot.
% Nt is the timepoint of the sequence. Ng is the number of cell lines (or
% spheroids) in group1/2/3. eg: F(1).feature is feature data for group1.
%
% 'NormFeature': a Nt*Ng matrix with the normalized feature data for G1/2/3.
% The feature value at first timepoint was normalized to the same among
% all cell lines (or spheroids) within each group.
% G1_norm(t) = G1_feature(t) – G1_feature(t0) + G1_initial .
% F(2).NormFeature is normalized feature data for group2.
%
% 'MEAN': a Nt*1 vector of the mean value for the feature you want to plot.
% Nt is the timepoint of the sequence for each group.
% eg: F(1).MEAN is average feature data for group1 at each timepoint.
%
% 'STD': a Nt*1 vector of the STD for the feature you want to plot.
% Nt is the timepoint of the sequence for each group.
% eg: F(1).STD is std of feature data for group1 at each timepoint.
%
% 'CI': a Nt*1 vector of 95% confidence interval for the feature.
% Nt is the timepoint of the sequence for each group.
% eg: F(1).CI is 95% confidence interval for group1 feature at each timepoint.
%
% 'X': a 1*Nt vector for timepoints, will be used in the bestfit function as input.
%
% 'Ng': number of spheroids in each group.
%
% 'YLabel':label for y axis, will be used in the bestfit function as input.
%
% 'FigName': name to save the figure.
%
% 'LegendLabel': 1*N or N*1 cell, legend label string for group 1/2/3.
% F(1).LegendLabel is the legend label string for group1.
%
% 'Linestyle': 1*N or N*1 cell, each cell is a string for group 1/2/3 linestyle.
%
% 'Marker': 1*N or N*1 cell, each cell is a string for group 1/2/3 markers.
%
% 'Width': LineWidth for all curves.
%
%
% H1legend: Ntotal*1 cell with legend for each individual spheroid
%
% Embedded function: shadedErrorBar
%
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% check input parameters, if not exist, pop-up an error
if ~exist('Summary_directory','var') || isempty(Summary_directory)
    Summary_directory = uigetdir([],'Please Choose the Summary Folder Containing Grouped Morphology Data');
end

% Group
if ~exist('Group','var') || isempty(Group)
    [Groupfile, ~] = uigetfile([Summary_directory filesep '*.mat'],...
        'Please Select Grouped Data');
    load([Summary_directory filesep Groupfile]);
end

% FieldName
if ~exist('FieldName','var') || isempty(FieldName)
    error('Please specify input3: FieldName, either "fill" or "unfill"');
end

% FeatureCol
if ~exist('FeatureCol','var') || isempty(FeatureCol)
    if strcmp('fill', FieldName)
        disp('Column Labels for "fill" spheroid are:');
        disp('Col 1      2         3          4            5             6          7           8            9          10            11               12           13      14         15         16');
        disp('x_pixel y_pixel Rcore_pixel Rinv_pixel Perimeter_pixel Area_pixel Intensity Intensity_STD Complexity Eccentricity Branch_Number Single_Cell_Number Rcore_um Rinv_um Perimeter_um Area_um2');
    else if strcmp('unfill', FieldName)
            disp('Column Labels for "unfill" spheroid are:');
            disp('Col 1                2       3         4          5           6              7           8          9          10           11          12         13      14        15         16');
            disp('Ncell_in_spheroid	x_pixel	y_pixel	Rcore_pixel Rinv_pixel Hole_Number Perimeter_pixel Area_pixel Intensity	Intensity_STD Complexity Eccentricity Rcore_um Rinv_um Perimeter_um	Area_um2');
        else
            error('input3:FieldName must be either "fill" or "unfill"');
        end
    end
    error('Please specify input4: Column Number for the Feature You Want to Plot');
end
% YLabel
if ~exist('YLabel','var') || isempty(YLabel)
    error('Please specify input5: YLabel, a string for your Y axis label');
end
% FigName
if ~exist('FigName','var') || isempty(FigName)
    error('Please specify input6: FigName, a string of the name to save your figure');
end
% TimeInterval
if ~exist('TimeInterval','var') || isempty(TimeInterval)
    error('Please specify input7: TimeInterval, a number in minutes');
end
% LegendLabel
if ~exist('LegendLabel','var') || isempty(LegendLabel)
    error('Please specify input8: LegendLabel, a 1*N or N*1 cell, with N strings for group1/2/3... legends');
end

% Linestyle
if ~exist('Linestyle','var') || isempty(Linestyle)
    error('Please specify input9: Linestyle, a 1*N or N*1 cell, with N strings for group1/2/3... line style');
end

% Marker
if ~exist('Marker','var') || isempty(Marker)
    error('Please specify input10: Marker, a 1*N or N*1 cell, with N strings for group1/2/3... marker symbols');
end

% fontsize
if ~exist('fontsize','var') || isempty(fontsize)
    disp('Warning: Please specify input11: fontsize, 3*1 vector for font size for [figures,fig1 legend,fig2/3 legend], otherwise it will be [14 10 14]');
    fontsize = [14 10 14];
end

% Width
if ~exist('Width','var') || isempty(Width)
    disp('Warning: Please specify input12: line width for all curves, otherwise it will be 2');
    Width = 2;
end

%% Initialization for saving output
N = length(Group);  % number of groups
Nt = zeros(1, N);   % number of time point matrix for all groups
feature = cell(1, N); % un-normalized original feature data for all groups
norm = cell(1, N); % normalized feature data for all groups
Mean = cell(1, N); % mean feature data for all groups
STD = cell(1, N); % STD of feature data for all groups
CI = cell(1, N); % 95% Confidence interval of feature data for all groups
X = cell(1, N);  % time points for x-axis
Ns = cell(1, N); % number of spheroids in each group
%% plot normalized curve for each spheroid in each group
Ntotal = sum([Group(:).Ng]);
H1 = figure; hold on
colors1 = hsv(Ntotal);
H1legend = cell(Ntotal,1);
index = 1;
%% extract individual feature's data from each spheroids matrix
for j = 1:N
    G = Group(j).(FieldName); % fill or unfill data for each group G
    Ng = length(G);  % number of spheroids in each group
    tp = size(G{1},1); % time points for x axis
    Nt(j) = tp;
    x = 0:TimeInterval:(tp-1)*TimeInterval; % x interval for plotting
    X{j} = x;
    Ns{j} = Ng;
    %% group feature for each group, put this feature data for each
    % spheroid in group j in a tp*Ng matrix called Feature
    Feature = zeros(tp, Ng);
    for i = 1:Ng
        Feature(:,i) = G{i}(:,FeatureCol);
    end
    feature{j} = Feature;
    %% calculate mean, STD and 95% confidence interval for each group
    G_mean = mean(Feature,2); % mean for each group
    G_std = std(Feature,0,2); % std for each group
    G_ci = abs(norminv(0.025).*(G_std/sqrt(Ng))); % 95% confidence interval
    Mean{j} = G_mean;
    STD{j} = G_std;
    CI{j} = G_ci;
    %% normalization for each group, using the following function:
    % Norm(t) = Feature(t) – Feature(t0) + Initial
    % Initial = G_mean(1),which is initial value for G1/G2/G3... feature
    Norm = zeros(tp, Ng);
    for k = 1:Ng
        Norm(:,k) = Feature(:,k)- Feature(1,k) + G_mean(1);
        plot(x, Norm(:,k),'Color',colors1(index, :), ...
            'LineStyle',Linestyle{j},'Marker',Marker{j},'LineWidth',Width);
        H1legend(index) = {[LegendLabel{j} num2str(k)]};
        index = index + 1;
    end
    norm{j} = Norm;
end
xlabel('Time (Minutes)');
ylabel(YLabel);
xlimit = max(Nt) * TimeInterval;  % maximum x interval
xlim([0 xlimit]); %xlim([0 1000]);
legend(H1legend, 'Location','best', 'FontSize',fontsize(2));
set(gca, 'FontSize', fontsize(1));
hold off
savefig(H1, [Summary_directory filesep FigName '.fig']);
saveas(H1, [Summary_directory filesep FigName '.tiff'], 'tiffn');
saveas(H1, [Summary_directory filesep FigName '.eps'], 'epsc');

%% plot shaded STD with mean curve for each group
H2 = figure; hold on
colors2 = {'r','g','b','y','m','c','k'};
legend2_hd = [];
for j = 1:N
    x = X{j};
    if strcmp(Marker{j}, 'none')
        mk = [];
    else mk = Marker{j};
    end
    cindex = mod(j-1,length(colors2))+1; % repeated cycle for colors2 index
    setting = [colors2{cindex} Linestyle{j} mk];
    h2 = shadedErrorBar(x, Mean{j}, STD{j},...
        {setting, 'markerfacecolor',colors2{cindex}}, 1);
    legend2_hd = [legend2_hd h2.mainLine];
end
xlabel('Time (Minutes)');
ylabel(['Mean ' YLabel]);
xlim([0 xlimit]);  %xlim([0 1000]);
title('STD');
legend(legend2_hd, LegendLabel, 'Location','best', 'FontSize',fontsize(3));
set(gca, 'FontSize', fontsize(1));
set(legend2_hd, 'LineWidth',Width);
hold off
% save figure
savefig(H2, [Summary_directory filesep FigName '-STD.fig']);
saveas(H2, [Summary_directory filesep FigName '-STD.tiff'], 'tiffn');
saveas(H2, [Summary_directory filesep FigName '-STD.eps'], 'epsc');
%% plot shaded 95% confidence interval with mean curve for each group
H3 = figure; hold on
colors2 = {'r','g','b','y','m','c','k'};
legend3_hd = [];
for j = 1:N
    x = X{j};
    if strcmp(Marker{j}, 'none')
        mk = [];
    else mk = Marker{j};
    end
    cindex = mod(j-1,length(colors2))+1; % repeated cycle for colors2 index
    setting = [colors2{cindex} Linestyle{j} mk];
    h3 = shadedErrorBar(x, Mean{j}, CI{j}, ...
        {setting, 'markerfacecolor',colors2{cindex}}, 1);
    legend3_hd = [legend3_hd h3.mainLine];
end
xlabel('Time (Minutes)');
ylabel(['Mean ' YLabel]);
xlim([0 xlimit]);  %xlim([0 1000]);
title('95% Confidence Interval');
legend(legend3_hd, LegendLabel, 'Location','best', 'FontSize',fontsize(3));
set(gca, 'FontSize', fontsize(1));
set(legend3_hd, 'LineWidth',Width);
hold off
% save figure
savefig(H3, [Summary_directory filesep FigName '-95CI.fig']);
saveas(H3, [Summary_directory filesep FigName '-95CI.tiff'], 'tiffn');
saveas(H3, [Summary_directory filesep FigName '-95CI.eps'], 'epsc');

%% save outputs
F = struct('Field',FieldName, 'FeatureCol',FeatureCol, 'feature',feature,...
    'NormFeature',norm, 'Mean',Mean, 'STD',STD, 'CI',CI, 'X',X, 'Ng',Ns,...
    'YLabel',YLabel, 'FigName',FigName, 'LegendLabel',LegendLabel, ...
    'Linestyle',Linestyle, 'Marker',Marker, 'Width',Width);
save([Summary_directory filesep FigName '.mat'], 'F', 'H1legend');