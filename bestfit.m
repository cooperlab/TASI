function [FitParameters,T_rank,best,order,FitMeanParameters,TMean_results]...
    = bestfit(Summary_directory,F,H1legend, fitline,fs)
% This function is to plot 4 figures for model fitting of dynamic
% morphology features and save them under Summary_directory.
% Figure1/2/3 are fitting for all individual spheroids using 3 models:
% poly1/poly2/exp1. Then perform student's t-test (or ANOVA test) 
% of the fitting parameters between 2 groups (or more than 2 groups)
% to calculate significant difference. 
% Finally select best fitting models based on the following criteria:
% score = adjusted_Rsquare/(mean_p_value_ttest)
% only use the pvalues of those fitting parameters measuring the slope, 
% not intercept among groups, and if the pvalue <=0.05, convert it to 0.05 
% for calculating the mean_p_value. Highest score is the best model.
% The outputs of fitting parameters, goodness of fit and best models 
% (for Fig1/2/3) will be saved as [FigName-model.mat];
% 
% Figure4 is the fitting of all 3 models for the mean curve+shaded STD for 
% each group. The outputs of fitting parameters, goodness of fit (for Fig4)
% will be saved as [FigName-mean-models.mat] under Summary_directory. 
%
% Inputs:
% Summary_directory: Input and Output folder for loading plotted features'
% data and saving output figures and data. If not specified, a window will
% pop up and let you select the summary directory.
%
% F: 1*N structure with 15 fields: 'Field', 'FeatureCol', 'feature',
% 'NormFeature', 'Mean', 'STD', 'CI', 'X', 'Ng', 'YLabel', 'FigName',
% 'LegendLabel', 'Linestyle', 'Marker', 'Width'.
% F(1).field is group1 data, F(2).field is group2 data, F(3).Field is
% group3 data, etc.
% If not specified, a window will pop up and let you select "Rinv.mat" or
% any "feature name.mat" you want to fit models.
% Detailed explaination for each field is the following:
% 'Field': 'fill' or 'unfill' for the feature plotted.
% 'FeatureCol': column of feature from morphology data.
% 'feature': a Nt*Ng matrix for the feature you want to plot.
% Nt is the timepoint of the sequence. Ng is the number of cell lines (or
% spheroids) in group1/2/3. eg: F(1).feature is feature data for group1.
% 'NormFeature': a Nt*Ng matrix with the normalized feature data for G1/2/3.
% The feature value at first timepoint was normalized to the same among
% all cell lines (or spheroids) within each group.
% G1_norm(t) = G1_feature(t) – G1_feature(t0) + G1_initial .
% F(2).NormFeature is normalized feature data for group2.
% 'MEAN': a Nt*1 vector of the mean value for the feature you want to plot.
% Nt is the timepoint of the sequence for each group.
% eg: F(1).MEAN is average feature data for group1 at each timepoint.
% 'STD': a Nt*1 vector of the STD for the feature you want to plot.
% Nt is the timepoint of the sequence for each group.
% eg: F(1).STD is std of feature data for group1 at each timepoint.
% 'CI': a Nt*1 vector of 95% confidence interval for the feature.
% Nt is the timepoint of the sequence for each group.
% eg: F(1).CI is 95% confidence interval for group1 feature at each timepoint.
% 'X': a 1*Nt vector for timepoints, will be used in the bestfit function as input.
% 'Ng': number of spheroids in each group.
% 'YLabel':label for y axis, will be used in the bestfit function as input.
% 'FigName': name to save the figure.
% 'LegendLabel': 1*N or N*1 cell, legend label string for group 1/2/3.
% F(1).LegendLabel is the legend label string for group1.
% 'Linestyle': 1*N or N*1 cell, each cell is a string for group 1/2/3 linestyle.
% 'Marker': 1*N or N*1 cell, each cell is a string for group 1/2/3 markers.
% 'Width': LineWidth for all curves.
%
%
% fitline: 1*N or N*1 cell, each cell is a string for fitting models'
% linestyle. eg. fitline = {'-.', '-.', '-.'}.
%
% fs: a 1*3 vector for fontsize for [all figures' fontsize; ...
% figure1/2/3 legends' fontsize; figure4 legend's fontsize]. 
% If not specified, the default = [14 8 14].
%
% Outputs:
% for Figure1/2/3:
% FitParameters: 1*1 structure with 3 fields, 'poly1', 'poly2', 'exp1'.
% Each field is a Ntotal * 5 table containing the fitting parameters
% (p1,p2,p3) and adjusted goodness of fit of each model for each spheroids. 
% Each row is a spheroid, each column is a fitting parameter. 
% Ntotal is the total number of spheroids.
% This output will also be saved in 3 csv files called
% 'FitParameters-poly1','FitParameters-poly2' and 'FitParameters-exp1.csv'.
%
% T_rank: 3*6 table, each row is a model. Column A is the score, the
% higher the better fitting for the model. Column 2 and 3 are mean adjusted
% Rsquare and mean Rsquare for all spheroids. Column 4,5,6 are the
% ttest (or anova) pvalues for fitting parameters p1,p2,p3 among groups.
% The order of each row (or each model) is ranked from best to worst.
% Tis output will also be saved as a csv file: 'T_rank.csv'.
%
% best: char, the best fitting model with the highest score.
% 
% order: 3*1 vector, the rank order number for 3 models (best to worst).
%
% for Figure4:
% FitMeanParameters: N*1 cell for each group.
% Each cell is a 3*5 matrix containing the fitting parameters for the mean
% Normal feature in each group. Row Labels are: 'poly1', 'poly2', 'exp1'.
% Column Labels are the fitting parameters (p1,p2,p3), adjusted Rsquare,
% and Rsquare for each model. 
% This output will also be saved in 3 csv files called
% 'FitParameters-Mean-[Group1 Name]','FitParameters-Mean-[Group2 Name]' 
% and 'FitParameters-Mean-[Group3 Name].csv'.
%
% TMean_results: 3*3 table, adjusted Rsquare for each model for each group.
% Each row is a group and each column is a model. 
% Tis output will also be saved as a csv file: 'Mean_gof.csv'.
%
% This function will also automatically generate 6 ANOVA table figures and
% 6 ANOVA box plots to show the statistical results, if N_group > 2.
% You can mannually save those figures if you need them.
%
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>


%% check input parameters, if not exist, pop-up an error
if ~exist('Summary_directory','var') || isempty(Summary_directory)
    Summary_directory = uigetdir([],'Please Choose the Summary Folder Containing Feature Figure');
end

% load F if non-exist
if ~exist('F','var') || isempty(F)
    [Ffile, ~] = uigetfile([Summary_directory filesep '*.mat'],...
        'Please Select the Feature Data For Curve Fitting, eg:"Rinv.mat"');
    load([Summary_directory filesep Ffile]);
end

% fitline style
if ~exist('fitline','var') || isempty(fitline)
    error('Please specify input3: fitline style, a 1*2 or 2*1 cell, with 2 strings for group1/2 MODEL FITTING line style, otherwise it will be "-."');
end

% fs
if ~exist('fs','var') || isempty(fs)
    disp('Warning: Please specify input4: fs, a 1*3 vector for font size, otherwise it will be [14 8 10]');
    fs = [14 8 12];
end
%% Initialization
% 3 models
fitType = {'poly1', 'poly2', 'exp1'};
FitParameters = struct(fitType{1},[], fitType{2},[], fitType{3},[]);
N = length(F);
% load plot setting parameters
Width = F.Width;
YLabel = F.YLabel;
FigName = F.FigName;
Ntotal = sum([F(:).Ng]);
colors1 = hsv(Ntotal);
% initialization for saving outputs
Glabel = zeros(Ntotal, 1);
fit_summary = zeros(length(fitType), 6);
%% model fitting and plotting for each individual spheroid
for type = 1:length(fitType)
    model = char(fitType(type));
    fit_results = zeros(Ntotal, 5);
    H1 = figure; hold on
    index = 1;
    % plot loop for each group j
    for j = 1:N
        x = [F(j).X]';
        Norm = F(j).NormFeature;
        Ng = size(Norm,2);
        % plot loop for each spheroid i in each group j
        for i = 1:Ng
            y = Norm(:,i);
            % fitting models
            if type == 1  %for poly1 model, parameters = [p1,p2]
                [f,gof] = fit(x,y,model);
                fit_results(index,:) = [f.p1 f.p2 NaN ...
                    gof.adjrsquare gof.rsquare];
            else if type == 2  %for poly2 model, parameters = [p1,p2,p3]
                    [f,gof] = fit(x,y,model);
                    fit_results(index,:) = [f.p1  f.p2  f.p3 ...
                        gof.adjrsquare gof.rsquare];
                else  %for exp model, parameters = [a,b]
                    [f,gof] = fit(x,y,model,'Exclude',x<=0); %exclude nonpositive points
                    fit_results(index,:) = [f.a f.b NaN ...
                        gof.adjrsquare gof.rsquare];
                end
            end
            % plot setting for each individual curve
            if strcmp(F(j).Marker, 'none')
                mk = [];
            else mk = F(j).Marker;
            end
            style = [F(j).Linestyle mk];
            c = colors1(index,:);
            h = plot(f,x,y,style);  % model and data
            set(h, 'color', c, 'LineWidth', Width);
            set(h(2), 'LineStyle',fitline{j}); %model line style
            % legend for each spheroid and its model
            legend1(2*index-1:2*index) = {[F(j).LegendLabel num2str(i)],...
                [F(j).LegendLabel num2str(i) ' ' model]};
            % add group labels and then increase index
            Glabel(index) = j;
            index = index + 1;
        end
    end
    % plot setting for fontsize, x y labels
    legend(legend1,'Location','best', 'FontSize', fs(2));
    xlabel('Time (minutes)'); ylabel(YLabel); title(model);
    set(gca, 'FontSize', fs(1));
    hold off
    savefig(H1, [Summary_directory filesep FigName '-' model '.fig']);
    saveas(H1,[Summary_directory filesep FigName '-' model '.tiff'],'tiffn');
    saveas(H1,[Summary_directory filesep FigName '-' model '.eps'],'epsc');
    % save fit_results to TFit table
    TFit = array2table(fit_results, 'RowNames',H1legend, ...
        'VariableNames',{'p1','p2','p3','Adjust_Rs','Rs'});
    FitParameters.(model) = TFit;
    writetable(TFit, [Summary_directory filesep ...
        'FitParameters-' model '.csv'], 'WriteRowNames',true);
    % calculate mean gof for each model
    gof1 = mean(fit_results(:,4));
    gof2 = mean(fit_results(:,5));
    %% statistical ttest (for 2 groups) or anova test (for >2 groups)
    % ttest for 2 groups
    if N == 2
        G1index = find(Glabel == 1);
        G2index = find(Glabel == 2);
        % equal variance
        if length(G1index) == length(G2index)
            [~,pv1] = ttest2(fit_results(G1index,1),fit_results(G2index,1));
            [~,pv2] = ttest2(fit_results(G1index,2),fit_results(G2index,2));
            % other than poly2 model, pv3 = NaN
            if type == 2
            [~,pv3] = ttest2(fit_results(G1index,3),fit_results(G2index,3));
            else pv3 = NaN;
            end
            % unequal variance
        else
            [~,pv1] = ttest2(fit_results(G1index,1),fit_results(G2index,1),...
                'Vartype','unequal');
            [~,pv2] = ttest2(fit_results(G1index,2),fit_results(G2index,2),...
                'Vartype','unequal');
            % other than poly2 model, pv3 = NaN
            if type == 2
            [~,pv3] = ttest2(fit_results(G1index,3),fit_results(G2index,3),...
                'Vartype','unequal');
            else pv3 = NaN;
            end
        end
    else if N > 2
            pv1 = anova1(fit_results(:,1), Glabel);
            title(['Fitting Parameter1 for ' model]);
            pv2 = anova1(fit_results(:,2), Glabel);
            title(['Fitting Parameter2 for ' model]);
            pv3 = anova1(fit_results(:,3), Glabel, 'off');
        end
    end
    %% calculate the score for each model
    % convert all p<=0.05 to p=0.05 for rank selection
    pv = [pv1 pv2 pv3];
    p_temp = pv;
    p_temp(pv<=0.05) = 0.05;
    if type == 1  %for poly1, only count parameter p1(slope) as pvalue
        score = gof1/p_temp(1);
    else score = gof1/((p_temp(1)+p_temp(2))/2); %for poly2,count p1,p2; exp:a,b
    end
    fit_summary(type,:) = [score gof1 gof2 pv];        
end
% convert fit_summary to table T_results
T_results = array2table(fit_summary, 'RowNames',fitType, ...
    'VariableNames',{'score','Adjust_Rs','Rs','pv1','pv2','pv3'});
% select best model
[~, best_index] = max(fit_summary(:,1));
best = char(fitType(best_index));
[~, order] = sort(fit_summary(:,1), 'descend');
T_rank = T_results(order,:);
% save outputs results
save([Summary_directory filesep FigName '-model.mat'],...
    'FitParameters', 'T_rank', 'best', 'order');
writetable(T_rank, [Summary_directory filesep 'T_rank.csv'],...
    'WriteRowNames',true);

%% Model fitting for mean curve for each group
% initialize Adjusted Rsquare output saving
Mean_gof = zeros(N, length(fitType));
FitMeanParameters = cell(N,1);
colors2 = {'ro','go','bo','yo','mo','co','ko'};
H2 = figure; hold on
hd = [];
for j = 1:N
    fitmean_results = zeros(length(fitType), 5);
    x = [F(j).X]';
    Mean = F(j).Mean;
    STD = F(j).STD;
    % plot mean curve data for each group in color dots first
    cindex = mod(j-1,length(colors2))+1; % repeated cycle for colors2 index
    h2 = shadedErrorBar(x, Mean, STD, colors2{cindex}, 1);
    hd = [hd h2.mainLine];
    % then plot each modeling on each mean curve
    for type = 1: length(fitType)
        model = char(fitType(type));
        % model fitting for mean
        if type == 1  %for poly1 model, parameters = [p1,p2]
            [f,gof] = fit(x,Mean,model);
            yfit1 = f.p1 * x + f.p2;
            m1 = plot(x,yfit1,'k-.','LineWidth',Width);
            fitmean_results(type,:) = [f.p1 f.p2 NaN ...
                gof.adjrsquare gof.rsquare];
        else if type == 2  %for poly2 model, parameters = [p1,p2,p3]
                [f,gof] = fit(x,Mean,model);
                yfit2 = f.p1 * x.^2 + f.p2 * x + f.p3;
                m2 = plot(x,yfit2,'k:','LineWidth',Width);
                fitmean_results(type,:) = [f.p1  f.p2  f.p3 ...
                    gof.adjrsquare gof.rsquare];
            else  %for exp model, parameters = [a,b]
                [f,gof] = fit(x,Mean,model,'Exclude',x<=0); %exclude nonpositive points
                yfit3 = f.a * exp(f.b * x);
                m3 = plot(x,yfit3,'k--','LineWidth',Width);
                fitmean_results(type,:) = [f.a f.b NaN ...
                    gof.adjrsquare gof.rsquare];
            end
        end
        Mean_gof(j,type) = gof.adjrsquare;
    end
    TFitmean = array2table(fitmean_results, 'RowNames',fitType,...
        'VariableNames',{'p1','p2','p3','Adjust_Rs','Rs'});
    FitMeanParameters{j} = fitmean_results;
    writetable(TFitmean, [Summary_directory filesep ...
      'FitParameters-Mean-' F(j).LegendLabel '.csv'],'WriteRowNames',true);
end
legend2_hd = [h2.mainLine m1 m2 m3];
xlabel('Time (Minutes)');
ylabel(['Models of Mean ' YLabel]);
legend(legend2_hd, ['Data',fitType], 'Location','best', 'FontSize',fs(3));
set(gca, 'FontSize', fs(1));
set(hd, 'LineWidth',1);
hold off
savefig(H2, [Summary_directory filesep FigName '-mean-models.fig']);
saveas(H2,[Summary_directory filesep FigName '-mean-models.tiff'],'tiffn');
saveas(H2,[Summary_directory filesep FigName '-mean-models.eps'],'epsc');
% save outputs for mean model fitting
TMean_results = array2table(Mean_gof, 'RowNames',{F(:).LegendLabel},...
    'VariableNames',fitType);
writetable(TMean_results, [Summary_directory filesep 'Mean_gof.csv'],...
    'WriteRowNames',true);
save([Summary_directory filesep FigName '-mean-models.mat'],...
    'FitMeanParameters', 'TMean_results');
end