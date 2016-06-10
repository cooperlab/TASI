%% Example1: for 3 groups comparison
% group morphology feature data into G1, G2 and G3
[Group, Summary_directory] = groupfeatures(3,[3 3 6], []);
% plot timelapse curve for each feature
[F, H1legend] = plotmorphology(Summary_directory,...
    Group,'fill',7,'Intensity of Filled Spheroids','Intensity', 10,...
    {'Leader','Follower','Parental'},{'-','--',':'},{'x','none','none'},...
    [],[]);
% model fitting and statistics
[FitParameters, T_rank, best, order, FitMeanParameters, TMean_results]...
    = bestfit(Summary_directory,F,H1legend,{'-.','-.','-.'},[]);

%% Example2: for 2 groups comparison
[Group, Summary_directory] = groupfeatures(2,[3 3], '2group');
[F, H1legend] = plotmorphology(Summary_directory, Group, 'fill', 14, ...
    'Invasive Radius of Filled Spheroids (\mum)', 'Rinv', 10, ...
    {'Leader','Follower'}, {'-','--'}, {'x','none'}, [], []);
[FitParameters, T_rank, best, order, FitMeanParameters, TMean_results]...
    = bestfit(Summary_directory,F,H1legend,{'-.','-.'},[]);