function [Group, Summary_directory] = groupfeatures(N, Ngi, GroupName)
% This function is to group individual spheroid's "morphology.mat" data 
% into N groups: g1, g2, ..., gN. 
%
% Inputs:
% N: interger, number of groups.
%
% Ngi: a 1*N or N*1 vector, [Ng1 Ng2 Ng3 ...]. eg. Ngi = [3 3 6]; 
% Ng1/Ng2/Ng3 is number of spheroids in group1/2/3.
% They can be same or different integers.
%
% GroupName: a string, the name for saving 'Group.mat' output structure.
% Default is 'Group'.
%
% Outputs:
% Group: N*1 structure with 3 fields: 'fill', 'unfill' and 'Ng'.
% row1/2/3 is the filled or unfilled morphology data for each group.
% eg: for group1, I have 3 spheroids, then Group(1).fill is a 3*1 cell,
% each cell is the filled morphology matrix for each spheroid in group1.
% The morphology matrix is a matrix with size [timepoint, feature number].
% Each column is a feature, each row is a timepoint.
% Similarly format for 'unfill' field. 
% 'Ng': number of spheroids in each group.
%
% Summary_directory: Output folder for saving Group structure data.
%
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% check input parameters, if not exist, pop-up an error  
if ~exist('N','var') || isempty(N)
    error('Please specify input1: N, integer, number of groups');
end

if ~exist('Ngi','var') || isempty(Ngi)
    error('Please specify input2: Ngi, 1*N or N*1 vector, number of spheroids in each group');
end

if length(Ngi) ~= N
    error(['Ngi must be a 1*' num2str(N) ' or ' num2str(N) '*1 vector!']);
end

if ~exist('GroupName','var') || isempty(GroupName)
    disp('Warning: Please specify input4: GroupName-a string, the name for saving "Group.mat" output; Otherwise it will be saved as"Group.mat"');
    GroupName = 'Group';
end
%% choose Input/Output folders and initialize output cell arrays
% choose the base folder
OutputBase = uigetdir([],'Please Choose the Base Output Folder');
% choose the summary/result/statistic folder to save the data
Summary_directory = uigetdir(OutputBase, 'Please Choose the Summary Folder to Save Statistical Results');
% use cell array to preallocate data, 20 times faster than multidimensional
% array preallocating.
Group_fill = cell(N,1);
Group_unfill = cell(N,1);
Ns = cell(N,1);
for j = 1:N
    Ng = Ngi(j);
    G_fill = cell(Ng,1);
    G_unfill = cell(Ng,1);
    %% load data and grouped them into g1, g2, ....
    % choose the morphology data for gj
    for i = 1:Ng
        [filename, filepath] = uigetfile([OutputBase filesep '.mat'],...
            ['Please Select the morphology.mat file for Group' num2str(j)...
            '-#'  num2str(i) ' CellLines/Spheroids']);
        load([filepath filename]);
        G_fill{i,1} = S_fill;
        G_unfill{i,1} = S_unfill;
    end
    Group_fill{j} = G_fill;
    Group_unfill{j} = G_unfill;
    Ns{j} = Ng;
end
%% save the output cell arrays into structures with 2 fields
Group = struct('fill',Group_fill, 'unfill',Group_unfill, 'Ng',Ns);
save([Summary_directory filesep GroupName '.mat'],'Group');