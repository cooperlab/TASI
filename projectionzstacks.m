function projection = projectionzstacks(Z,n_channel,channel,saveimage,Input_directory,Output_directory)
% this function is to import timelapse and z-stacked images (image format:
% t1_Z1_channel1,t1_Z1_channel2,t1_Z2_channel1,t1_Z2_channel2,
% t2_Z1_channel1,t2_Z1_channel2,t2_Z2_channel1,t2_Z2_channel2,...)
% load only one channel of images, projection all z-stacks at each 
% timepoint, and then save all the projected images as outputs.
% Projection is to calculate the STD of intensities in z direction.
% Background has a low STD, while cells or spheroids have a high STD.
% In this way we can enhance the contrast of our image, especially the 
% dimmer branches of our spheroid.
%
% inputs:
% Z: number of Z slices, eg. Z = 30.
%
% n_channel: number of channels, eg. n_channel = 1. 
%
% channel: indicate which channel you want to process, eg. channel = 2,
% means the second channel. If you named your channel as ch0, ch1, ch2, 
% channel = 2 processing for ch1 images.
%
% saveimage: 0 or 1 scalar, whether to save individual projected images 
% under Output_directory. 1 is save, 0 is not save. Default is 1.
% 
% Input_directory: input folder containing original images
% Output_directory: output folder to save results and projected images
%
% outputs:
% proj: M*N*T matrix for projected image sequences, 
% (M*N:image size; T:timepoint).
% 
% The function will also save all the projected images as 'proj00##.tif' 
% under Output_directory (when saveimage=1), ## is timepoints. 
% 
% Embedded functions needed: projstack
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% check input parameters, if not exist, pop-up an error  
if ~exist('Z','var') || isempty(Z)
    error('Please specify input1: number of Z slices');
end

if ~exist('n_channel','var') || isempty(n_channel)
    error('Please specify input2: number of channels');
end

if ~exist('channel','var') || isempty(channel)
    error('Please specify input3: which channel to analyze');
end

if ~exist('saveimage','var') || isempty(saveimage)
    saveimage = 1;
    display('Warning: Projected images will be saved under Output_directory. This will slow down the processes.');
elseif exist('saveimage','var') && (saveimage~=1) && (saveimage~=0)
    error('Input4 "saveimage" must be 0 or 1');
end

if ~exist('Input_directory','var') || isempty(Input_directory)
    Input_directory = uigetdir([],'Please Choose the Input Folder Containing Original Images');
end

if ~exist('Output_directory','var') || isempty(Output_directory)
    Output_directory = uigetdir([],'Please Choose the Output Folder for Saving Projection Results');
end

%% load raw images (z first, then time domain)
IList = dir([Input_directory filesep '*.tif']);
% total number of images
frame_number = length(IList);
% number of timepoints
timepoints = frame_number/Z/n_channel;
% index number for the images under the channel you want to analyze
vector = (channel:n_channel:frame_number)'; 
% 
breakpoint = (0:Z:frame_number/n_channel)';
% load first image I0 for size information
I0 = imread([Input_directory filesep IList(1).name]);
% preallocate for output matrix
stack = zeros(size(I0,1),size(I0,2),Z);
projection = zeros(size(I0,1),size(I0,2),timepoints);

%% stack in z direction
tic
for i = 2:length(breakpoint);
    index = vector((breakpoint(i-1)+1):breakpoint(i));
    % at each timepoint, convert all z slices into 3d matrix "stack"
    for j = 1:Z
        % load each z slice at timepoint i-1
        image_name = IList(index(j)).name;
        image = imread([Input_directory filesep image_name]);
        % convert to grayscale images and save in 3d matrix
        stack(:,:,j) = mat2gray(image);
    end
    % perform 3d STD projection in Z direction
    pj = projstack(stack);  
    pj = pj.std;
    projection(:,:,i-1) = pj;
    % save projected image for each timepoint when save = 1
    if saveimage == 1
    proj_name = [Output_directory filesep 'proj' num2str(i-1,'%04g') '.tif'];
    imwrite(pj, proj_name, 'tiff', 'Compression', 'none');
    end
end
% save output 3d projection matrix under Output_directory
save([Output_directory filesep 'projection.mat'],'projection');  

% display function elapse time
disp('Projection:');
toc