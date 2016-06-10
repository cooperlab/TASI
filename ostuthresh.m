function []=ostuthresh(Output_directory, graphcuts, thresh, area, se)
% this function is to threshold the graphcuts results into binary masks,
% and save all the binary masks under the path folder. Then generate a text
% file "thresh_parameters.txt" saving all the input parameter settings.
%
% inputs:
% Output_directory: output folder for saving threshold binary images.
%
% graphcuts: 3d matrix (x,y,t) or (x,y,z) from graphcuts segmentation 
% or other segmentation methods.
% 
% thresh: thresholding value for im2bw, range[0,1].
% if not specified, use Ostu graythresh method.
% 
% area: threshold area for background noises, remove all objects smaller
% than area value, default = 10 pixel.
%
% se: disk structure size for close gaps. Se must be integers, possible 
% ranges for se: [3-8], the larger the gaps on the boundary, se should be 
% larger. Default = 0, no close gaps.
%
% Outputs: NA
% This function will automatically save all the binary images under
% Output_directory as 'mask00##.tif', ## is timepoint.
% It will also automatically generate a 'thresh_parameters.txt' file to 
% save all the thresholding parameters (thresh, area, se) for future
% refining/adjustment of the thresholding.
%
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% check input parameters, if not exist, pop-up an error or a warning
if ~exist('Output_directory','var') || isempty(Output_directory)
    Output_directory = uigetdir([],'Please Choose the Output Folder for Saving Binary Images');
end

if ~exist('graphcuts','var') || isempty(graphcuts)
    graphcuts_filename = uigetfile([Output_directory filesep '*.mat'],...
        'Please Select the Input 3D Graphcuts Matrix for Thresholding');
    load([Output_directory filesep graphcuts_filename]);
end

if ~exist('thresh','var') || isempty(thresh)
    disp('Warning: Please specify input3:thresh value, otherwise will use Ostu graythresh method');
end

if ~exist('area','var') || isempty(area)
    disp('Warning: Please specify input4: background area thershold(pixel), otherwise area=10');
    area = 10;
end

if ~exist('se','var') || isempty(se)
    disp('Warning: Please specify input5: must be integers, gap close size, otherwise no gap close');
    se = 0;
end

%% postprocessing 
tic
for i = 1:size(graphcuts,3)
    % load each timepoint image as mask
    mask = graphcuts(:,:,i);
    % convert from single to double (grayscale) image
    maskd = mat2gray(mask);
    % if exist thresh value, use fixed thresholding; if not, use Ostu
    % thresh
    if ~exist('thresh','var') || isempty(thresh)
        bw = im2bw(maskd, graythresh(maskd));
    else bw = im2bw(maskd, thresh);
    end
    % remove small background noisies, which is less than area1
    bw1 = bwareaopen(bw, area);
    %% imclose for follower spheroids or dimmer inside images
    if se ~= 0
        bw2 = imclose(bw1, strel('disk',se));
        bw3 = imfill(bw2, 'holes');
        bw1 = bw3;
    end
    % delete foreground pixels that touched the border of the image
    bw_shrink = bw1(2:end-1, 2:end-1);
    final_mask = padarray(bw_shrink, [1 1]);
    % save the threshod images as 'mask0000.tif', 'mask0001.tif', ...
    filename = [Output_directory filesep 'mask' num2str(i,'%04g') '.tif'];
    imwrite(final_mask, filename, 'tiff', 'Compression', 'None');
end

%% write all the parameters to a .txt file
% Open or create new text file for reading and writing. Discard existing contents, if any.
fid = fopen([Output_directory filesep 'thresh_parameters.txt'],'wt+');
% save the parameters value for rerun or refine analysis
if ~exist('thresh','var') || isempty(thresh)
    fprintf(fid, 'Ostu Graythresh;\n');
else fprintf(fid, 'thresh = %f;\n', thresh);
end
fprintf(fid,'area = %d;\n', area);
if se == 0
    fprintf(fid, 'No Gap Close;\n');
else fprintf(fid, 'se = %d;\n', se);
end
fclose(fid);

% display function elapse time
disp('Thresholding:');
toc