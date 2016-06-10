function [proj_3dgauss,proj_adjust] = preprocess(Output_directory,projection,sigma,adjust,saveimage)
% this function is to perform 3d gaussian filtering and then enhance
% contrast of the filtered images.
% Inputs:
% Output_directory: output folder to save results and images. Don't change
% it. 
% 
% projection: 3d M*N*T matrix, which is the output of projectionzstacks
% function. Don't change it.
% 
% sigma: 1*3 vector with sigma values in x, y, t(or z) direction for 
% gaussian filter. The larger the sigma, more smooth but blurry filtering 
% in that direction. Default is [2.5, 2.5, 0.5].
% 
% adjust: 0 or 1 scalar, whether to perform contrast enhancement on 
% filtered images, 1 is adjust, 0 is not adjust. Default is 1.
% 
% saveimage: 0 or 1 scalar, whether to save individual filtered images or/and
% adjusted images under the Output_directory. 1 is save, 0 is not save.
% Default is 1.
%
% Outputs:
% proj_3dgauss: 3d M*N*T gaussian filtered matrix, [M, N] is the image
% size. T is timepoint.
% 
% proj_adjust: when adjust = 1, proj_adjust is a 3d M*N*T adjust matrix, 
% which is gaussian filtered and then enhanced contrast, [M, N] is the 
% image size. T is timepoint. When adjust = 0, proj_adjust is empty.
%
% The function will also save all the filtered or/and adjusted images as
% 'gaussian00##.tif' or 'gadj00##.tif' under Output_directory
% (when saveimage=1),## is timepoints. 
% It will also automatically generate a 'preprocess_parameters.txt' file to 
% save all the gaussian filtering parameters (sigma, filter method, size)  
% for future refining/adjustment of the gaussian filtering.
% 
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% check input parameters, if not exist, pop-up an error or use default 
if ~exist('Output_directory','var') || isempty(Output_directory)
    Output_directory = uigetdir([],'Please Choose the Output Folder for Saving Filtered Results');
end

if ~exist('projection','var') || isempty(projection)
    projection_filename = uigetfile([Output_directory filesep '*.mat'],...
        'Please Select the Input 3D Projected Matrix');
    load([Output_directory filesep projection_filename]);
end

if ~exist('sigma','var') || isempty(sigma)
    sigma = [2.5, 2.5, 0.5]; % default sigma for 3d gaussian filter
    display('Warning: Please specify input1: sigma, otherwise sigma = [2.5, 2.5, 0.5];');
end

if ~exist('adjust','var') || isempty(adjust)
    adjust = 1;
    display('Warning: will enhance contrast of 3d gaussian filtered images. This will slow down the processes.');
elseif exist('adjust','var') && (adjust~=1) && (adjust~=0)
    error('Input2 "adjust" must be 0 or 1');
end

if ~exist('saveimage','var') || isempty(saveimage)
    saveimage = 1;
    display('Warning: Projected images will be saved under Output_directory. This will slow down the processes.');
elseif exist('saveimage','var') && (saveimage~=1) && (saveimage~=0)
    error('Input3 "saveimage" must be 0 or 1');
end

%% perform 3d gaussian filter to smooth both in xy and time domain
% check the Matlab version, if later than R2015a, can use imgaussfilt3
% function, if not, have to use the simple 3d_gaussian_filter G
tic
v = ver('MATLAB');
version = str2double(v(1).Version);
if version >= 8.5
    proj_3dgauss = imgaussfilt3(projection, sigma);
else
    % if Matlab version earlier than R2015a, use the simple
    % 3d_gaussian_filter G
    siz = [2, 2, 4]; % default filter size, can be changed if needed
    % modified according to fspecial3.m function on FileExchange, url:
    % http://www.mathworks.com/matlabcentral/fileexchange/21130-dti-and-fiber-tracking/content/fspecial3.m
    [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
    G = exp(-(x.*x/2/sigma(1)^2 + y.*y/2/sigma(2)^2 + z.*z/2/sigma(3)^2));
    G = G/sum(G(:));
    proj_3dgauss = imfilter(projection, G);
end
% save the 3d gaussian filtered matrix for graphcuts processing
save([Output_directory filesep 'proj_3dgauss.mat'],'proj_3dgauss');
disp('3D Gaussian Filter:');
toc

%% optional: enhance contrast on 3d gaussian filtered images
% for follower spheroid, which are very dim both on inside core and outside
% boundaries, need to enhanced the contrast for better segmentations
% when adjust=1, enhance contrast
tic
if adjust == 1
    proj_adjust = zeros(size(proj_3dgauss));
    for i = 1:size(proj_3dgauss,3)
        gaussian = proj_3dgauss(:,:,i);
        gaussian_adj = imadjust(gaussian);
        proj_adjust(:,:,i) = gaussian_adj;
        if saveimage == 1
            gaussian_name = [Output_directory filesep 'gaussian' ...
                num2str(i, '%04g') '.tif'];
            imwrite(gaussian, gaussian_name, 'tiff','Compression','none');
            adjust_name = [Output_directory filesep 'gadj' ...
                num2str(i, '%04g') '.tif'];
            imwrite(gaussian_adj, adjust_name,'tiff','Compression','none');
        end
    end
    % save the enhanced contrast 3d matrix for segmentation inputs
    save([Output_directory filesep 'proj_adjust.mat'], 'proj_adjust');
    % if adjust=0, but saveimage=1, only save 3d gaussian filtered image
elseif saveimage == 1
    for i = 1:size(proj_3dgauss,3)
        gaussian = proj_3dgauss(:,:,i);
        gaussian_name = [Output_directory filesep 'gaussian' ...
            num2str(i, '%04g') '.tif'];
        imwrite(gaussian, gaussian_name, 'tiff','Compression','none');
    end
    proj_adjust = [];
else proj_adjust = [];
end
disp('Adjust (or Enhance Contrast) of 3d gaussian filtered image:');
toc

%% write all the parameters to a .txt file
% Open or create new text file for reading and writing. Discard existing contents, if any.
fid = fopen([Output_directory filesep 'preprocess_parameters.txt'],'wt+');
% save the parameters value for rerun or refine segmentation
fprintf(fid,'sigma = [%f, %f, %f];\n', sigma);
if version >= 8.5
    fprintf(fid,'use imgaussfilt3 function;\n');
else fprintf(fid,'use simple filter G;\n');
    fprintf(fid,'size = [%d, %d, %d];\n', siz);
end
fclose(fid);