function drawboundaries(Output_directory, enhance,movie_name,Frame_rate,Linewidth,color)
% draw boundaries from masks on original images
% example: drawboundaries(1,[],'Mg.avi',7); %enhance contrast for dim image.
% example: drawboundaries(0,[],'Mr.avi',7); %don't enhance for bright image.
%
% Inputs:
% Output_directory: output folder for saving output movie.
% 
% enhance: whether to enhance the contrast of the original images for
% better visualization, but will slow down the speed. It must be either 1
% or 0. 1 is to enhance, 0 is not. Default is 0, not enhance.
%
% movie_name: 'xxxxx' any name you want to call the output movie; If not
% specified, default is "movie".
%
% Frame_rate: scalar or integer. fps when playing the movie; default = 10.
%
% Linewidth: scaler or integer. pixels for the line width drawing boundaries;
%
% color: must be in format [x x x], where x is any number between 0 and 1.
% For example, [0 0 1] is blue; [1 1 1] is white(default); 
% [1 0 0] is red; [0 1 0] is green, etc.
%
% Outputs: NA
% This function will automatically create an output movie called 
% [movie_name '.avi'] file contains all original images embedded with mask 
% boundaries.
%
% Embedded function or code: embedboundaries.m
%
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% parameters for drawing boundaries
if ~exist('Output_directory','var') || isempty(Output_directory)
    Output_directory = uigetdir([],'Please Choose the Output Folder for Saving Boundary Movie');
end

if ~exist('enhance','var') || isempty(enhance)
    disp('Warning:Please choose whether to enhance contrast or not, otherwise it will not be enhanced');
    enhance = 0; % default: not enhance the contrast of original images.
end

if ~exist('movie_name','var') || isempty(movie_name)
    disp('Warning:Please name the output movie, otherwise it will be called movie');
    movie_name = 'movie'; % default name for output movie
end

if ~exist('Frame_rate','var') || isempty(Frame_rate)
    disp('Warning:Please setup movie play frame rate, otherwise it will be 10 fps');
    Frame_rate = 10; % frame rate (frame per second) for playing the movie
end

if ~exist('Linewidth','var') || isempty(Linewidth)
    disp('Warning:Please setup Linewidth, otherwise it will be 2');
    Linewidth = 2;
end

if ~exist('color','var') || isempty(color)
    disp('Warning:Please choose boundary color, otherwise it will be white')
    color = [1 1 1];  % white boundaries, or blue:[0 0 1]
end

se = strel('disk',Linewidth);  % LineWidth = 2

%% load all input images
% load projected images, if your projected image is not named as 
% 'proj0001.tif', 'proj0002.tif', change the following line according to 
% the name of your image. eg. if your image is called 'image0001.tif', etc,
% change the following line to:
% IList = dir([Output_directory filesep 'image*.tif']);
IList = dir([Output_directory filesep 'proj*.tif']);
number_of_frames = length(IList);
% load masks, similar as load projected images, if your mask is not named
% as 'mask0001.tif', etc, change the name accordingly. eg. if your mask is
% called 'binary0001.tif', etc. change the following line to:
% mList = dir([Output_directory filesep 'binary*.tif']);
mList = dir([Output_directory filesep 'mask*.tif']);
% create an ".avi" file called movie_name and save it under
% Output_directory
writerObj = VideoWriter(fullfile(Output_directory, [movie_name '.avi']));
writerObj.FrameRate = Frame_rate;   % set up the movie play rate is 10 fps
open(writerObj);                    % open the ".avi" file for writing

%% merge green and red and then draw mask boundaries on merged images
for frame_number = 1:number_of_frames
    % read green image
    original_filename = IList(frame_number).name;
    original = imread([Output_directory filesep original_filename]);
    % read mask image
    mask_filename = mList(frame_number).name;
    mask = imread([Output_directory filesep mask_filename]);
    [M, N] = size(mask);
    
    % save boundary matrix
    B = bwperim(mask);
    % make the boundary LineWidth = 2 pixels
    Bw = imdilate(B,se);
    % enhance the contrast
    if enhance == 1
    limits = stretchlim(original, 0.01); %return the bottom and top 1% of all pixels
    original_adjusted = imadjust(original, limits, []); %limits are scaled between 0 and 100% of all pixels
    original = original_adjusted;
    end
    % create merged green and red original image
    if(isa(original, 'uint16'))
    Merge = uint16(zeros(M,N,3));
    Merge(:,:,1) = original; 
    Merge(:,:,2) = original;
    elseif(isa(original, 'uint8'))
        Merge = uint8(zeros(M,N,3));
        Merge(:,:,1) = original;
        Merge(:,:,2) = original;
    else
        error('Original Image should be 16bit or 8bit');
    end
    % embed boundaries in white color
    Mg = embedboundaries(Merge, Bw, color);
    % save the image into video
    if (isa(original, 'uint16'))
    Mg_uint8 = im2uint8(Mg);
    writeVideo(writerObj,Mg_uint8);
    else
        writeVideo(writerObj,Mg);
    end
end
close(writerObj);