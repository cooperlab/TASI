function [T_unfill,S_unfill,T_fill,S_fill,single] = featureextraction...
    (ZSliceName, Resolution, SingleCellArea, smooth, plotsetting,...
    Input_directory, Output_directory)
% This function is to smooth the binary mask of a spheroid, then extract 16
% morphology features of the unfilled spheriod, filled spheroid and single
% cells isolated from the spheroid. All the output results (matrix, text
% files and images) will be saved under Output_directory folder.
%
% Inputs:
% ZSliceName: a string, part of the name of your original fluorescent 
% images. e.g. If you want to select Z5 plane and channel 00, then 
% ZSliceName = 'z5_ch00'. If you want to select Z1 and channel 1, then 
% ZSliceName = 'z1_ch1', according to the name of your original images.
%
% Resolution: resolution of your camera, unit: um/pixel
% 
% SingleCellArea: threshold value to remove debris smaller than
% SingleCellArea, unit: um^2; default is 200 um^2, a typical HeLa cell on a
% confluent 2d substrate is ~10um in radius, so if in a compressed 3d 
% spheroid, use a threshold radius as 8um, then area = pi*8^2 = 201 um^2.
%
% smooth: whether have Ismooth image sequences or not. default is 0, which
% means haven't smooth the binary images before, so the function will
% perform the smoothing using smoothboundaries function. If smooth is 1,
% then load the already saved smooth images.
%
% plotsetting: a 3*1 or 1*3 vector, [width szB szC], which are the setting 
% for plot merge images. Default is [2 5 5].
% width is the line width for 2 radius circles, default = 2;
% szB is the marker size for branch points, default = 5; 
% szC is the marker size for centroids, default = 5.
%
% Input_directory: input folder containing original images.
%
% Output_directory: output folder to save morphology results and images.
%
% Outputs:
% T_unfill: T*16 table, containing morphology features of unfilled
% spheroid. T is timepoints of the spheroid movie. 16 features labeled on
% top of each column, are: {'Ncell_in_spheroid', 'x_pixel', 'y_pixel',
% 'Rcore_pixel', 'Rinv_pixel', 'Hole_Number', 'Perimeter_pixel', 
% 'Area_pixel', 'Intensity', 'Intensity_STD', 'Complexity', 
% 'Eccentricity', 'Rcore_um', 'Rinv_um', 'Perimeter_um','Area_um2'}.
%
% S_unfill: T*16 matrix, same data as T_unfill, just without feature
% labels. Save for future data comparison and statistical tests.
%
% T_fill: T*16 table, containing morphology features of filled
% spheroid. T is timepoints of the spheroid movie. 16 features labeled on
% top of each column, are: {'x_pixel', 'y_pixel', 'Rcore_pixel', 
% 'Rinv_pixel', 'Perimeter_pixel', 'Area_pixel', 'Intensity', 
% 'Intensity_STD', 'Complexity', 'Eccentricity', 'Branch_Number', 
% 'Single_Cell_Number', 'Rcore_um', 'Rinv_um', 'Perimeter_um', 'Area_um2'}.
%
% S_fill: T*16 matrix, same data as T_fill, just without feature
% labels. Save for future data comparison and statistical tests.
%
% single: T*1 struct with four fields: area, centroids, intensity, and isd.
% Each field is a T*1 cell array, contains the area/x,y coordinates/
% mean intensity of the original fluorescent image at ZSliceName/
% std of intensity from projected image at all Z slices,
% for each single cell at timepoint T.
%
% This function will also automatically save all the smoothed and merge 
% images under Output_directory. Smoothed images are binary masks with
% smoother boundaries. Merge images are binary masks with invasive radius,
% core radius, branch points and centroids labeled in different color, as a
% reference to check the morphology quantification. All the results matrix
% and structure will be saved into a "morphology.mat" file. And the tables
% will be saved both in excel and text format.
%
% Embedded functions: smoothboundaries, embedboundaries, enlargepoint, lowb
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% check the inputs 
if ~exist('ZSliceName','var') || isempty(ZSliceName)
    error('Please specify input2: ZSliceName for original fluorescent images');
end

if ~exist('Resolution','var') || isempty(Resolution)
    error('Please specify input3: Resolution for the image (um/pixel)');
end

if ~exist('SingleCellArea','var') || isempty(SingleCellArea)
    display('Warning: Please specify input4, SingleCellArea in um^2, otherwise it will be 200 um^2');
    SingleCellArea = 200;
end

if ~exist('smooth','var') || isempty(smooth)
    display('Warning: Please specify input5, smooth or not, otherwise it will be 0 and smooth the image');
    smooth = 0;
end

if ~exist('plotsetting','var') || isempty(plotsetting)
    display('Warning: Please specify input6, plotsetting, otherwise it will be [2 5 5]');
    plotsetting = [2 5 5];
end

if ~exist('Input_directory','var') || isempty(Input_directory)
    Input_directory = uigetdir([],'Please Choose the Input Folder Containing Original Images');
end

if ~exist('Output_directory','var') || isempty(Output_directory)
    Output_directory = uigetdir([],'Please Choose the Output Folder for Saving Morphology Results');
end
%% main codes
% calculate threshold area and perimeter for single cell (pixel)
area1 = floor(SingleCellArea/Resolution^2);  % area
p1 = ceil(2 * pi * sqrt(area1/pi));  % perimeter

% load segmented masks
% if your segmented binary image is not named 'mask---.tif', change it
% according to your binary image's name. eg. if your segmented images
% are named as 'binary01.tif','binary02.tif',etc.
% change line 117 to:
% mask_list = dir([Output_directory filesep 'binary*.tif']);
mask_list = dir([Output_directory filesep 'mask*.tif']);
% load original 3d_projection images for intensity detection, similarly,
% change 'proj*.tif' in the following line according to the name of your 
% original grayscale/RGB image sequences.
projection_list = dir([Output_directory filesep 'proj*.tif']);

% load original fluorescent images at one z slice for intensity detection, 
% similarly, change the following line according to the name of your 
% original grayscale/RGB image sequences.
original_list = dir([Input_directory filesep '*' ZSliceName '.tif']);

% number of timepoints/frames/slices
Nt = length(mask_list);
% create pre-allocate matrix for saving all extracted features
S_unfill = zeros(Nt, 12);
S_fill = zeros(Nt, 12);
single_area = cell(Nt,1);
single_centroids = cell(Nt,1);
single_intensity = cell(Nt,1);
single_isd = cell(Nt,1);
mask1 = imread([Output_directory filesep mask_list(1).name]);
[X, Y]=meshgrid(1:size(mask1,1), 1:size(mask1,2));
%%
tic
for t = 1: Nt
    %% load binary segmented mask, projected and original z5 images
    mask = imread([Output_directory filesep mask_list(t).name]);
    projection = imread([Output_directory filesep projection_list(t).name]);
    original = imread([Input_directory filesep original_list(t).name]);
    %% postprocessing on mask, delete small debris and holes, smooth boundary
    % delete any noises or debris smaller than SingleCellArea threshold,
    % default area is 200 um^2
    I1 = bwareaopen(mask, area1); %figure,imshow(I1);title('I1');
    % smooth boundary of I1, fill all holes smaller than 45 pixel (~9um
    % radius circle)
    Ismooth_filename = [Output_directory filesep 'smooth' ...
        num2str(t, '%04g') '.tif'];
    % smooth =1, means already have Ismooth images, can just load them.
    if smooth == 1;
        Ismooth = imread(Ismooth_filename);
    else
        Ismooth = smoothboundaries(I1,p1); %figure,imshow(Ismooth);title('Ismooth');
        % remove any left small objects after smoothing, smaller than half size
        % of the SingleCellArea
        Ismooth = bwareaopen(Ismooth,round(area1/2)); %figure,imshow(Ismooth);title('Ismooth');
        % save Ismooth sequences in morphology folder, named with
        % 'smooth0001.tif', 'smooth0002.tif', ..., etc.
        imwrite(Ismooth, Ismooth_filename, 'tiff', 'Compression', 'None');
    end
    
    %% Classify Ismooth objects into 3 categories:
    % Isingle (only isolated single cells);
    % Ispheroid (only spheroid, with inside holes);
    % Ispheroid_fill (filled Ispheroid, no inside holes).
    
    % delete all isolated small parts except the spheroid
    a_smooth = regionprops(Ismooth, 'Area');
    thresh_area = max(cat(1,a_smooth.Area))-1;
    Ispheroid = bwareaopen(Ismooth, thresh_area); 
    %figure,imshow(Ispheroid);title('Ispheroid');
    % fill holes inside Ispheroid
    Ispheroid_fill = imfill(Ispheroid,'holes'); 
    %figure,imshow(Ispheroid_fill);title('Ispheroid fill');
    % create isolated single cell mask
    Isingle = imabsdiff(Ismooth, Ispheroid);  
    %figure,imshow(Isingle);title('Isingle');
    
    %% extract features for Ispheroid (unfilled)
    % core radius: the largest circle inside the filled spheroid. 
    % find the boundary of the filled spheroid
    B = bwboundaries(Ispheroid);
    % calculate the length for each cell in B
    length_B = cellfun(@numel,B)/2;
    % the maximum length_B is the spheroid boundary, bd_idx: index for the
    % spheroid boundary, p: perimeter for the spheroid_fill (unit:pixel).
    [p_fill, bd_idx] = max(length_B);
    % save all x,y coordinates on the spheroid boundary in boundary, column
    % 1: x, column2: y
    boundary = cell2mat(B(bd_idx));
    boundary = fliplr(boundary);
    % find the centroid of the spheroid
    c = regionprops(Ispheroid, 'Centroid');
    center = c.Centroid;
    % define core radius and invasive radius 
    % core radius Rcore = distance(nearest neighbour point on boundary, center)
    % invasive radius Rinv = distance(furthest neighbour point on boundary, center) 
    [~, R] = knnsearch(boundary, center, 'k', p_fill);
    Rcore = R(1);
    Rinv = R(end);
    % number of holes inside spheroid
    Nhole = numel(B) - 1;
    % perimeter of unfilled spheroid, including holes
    p_unfill = sum(length_B);
    % area of unfilled spheroid
    a_unfill = regionprops(Ispheroid, 'Area');
    area_unfill = a_unfill.Area;
    % intensity of unfilled spheroid from original
    i_unfill = regionprops(Ispheroid, original, 'MeanIntensity');
    intensity_unfill = i_unfill.MeanIntensity;
    % intensity of unfilled spheroid from projection
    isd_spheroid = regionprops(Ispheroid, projection, 'MeanIntensity');
    isd_unfill = isd_spheroid.MeanIntensity;
    % complexity of unfilled spheroid
    complexity_unfill = p_unfill^2/(4 * pi * area_unfill);
    % eccentricity of unfilled spheroid
    e_unfill = regionprops(Ispheroid, 'Eccentricity');
    eccentricity_unfill = e_unfill.Eccentricity;
    % approximate cell number in spheroid = area_unfill/area1;
    Nunfill = area_unfill/area1;
    % save all features in spheroid_unfill matrix
    S_unfill(t,:) = [Nunfill center Rcore Rinv Nhole p_unfill area_unfill...
        intensity_unfill isd_unfill complexity_unfill eccentricity_unfill];
    
    %% extract features for Ispheroid_filled
    % find the centroid of the filled spheroid
    c_fill = regionprops(Ispheroid_fill, 'Centroid');
    center_fill = c_fill.Centroid;
    % define core radius and invasive radius for filled spheroid
    % core radius Rcore = distance(nearest neighbour point on boundary, center_fill)
    % invasive radius Rinv = distance(furthest neighbour point on boundary, center_fill) 
    [~, R_fill] = knnsearch(boundary, center_fill, 'k', p_fill);
    Rcore_fill = R_fill(1);
    Rinv_fill = R_fill(end);
    % area of filled spheroid
    a_fill = regionprops(Ispheroid_fill, 'Area');
    area_fill = a_fill.Area;
    % intensity of filled spheroid from original
    i_fill = regionprops(Ispheroid_fill, original, 'MeanIntensity');
    intensity_fill = i_fill.MeanIntensity;
    % intensity of filled spheroid from projection
    isd_spheroid_fill = regionprops(Ispheroid_fill, projection, 'MeanIntensity');
    isd_fill = isd_spheroid_fill.MeanIntensity;
    % complexity of filled spheroid
    complexity_fill = p_fill^2/(4 * pi * area_fill);
    % eccentricity of filled spheroid
    e_fill = regionprops(Ispheroid_fill, 'Eccentricity');
    eccentricity_fill = e_fill.Eccentricity;
    % use thin to count the effective branch points of filled spheroid
    skel = bwmorph(Ispheroid_fill, 'thin', Inf);
    ep = bwmorph(skel, 'endpoints');
    % save all features in spheroid_unfill matrix
    S_fill(t,1:10) = [center_fill Rcore_fill Rinv_fill p_fill area_fill...
        intensity_fill isd_fill complexity_fill eccentricity_fill];
    
    %% extract features for single cells
    % if no single cell, save as an empty cell array for each feature
    if isequal(Ismooth, Ispheroid)
        single_area{t} = [];
        single_centroids{t} = [];
        single_intensity{t} = [];
        single_isd{t} = [];
        Nsingle = 0;
    else
        %% save each feature in a separate cell array
        % single cell area
        a_single = regionprops(Isingle, 'Area');
        single_area{t} = cat(1, a_single.Area);
        % single cell centroid
        c_single = regionprops(Isingle, 'Centroid');
        single_centroids{t} = cat(1, c_single.Centroid);
        % single cell intensity from original zslice fluorescent image
        int_single = regionprops(Isingle, original, 'MeanIntensity');
        single_intensity{t} = cat(1, int_single.MeanIntensity);
        % single cell intensity STD from projection image
        isd_single = regionprops(Isingle, projection, 'MeanIntensity');
        single_isd{t} = cat(1, isd_single.MeanIntensity);
        % single cell number
        Nsingle = length(a_single);
    end
    S_fill(t, 12) = Nsingle;
    
    %% plot the merged image with every features
    % load plotsetting for line width, marker size
    width = plotsetting(1);  % line width for 2 radius circles
    szB = plotsetting(2);    % marker size for branch points
    szC = plotsetting(3);    % marker size for centroids
    % green circle: core boundaries of filled spheroid
    core_circle = (Y-center_fill(2)).^2 + (X-center_fill(1)).^2 <= Rcore_fill^2;
    corebd = imdilate(bwperim(core_circle), strel('disk',width));
    merge1 = embedboundaries(mat2gray(Ispheroid),corebd, [0 1 0]);
    % red circle: invasive boundaries of filled spheroid
    inv_circle = (Y-center_fill(2)).^2 + (X-center_fill(1)).^2 <= Rinv_fill^2;
    invbd = imdilate(bwperim(inv_circle), strel('disk',width));
    merge2 = embedboundaries(merge1,invbd, [1 0 0]);
    % remove unreal branch endpoints, if an endpoint is within the core circle,
    % remove it.
    [epx, epy] = find(ep);
    position_idx = zeros(length(epx),1);
    for i = 1:length(epx)
        position_idx(i) = core_circle(epx(i),epy(i));
    end
    epx_real = epx(~position_idx);
    epy_real = epy(~position_idx);
    Nbranch = length(epx_real);
    S_fill(t,11) = Nbranch;
    if Nbranch ~= 0
        % make the branch points larger to visualize
        [epxlarge, epylarge] = enlargepoint(epx_real, epy_real, [szB szB]);
        % exclude out of field points
        outside_idx = find(epxlarge < 1 | epxlarge > size(mask,1) | ...
            epylarge < 1 | epylarge > size(mask,2));
        epxlarge(outside_idx) = [];
        epylarge(outside_idx) = [];
        % add real branch endpoints onto merge image as cyan dots
        for i = 1:length(epxlarge)
            merge2(epxlarge(i), epylarge(i), 1) = 0;
            merge2(epxlarge(i), epylarge(i), 2) = 1;
            merge2(epxlarge(i), epylarge(i), 3) = 1;
        end
    end
        
    % add centroids onto merge image as blue dots
    %center(1) = x, center(2)=y, but for saving image, flip x and y is the correct format
    % if using imshow(Ispheroid); hold on; plot(center(1),center(2),'rx') to
    % display, then center(1),center(2) is the correct format.
    [ylarge, xlarge] = enlargepoint(round(center_fill(2)), ...
        round(center_fill(1)), [szC szC]);
    merge2(ylarge, xlarge, 1) = 0;
    merge2(ylarge, xlarge, 2) = 0;
    merge2(ylarge, xlarge, 3) = 1; %figure,imshow(merge);
    
    merge_name = [Output_directory filesep 'merge' num2str(t,'%04g') '.tif'];
    imwrite(merge2, merge_name, 'tiff', 'Compression', 'None');
end
%% convert some features in pixel to um
% radius and perimeter
S_fill(:, 13:15) = S_fill(:, 3:5) .* Resolution;
S_unfill(:, 13:15) = S_unfill(:, [4 5 7]) .* Resolution;
% area
S_fill(:, 16) = S_fill(:, 6) .* (Resolution^2);
S_unfill(:, 16) = S_unfill(:, 8) .* (Resolution^2);
% convert results from matrix to table with variable names
T_fill = array2table(S_fill, 'VariableNames',{'x_pixel','y_pixel',...
    'Rcore_pixel','Rinv_pixel','Perimeter_pixel','Area_pixel',...
    'Intensity','Intensity_STD','Complexity','Eccentricity',...
    'Branch_Number','Single_Cell_Number','Rcore_um','Rinv_um',...
    'Perimeter_um','Area_um2'});
T_unfill = array2table(S_unfill, 'VariableNames',{'Ncell_in_spheroid',...
    'x_pixel','y_pixel','Rcore_pixel','Rinv_pixel','Hole_Number',...
    'Perimeter_pixel','Area_pixel','Intensity','Intensity_STD',...
    'Complexity','Eccentricity','Rcore_um','Rinv_um',...
    'Perimeter_um','Area_um2'});

%% convert cell arrays of single_cell into structure
single = struct('area',single_area,'centroids',single_centroids,...
    'intensity',single_intensity,'isd',single_isd);
% if you want to convert back from structure to cell array, use the
% following 2 line codes:
%test = struct2cell(single);
%test_area = test(1,:)';  %test_area = single_area;
%% save the results
save([Output_directory filesep 'morphology.mat'], 'S_unfill', 'S_fill',...
        'T_unfill', 'T_fill', 'single');
writetable(T_unfill, [Output_directory filesep 'unfilled_spheroid.csv']);
writetable(T_fill, [Output_directory filesep 'filled_spheroid.csv']);

disp('Morphology Feature Extraction:');
toc