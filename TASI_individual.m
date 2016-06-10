%% Data Format: Provide the input and output directory, 
% see "Data Format" part of the instruction document.
Input_directory = uigetdir([],'Please Choose the Input Folder Containing Original Images');
Output_directory = uigetdir([],'Please Choose the Output Folder for Saving Cell Line Results');

%% Projection: import the original timelapse,Z-stacks images and perform 
% 3D projection to enhance the contrast for easier segmentation.
% see "3D Projection" part of the instruction document.
projection = projectionzstacks(7,2,1,[],Input_directory,Output_directory);

%% Preprocess: 3D gaussian filter and/or enhancing contrast
% 3d gaussian filter on 3d image matrix to smooth in time domain and remove
% noisies on background. Then use imadjust to enhance contrast if needed.
% see "Preprocess" part of the instruction document.
[proj_3dgauss,proj_adjust] = preprocess(Output_directory,projection,[2.5,2.5,0.5],1,1);

%% 3D Graphcuts: segmentation
% see "3D Graphcuts" part of the instruction document.
graphcuts = graphcuts3d(Output_directory, proj_adjust, 0.1, 0, 0.7);
% optional: display one graphcuts image in the middle of time sequence, to
% check the segmentation results and help to determine the threshold
% settings
t = size(graphcuts,3);
i_test = round(t/2);
I_test = imread([Output_directory filesep 'proj' num2str(i_test, '%04g') '.tif']);
figure,subplot(121),imshow(imadjust(I_test)); title('projected image');
subplot(122),imshow(graphcuts(:,:,i_test),[]); title('graphcuts image');
linkaxes;

%% Ostuthresh: postprocessing using ostu threshold, convert 3d matrix to 
% binary image sequences, see "Ostuthresh" part of the instruction document.
ostuthresh(Output_directory, graphcuts, [], 100, 1);
% draw boundaries of the segmented images for visualization
drawboundaries(Output_directory, 1,'boundary.avi',7);

%% morphology feature extraction
[T_unfill, S_unfill, T_fill, S_fill, single] = featureextraction...
    ('z3_ch00',1.13,200,0,[],Input_directory,Output_directory);
% resolution = 1.13 for leader and follower; 2.27 for parental.