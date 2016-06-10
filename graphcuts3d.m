function graphcuts = graphcuts3d(Output_directory, proj_adjust, alpha, u1, u2)
%
%   Function CMF3D
%
%   The matlab function to show how to use the functions CMF3D_Mex 
%
%   Before using the functions CMF3D_mex, you should compile it as follows:
%       >> mex CMF3D_mex.c
%
%   After compilation, you can define all the parameters (penalty, C_s, C_t, para) as follows: 
%   
%        - penalty: point to the edge-weight penalty parameters to
%                   total-variation function.
% 
%          For the case without incorporating image-edge weights, 
%          penalty is given by the constant everywhere. For the case 
%          with image-edge weights, penalty is given by the pixelwise 
%          weight function:
% 
%          for example, penalty(x) = b/(1 + a*| grad f(x)|) where b,a > 0 .
% 
%        - C_s: point to the capacities of source flows ps
% 
%        - C_t: point to the capacities of sink flows pt
% 
%        - para: a sequence of parameters for the algorithm
%             para[0,1,2]: rows, cols, heights of the given image
%             para[3]: the maximum iteration number
%             para[4]: the error bound for convergence
%             para[5]: cc for the step-size of augmented Lagrangian method
%             para[6]: the step-size for the graident-projection step to the
%                    total-variation function. Its optimal range is [0.1, 0.17].
% 
%
%       Example:
% 
%             >> [u, erriter, i, timet] = CMF3D_mex(single(penalty), single(Cs), single(Ct), single(para));
%
%             >> us = max(u, beta);  % where beta in (0,1)
%
%             >> figure, loglog(erriter,'DisplayName','erriterN');figure(gcf)
%
%             >> isosurface(u,0.5), axis([1 rows 1 cols 1 heights]), daspect([1 1 1]);
%
%
%   Please email Jing Yuan (cn.yuanjing@gmail.com) for any questions, 
%   suggestions and bug reports
%
%   The Software is provided "as is", without warranty of any kind.
%
%               Version 1.0
%   https://sites.google.com/site/wwwjingyuan/       
%
%   Copyright 2011 Jing Yuan (cn.yuanjing@gmail.com)   
%
%%%%Addition by Yue Hou 2016 <lotushouyue@hotmail.com>
% Embedded codes: CMF3D_mex.c
% Inputs:
% Output_directory: folder for saving outputs. Don't change it.
%
% proj_adjust (or proj_3dgauss if you don't have proj_adjust): 
% input 3d image matrix in (x,y,t) or (x,y,z) for segmentation.
% proj_adjust (or proj_3dgauss) are outputs from function preprocess. 
% Don't change it.
%
% alpha: scalar within [0,1]. Constant for penalty, smaller the better. 
% Default is 0.2.
%
% u1 and u2: lower and upper limit for ulab (or u(x) in the Readme file),
% u(x)= [0, 1]. Smaller the better segmentation. Default is u1=0.2, u2=0.7.
% 

% Outputs: 
% graphcuts: grayscale 3d image matrix after graphcuts segmentation, single
% format. M*N*T, [M,N] is image size, T is time point.
%
% The function will also automatically generate a
% 'graphcuts_parameters.txt' file to save all the alpha, u1, u2 for future
% refining/adjustment of the graphcuts segmentation.

%% check input parameters, if not exist, pop-up an error or a warning 
if ~exist('Output_directory','var') || isempty(Output_directory)
    Output_directory = uigetdir([],'Please Choose the Output Folder for Saving Graphcuts Results');
end

if ~exist('proj_adjust','var') || isempty(proj_adjust)
    proj_adjust_filename = uigetfile([Output_directory filesep '*.mat'],...
        'Please Select the Input 3D Matrix, either 3d gaussian filtered or adjusted results');
    load([Output_directory filesep proj_adjust_filename]);
end

if ~exist('alpha','var') || isempty(alpha)
    disp('Warning: Please specify alpha for penalty constant, otherwise it will be 0.2');
    alpha = 0.2;
end

if ~exist('u1','var') || isempty(u1)
    disp('Warning: Please specify u1 for lower limit of ulab, otherwise it will be 0.2');
    u1 = 0.2;
end

if ~exist('u2','var') || isempty(u2)
    disp('Warning: Please specify u2 for upper limit of ulab, otherwise it will be 0.7');
    u2 = 0.7;
end

%% 3D graphcuts
tic
% graphcuts variables
[rows,cols,heights] = size(proj_adjust);

varParas = [rows; cols; heights; 300; 5e-4; 0.2; 0.11];
%                para 0,1,2 - rows, cols, heights of the given image.
%                para 3 - the maximum number of iterations, default=200~300.
%                para 4 - the error bound for convergence, default=5e-4.
%                para 5 - cc for the step-size of augmented Lagrangian
%                method, default=0.2~0.35.
%                para 6 - the step-size for the graident-projection of p,
%                default=0.11.

penalty = alpha*ones(rows, cols, heights);

% build up the priori L_2 data terms
fCs = abs(proj_adjust - u1);
fCt = abs(proj_adjust - u2);

% ----------------------------------------------------------------------
%  Use the function CMF3D_mex to run the algorithm on CPU
% ----------------------------------------------------------------------

[graphcuts, ~, ~, ~] = CMF3D_mex(single(penalty), single(fCs), single(fCt), single(varParas));
% save the output matrix under graphcuts_path
save([Output_directory filesep 'graphcuts.mat'],'graphcuts');

%% write all the parameters to a .txt file
% Open or create new text file for reading and writing. Discard existing contents, if any.
fid = fopen([Output_directory filesep 'graphcuts_parameters.txt'],'wt+');
% save the parameters value for rerun or refine segmentation
fprintf(fid,'alpha = %f;\n', alpha);
fprintf(fid,'u1 = %f;\n', u1);
fprintf(fid,'u2 = %f;\n', u2);
fclose(fid);

% display function elapse time
disp('3D Graphcuts:');
toc