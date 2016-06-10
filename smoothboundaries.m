function Ismooth = smoothboundaries(I, ThreshLength)
% use lowb function to smooth the boundaries of image I, including all 
% the inside holes larger than ThreshLength, and then returned the 
% smoothed binary image as Ismooth
%
% Inputs:
% I: M*N binary (or logical) image, objects is white and background is
% black.
%
% ThreshLength: threshold for eliminating small holes inside main object,
% if the perimeter of a hole is smaller than ThreshLength, it will be 
% filled rather than be smoothed. Only the large holes will be smoothed.
% Minimum ThreshLength = 10 pixel. 
%
% Outputs:
% Ismooth: M*N binary image after smoothing the boundaries (including
% holes).
%
% Embedded functions: lowb.m
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

%% check the formats of inputs
% check image I, must exist and be logical (binary)
if ~exist('I','var') || isempty(I)
    error('Please specify input1: binary M*N image for smoothing');
else if ~isa(I, 'logical')
        error('Input image I must be binary (logical)');
    end
end
% check the ThreshLength, make sure it exists and >= 10 pixel.
if ~exist('ThreshLength','var')||isempty(ThreshLength)
    disp('Warning: Please specify ThreshLength, otherwise ThreshLength=10')
    ThreshLength = 10; %default is 10 pixel
else if ThreshLength < 10
        disp('Warning: Minimum ThresLength should be 10 pixel! Set ThreshLength=10')
        ThreshLength = 10; %default is 10 pixel
    end
end

%% smoothing boundary
B = bwboundaries(I);   % trace boundaries for all objects, including holes
len_B = cellfun(@length, B);
len_B(:,2) = 1:length(B);
len_B = sortrows(len_B, -1);
Bsort = B(len_B(:,2));
% exclude small holes shorter than 10 pixels for perimeter, because lowB
% has a frequency requirement for longer than 10 pixels vector
sh_ind = len_B(:,1) <= ThreshLength;
Bsort(sh_ind) = [];

% plot the smoothed outside boundaries and fill all the holes
max_B = Bsort{1};  % the longest boundary is the cell boundary
% lowB is to smooth the sharp curvature of the boundary, see fig_lowB
b(:,2) = lowb(max_B(:,2));
b(:,1) = lowb(max_B(:,1));
if ceil(max(b(:))) > max(size(I))
    indh = find(ceil(b(:)) > max(size(I)));
    b(indh) = max(size(I));
end
if floor(min(b(:))) < 1
    indl = find(floor(b(:)) < 1);
    b(indl) = 1;
end
img1 = false(size(I));
img1(sub2ind(size(I), round(b(:,1)), round(b(:,2)))) = true;
img1(sub2ind(size(I), floor(b(:,1)), floor(b(:,2)))) = true;
img1(sub2ind(size(I), ceil(b(:,1)), ceil(b(:,2)))) = true;
img1 = imclose(img1,strel('disk',3));
s1 = imfill(img1,'holes');

% plot the smoothed inside holes 
img2 = false(size(I));
for i = 2:length(Bsort)
    holebd = [];
    holebd(:,2) = lowb(Bsort{i}(:,2));
    holebd(:,1) = lowb(Bsort{i}(:,1));
    if ceil(max(holebd(:))) > max(size(I))
        indh = find(ceil(holebd(:)) > max(size(I)));
        holebd(indh) = max(size(I));
    end
    if floor(min(holebd(:))) < 1
        indl = find(floor(holebd(:)) < 1);
        holebd(indl) = 1;
    end
    img2(sub2ind(size(I), round(holebd(:,1)), round(holebd(:,2)))) = true;
    img2(sub2ind(size(I), floor(holebd(:,1)), floor(holebd(:,2)))) = true;
    img2(sub2ind(size(I), ceil(holebd(:,1)), ceil(holebd(:,2)))) = true;
    img2 = imclose(img2,strel('disk',2));
    %s = imfill(img,'holes');
end
s2 = imfill(img2,'holes'); 
Ismooth = s1 - s2; 
Ismooth(Ismooth<0)=1;
Ismooth = im2bw(Ismooth,0.5);