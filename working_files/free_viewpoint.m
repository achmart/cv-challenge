function [output_image]  = free_viewpoint(image1, image2, p, img_src, debug)
% This function generates an image from a virtual viewpoint between two
% real images. The output image has the same size as the input images.

if debug == 1
    do_plot = true;
else
    do_plot = false;
end

%% Return original images for p = {0, 1}
warning('off', 'Images:initSize:adjustingMag');

if (p == 0)
    output_image = image1;
    return;
end

if (p == 1)
    output_image = image2;
    return;
end

%% Determine correct K matrix

% Calibration matrix obtained using Matlab Camera calibrator on Schach_1
K1 = [2907.774829899128 0 1501.912926644449;
       0.0 2907.741631813634 990.5798681074887;
       0 0 1];
   
% Calibration matrix obtained using Matlab Camera calibrator on Schach_2
K2 = [3940.92356754738 0 1507.04990870178;
      0 3942.81592663649 988.776774304098;
      0 0 1];

K = [];

if(img_src == 1)
    K = K1;
elseif(img_src == 2)
    K = K2;
end
clear K1 K2

%% Convert images to grayscale

disp('Converting images to grayscale...');
tic
I_L = rgb2gray(image1);
I_R = rgb2gray(image2);
toc
disp('Finished loading images I_L & I_R');

%% Extract features using the Harris Detector

min_dist = 50;          % default: 20
N = 5;                 % default: 5
tile_size = [100,100];  % default: [200, 200]

tic
disp('Harris detector running...');
[merkmale_L,winkel_L] = harris_detektor(I_L,'segment_length',15,'k',0.06,...
    'tau',2e6,'min_dist',min_dist,'N',N,'tile_size',tile_size,'do_plot',do_plot);
[merkmale_R,winkel_R] = harris_detektor(I_R,'segment_length',15,'k',0.06,...
    'tau',2e6,'min_dist',min_dist,'N',N,'tile_size',tile_size,'do_plot',do_plot);
disp('Finished Feature Extraction for images I_L & I_R.');
toc
clear min_dist N tile_size

%% Find stereo correspondences in both images 
tic
disp('Finding correspondences...')
Korr = punkt_korrespondenzen(I_L,I_R,merkmale_L,merkmale_R,...
    winkel_L,winkel_R,false,'window_length',25,'min_corr', 0.73,'do_plot',do_plot);
toc
clear merkmale_L merkmale_R winkel_L winkel_R
% default window_length: 25

%%  Find robust Correspondence-pairs using the RANSAC-Algorithm

tic
disp('Ransac...');
Korr_robust = F_ransac(Korr,'tolerance', 5, 'epsilon',0.75,'p',0.8); % 40-200 tol seems ok
%Korr_robust = Korr;
toc
clear Korr

%% Show the robust Correspondence-pairs
if debug == 1
    figure
    imshow(I_L)
    alpha 0.5
    hold on
    imshow(I_R)
    alpha 0.5
    plot(Korr_robust(1,:), Korr_robust(2,:),'or')
    plot(Korr_robust(3,:), Korr_robust(4,:),'xg')

    for i=1:size(Korr_robust,2)
        plot(Korr_robust([1, 3],i), Korr_robust([2,4],i),'-')
    end

    hold off
    title('Correspondences selected by RanSaC');
end
clear I_L I_R


%% Determine E using the 8-point algorithm

disp('Achtpunktalgorithmus...');
E = achtpunktalgorithmus(Korr_robust, K);
disp('Successfully determined the Euclidian transformation E.');
disp(E);

%% Get rotation and translation using E
[T1, R1, T2, R2] = TR_aus_E(E);

%[T, R, lambda, P1, ~, ~] = rekonstruktion(T1, T2, R1, R2, Korr_robust, K);
[T, R, ~, ~, ~, ~] = rekonstruktion(T1, T2, R1, R2, Korr_robust, K, do_plot);
clear Korr_robust E T1 R1 T2 R2
if debug == 1
    title('Your Solution')
end

disp('Rotation matrix R:')
disp(R)

disp('Translation vector T:')
disp(T)

%% Use E and percentage p to interpolate the E for free-viewpoint
% Code has been changed such that Quaternion.m is used instead

% Interpolate the Rotation
[R_int,T_int] = interpol_R_T(R,T,p);


%% Rectify Stereo Images

disp('Rectifying Stereo Images...');
tic
[I1Rect,I2Rect,~,~, Transform_int] = RectifyImages(image1,image2,K,R,T,...
    R_int,T_int, do_plot);
toc
clear image1 image2 T R K T_int R_int


%% Compute disparity
% Warning: Using Matlab CV Toolbox is not allowed!

disp('Computing disparity...');
min_disp = 600;
max_disp = 1200;

if(img_src == 1)
    min_disp = 1100;
    max_disp = 2400;
elseif(img_src == 2)
    min_disp = 600;
    max_disp = 1200;
end

% Correct if difference is not divisible by 16 (required by Toolbox func.)
if(mod(max_disp-min_disp,16) ~= 0)
    max_disp = ceil((max_disp-min_disp) / 16) * 16 + min_disp;
end

dRange = [min_disp max_disp];

sc = 4; % scale factor
w_radius = 2; % 4 -> 9x9 // 3 -> 7x7 // 2 -> 5x5 // 1 -> 3x3

I1small = uint8(imageresize(I1Rect,1/sc,1/sc,'bicubic'));
I2small = uint8(imageresize(I2Rect,1/sc,1/sc,'bicubic'));

if debug == 1
    figure
    imshow(I1small)
    title('I1 small');

    figure
    imshow(I2small)
    title('I2 small');
end

[d12,~,~] = stereo_sad(I1small,I2small,min_disp/sc,max_disp/sc,w_radius);

% http://answers.opencv.org/question/9026/in-stereo-matching-is-the-disparity-always-to-the-left/
% Flip the input images
a = fliplr(I2small);
b = fliplr(I1small);
[d21,~,~] = stereo_sad(a,b,min_disp/sc,max_disp/sc,w_radius);

% Correct disparity values
d12 = d12*sc;
d21 = d21*sc;

% Show raw images
if debug == 1
    figure
    imshow(imageresize(d12,sc,sc,'bicubic'),dRange);
    title('RAW Disparity Map (left) d12');
    colorbar

    figure
    imshow(imageresize(fliplr(d21),sc,sc,'bicubic'),dRange);
    title('RAW Disparity Map (right) d21');
    colorbar
end


%% Post-processing disparity maps
disp('Post-processing disparity maps...');
tic
lambda = 2.0;
alpha2 = 1.2;

if(img_src == 1)
    lambda = 0.2;
    alpha2 = 0.2;
elseif(img_src == 2)
    lambda = 0.4;
    alpha2 = 1.2;
end
d12_min = min(min(d12));
d12_max = max(max(d12));
d21_min = min(min(d21));
d21_max = max(max(d21));
d12 = wlsFilter(d12,lambda,alpha2);
d21 = wlsFilter(d21,lambda,alpha2);
d12 = rescale(d12,d12_min,d12_max);
d21 = rescale(d21,d21_min,d21_max);
toc

% Flip the second result
d21 = fliplr(d21);

% Resize to original size
d12 = imageresize(d12,sc,sc,'bicubic');
d21 = imageresize(d21,sc,sc,'bicubic');

if debug == 1
    figure
    imshow(d12,dRange);
    title('FILTERED Disparity Map (left) d12');
    colorbar

    figure
    imshow(d21,dRange);
    title('FILTERED Disparity Map (right) d21');
    colorbar
end


%% Synthesize new view

disp('Synthesizing new view...');
tic
[synth_L,synth_R] = synthesizeView(I1Rect,I2Rect,d12,d21,p);
toc
clear I1Rect I2Rect d12 d21

if debug == 1
    figure
    imshow(uint8(synth_L));
    title('Synth L');

    figure
    imshow(uint8(synth_R));
    title('Synth R');
end

%% Post-process both images

% disp('Post-processing 1 of the 2 shifted views...');
% tic
% [synth_L,~] = postprocess(synth_L,0,0,4,false,0,0); 
% [synth_R,~] = postprocess(synth_R,0,0,4,false,0,0);
% toc
% 
% figure
% imshow(uint8(synth_L));
% title('Post-processed Synth L');
% 
% figure
% imshow(uint8(synth_R));
% title('Post-processed Synth R');

%% Merge views

disp('Merging views...')
tic
synth_image = mergeImages(synth_L,synth_R,p);
toc
clear synth_L synth_R

if debug == 1
    figure
    imshow(uint8(synth_image));
    title('Synthesized image');
end

%% Post-process merged views

% find middle pixels
[c1,~,~] = size(synth_image);

% Approach 2: Find first colored pixel in middle image
cx = floor(c1/2);
middle = synth_image(cx,:,:);
middle = sum(middle,3);
min_y = find(middle,1,'first');

%% Post-process before warp does work

disp('post-processing synth _ image');
tic
scl = 2;
synth_small = imageresize(synth_image,1/scl,1/scl,'nearest');

figure
imshow(uint8(synth_small));
title('synth small');

[synth_small,~] = postprocess(synth_small,0,0,4,false,0,0);
synth_image = imageresize(synth_small,scl,scl,'bicubic');

%synth_image = postprocess(synth_image,0,0,4,false,0,0);
toc

%% De-Rectification
% Comupte the final derectified image I3 from I3' using inverse
% Homography inv(H3)

disp('Image De-Rectification...');
tic
parfor c = 1:size(synth_image, 3)
    [I3(:,:,c),~,~] = imagewarp(synth_image(:,:,c),inv(Transform_int), 'bilinear','valid');
end
toc
if debug == 1
    figure
    imshow(I3);
end
clear synth_image

%% Cropping to 3000x2000 image

disp('Cropping image...');

out = I3(cx-999:cx+1000, min_y:min_y+2999,:);
clear I3

if debug == 1
    figure
    imshow(out);
    title('Cropped image');
end

%% Post-process after warping does not work!
% (prob. because of imwarp / datatypes)
out = postprocess(out,0,0,4,false,2,1.2);

%%
output_image = uint8(out);
end