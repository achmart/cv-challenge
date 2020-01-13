function [out] = denoise_disparity(img,width2,sigma2,k,maxIter,lambda,range)
%DENOISE_DISPARITY Summary of this function goes here
%   Detailed explanation goes here
% https://www.utdallas.edu/~cxc123730/SPIE2016.pdf

img_min = range(1);
img_max = range(2);

% Reset values that don't make sense
img(img < img_min) = img_min;
img(img > img_max) = img_max;

%img = double(img);
img = double(img);
%img = double(rescale(img,0,255));

%% Denoising

% 1) apply median filter
%img = medfilt2(img,[3 3]);
%img = medfilt2(img,[9 9]);

% 2) Apply gaussian blur
%size = 5;
%W1 = gauss_filter(gauss_width1,sigma1);



%img = conv2(img,W1,'same');


%% Hole filling

% 1) Morpholocial closing
%se = strel('disk',3);
%img = imclose(img,se);

%figure
%imshow(uint8(img));
%title('After denoising before hole-filling');

min(img(:))
max(img(:))

% 2) Hole search and filling

mask = zeros(k,k);
% dynamic threshold that increases with each iteration
threshold = 1;

for it=1:maxIter    
    %W = gauss_filter(gauss_width,sigma);
    %img = conv2(img,W,'same');
    for i=k+1:size(img,1)-k
        for j=k+1:size(img,2)-k
            % If pixel is zero
            if(img(i,j) < threshold)
                % Get neighbor pixels
                mask = img((i-k):(i+k),(j-k):(j+k));
                % Search non-zero pixels in mask
                if(max(max(mask)) < threshold)
                    img(i,j) = 0;
                else
                    % Search the maximal value m of all elements in mask
                    % replace zero pixels with m
                    %m = ones(2*k+1,2*k+1) .* max(mask(:));
                    %m = m .* W_mask;
                    m = max(mask(:));
                    img((i-k):(i+k),(j-k):(j+k)) = m;
                end
            end
        end
    end
end

%W2 = gauss_filter(width2,sigma2);
%img = conv2(img,W2,'same');
%img = wlsFilter(img,lambda,1.2);
out = rescale(img,img_min,img_max);
end