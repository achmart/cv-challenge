function [out1,out2] = postprocess(img,iter,medianSize,method,enableWLS,lambda,alpha)
%POSTPROCESS Summary of this function goes here
% https://de.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
% img:          input image 3 channel MxNx3
% medianIter:   number of median iterations
% medianSize:   Scalar for building the MxM median window
% method: 1-5

%% WARNING: medfilt2 needs Image Processing Toolbox (should not be used!)
%

%% Code

% 2) Convert missing values to NaN
img = double(img);

for i=1:size(img,1)
    for j=1:size(img,2)
        if(sum(img(i,j,:)) < 1.0)
            img(i,j,:) = [NaN, NaN, NaN];
        end
    end
end

% 3) Median Filtering
if(iter > 0)
    med_range = [medianSize medianSize];
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);
    %parfor i=1:iter
     %   R = medfilt2(R,med_range);
      %  G = medfilt2(G,med_range);
       % B = medfilt2(B,med_range);
    %end
    img = cat(3,R,G,B);
end

% 4) Call inpaint_nans on all RGB color channels
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);
clear img
R = inpaint_nans(R,method);
G = inpaint_nans(G,method);
B = inpaint_nans(B,method);
out1 = uint8(cat(3,R,G,B));

% 5) Call wlsFilter to smoothen out more
if(enableWLS)
    R = wlsFilter(R,lambda,alpha);
    G = wlsFilter(G,lambda,alpha);
    B = wlsFilter(B,lambda,alpha);
end
out2 = uint8(cat(3,R,G,B));
end

