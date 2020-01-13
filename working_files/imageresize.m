function img_out=imageresize(img_in, rowScale, colScale,method)
%IMAGERESIZE Resize image 
%   IMAGERESIZE resize an image to any scale. 
%   This is a simple implementation of IMRESIZE.

%   Structure of function from:
%   https://de.mathworks.com/matlabcentral/fileexchange/30787-image-resize

if nargin<3
    colScale = rowScale;
end

%% set axis scaling
xscale=colScale;
yscale=rowScale;

sz=size(img_in);
numchans=size(img_in,3);
yy=linspace(1,sz(1),ceil(sz(1)*yscale));
xx=linspace(1,sz(2),ceil(sz(2)*xscale));

img_out=zeros([ceil(sz(1:2).*[yscale xscale]) numchans]);

% Interpolate image with new size
parfor c=1:numchans
    img_out(:,:,c)=interp2(double(img_in(:,:,c)),xx',yy,method);
    %img_out(:,:,c)=interp2(double(img_in(:,:,c)),xx',yy,'bicubic');
end
%%%%%%%%%%%%%%

%img_out = uint8(img_out);

% imshow(img_in); 
% title('original image');
% figure;imshow(img_out); 
% title('resized image');
% truesize(size(img_in)); % omit this to display as it is

    