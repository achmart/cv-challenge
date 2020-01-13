function synth_image = mergeImages(I_L,I_R,p)
%MERGEIMAGES Summary of this function goes here
%   Detailed explanation goes here

%holes1 = I_L == 0;
%holes2 = I_R == 0;

if(p <= 0.5)
    synth_image = I_L;
    %synth_image(holes1) = I_R(holes1);
else
    synth_image = I_R;
    %synth_image(holes2) = I_L(holes2);
end

end