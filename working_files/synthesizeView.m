function [synth_L,synth_R] = synthesizeView(rect_L,rect_R,disp12,disp21,p)
%SYNTHESIZEVIEW Summary of this function goes here
%   Detailed explanation goes here

[rows,cols,~] = size(rect_L);
synth_L = zeros(rows,cols,3);
synth_R = zeros(rows,cols,3);

for i=1:rows
    for j=1:cols
       % calculate new column (x-coordinate) for current pixel
       newColLeft = round(j - p * disp12(i,j));
       newColRight = round(j - (p-1) * disp21(i,j));
       
       % Pixel warping
       if(newColLeft >= 1 && newColLeft <= cols)
          synth_L(i,newColLeft,:) = rect_L(i,j,:); 
       end
       if(newColRight >= 1 && newColRight <= cols)
          synth_R(i,newColRight,:) = rect_R(i,j,:); 
       end
    end 
end


end

