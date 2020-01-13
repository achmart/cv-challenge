function [I_moved] = immove(Image, translation)
%IMMOVE Summary of this function goes here
%   Detailed explanation goes here
    x = translation(1);
    y = translation(2);
    sz = size(Image);
    I_moved = zeros(sz);
    I_moved(1+y:sz(1),1+x:sz(2)) = Image(1:(sz(1)-y),1:(sz(2)-x));
end

