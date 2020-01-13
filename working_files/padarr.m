function [padded] = padarr(array, dims, method)
%PADARR Summary of this function goes here
%   Detailed explanation goes here
    
% Pad arrays
    if isvector(array) == 1
        % Check if vector is vertical
        if size(array,1) > 1
            if strcmp(method,'pre')
                padded = [zeros(dims,1);array];
            elseif strcmp(method,'post')
                padded = [array;zeros(dims,1)];
            end
        else
            if strcmp(method,'pre')
                padded = [zeros(dims,1),array];
            elseif strcmp(method,'post')
                padded = [array,zeros(dims,1)];
            end
        end
    % Pad matrices
    else
        if strcmp(method,'pre')
            padded = [zeros(dims(1),size(array,2));array];
            padded = [zeros(size(padded,1),dims(2)),padded];
        elseif strcmp(method,'post')
            padded = [array;zeros(dims(1),size(array,2))];
            padded = [padded,zeros(size(padded,1),dims(2))];
        end
    end

end

