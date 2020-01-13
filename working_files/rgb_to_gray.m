function gray_image = rgb_to_gray(input_image)
    % Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
    % das Bild bereits in Graustufen vorliegt, soll es direkt zurueckgegeben werden.
    dims = size(input_image);
    if numel(dims) == 2
        gray_image = input_image;
    else
        m = dims(1);
        n = dims(2);
        rgb_vec = [0.299, 0.587, 0.114];
        gray_image = zeros(m,n);
        for i = 1:m
            for j = 1:n           
                vec = double(squeeze(input_image(i,j,1:3)));
                gray_image(i,j) = dot(rgb_vec, vec);
            end
        end
        gray_image = uint8(gray_image);
    end
end