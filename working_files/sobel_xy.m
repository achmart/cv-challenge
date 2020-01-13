function [Fx, Fy] = sobel_xy(input_image)
    % In dieser Funktion soll das Sobel-Filter implementiert werden, welches
    % ein Graustufenbild einliest und den Bildgradienten in x- sowie in
    % y-Richtung zurueckgibt.
    sobelX = [1 0 -1; 2 0 -2; 1 0 -1];
    sobelY = [1 2 1; 0 0 0; -1 -2 -1];
    Fx = conv2(input_image, sobelX, 'same');
    Fy = conv2(input_image, sobelY, 'same');    
end
